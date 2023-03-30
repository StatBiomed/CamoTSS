import numpy as np
import pandas as pd
import os
import pickle
from functools import reduce
import anndata as ad
import multiprocessing 
import pysam
from sklearn.cluster import AgglomerativeClustering
from scipy.optimize import linear_sum_assignment
import time
import random
import pickle
import statistics
import editdistance
import warnings
from pathlib import Path
import sys
from brie.utils import fetch_reads



warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=Warning)







global CACHE_CHROM
global CACHE_SAMFILE
CACHE_CHROM = None
CACHE_SAMFILE = None

def check_pysam_chrom(samFile,chrom=None):
    """
    check samFile is a file name or a pysam object, and if chrom format.
    """
    global CACHE_CHROM
    global CACHE_SAMFILE

    if CACHE_CHROM is not None:
        if (samFile==CACHE_SAMFILE) and (chrom==CACHE_CHROM):
            return CACHE_SAMFILE,CACHE_CHROM

    
    if type(samFile)==str or type(samFile)==np.str_:
        ftype=samFile.split(".")[-1]
        if ftype != "bam" and ftype!="sam" and ftype!="cram":
            print("Error: file type need suffix of bam, sam or cram")
            sys.exit(1)
        if ftype == "cram":
            samFile=pysam.AlignmentFile(samFile,'rc')
        elif ftype=="bam":
            samFile=pysam.AlignmentFile(samFile,'rb')
        else:
            samFile=pysam.AlignmentFile(samFile,'r')

    if chrom is not None:
        if chrom not in samFile.references:
            if chrom.startswith("chr"):
                chrom=chrom.split('chr')[1]
            else:
                chrom="chr"+chrom
        if chrom not in samFile.references:
            print("Can't find reference %s in samFile"%chrom)
            return samFile,None

    CACHE_CHROM=chrom
    CACHE_SAMFILE=samFile
    return samFile,chrom



def get_fastq_file(fastqFilePath):
    fastqFile=pysam.FastaFile(fastqFilePath)
    return fastqFile




class get_TSS_count():



    def __init__(self,generefPath,tssrefPath,bamfilePath,fastqFilePath,outdir,cellBarcodePath,nproc,minCount,maxReadCount,clusterDistance,InnerDistance,windowSize,minCTSSCount,minFC):
        self.generefdf=pd.read_csv(generefPath,delimiter='\t')
        self.generefdf['len']=self.generefdf['End']-self.generefdf['Start']
        self.tssrefdf=pd.read_csv(tssrefPath,delimiter='\t')
        self.bamfilePath=bamfilePath
        self.count_out_dir=str(outdir)+'/count/'
        if not os.path.exists(self.count_out_dir):
            os.mkdir(self.count_out_dir)

        self.ctss_out_dir=str(outdir)+'/CTSS/'
        if not os.path.exists(self.ctss_out_dir):
            os.mkdir(self.ctss_out_dir)


        self.minCount=minCount
        self.cellBarcode=pd.read_csv(cellBarcodePath,delimiter='\t')['cell_id'].values
        self.nproc=nproc
        self.maxReadCount=maxReadCount
        self.clusterDistance=clusterDistance
        self.fastqFilePath=fastqFilePath
        self.InnerDistance=InnerDistance
        self.windowSize=windowSize
        self.minCTSSCount=minCTSSCount
        self.minFC=minFC

        

    def _getreads(self,bamfilePath,fastqFilePath,geneid):
        #fetch reads1 in gene 
        samFile, _chrom = check_pysam_chrom(bamfilePath, str(self.generefdf.loc[geneid]['Chromosome']))
        reads = fetch_reads(samFile, _chrom,  self.generefdf.loc[geneid]['Start'] , self.generefdf.loc[geneid]['End'],  trimLen_max=100)
        reads1_umi = reads["reads1"]



        #select according to GX tag and CB (filter according to user owned cell)
        reads1_umi=[r for r in reads1_umi if r.get_tag('GX')==geneid]
        reads1_umi=[r for r in reads1_umi if r.get_tag('CB') in self.cellBarcode]

        #filter strand invasion
        fastqFile=get_fastq_file(fastqFilePath)
        reads1_umi=[r for r in reads1_umi if editdistance.eval(fastqFile.fetch(start=r.reference_start-14, end=r.reference_start-1, region='chr'+str(self.generefdf.loc[geneid]['Chromosome'])),'TTTCTTATATGGG') >3 ]



        #filter according to the cateria of SCAFE
        if self.generefdf.loc[geneid]['Strand']=='+':
            reads1_umi=[r for r in reads1_umi if r.is_reverse==False]
            reads1_umi=[r for r in reads1_umi if editdistance.eval(r.query_sequence[9:14],'ATGGG')<=4]
            reads1_umi=[r for r in reads1_umi if len(r.cigartuples)>=2]
            reads1_umi=[r for r in reads1_umi if (r.cigartuples[0][0]==4)&(r.cigartuples[0][1]>6)&(r.cigartuples[0][1]<20)&(r.cigartuples[1][0]==0)&(r.cigartuples[1][1]>5)]
        
        elif self.generefdf.loc[geneid]['Strand']=='-':
            reads1_umi=[r for r in reads1_umi if r.is_reverse==True]
            reads1_umi=[r for r in reads1_umi if editdistance.eval(r.query_sequence[-13:-8],'CCCAT')<=4]
            reads1_umi=[r for r in reads1_umi if len(r.cigartuples)>=2]
            reads1_umi=[r for r in reads1_umi if (r.cigartuples[0][0]==0)&(r.cigartuples[0][1]>5)&(r.cigartuples[1][0]==4)&(r.cigartuples[1][1]>6)&(r.cigartuples[1][1]<20)]



        #store start of reads and CB 
        reads_info=[]
        reads_info=[(r.positions[0],r.get_tag('CB'),r.cigarstring) for r in reads1_umi]
        #,r.query_sequence[:]
        

        return reads_info
    

        
    def _get_gene_reads(self):

        pool = multiprocessing.Pool(processes=self.nproc)

        self.generefdf.set_index('gene_id',inplace=True)
        bamfilePath=self.bamfilePath
        fastqFilePath=self.fastqFilePath

        readinfodict={}
        results=[]

        #get reads because pysam object cannot be used for multiprocessing so inputting bam file path 
        for i in self.generefdf.index:
            results.append(pool.apply_async(self._getreads,(bamfilePath,fastqFilePath,i)))
        pool.close()
        pool.join()
        results=[res.get() for res in results]

        for geneid,resls in zip(self.generefdf.index,results):
            readinfodict[geneid]=resls  


        #delete gene whose reads length is larger than maxReadCount
        for i in list(readinfodict.keys()):
            if len(readinfodict[i])>self.maxReadCount:
                readinfodict[i]=random.sample(readinfodict[i],self.maxReadCount)
            if len(readinfodict[i])<2:
                del readinfodict[i] 

        #print('hello,we finish get readinfodict')
        #store reads fetched
        outfilename=self.count_out_dir+'fetch_reads.pkl'
        with open(outfilename,'wb') as f:
            pickle.dump(readinfodict,f)


        return readinfodict




    def _do_clustering(self,dictcontentls):
        #geneid=success[0]
        readinfo=dictcontentls

        # do hierarchical cluster
        clusterModel = AgglomerativeClustering(n_clusters=None,linkage='average',distance_threshold=self.InnerDistance)

        posiarray=np.array([t[0] for t in readinfo]).reshape(-1,1)

        #print(posiarray.shape)

        CBarray=np.array([t[1] for t in readinfo]).reshape(-1,1)
        cigartuplearray=np.array([t[2] for t in readinfo]).reshape(-1,1)
        #seqarray=np.array([t[3] for t in readinfodict[geneid]]).reshape(-1,1)  #have more opportunity that this step has question
        clusterModel=clusterModel.fit(posiarray)

        #print('finish clustering fit')

        labels=clusterModel.labels_
        label,count=np.unique(labels,return_counts=True)

        #print('finish label unique')
        selectlabel=label[count>=self.minCount]
        selectcount=count[count>=self.minCount]
        #finalcount=list(selectcount[np.argsort(selectcount)[::-1]])
        finallabel=list(selectlabel[np.argsort(selectcount)[::-1]])

        #print(finallabel)

        #after adding 
        #numlabel=len(finallabel)    

        altTSSls=[]
        #if len(finallabel)>=2:
        for i in range(0,len(finallabel)):
            altTSSls.append([posiarray[labels==finallabel[i]],CBarray[labels==finallabel[i]],cigartuplearray[labels==finallabel[i]]])
        
        #print(altTSSls)
                       
        return altTSSls




    def _do_hierarchial_cluster(self):
        start_time=time.time()

        pool = multiprocessing.Pool(processes=self.nproc)
        readinfodict=self._get_gene_reads() 
        #print(len(readinfodict))

        altTSSdict={}
        altTSSls=[]
        dictcontentls=[]
        readls=list(readinfodict.keys())
        #print(len(readls))
        #print('unique gene id %i'%(len(set(readls))))
        for i in readls:
            dictcontentls.append(readinfodict[i])

        #print(inputpar[0])
        #print(len(dictcontentls))

        #print(len(inputpar))



        with multiprocessing.Pool(self.nproc) as pool:
            altTSSls=pool.map_async(self._do_clustering,dictcontentls).get()



        for geneidSec, reslsSec in zip(readls,altTSSls):
            altTSSdict[geneidSec]=reslsSec
        altTSSdict={k: v for k, v in altTSSdict.items() if v}

        # tss_output=self.count_out_dir+'gene_with_onecluster.pkl'
        # with open(tss_output,'wb') as f:
        #     pickle.dump(altTSSdict,f)

        print('do clustering Time elapsed',int(time.time()-start_time),'seconds.')

        return altTSSdict


    def _filter_false_positive(self):

        altTSSdict=self._do_hierarchial_cluster()      
        #get testX
        ## get RNA-seq X
        #make a new dictionary
        clusterdict={}
        for i in altTSSdict.keys():
            for j in range(0,len(altTSSdict[i])):
                #print(altTSSdict[i][j])
                startpos=np.min(altTSSdict[i][j][0])
                stoppos=np.max(altTSSdict[i][j][0])
                clustername=str(i)+'*'+str(startpos)+'_'+str(stoppos)
                

                count=len(altTSSdict[i][j][0])
                std=statistics.stdev(altTSSdict[i][j][0].flatten())
                summit_count=np.max(np.unique(altTSSdict[i][j][0].flatten(),return_counts=True)[1])
                unencoded_G_percent=sum([('14S' in ele)or('15S' in ele)or('16S' in ele) for ele in altTSSdict[i][j][2].flatten()])/count
                
                #summit position
                tempposi,tempposicount=np.unique(altTSSdict[i][j][0].flatten(),return_counts=True)
                maxpos=np.argmax(tempposicount)
                summitpos=tempposi[maxpos]   

                clusterdict[clustername]=(count,std,summit_count,unencoded_G_percent,j,i,summitpos)    
                #summitpos,altTSSdict[i][j][0],altTSSdict[i][j][1],altTSSdict[i][j][2]
                

        # cluster_output=self.count_out_dir+'before_filter_cluster.pkl'
        # with open(cluster_output,'wb') as f:
        #     pickle.dump(clusterdict,f)

        
        fourfeaturedf=pd.DataFrame(clusterdict).T 
        fourfeature_output=self.count_out_dir+'fourFeature.csv'
        fourfeaturedf.to_csv(fourfeature_output)

        print('one_gene_with_two_TSS_fourfeature : %i'%(len(fourfeaturedf)))
        test_X=fourfeaturedf.iloc[:,0:4]
        print('hello')

        print(os.path.abspath(__file__))

        print(os.path.dirname(os.path.abspath(__file__)))

        print(Path(os.path.dirname(os.path.abspath(__file__))))

        print(Path(os.path.dirname(os.path.abspath(__file__))).parents[1])

        pathstr=str(Path(os.path.dirname(os.path.abspath(__file__))).parents[0])+'/model/logistic_4feature_model.sav'
        loaded_model = pickle.load(open(pathstr, 'rb'))
        test_Y=loaded_model.predict(test_X.values)

        #do filtering, the result of this step should be output as final h5ad file display at single cell level. 
        afterfiltereddf=fourfeaturedf[test_Y==1]
        afterfiltereddf.columns=['all_counts','std','summit_count','unencoded_G_percent','TSS.no','gene_id','summit_position']

        # afterfilter_output=self.count_out_dir+'afterfiltered.csv'
        # afterfiltereddf.to_csv(afterfilter_output)


        allgeneID=afterfiltereddf['gene_id'].unique()
        keepdict={}
        for i in allgeneID:
            selectgeneiddf=afterfiltereddf[afterfiltereddf['gene_id']==i]
            keeptranscriptls=[]
            for j in selectgeneiddf.index:
                index=afterfiltereddf.loc[j]['TSS.no']
                keeptranscriptls.append(altTSSdict[i][index])
            keepdict[i]=keeptranscriptls



        tss_output=self.count_out_dir+'keepdict.pkl'
        with open(tss_output,'wb') as f:
            pickle.dump(keepdict,f)

        
        return keepdict


    def _do_anno_and_filter(self,inputpar):
        geneid=inputpar[0]
        altTSSitemdict=inputpar[1]
        temprefdf=self.tssrefdf[self.tssrefdf['gene_id']==geneid]

        #print(geneid)
        # print(altTSSdict)


        #use Hungarian algorithm to assign cluster to corresponding transcript
        cost_mtx=np.zeros((len(altTSSitemdict),temprefdf.shape[0]))
        for i in range(len(altTSSitemdict)):
            for j in range(temprefdf.shape[0]):
                cluster_val=altTSSitemdict[i][0]

                #this cost matrix should be corrected
                position,count=np.unique(cluster_val,return_counts=True)
                mode_position=position[np.argmax(count)]
                cost_mtx[i,j]=np.absolute(np.sum(mode_position-temprefdf.iloc[j,5]))
        row_ind, col_ind = linear_sum_assignment(cost_mtx)
        transcriptls=list(temprefdf.iloc[col_ind,:]['transcript_id'])

        # print(row_ind)
        # print(col_ind)

        #do quality control
        tssls=list(temprefdf.iloc[col_ind,:]['TSS'])
        #print(tssls)

        transcriptdict={}
        for i in range(0,len(tssls)):
            if (tssls[i]>=np.min(altTSSitemdict[i][0])) & (tssls[i]<=np.max(altTSSitemdict[i][0])):
                name1=str(geneid)+'_'+str(transcriptls[i])
                transcriptdict[name1]=(altTSSitemdict[row_ind[i]][0],altTSSitemdict[row_ind[i]][1],altTSSitemdict[row_ind[i]][2])
            else:
                newname1=str(geneid)+'_newTSS'
                transcriptdict[newname1]=(altTSSitemdict[row_ind[i]][0],altTSSitemdict[row_ind[i]][1],altTSSitemdict[row_ind[i]][2])
        #print(transcriptdict)



        return transcriptdict



    def _TSS_annotation(self):
        start_time=time.time()

        keepdict=self._filter_false_positive()

        keepIDls=list(keepdict.keys())
        
        inputpar=[]
        for i in keepIDls:
            inputpar.append((i,keepdict[i]))

        pool = multiprocessing.Pool(processes=self.nproc)
        with multiprocessing.Pool(self.nproc) as pool:
            #transcriptdictls=pool.map_async(self.filter_false_positive,inputpar).get()
            transcriptdictls=pool.map_async(self._do_anno_and_filter,inputpar).get()


        # tss_output=self.count_out_dir+'temp_tss.pkl'
        # with open(tss_output,'wb') as f:
        #     pickle.dump(transcriptdictls,f)

        extendls=[]
        for d in transcriptdictls:
            extendls.extend(list(d.items()))


        
        ### organize the output result
        d={'transcript_id':[transcript[0] for transcript in extendls],'TSS_start':[np.min(transcript[1][0]) for transcript in extendls],
        'TSS_end':[np.max(transcript[1][0]) for transcript in extendls]}

        regiondf=pd.DataFrame(d)
        print('do annotation Time elapsed',int(time.time()-start_time),'seconds.')
        # print(extendls)
        # print(regiondf)

        return extendls,regiondf




    def produce_sclevel(self):
        ctime=time.time()
        extendls,regiondf=self._TSS_annotation()
        #transcriptdfls=[]

        cellIDls=[]
        for i in range(0,len(extendls)):
            cellID=np.unique(extendls[i][1][1])
            cellIDls.append(list(cellID))
        cellIDset = set([item for sublist in cellIDls for item in sublist])
        finaldf=pd.DataFrame(index=cellIDset)



        for i in range(0,len(extendls)):
            transcriptid=extendls[i][0]       
            cellID,count=np.unique(extendls[i][1][1],return_counts=True)
            transcriptdf=pd.DataFrame({'cell_id':cellID,transcriptid:count})
            transcriptdf.set_index('cell_id',inplace=True)
            finaldf[transcriptid]=finaldf.index.map(transcriptdf[transcriptid])


        finaldf.fillna(0,inplace=True)
        adata=ad.AnnData(finaldf)
        vardf=pd.DataFrame(adata.var.copy())
        vardf.reset_index(inplace=True)
        vardf.columns=['transcript_id']
        vardf=vardf.join(regiondf.set_index('transcript_id'), on='transcript_id')
        vardf['gene_id']=vardf['transcript_id'].str.split('_',expand=True)[0]
        vardf=vardf.merge(self.generefdf,on='gene_id')
        vardf.set_index('transcript_id',drop=True,inplace=True)

        adata.var=vardf
        sc_output_h5ad=self.count_out_dir+'scTSS_count_all.h5ad'
        adata.write(sc_output_h5ad)

        #filter according to user' defined distance
        newdf=adata.var.copy()
        newdf.reset_index(inplace=True)
        selectedf=newdf[newdf.duplicated('gene_id',keep=False)]  #get data frame which includes two transcript for one gene
        geneID=selectedf['gene_id'].unique()

        keepdfls=[]
        for i in geneID:
            tempdf=selectedf[selectedf['gene_id']==i]

            tempdf=tempdf.sort_values('transcript_id',ascending=False)
            tempdf['diff']=tempdf['TSS_start'].diff()
            keepdf=tempdf[tempdf['diff'].isna()|tempdf['diff'].abs().ge(self.clusterDistance)]    #want to get TSS whose cluster distance is more than user defined.
            #keepdf=keepdf.iloc[:2,:]
            keepdfls.append(keepdf) 

        print(keepdfls)


        allkeepdf=reduce(lambda x,y:pd.concat([x,y]),keepdfls)
        finaltwodf=allkeepdf[allkeepdf.duplicated('gene_id',keep=False)] 
        finaltwoadata=adata[:,adata.var.index.isin(finaltwodf['transcript_id'])]  

        sc_output_h5ad=self.count_out_dir+'scTSS_count_two.h5ad'
        finaltwoadata.write(sc_output_h5ad)


        print('produce h5ad Time elapsed',int(time.time()-ctime),'seconds.')


        return adata






    def window_sliding(self,genereads,TSS_start,TSS_end,strand):

        leftIndex=0

        # do filtering; drop reads which does not include unencoded G
        filterls=[]
        for i in genereads:
            if ('14S' in i[2]) or ('15S' in i[2]) or ('16S' in i[2]):
                filterls.append(i)


        #calculate the TSS position and corresponding counts
        promoterTSS=[]
        for read in filterls:
            tss=read[0]
            if (tss>=TSS_start)&(tss<=TSS_end):
                promoterTSS.append(tss)
        TSS,count=np.unique(promoterTSS,return_counts=True)

        nonzeroarray=np.asarray((TSS, count)).T


        if strand=='+':
            sortfinalarray=nonzeroarray[nonzeroarray[:, 0].argsort()]
            TSS=sortfinalarray.T[0]
            count=sortfinalarray.T[1]
        elif strand=='-':
            sortfinalarray=nonzeroarray[nonzeroarray[:, 0].argsort()[::-1]]
            TSS=sortfinalarray.T[0]
            count=sortfinalarray.T[1]


        #do something with sliding windows algorithm   
        storels=[]
        for i in range(len(TSS) - self.windowSize + 1):
            onewindow=TSS[i: i + self.windowSize]
            correspondingcount=count[i: i + self.windowSize]
            middlecount=correspondingcount[leftIndex]
            foldchange=(middlecount+1)/(sum(correspondingcount)/len(correspondingcount)+1)
            storels.append([onewindow[leftIndex],correspondingcount[leftIndex],foldchange])
            
        foldchangels=[i[2] for i in storels]
        sortindex=sorted(range(len(foldchangels)), key=lambda k: foldchangels[k],reverse=True)
        allsortls=[storels[i] for i in sortindex]

        return allsortls



    def _get_CTSS(self,fetchadata):
        oneclusterfilePath=self.count_out_dir+'afterfiltered.csv'
        alloneclusterdf=pd.read_csv(oneclusterfilePath)
        alloneclusterdf['gene_id']=alloneclusterdf['Unnamed: 0'].str.split('*',expand=True)[0]
        alloneclusterdf['TSS_start']=alloneclusterdf['Unnamed: 0'].str.split('*',expand=True)[1].str.split('_',expand=True)[0].astype('float')
        alloneclusterdf['TSS_end']=alloneclusterdf['Unnamed: 0'].str.split('_',expand=True)[1].astype('float')

        self.generefdf.reset_index(inplace=True)
        

        print(self.generefdf)
        stranddf=self.generefdf[['Strand','gene_id']]
        alloneclusterdf=alloneclusterdf.merge(stranddf,on='gene_id')






        pool = multiprocessing.Pool(processes=self.nproc)
        allsortls=[]

        for i in range(0,len(alloneclusterdf)):
            
            geneID=alloneclusterdf['gene_id'][i]
            genereads=fetchadata[geneID]
            TSS_start=alloneclusterdf['TSS_start'][i]
            TSS_end=alloneclusterdf['TSS_end'][i]
            strand=alloneclusterdf['Strand'][i]
            allsortls.append(pool.apply_async(self.window_sliding,(genereads,TSS_start,TSS_end,strand)))

        pool.close()
        pool.join()
        results=[res.get() for res in allsortls]

        allsortfddict={}
        for geneid,resls in zip(alloneclusterdf['Unnamed: 0'],results):
            allsortfddict[geneid]=resls  

        
        ctssOutPath=self.ctss_out_dir+'CTSS_foldchange.pkl'
        with open(ctssOutPath,'wb') as f:
            pickle.dump(allsortfddict,f)

        return allsortfddict


    def pickCTSS(self,ctssls):

        keepCTSS=[]
        for ele in ctssls:
            if (ele[1]>self.minCTSSCount)&(ele[2]>self.minFC):
                keepCTSS.append(ele)
        return keepCTSS
    





    def produce_CTSS_adata(self):
        ctime=time.time()

        readspath=self.count_out_dir+'fetch_reads.pkl'
        with open(readspath,'rb') as f:
            fetchadata=pickle.load(f)

        allsortfddict=self._get_CTSS(fetchadata)
        keepdict={}
        for ctssid in allsortfddict.keys():
            keepdict[ctssid]=self.pickCTSS(allsortfddict[ctssid])
        
        #print(keepdict)


        #get the cellID meeting our requirement
        cellIDdict={}
        for i in keepdict.keys():
            
            for j in keepdict[i]:
                geneid=i.split('*')[0]
                newid=i+'#'+str(j[0])+'@'+str(j[1])+'$'+str(j[2])
                cellIDls=[]
                for ele in fetchadata[geneid]:
                    if j[0]==ele[0]:
                        cellIDls.append(ele[1])
                cellIDdict[newid]=cellIDls

        #print(len(cellIDdict))


        #create a big matrix including cell ID
        cellidls=list(cellIDdict.values())
        cellidset = set([item for sublist in cellidls for item in sublist])
        ctssfinaldf=pd.DataFrame(index=cellidset)



        
        for clusterID in cellIDdict.keys():
            cellID,count=np.unique(cellIDdict[clusterID],return_counts=True)
            CTSSdf=pd.DataFrame({'cell_id':cellID,clusterID:count})
            CTSSdf.set_index('cell_id',inplace=True)
            ctssfinaldf[clusterID]=ctssfinaldf.index.map(CTSSdf[clusterID])


        ctssfinaldf.fillna(0,inplace=True)
        ctssadata=ad.AnnData(ctssfinaldf)

        ctssvardf=pd.DataFrame(ctssadata.var.copy())
        ctssvardf.reset_index(inplace=True)
        ctssvardf.columns=['clusterID']
        ctssvardf['gene_id']=ctssvardf['clusterID'].str.split('*',expand=True)[0]
        ctssvardf['CTSS']=ctssvardf['clusterID'].str.split('#',expand=True)[1].str.split('@',expand=True)[0]
        ctssvardf['counts_dropped_UnencodedG']=ctssvardf['clusterID'].str.split('@',expand=True)[1].str.split('$',expand=True)[0]
        ctssvardf['fold_change']=ctssvardf['clusterID'].str.split('$',expand=True)[1]


        ctssvardf=ctssvardf.merge(self.generefdf,on='gene_id')
        ctssvardf.set_index('clusterID',drop=True,inplace=True)
        ctssadata.var=ctssvardf

        ctss_output_h5ad=self.ctss_out_dir+'all_ctss.h5ad'
        ctssadata.write(ctss_output_h5ad)


        twoctssselect=ctssadata.var[ctssadata.var.duplicated('gene_id',keep=False)].index
        twoctssadata=ctssadata[:,twoctssselect]

        sc_output_h5ad=self.ctss_out_dir+'all_ctss_two.h5ad'
        twoctssadata.write(sc_output_h5ad)

        print('produce CTSS h5ad Time elapsed',int(time.time()-ctime),'seconds.')


        return twoctssadata



