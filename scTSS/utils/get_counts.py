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
    def __init__(self,generefPath,tssrefPath,bamfilePath,fastqFilePath,outdir,cellBarcodePath,nproc,minCount=50,maxReadCount=10000,clusterDistance=300):
        self.generefdf=pd.read_csv(generefPath,delimiter='\t')
        self.generefdf['len']=self.generefdf['End']-self.generefdf['Start']
        self.tssrefdf=pd.read_csv(tssrefPath,delimiter='\t')
        self.bamfilePath=bamfilePath
        self.count_out_dir=str(outdir)+'/count/'
        if not os.path.exists(self.count_out_dir):
            os.mkdir(self.count_out_dir)
        self.minCount=minCount
        self.cellBarcode=pd.read_csv(cellBarcodePath,delimiter='\t')['cell_id'].values
        self.nproc=nproc
        self.maxReadCount=maxReadCount
        self.clusterDistance=clusterDistance
        self.fastqFilePath=fastqFilePath

        

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
        clusterModel = AgglomerativeClustering(n_clusters=None,linkage='average',distance_threshold=100)

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
        if len(finallabel)>=2:
            for i in range(0,len(finallabel)):
                altTSSls.append([posiarray[labels==finallabel[i]],CBarray[labels==finallabel[i]],cigartuplearray[labels==finallabel[i]]])
        
        #print(altTSSls)
                       
        return altTSSls




    def _do_hierarchial_cluster(self):
        start_time=time.time()

        pool = multiprocessing.Pool(processes=self.nproc)
        readinfodict=self._get_gene_reads() 
        print(len(readinfodict))

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

        print('finish multi-processing')
        #print("I need success")
        # print(altTSSls)
        # print(len(readls))


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

        pathstr=str(Path(os.path.dirname(os.path.abspath(__file__))).parents[1])+'/test/logistic_4feature_model.sav'
        loaded_model = pickle.load(open(pathstr, 'rb'))
        test_Y=loaded_model.predict(test_X.values)

        #do filtering
        afterfiltereddf=fourfeaturedf[test_Y==1]
        afterfilter_output=self.count_out_dir+'afterfiltered.csv'
        afterfiltereddf.to_csv(afterfilter_output)
        print('after_filter_false_positive_TSS_afterfiltereddf : %i'%(len(afterfiltereddf)))


        selectedf=afterfiltereddf[afterfiltereddf.duplicated(5,keep=False)] #select number of transcript is more than 2
        geneID=selectedf[5].unique()

        #make sure the distance of cluster is more than user's defination distance 
        keepdfls=[]
        for i in geneID:
            tempdf=selectedf[selectedf[5]==i]
            
            tempdf=tempdf.sort_values(0,ascending=False)
            tempdf['diff']=tempdf[6].diff()
            keepdf=tempdf[tempdf['diff'].isna()|tempdf['diff'].abs().ge(self.clusterDistance)]
            keepdf=keepdf.iloc[:2,:]
            keepdfls.append(keepdf) 

        allkeepdf=reduce(lambda x,y:pd.concat([x,y]),keepdfls)
        # allkeepdf.to_csv('/storage/yhhuang/users/ruiyan/15organ/SRR13075718_scTSS_out/count/allkeeped.csv')
        # print('after_filter_accordingtoDistance_afterfiltereddf : %i'%(len(allkeepdf)))

        tss_output=self.count_out_dir+'allkeep.csv'
        allkeepdf.to_csv(tss_output)

        allkeepdf=allkeepdf[allkeepdf.duplicated(5,keep=False)]
        # allkeepdf.to_csv('/storage/yhhuang/users/ruiyan/15organ/SRR13075718_scTSS_out/count/final_keeped.csv')
        # print('after_filtergene_toget_two_TSS_again_afterfiltereddf : %i'%(len(allkeepdf)))


        allgeneID=allkeepdf[5].unique()
        keepdict={}
        for i in allgeneID:
            selectgeneiddf=allkeepdf[allkeepdf[5]==i]
            twotranscriptls=[]
            for j in selectgeneiddf.index:
                index=allkeepdf.loc[j][4]
                twotranscriptls.append(altTSSdict[i][index])
            keepdict[i]=twotranscriptls



        tss_output=self.count_out_dir+'keepdict.pkl'
        with open(tss_output,'wb') as f:
            pickle.dump(keepdict,f)

        
        return keepdict


    def _do_anno_and_filter(self,inputpar):
        geneid=inputpar[0]
        altTSSdict=inputpar[1]
        temprefdf=self.tssrefdf[self.tssrefdf['gene_id']==geneid]

        #print(geneid)
        # print(altTSSdict)


        #use Hungarian algorithm to assign cluster to corresponding transcript
        cost_mtx=np.zeros((2,temprefdf.shape[0]))
        for i in range(2):
            for j in range(temprefdf.shape[0]):
                cluster_val=altTSSdict[geneid][i][0]
                quanp=cluster_val[cluster_val<np.percentile(cluster_val,10)]
                cost_mtx[i,j]=np.absolute(np.sum(quanp-temprefdf.iloc[j,5]))
        row_ind, col_ind = linear_sum_assignment(cost_mtx)
        transcriptls=list(temprefdf.iloc[col_ind,:]['transcript_id'])

        #do quality control
        tssls=list(temprefdf.iloc[col_ind,:]['TSS'])
        #print(tssls)


        transcriptdict={}
        if (tssls[0]>np.min(altTSSdict[geneid][0][0]-10)) & (tssls[0]<np.max(altTSSdict[geneid][0][0]+10)):
            name1=str(geneid)+'_'+str(transcriptls[0])
            transcriptdict[name1]=(altTSSdict[geneid][row_ind[0]][0],altTSSdict[geneid][row_ind[0]][1],altTSSdict[geneid][row_ind[0]][2]) 
        else:
            newname1=str(geneid)+'_newTSS_1'
            transcriptdict[newname1]=(altTSSdict[geneid][row_ind[0]][0],altTSSdict[geneid][row_ind[0]][1],altTSSdict[geneid][row_ind[0]][2])

        if (tssls[1]>np.min(altTSSdict[geneid][1][0]-10)) & (tssls[1]<np.max(altTSSdict[geneid][1][0]+10)):
            name2=str(geneid)+'_'+str(transcriptls[1])
            transcriptdict[name2]=(altTSSdict[geneid][row_ind[1]][0],altTSSdict[geneid][row_ind[1]][1],altTSSdict[geneid][row_ind[1]][2])  
        else:
            newname2=str(geneid)+'_newTSS_2'
            transcriptdict[newname2]=(altTSSdict[geneid][row_ind[1]][0],altTSSdict[geneid][row_ind[1]][1],altTSSdict[geneid][row_ind[1]][2])


        # cluster_output=self.count_out_dir+'do_annotation.pkl'
        # with open(cluster_output,'wb') as f:
        #     pickle.dump(transcriptdict,f)


        #print(transcriptdict)

        return transcriptdict



    def _TSS_annotation(self):
        start_time=time.time()

        keepdict=self._filter_false_positive()
        
        inputpar=[]
        for i in keepdict:
            inputpar.append((i,keepdict))

        #print(inputpar)

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
        sc_output_h5ad=self.count_out_dir+'sc_TSS_count.h5ad'
        adata.write(sc_output_h5ad)
        print('produce h5ad Time elapsed',int(time.time()-ctime),'seconds.')


        return adata