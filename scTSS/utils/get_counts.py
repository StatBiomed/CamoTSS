import brie
import numpy as np
import pandas as pd
from brie.utils import fetch_reads,load_samfile
import os
import pickle
from functools import reduce
import anndata as ad
import multiprocessing 
from pathos.multiprocessing import ProcessingPool 
import pathos.pools as pp
import pysam
from sklearn.cluster import AgglomerativeClustering
from scipy.optimize import linear_sum_assignment
import time
from sklearn.cluster import KMeans
from collections import defaultdict
from os import getpid
import random
import pickle
import numpy.lib.recfunctions as rfn
import numpy.ma as ma
from itertools import compress
import statistics
import editdistance
import warnings
from collections import ChainMap
from sklearn.preprocessing import MinMaxScaler
from sklearn import linear_model
from pathlib import Path
import sys
from sklearn.neighbors import NearestNeighbors




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
        self.tssrefdf['transcript_id']=self.tssrefdf['gene_id']+'_'+self.tssrefdf['transcript_id']
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
        clusterModel = AgglomerativeClustering(n_clusters=None,linkage='average',distance_threshold=10)
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

        # with open('/storage/yhhuang/users/ruiyan/iPSC/scRNAseq/scTSS_out/count/fetch_reads.pkl','rb') as f:
        #     readinfodict=pickle.load(f)


        readinfodict=self._get_gene_reads() 
        print(len(readinfodict))

        altTSSdict={}
        altTSSls=[]
        dictcontentls=[]
        readls=list(readinfodict.keys())
        #print(len(readls))
        print('unique gene id %i'%(len(set(readls))))
        for i in readls:
            dictcontentls.append(readinfodict[i])

        #print(inputpar[0])
        print(len(dictcontentls))

        #print(len(inputpar))



        with multiprocessing.Pool(self.nproc) as pool:
            altTSSls=pool.map_async(self._do_clustering,dictcontentls).get()

        #print('finish multi-processing')
        #print("I need success")
        # print(altTSSls)
        # print(len(readls))


        for geneidSec, reslsSec in zip(readls,altTSSls):
            altTSSdict[geneidSec]=reslsSec
        altTSSdict={k: v for k, v in altTSSdict.items() if v}

        tss_output=self.count_out_dir+'gene_with_onecluster.pkl'
        with open(tss_output,'wb') as f:
            pickle.dump(altTSSdict,f)

        print('do clustering Time elapsed',int(time.time()-start_time),'seconds.')

        return altTSSdict



    def filter_closer_TSS(self,anno_unanno_gtf_df,gene_id):
        specificGenedf=anno_unanno_gtf_df[anno_unanno_gtf_df['gene_id']==gene_id]
        #print(specificGenedf)
        specificGenedf.sort_values('TSS',inplace=True)
        specificGenedf['diff']=specificGenedf['TSS'].diff()
        specificGenedf.reset_index(inplace=True,drop=True)
        annoIndexls=list(specificGenedf[specificGenedf['transcript_id'].str.contains('ENST')].index.values)
        unannoIndexls=list(specificGenedf[specificGenedf['transcript_id'].str.contains('\*')].index.values)

        keepindex=[]
        for annoIndex in annoIndexls:
            if (specificGenedf.loc[annoIndex]['diff']>5)|(np.isnan(specificGenedf.loc[annoIndex]['diff'])==True):
                keepindex.append(annoIndex)

        for unannoIndex in unannoIndexls:
            try:
                if ((specificGenedf.loc[unannoIndex]['diff']>50)|(np.isnan(specificGenedf.loc[annoIndex]['diff'])==True))&(specificGenedf.loc[unannoIndex+1]['diff']>50):
                    keepindex.append(unannoIndex)
            except KeyError:
                if (specificGenedf.loc[unannoIndex]['diff']>50):
                    keepindex.append(unannoIndex)
        
        specificGenedf=specificGenedf.loc[keepindex]

                
        if len(specificGenedf)>=2:
            return specificGenedf


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
                clustername=str(i)+'_*'+str(startpos)+'_'+str(stoppos)
                

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
        fourfeaturedf.reset_index(inplace=True)
        fourfeaturedf.columns=['cluster_name','count','std','summit_count','unencodedG_percent','index','gene_id','summit_position']
        fourfeature_output=self.count_out_dir+'fourFeature.csv'
        fourfeaturedf.to_csv(fourfeature_output)

        geneID=fourfeaturedf['gene_id'].unique()

        keepdfls=[]
        for i in geneID:
            tempdf=fourfeaturedf[fourfeaturedf['gene_id']==i]
            
            tempdf=tempdf.sort_values('count',ascending=False)
            tempdf['diff']=tempdf['summit_position'].diff()
            keepdf=tempdf[tempdf['diff'].isna()|tempdf['diff'].abs().ge(self.clusterDistance)]
            keepdfls.append(keepdf) 
        filterDistancedf=reduce(lambda x,y:pd.concat([x,y]),keepdfls)



        #print('one_gene_with_two_TSS_fourfeature : %i'%(len(fourfeaturedf)))
        test_X=filterDistancedf.iloc[:,1:5]

        pathstr=str(Path(os.path.dirname(os.path.abspath(__file__))).parents[1])+'/test/logistic_4feature_model.sav'
        loaded_model = pickle.load(open(pathstr, 'rb'))
        test_Y=loaded_model.predict(test_X.values)

        #do filtering
        afterfiltereddf=filterDistancedf[test_Y==1] 

        #print('after_filter_false_positive_TSS_afterfiltereddf : %i'%(len(afterfiltereddf)))
        afterfiltereddf['cluster_start']=afterfiltereddf['cluster_name'].str.split('*',expand=True)[1].str.split('_',expand=True)[0]
        afterfiltereddf['cluster_stop']=afterfiltereddf['cluster_name'].str.split('*',expand=True)[1].str.split('_',expand=True)[1]
        afterfiltereddf=afterfiltereddf.merge(self.generefdf,on='gene_id')

        afterfiltereddf=afterfiltereddf[['cluster_name','gene_id','gene_name','Chromosome','Strand','summit_position']]
        afterfiltereddf.rename(columns={'cluster_name':'transcript_id','summit_position':'TSS'},inplace=True)


        anno_unanno_gtf_df=pd.concat([afterfiltereddf,self.tssrefdf],axis=0)
        anno_unanno_gtf_df['transcript_id']=anno_unanno_gtf_df['transcript_id']+'@'+anno_unanno_gtf_df['TSS'].astype('str')

        outfilename=self.count_out_dir+'anno_unanno_gtf.csv'
        anno_unanno_gtf_df.to_csv(outfilename)
        #print(anno_unanno_gtf_df)

        keepTSSls=[]
        for geneid in anno_unanno_gtf_df['gene_id'].unique():
            specificGenedf=self.filter_closer_TSS(anno_unanno_gtf_df,geneid)
            keepTSSls.append(specificGenedf)

        filterTSSdf=reduce(lambda x,y:pd.concat([x,y],axis=0) ,keepTSSls)
        filterTSSdf.dropna(subset=['Chromosome'],inplace=True)

        return filterTSSdf



    def distribution(self,geneID,readinfo,filterTSSdf):

        
        singletranscriptdict={}
        try:
            generefdf=filterTSSdf[filterTSSdf['gene_id']==geneID]
            generefdf.reset_index(inplace=True,drop=True)
            trainX=generefdf['TSS'].values.reshape(-1,1)
            posiarray=np.array([t[0] for t in readinfo]).reshape(-1,1)
            CBarray=np.array([t[1] for t in readinfo]).reshape(-1,1)
            nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(trainX)
            distances, indices = nbrs.kneighbors(posiarray)
            select=indices[distances<=5]
            selectCB=CBarray[distances<=5]
            label=np.unique(select)
            for i in range(0,len(label)):
                keyid=generefdf.loc[label[i]]['transcript_id']
                if len(select[select==label[i]])>=self.minCount:
                    singletranscriptdict[keyid]=(len(select[select==label[i]]),selectCB[select==label[i]])

            #print(singletranscriptdict)
            if len(singletranscriptdict)>=2:
                return singletranscriptdict
            else:
                return {}
            
        except ValueError:
            return {}



    def get_transcriptdict(self,):

        start_time=time.time()
        pool = multiprocessing.Pool(processes=self.nproc)

        readsfilename=self.count_out_dir+'fetch_reads.pkl'
        with open(readsfilename,'rb') as f:
            readinfodict=pickle.load(f)

        filterTSSdf=self._filter_false_positive()
        transcriptdictls=[]

        for geneID in readinfodict.keys():
            transcriptdictls.append(pool.apply_async(self.distribution,(geneID,readinfodict[geneID],filterTSSdf)))
        pool.close()
        pool.join()

        transcriptdictls=[res.get() for res in transcriptdictls]

        #print(transcriptdictls)
        transcriptdict = {}
        for d in transcriptdictls:
            transcriptdict.update(d)


        print('do distribution',int(time.time()-start_time),'seconds.')

        #store reads fetched
        outfilename=self.count_out_dir+'transcript_store.pkl'
        with open(outfilename,'wb') as f:
            pickle.dump(transcriptdict,f)

        return transcriptdict



    def produce_sclevel(self):
        ctime=time.time()
        transcriptdict=self.get_transcriptdict()

        cellIDls=[]
        for transcriptID in transcriptdict.keys():
            cellID=np.unique(transcriptdict[transcriptID][1])
            cellIDls.append(list(cellID))
        cellIDset = set([item for sublist in cellIDls for item in sublist])
        finaldf=pd.DataFrame(index=cellIDset)


        for transcriptID in transcriptdict.keys():
            cellID,count=np.unique(transcriptdict[transcriptID][1],return_counts=True)
            transcriptdf=pd.DataFrame({'cell_id':cellID,transcriptID:count})
            transcriptdf.set_index('cell_id',inplace=True)
            finaldf[transcriptID]=finaldf.index.map(transcriptdf[transcriptID])
        finaldf.fillna(0,inplace=True)
        #print(finaldf)
        adata=ad.AnnData(finaldf) 
  
        
        vardf=pd.DataFrame(adata.var.copy())
        vardf.reset_index(inplace=True)

        vardf['gene_id']=vardf['index'].str.split('_',expand=True)[0]
        #print(vardf)
        vardf=vardf.merge(self.generefdf,on='gene_id')




        vardf.set_index('index',drop=True,inplace=True)
        adata.var=vardf
        sc_output_h5ad=self.count_out_dir+'sc_TSS_count.h5ad'
        adata.write(sc_output_h5ad)
        print('produce h5ad Time elapsed',int(time.time()-ctime),'seconds.')


        return adata

