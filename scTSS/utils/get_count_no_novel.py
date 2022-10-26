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




class get_old_TSS_count():
    def __init__(self,generefPath,tssrefPath,bamfilePath,fastqFilePath,outdir,cellBarcodePath,nproc,minCount=50,maxReadCount=50000,clusterDistance=200,psi=0.1):
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
        self.psi=psi
        self.maxReadCount=maxReadCount
        self.clusterDistance=clusterDistance
        self.fastqFilePath=fastqFilePath

        

    def _getreads(self,bamfilePath,fastqFilePath,geneid):
        #fetch reads1 in gene 
        samFile, _chrom = check_pysam_chrom(bamfilePath,'chr'+str(self.generefdf.loc[geneid]['Chromosome']))
        reads = fetch_reads(samFile, _chrom,  self.generefdf.loc[geneid]['Start'] , self.generefdf.loc[geneid]['End'],  trimLen_max=100)
        reads1_umi = reads["reads1"]



        #select according to GX tag and CB (filter according to user owned cell)
        reads1_umi=[r for r in reads1_umi if r.get_tag('GX')==geneid]
        reads1_umi=[r for r in reads1_umi if r.get_tag('CB') in self.cellBarcode]

        # #filter strand invasion
        # fastqFile=get_fastq_file(fastqFilePath)
        # reads1_umi=[r for r in reads1_umi if editdistance.eval(fastqFile.fetch(start=r.reference_start-14, end=r.reference_start-1, region=str(self.generefdf.loc[geneid]['Chromosome'])),'TTTCTTATATGGG') >3 ]


        # #filter according to the cateria of SCAFE
        # reads1_umi=[r for r in reads1_umi if editdistance.eval(r.query_sequence[9:14],'ATGGG')<=4]
        # reads1_umi=[r for r in reads1_umi if len(r.cigartuples)>=2]
        # reads1_umi=[r for r in reads1_umi if (r.cigartuples[0][0]==4)&(r.cigartuples[0][1]>6)&(r.cigartuples[0][1]<20)&(r.cigartuples[1][0]==0)&(r.cigartuples[1][1]>5)]


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

        # #store reads fetched
        # outfilename=self.count_out_dir+'fetch_reads.pkl'
        # with open(outfilename,'wb') as f:
        #     pickle.dump(readinfodict,f)

        return readinfodict
















        




    def _do_clustering(self,success):
        geneid=success[0]
        readinfodict=success[1]  
        altTSSls=[]

        # do hierarchical cluster
        clusterModel = AgglomerativeClustering(n_clusters=None,linkage='average',distance_threshold=100)
        posiarray=np.array([t[0] for t in readinfodict[geneid]]).reshape(-1,1)
        CBarray=np.array([t[1] for t in readinfodict[geneid]]).reshape(-1,1)
        cigartuplearray=np.array([t[2] for t in readinfodict[geneid]]).reshape(-1,1)
        #seqarray=np.array([t[3] for t in readinfodict[geneid]]).reshape(-1,1)  #have more opportunity that this step has question
        clusterModel=clusterModel.fit(posiarray)
        labels=clusterModel.labels_
        label,count=np.unique(labels,return_counts=True)
        selectlabel=label[count>=self.minCount]
        selectcount=count[count>=self.minCount]
        finalcount=list(selectcount[np.argsort(selectcount)[::-1]])
        finallabel=list(selectlabel[np.argsort(selectcount)[::-1]])

        if len(finallabel)>=2:
            altTSSls=[]
            for i in range(0,len(finallabel)):
                altTSSls.append([posiarray[labels==finallabel[i]],CBarray[labels==finallabel[i]],cigartuplearray[labels==finallabel[i]]])
       
        return altTSSls




    def _do_hierarchial_cluster(self):
        start_time=time.time()

        pool = multiprocessing.Pool(processes=self.nproc)
        readinfodict=self._get_gene_reads() 

        altTSSdict={}
        altTSSls=[]
        inputpar=[]
        readls=list(readinfodict.keys())
        for i in readls:
            inputpar.append((i,readinfodict))

        with multiprocessing.Pool(self.nproc) as pool:
            altTSSls=pool.map_async(self._do_clustering,inputpar).get()

        for geneidSec, reslsSec in zip(readls,altTSSls):
            altTSSdict[geneidSec]=reslsSec
        altTSSdict={k: v for k, v in altTSSdict.items() if v}

        tss_output=self.count_out_dir+'without_anno_temp_tss.pkl'
        with open(tss_output,'wb') as f:
            pickle.dump(altTSSdict,f)

        print('do clustering Time elapsed',int(time.time()-start_time),'seconds.')

        return altTSSdict



    def _do_anno_and_filter(self,inputpar):
        geneid=inputpar[0]
        altTSSdict=inputpar[1]
        temprefdf=self.tssrefdf[self.tssrefdf['gene_id']==geneid]


        #use Hungarian algorithm to assign cluster to corresponding transcript
        cost_mtx=np.zeros((len(altTSSdict[geneid]),temprefdf.shape[0]))
        for i in range(0,len(altTSSdict[geneid])):
            for j in range(temprefdf.shape[0]):
                #print(altTSSdict[geneid])
                cluster_val=altTSSdict[geneid][i][0]
                quanp=cluster_val[cluster_val<np.percentile(cluster_val,10)]
                cost_mtx[i,j]=np.absolute(np.sum(quanp-temprefdf.iloc[j,5]))
        row_ind, col_ind = linear_sum_assignment(cost_mtx)
        transcriptls=list(temprefdf.iloc[col_ind,:]['transcript_id'])

        #do quality control
        tssls=list(temprefdf.iloc[col_ind,:]['TSS'])


        transcriptdict={}
        for i in range(0,len(tssls)):
            if np.min([np.absolute(j-np.min(altTSSdict[geneid][i][0])) for j in tssls])<300:
                transcriptdict[transcriptls[i]]=(altTSSdict[geneid][row_ind[i]][0],altTSSdict[geneid][row_ind[i]][2])
            else:
                newname=str(geneid)+'_newTSS'+'_'+str(i)
                transcriptdict[newname]=(altTSSdict[geneid][row_ind[i]][0],altTSSdict[geneid][row_ind[i]][2])
        


        return transcriptdict



    def _TSS_annotation(self):
        start_time=time.time()

        altTSSdict=self._do_hierarchial_cluster()
        
        inputpar=[]
        for i in altTSSdict:
            inputpar.append((i,altTSSdict))

        pool = multiprocessing.Pool(processes=self.nproc)
        with multiprocessing.Pool(self.nproc) as pool:
            transcriptdictls=pool.map_async(self._do_anno_and_filter,inputpar).get()

        #print(transcriptdictls)
        # for i in transcriptdictls:
        #     print(i)
        #add some limitation to filter novel TSS
        # remainingtranscriptls=[]
        # if not any('newTSS' in i.keys() for i in transcriptdictls):
        #     remainingtranscriptls.append(i)

        # print('hello')
        # print(remainingtranscriptls)
        extendls=[]
        for d in transcriptdictls:
            extendls.extend(list(d.items()))

        print(extendls)

        # remainingtranscriptls=[]
        # for i in extendls:
        #     remainingtranscriptls.append({k:v for k,v in i.items() if not 'null' in k})


        
        d={'transcript_id':[transcript[0] for transcript in extendls],'TSS_start':[np.min(transcript[1][0]) for transcript in extendls],
        'TSS_end':[np.max(transcript[1][0]) for transcript in extendls]}


        regiondf=pd.DataFrame(d)

        return extendls,regiondf




    def get_transcriptls(self,eachtranscript):
        
        transcriptid=eachtranscript[0]

        #here, using numpy structure
        cellID,count=np.unique(eachtranscript[1][1],return_counts=True)
        datatype=np.dtype({'names':['cellID',transcriptid],'formats':['S32','i']})
        transnp=np.fromiter(zip(cellID,count),dtype=datatype)

        return transnp




    def produce_sclevel(self):
        ctime=time.time()
        extendls,regiondf=self._TSS_annotation()
        #transcriptdfls=[]
        pool = multiprocessing.Pool(processes=self.nproc)

        with multiprocessing.Pool(self.nproc) as pool:
            transcriptdfls=pool.map_async(self.get_transcriptls,extendls).get()



        finalnp=reduce(lambda x,y: rfn.join_by('cellID',x,y,jointype='outer',usemask=True), transcriptdfls )
        finalnp=finalnp.filled(0)
        unstrfinalnp=rfn.structured_to_unstructured(finalnp)


        #ctime=time.time()
        adata=ad.AnnData(unstrfinalnp[:,1:])
        adata.obs.index=unstrfinalnp[:,0]
        adata.var.index=finalnp.dtype.names[1:]
        vardf=pd.DataFrame(adata.var.copy())
        vardf.reset_index(inplace=True)
        vardf.columns=['transcript_id']
        vardf=vardf.join(regiondf.set_index('transcript_id'), on='transcript_id')
        vardf.set_index('transcript_id',drop=True,inplace=True)
        adata.var=vardf
        sc_output_h5ad=self.count_out_dir+'sc_TSS_count.h5ad'
        adata.write(sc_output_h5ad)
        print('produce h5ad Time elapsed',int(time.time()-ctime),'seconds.')


        return adata
