import numpy as np
import pandas as pd
from brie.utils import fetch_reads
import os
import pickle
from functools import reduce
import anndata as ad
import multiprocessing 
import pysam
import time
from os import getpid
import random
import pickle
import editdistance
import warnings
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




class get_old_TSS_count():
    def __init__(self,generefPath,filterTssPath,bamfilePath,fastqFilePath,outdir,cellBarcodePath,nproc,minCount,maxReadCount,clusterDistance):
        self.generefdf=pd.read_csv(generefPath,delimiter='\t')
        self.generefdf['len']=self.generefdf['End']-self.generefdf['Start']
        self.tssrefdf=pd.read_csv(filterTssPath,delimiter='\t')
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
        samFile, _chrom = check_pysam_chrom(bamfilePath,'chr'+str(self.generefdf.loc[geneid]['Chromosome']))
        reads = fetch_reads(samFile, _chrom,  self.generefdf.loc[geneid]['Start'] , self.generefdf.loc[geneid]['End'],  trimLen_max=100)
        reads1_umi = reads["reads1"]



        #select according to GX tag and CB (filter according to user owned cell)
        reads1_umi=[r for r in reads1_umi if r.get_tag('GX')==geneid]
        reads1_umi=[r for r in reads1_umi if r.get_tag('CB') in self.cellBarcode]

        # #filter strand invasion
        fastqFile=get_fastq_file(fastqFilePath)

        reads1_umi=[r for r in reads1_umi if editdistance.eval(fastqFile.fetch(start=r.reference_start-14, end=r.reference_start-1, region='chr'+str(self.generefdf.loc[geneid]['Chromosome'])),'TTTCTTATATGGG') >3 ]


        # #filter according to the cateria of SCAFE
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

        #store reads fetched
        outfilename=self.count_out_dir+'fetch_reads.pkl'
        with open(outfilename,'wb') as f:
            pickle.dump(readinfodict,f)

        return readinfodict



    def distribution(self,geneID,readinfo):
        singletranscriptdict={}
        try:
            generefdf=self.tssrefdf[self.tssrefdf['gene_id']==geneID]
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
                keyid=generefdf.loc[label[i]]['gene_id']+'_'+generefdf.loc[label[i]]['transcript_id']+'@'+generefdf.loc[label[i]]['TSS'].astype('str')
                if len(select[select==label[i]])>=self.minCount:
                    singletranscriptdict[keyid]=(len(select[select==label[i]]),selectCB[select==label[i]])

            if len(singletranscriptdict)>=2:
                return singletranscriptdict
            else:
                return {}
            
        except ValueError:
            return {}
            



    def get_transcriptdict(self,):

        start_time=time.time()
        pool = multiprocessing.Pool(processes=self.nproc)
        readinfodict=self._get_gene_reads()
        transcriptdictls=[]

        for geneID in readinfodict.keys():
            transcriptdictls.append(pool.apply_async(self.distribution,(geneID,readinfodict[geneID])))
        pool.close()
        pool.join()

        transcriptdictls=[res.get() for res in transcriptdictls]

        #print(transcriptdictls)
        transcriptdict = {}
        for d in transcriptdictls:
            transcriptdict.update(d)


        

        print('do annotation',int(time.time()-start_time),'seconds.')


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
        adata=ad.AnnData(finaldf) 

        print('we get adata')    
        
        vardf=pd.DataFrame(adata.var.copy())
        vardf.reset_index(inplace=True)
        vardf['gene_id']=vardf['index'].str.split('_',expand=True)[0]
        vardf=vardf.merge(self.generefdf,on='gene_id')
        vardf.set_index('index',drop=True,inplace=True)
        adata.var=vardf
        sc_output_h5ad=self.count_out_dir+'sc_TSS_count.h5ad'
        adata.write(sc_output_h5ad)
        print('produce h5ad Time elapsed',int(time.time()-ctime),'seconds.')


        return adata


