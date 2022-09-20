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




class get_TSS_count():
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

        #filter strand invasion
        fastqFile=get_fastq_file(fastqFilePath)
        reads1_umi=[r for r in reads1_umi if editdistance.eval(fastqFile.fetch(start=r.reference_start-14, end=r.reference_start-1, region=str(self.generefdf.loc[geneid]['Chromosome'])),'TTTCTTATATGGG') >3 ]


        #filter according to the cateria of SCAFE
        reads1_umi=[r for r in reads1_umi if editdistance.eval(r.query_sequence[9:14],'ATGGG')<=4]
        reads1_umi=[r for r in reads1_umi if len(r.cigartuples)>=2]
        reads1_umi=[r for r in reads1_umi if (r.cigartuples[0][0]==4)&(r.cigartuples[0][1]>6)&(r.cigartuples[0][1]<20)&(r.cigartuples[1][0]==0)&(r.cigartuples[1][1]>5)]


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
                    
        # if there are more than 3 cluster then select which can assign to annotation position cluster
        if len(finallabel)>=3:
            selectannlabel=[]
            enstdf=self.tssrefdf[self.tssrefdf['gene_id']==geneid]
            enstls=list(enstdf['TSS'])
            for i in range(0,len(finallabel)):
                if any(tss>np.min(posiarray[labels==finallabel[i]])-300 and tss<np.max(posiarray[labels==finallabel[i]])+300 for tss in enstls):
                    selectannlabel.append(finallabel[i])
            if len(selectannlabel)>=2:
                finallabel=selectannlabel 
    

        #require two cluster's distance should reach the requirement.
        if len(finallabel)>=2:
            i=1
            psi=len(posiarray[labels==finallabel[i]])/(len(posiarray[labels==finallabel[i]])+len(posiarray[labels==finallabel[0]]))
            try:
                while (np.abs(np.min(posiarray[labels==finallabel[0]])-np.min(posiarray[labels==finallabel[i]]))<self.clusterDistance) or (psi<self.psi):
                    i=i+1
                    psi=len(posiarray[labels==finallabel[i]])/(len(posiarray[labels==finallabel[i]])+len(posiarray[labels==finallabel[0]]))
                altTSSls=[[posiarray[labels==finallabel[0]],CBarray[labels==finallabel[0]],cigartuplearray[labels==finallabel[0]]]
                #,seqarray[labels==finallabel[0]]
                            ,[posiarray[labels==finallabel[i]],CBarray[labels==finallabel[i]],cigartuplearray[labels==finallabel[i]]]]
                            #,seqarray[labels==finallabel[i]]
            except IndexError:
                altTSSls=[]  
                
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


        transcriptdict={}
        if np.absolute(tssls[0]-np.min(altTSSdict[geneid][0][0]))<300:
            transcriptdict[transcriptls[0]]=(altTSSdict[geneid][row_ind[0]][0],altTSSdict[geneid][row_ind[0]][1],altTSSdict[geneid][row_ind[0]][2])
        else:
            newname1=str(geneid)+'_newTSS_1'
            transcriptdict[newname1]=(altTSSdict[geneid][row_ind[0]][0],altTSSdict[geneid][row_ind[0]][1],altTSSdict[geneid][row_ind[0]][2])

        if np.absolute(tssls[1]-np.min(altTSSdict[geneid][1][0]))<300:
            transcriptdict[transcriptls[1]]=(altTSSdict[geneid][row_ind[1]][0],altTSSdict[geneid][row_ind[1]][1],altTSSdict[geneid][row_ind[1]][2])  
        else:
            newname2=str(geneid)+'_newTSS_2'
            transcriptdict[newname2]=(altTSSdict[geneid][row_ind[1]][0],altTSSdict[geneid][row_ind[1]][1],altTSSdict[geneid][row_ind[1]][2])

        #print(transcriptdict)

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


        tss_output=self.count_out_dir+'temp_tss.pkl'
        with open(tss_output,'wb') as f:
            pickle.dump(transcriptdictls,f)


        # add a classifier

        ## extract some features from reads 
        final_map = ChainMap(*transcriptdictls)
        final_map

        ### unencoded G percentage 
        gpercentls=[]
        for i in final_map.keys():
            newls=final_map[i][2]
            secondls=[j[0] for j in newls]
            per=[('14S' in ele)or('15S' in ele)or('16S' in ele) for ele in secondls]
            gpercentls.append(sum(per)/len(secondls))

        gpercentarray=np.array(gpercentls).reshape(-1,1)
        gpercentarray

        scaler = MinMaxScaler()
        scaler.fit(gpercentarray)
        scale_gpercentls=list(scaler.transform(gpercentarray).flatten())
        scale_gpercentls

        ### summit umi count 
        intsmmit_countls=[np.max(np.unique(final_map[i][0].flatten(),return_counts=True)[1]) for i in final_map.keys()]
        intsmmit_countlog2array=np.log2(np.array(intsmmit_countls)).reshape(-1,1)
        scaler = MinMaxScaler()
        scaler.fit(intsmmit_countlog2array)
        logscale_summitls=list(scaler.transform(intsmmit_countlog2array).flatten())


        ### flank count
        flankcountls=[]
        for i in final_map.keys():
            posls=list(final_map[i][0].flatten())
            summitpos=max(set(posls), key = posls.count)
            pos=np.unique(final_map[i][0].flatten(),return_counts=True)[0]
            count=np.unique(final_map[i][0].flatten(),return_counts=True)[1]
            flank_umi=sum(count[(pos>summitpos-75)&(pos<summitpos+75)])
            flankcountls.append(flank_umi)
        flankcountlog2array=np.log2(np.array(flankcountls)).reshape(-1,1)
        scaler = MinMaxScaler()
        scaler.fit(flankcountlog2array)
        logscale_flankls=list(scaler.transform(flankcountlog2array).flatten())


        ### cluster count 
        clustercountls=[len(final_map[i][0]) for i in final_map.keys()]
        clustercountlog2array=np.log2(np.array(clustercountls)).reshape(-1,1)
        scaler = MinMaxScaler()
        scaler.fit(clustercountlog2array)
        logscale_clusterls=list(scaler.transform(clustercountlog2array).flatten())


        ### create test dataset
        testdf=pd.DataFrame({'unencodeG':scale_gpercentls,'summit_umi_count':logscale_summitls,'cluster_count':logscale_clusterls})
        pathstr=str(Path(os.path.dirname(os.path.abspath(__file__))).parents[1])+'/test/train_fromSCAFE.csv'


        #print(pathstr)

        ### import train dataset
        traindf=pd.read_csv(pathstr)


        ### train model 
        oneX=traindf.iloc[:,1:].values
        oneY=traindf.iloc[:,0].values
        anotherX=testdf.values

        LogisticRegression = linear_model.LogisticRegression(solver='lbfgs')
        LogisticRegression.fit(oneX,oneY)
        anotherpredictY=LogisticRegression.predict(anotherX)

        ### filter according to the result of model
        finalresultdf=pd.DataFrame({'transcript_id':final_map.keys(),'label':anotherpredictY})
        dropdf=finalresultdf[finalresultdf['label']==0]

        newdict={}
        for ele in transcriptdictls:
            newdict["+".join(list(ele.keys()))]=[]
        newdict

        remainingls=[]
        for i in newdict.keys():
            if ( i.split('+')[0] not in list(dropdf['transcript_id'])) and (i.split('+')[1] not in list(dropdf['transcript_id'])):
                remainingls.append(i)

        newnewls=[]
        for ele in transcriptdictls:
            if "+".join(list(ele.keys())) in remainingls:
                newnewls.append(ele) 

        extendls=[]
        for d in newnewls:
            extendls.extend(list(d.items()))

        
        ### organize the output result
        d={'transcript_id':[transcript[0] for transcript in extendls],'TSS_start':[np.min(transcript[1][0]) for transcript in extendls],
        'TSS_end':[np.max(transcript[1][0]) for transcript in extendls]}

        regiondf=pd.DataFrame(d)
        print('do annotation Time elapsed',int(time.time()-start_time),'seconds.')

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
