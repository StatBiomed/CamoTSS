import brie
import numpy as np
import pandas as pd
from brie.utils import fetch_reads
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




class get_TSS_count():
    def __init__(self,generefPath,tssrefPath,bamfilePath,outdir,cellBarcodePath,nproc,minCount=50,maxReadCount=50000,clusterDistance=500,psi=0.1):
        self.generefdf=pd.read_csv(generefPath,delimiter='\t')

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
        self.psi=psi


        

    def _getreads(self,bamfilePath,geneid):

        #print('Get_reads the pid is %s, gene_id=%s' % (getpid(), geneid))
        #fetch reads1 in gene 
        samFile, _chrom = check_pysam_chrom(bamfilePath,'chr'+str(self.generefdf.loc[geneid]['Chromosome']))
        reads = fetch_reads(samFile, _chrom,  self.generefdf.loc[geneid]['Start'] , self.generefdf.loc[geneid]['End'],  trimLen_max=100)
        reads1_umi = reads["reads1"]

        #do filtering，selecting reads1 with cell barcode and umi+barcode
        reads1_umi = [r for r in reads1_umi if r.has_tag('CB')]  
        reads1_umi = [r for r in reads1_umi if r.has_tag('UR')]

        #filtering, unique umi and cell barcode
        _umi = [r.get_tag('CB') + r.get_tag('UR') for r in reads1_umi]
        _umi_uniq, _umi_idx = np.unique(_umi, return_index=True) 
        reads1_umi = [reads1_umi[x] for x in _umi_idx]
        #print(self.cellBarcode)
        #print(reads1_umi.get_tag('CB'))

        #filter according to the CB provided by users
        reads1_umi=[r for r in reads1_umi if r.get_tag('CB') in self.cellBarcode]


        reads_info=[]
        #store start of reads,CB and UR
        reads_info=[(r.positions[0],r.get_tag('CB')) for r in reads1_umi]
        #print(len(reads_info))
        #print(type(reads_info))
        return reads_info
    

        
    def _get_gene_reads(self):
        ctime=time.time()
        pool = multiprocessing.Pool(processes=self.nproc)
        #print(self.nproc)
        self.generefdf.set_index('gene_id',inplace=True)
        bamfilePath=self.bamfilePath
        readinfodict={}
        results=[]


        for i in self.generefdf.index:
            #readinfodict[i]=self._getreads(bamfilePath,i)
            results.append(pool.apply_async(self._getreads,(bamfilePath,i)))
        pool.close()
        pool.join()
        results=[res.get() for res in results]


        for geneid,resls in zip(self.generefdf.index,results):
            readinfodict[geneid]=resls  


        #delete gene whose reads length is smaller than 
        for i in list(readinfodict.keys()):
            if len(readinfodict[i])>self.maxReadCount:
                readinfodict[i]=random.sample(readinfodict[i],self.maxReadCount)
            if len(readinfodict[i])<2:
                del readinfodict[i] 



        # lenls=[]
        # for i in list(readinfodict.keys()):
        #     lenls.append(len(readinfodict[i]))
        # print(lenls)
        # lendf=pd.DataFrame([readinfodict.keys(),lenls])
        # lendf.to_csv('/storage/yhhuang/users/ruiyan/15organ/SRR13075718_scTSS_out/length.csv')  

        print('Get reads Time elapsed',int(time.time()-ctime),'seconds.')

        outfilename=self.count_out_dir+'fetch_reads.pkl'

        with open(outfilename,'wb') as f:
            pickle.dump(readinfodict,f)


        return readinfodict




    def _do_clustering(self,success):

        geneid=success[0]
        readinfodict=success[1]  

        # genlen=self.generefdf[self.generefdf.index==geneid]['End']-self.generefdf[self.generefdf.index==geneid]['Start']
        # print(genlen)

        altTSSls=[]

        clusterModel = AgglomerativeClustering(n_clusters=None,linkage='average',distance_threshold=100)
        posiarray=np.array([t[0] for t in readinfodict[geneid]]).reshape(-1,1)
        CBarray=np.array([t[1] for t in readinfodict[geneid]]).reshape(-1,1)
        clusterModel=clusterModel.fit(posiarray)
        labels=clusterModel.labels_
        label,count=np.unique(labels,return_counts=True)
        selectlabel=label[count>=self.minCount]
        selectcount=count[count>=self.minCount]
        finalcount=selectcount[np.argsort(selectcount)[::-1]]
        finallabel=selectlabel[np.argsort(selectcount)[::-1]]

        posi=[t[0] for t in readinfodict[geneid]]


        if len(selectlabel)>=2:
            i=1
            psi=len(posiarray[labels==finallabel[i]])/(len(posiarray[labels==finallabel[i]])+len(posiarray[labels==finallabel[0]]))
            try:

                while (np.abs(np.min(posiarray[labels==finallabel[0]])-np.min(posiarray[labels==finallabel[i]]))<self.clusterDistance) or (psi<self.psi):
                    i=i+1
                    psi=len(posiarray[labels==finallabel[i]])/(len(posiarray[labels==finallabel[i]])+len(posiarray[labels==finallabel[0]]))

                altTSSls=[[posiarray[labels==finallabel[0]],CBarray[labels==finallabel[0]]],[posiarray[labels==finallabel[i]],CBarray[labels==finallabel[i]]]]
            except IndexError:
                altTSSls=[]
                

            

            #print('start')
            # # print(np.min(altTSSls[0][0]))
            # # print(np.max(altTSSls[0][0]))
            # posls=[posiarray[labels==selectlabel[0]],posiarray[labels==selectlabel[1]]]
            # print(posls)
            # exit(0)


            # distance1=np.max(altTSSls[0][0])-np.min(altTSSls[0][0])
            # distance2=np.max(altTSSls[1][0])-np.min(altTSSls[1][0])
            # if np.max(altTSSls[0][0])<np.min(altTSSls[1][0]):
            #     cluster_distance=np.min(altTSSls[1][0])-np.max(altTSSls[0][0])
            # elif np.max(altTSSls[1][0])<np.min(altTSSls[0][0]):
            #     cluster_distance=np.min(altTSSls[0][0])-np.min(altTSSls[1][0])
            #print(distance)

            # with open('/storage/yhhuang/users/ruiyan/15organ/SRR13075718_scTSS_out/distance_average.txt','a') as f:
            #     f.write('gene\t%s\tdistance1\t%i\tdistance2\t%i\tcluster_distance\t%i\n'%(geneid,distance1,distance2,cluster_distance))


            # cluster1=altTSSls[0][0]
            # cluster2=altTSSls[1][0]

            # with open('/storage/yhhuang/users/ruiyan/15organ/SRR13075718_scTSS_out/altTSSnp.txt','a') as f:
            #     f.write('%\t%i\tdistance2\t%i\tcluster_distance\t%i\n'%(cluster1,cluster2))

            


            # print(np.max(altTSSls[0][0])-np.min(altTSSls[0][0]))
            # print(np.max(altTSSls[1][0])-np.min(altTSSls[1][0]))
            # print(np.min(altTSSls[1][0])-np.max(altTSSls[0][0]))
            # print(np.min(altTSSls[0][0])-np.max(altTSSls[1][0]))
            #print('end')



            # with open('/storage/yhhuang/users/ruiyan/15organ/SRR13075718_scTSS_out/length.txt','a') as f:
            #     len1=len(labels[labels==selectlabel[0]])
            #     len2=len(labels[labels==selectlabel[1]])
            #     f.write('%i\n%i'%(len1,len2))
        return altTSSls




    def _do_hierarchial_cluster(self):
        ctime=time.time()
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
        #print(len(altTSSdict))
        altTSSdict={k: v for k, v in altTSSdict.items() if v}
        #print(altTSSdict)

        # with open('/storage/yhhuang/users/ruiyan/15organ/SRR13075718_scTSS_out/test/altTSSnp_average.pkl','wb') as f:
        #     pickle.dump(altTSSdict,f)


        print('Do clustering Time elapsed',int(time.time()-ctime),'seconds.')
        return altTSSdict



    def _do_anno_and_filter(self,inputpar):
        geneid=inputpar[0]
        altTSSdict=inputpar[1]
        # print(len(altTSSdict))


        #print('Start annotation the pid is %s, gene_id=%s' % (getpid(), geneid))

        # with open('/storage/yhhuang/users/ruiyan/15organ/SRR13075718_scTSS_out/annotation_gene.txt','a') as f:
        #     f.write('%s\n'%geneid)

        
        temprefdf=self.tssrefdf[self.tssrefdf['gene_id']==geneid]
        # print(temprefdf)
        # exit(0)

        cost_mtx=np.zeros((2,temprefdf.shape[0]))
        for i in range(2):
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
        if np.absolute(tssls[0]-np.min(altTSSdict[geneid][0][0]))<100:
            transcriptdict[transcriptls[0]]=(altTSSdict[geneid][row_ind[0]][0],altTSSdict[geneid][row_ind[0]][1])
        else:
            newname1=str(geneid)+'_newTSS_1'
            transcriptdict[newname1]=(altTSSdict[geneid][row_ind[0]][0],altTSSdict[geneid][row_ind[0]][1])

        if np.absolute(tssls[1]-np.min(altTSSdict[geneid][1][0]))<100:
            transcriptdict[transcriptls[1]]=(altTSSdict[geneid][row_ind[1]][0],altTSSdict[geneid][row_ind[1]][1])  
        else:
            newname2=str(geneid)+'_newTSS_2'
            transcriptdict[newname2]=(altTSSdict[geneid][row_ind[1]][0],altTSSdict[geneid][row_ind[1]][1])


        # with open('/storage/yhhuang/users/ruiyan/15organ/SRR13075718_scTSS_out/annotation_gene.txt','a') as f:
        #     f.write('%s\n'%geneid)
        

        #print('Finish annotation the pid is %s, gene_id=%s' % (getpid(), geneid))        
        return transcriptdict



    def _TSS_annotation(self):
        ctime=time.time()
        altTSSdict=self._do_hierarchial_cluster()
        
        inputpar=[]
        for i in altTSSdict:
            inputpar.append((i,altTSSdict))

        pool = multiprocessing.Pool(processes=self.nproc)
        with multiprocessing.Pool(self.nproc) as pool:
            transcriptdictls=pool.map_async(self._do_anno_and_filter,inputpar).get()

        # super_dict = {}
        # for d in transcriptdictls:
        #     for k, v in d.items():  
        #         super_dict.setdefault(k, []).append(v)

        extendls=[]
        for d in transcriptdictls:
            extendls.extend(list(d.items()))
        



        d={'transcript_id':[transcript[0] for transcript in extendls],'TSS_start':[np.min(transcript[1][0]) for transcript in extendls],
        'TSS_end':[np.max(transcript[1][0]) for transcript in extendls]}


        regiondf=pd.DataFrame(d)

        # print(regiondf)

        # exit(0)
        #notice, a dictionary contain a list. so there are only one element. then we can select cellbarcode or TSS
        print('TSS annotation Time elapsed',int(time.time()-ctime),'seconds.')
        return extendls,regiondf




    def get_transcriptls(self,eachtranscript):
        
        transcriptid=eachtranscript[0]
        #print('Get_reads the pid is %s, gene_id=%s' % (getpid(), transcriptid))
        #print(eachtranscript[0])
        #transcriptdf=pd.DataFrame(pd.Series(eachtranscript[1][1][:,0]).value_counts(),columns=[eachtranscript[0]])
        #print(transcriptdf)

        with open('/storage/yhhuang/users/ruiyan/15organ/SRR13075718_scTSS_out/anno_transcript.txt','a') as f:
            f.write('%s\n'%transcriptid)


        cellID,count=np.unique(eachtranscript[1][1],return_counts=True)
        transcriptdf=pd.DataFrame({'cell_id':cellID,transcriptid:count})
        transcriptdf.set_index('cell_id',inplace=True)

        with open('/storage/yhhuang/users/ruiyan/15organ/SRR13075718_scTSS_out/anno_transcript.txt','a') as f:
            f.write('%s\n'%transcriptid)



        #print('finish the pid is %s, gene_id=%s' % (getpid(), transcriptid))
        return transcriptdf




    def produce_sclevel(self):
        ctime=time.time()
        extendls,regiondf=self._TSS_annotation()
        #transcriptdfls=[]
        pool = multiprocessing.Pool(processes=self.nproc)

        with multiprocessing.Pool(self.nproc) as pool:
            transcriptdfls=pool.map_async(self.get_transcriptls,extendls).get()
        print('multi process Time elapsed',int(time.time()-ctime),'seconds.')

        ctime=time.time()

        # for i in range(0,len(extendls)):
        #     transcriptid=extendls[i][0]       
        #     #print(transcriptid)
        #     cellID,count=np.unique(extendls[i][1][1],return_counts=True)
        #     transcriptdf=pd.DataFrame({'cell_id':cellID,transcriptid:count})
        #     transcriptdf.set_index('cell_id',inplace=True)
        #     transcriptdfls.append(transcriptdf)
        #print(len(transcriptdfls))


        finaldf=reduce(lambda x,y: x.join(y,how='outer'), transcriptdfls )
        print('reduce union Time elapsed',int(time.time()-ctime),'seconds.')

        ctime=time.time()

        finaldf.fillna(0,inplace=True)
        adata=ad.AnnData(finaldf)  
        vardf=pd.DataFrame(adata.var.copy())
        vardf.reset_index(inplace=True)
        vardf.columns=['transcript_id']
        #vardf=pd.merge(vardf,self.tssrefdf,on='transcript_id')
        vardf=vardf.join(regiondf.set_index('transcript_id'), on='transcript_id')

        adata.var=vardf

        sc_output_h5ad=self.count_out_dir+'sc_TSS_count.h5ad'
        adata.write(sc_output_h5ad)
        print('produce h5ad Time elapsed',int(time.time()-ctime),'seconds.')
        return adata
