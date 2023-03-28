import pandas as pd
import scanpy as sc
import pickle
import numpy as np
import multiprocessing 



def window_sliding(self,genereads,TSS_start,TSS_end):

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


    #do something with sliding windows algorithm   
    storels=[]
    for i in range(len(TSS) - self.window_size + 1):
        onewindow=TSS[i: i + self.window_size]
        correspondingcount=count[i: i + self.window_size]
        middlecount=correspondingcount[leftIndex]
        foldchange=(middlecount+1)/(sum(correspondingcount)/len(correspondingcount)+1)
        storels.append([onewindow[leftIndex],correspondingcount[leftIndex],foldchange])
        
    foldchangels=[i[2] for i in storels]
    sortindex=sorted(range(len(foldchangels)), key=lambda k: foldchangels[k],reverse=True)
    allsortls=[storels[i] for i in sortindex]

    return allsortls




def get_CTSS(self,windowSize,minCTSSCount):
    oneclusterfilePath=self.count_out_dir+'afterfiltered.csv'
    alloneclusterdf=pd.read_csv(oneclusterfilePath)
    alloneclusterdf['gene_id']=alloneclusterdf['Unnamed: 0'].str.split('*',expand=True)[0]
    alloneclusterdf['TSS_start']=alloneclusterdf['Unnamed: 0'].str.split('*',expand=True)[1].str.split('_',expand=True)[0]
    alloneclusterdf['TSS_end']=alloneclusterdf['Unnamed: 0'].str.split('_',expand=True)[1]
    alloneclusterdf

    readspath=self.count_out_dir+'fetch_reads.pkl'
    with open(readspath,'rb') as f:
        fetchadata=pickle.load(f)



    pool = multiprocessing.Pool(processes=self.nproc)
    allsortls=[]

    for i in range(0,len(alloneclusterdf)):
        
        geneID=alloneclusterdf['gene_id'][i]
        genereads=fetchadata[geneID]
        TSS_start=alloneclusterdf['TSS_start'][i]
        TSS_end=alloneclusterdf['TSS_end'][i]
        allsortls.append(pool.apply_async(self.window_sliding,(genereads,TSS_start,TSS_end)))

    pool.close()
    pool.join()
    results=[res.get() for res in allsortls]

    allsortfddict={}
    for geneid,resls in zip(self.alloneclusterdf['Unnamed: 0'],results):
        allsortfddict[geneid]=resls  

    
    ctssOutPath=self.ctss_out_dir+'CTSS_foldchange.pkl'
    with open(ctssOutPath,'wb') as f:
        pickle.dump(allsortfddict,f)

    return allsortfddict

