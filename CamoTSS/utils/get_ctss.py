import numpy as np
import pandas as pd
import os
import pickle
import anndata as ad
import multiprocessing 
import time
import pickle
import warnings
import random



warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=Warning)


class get_CTSS_count():


    def __init__(self,out_dir,minCTSSCount,minFC,n_proc,windowSize):
        self.out_dir=out_dir
        generefPath=os.path.join(out_dir,'ref_file','ref_gene.tsv')
        fetched_reads_path=os.path.join(out_dir,'count','fetch_reads.pkl')
        afterfiltered_path=os.path.join(out_dir,'count','afterfiltered.csv')

        self.generefdf=pd.read_csv(generefPath,delimiter='\t')
        with open(fetched_reads_path,'rb') as f:
            self.fetchadata=pickle.load(f)
        self.alloneclusterdf=pd.read_csv(afterfiltered_path)

        self.minCTSSCount=minCTSSCount
        self.minFC=minFC
        self.n_proc=n_proc
        self.windowSize=windowSize



    def window_sliding(self,genereads,TSS_start,TSS_end,strand):


        leftIndex=0
        filterls=[i for i in genereads if ('14S' in i[2]) or ('15S' in i[2]) or ('16S' in i[2])  ]
        #print("the pid of current process：",os.getpid())
        
        #calculate the TSS position and corresponding counts
        promoterTSS = [read[0] for read in filterls if (read[0] >= TSS_start) and (read[0] <= TSS_end)]
        TSS,count=np.unique(promoterTSS,return_counts=True)



        originaldf=pd.DataFrame({'count':count},index=TSS)
        desireddf=pd.DataFrame({'count':0},index=[i for i in range(int(TSS_start),int(TSS_end))])
        desireddf['count']=desireddf.index.map(originaldf['count']).fillna(0)
        desireddf.reset_index(inplace=True)
        nonzeroarray=np.array(desireddf)


        if strand=='+':
            sortfinalarray=nonzeroarray[nonzeroarray[:, 0].argsort()]
        elif strand=='-':
            sortfinalarray=nonzeroarray[nonzeroarray[:, 0].argsort()[::-1]]
        
        TSS=sortfinalarray.T[0]
        count=sortfinalarray.T[1]


        #do something with sliding windows algorithm   
        storels=[]
        for i in range(len(TSS) - self.windowSize + 1):
            #print(i)
            onewindow=TSS[i: i + self.windowSize]
            correspondingcount=count[i: i + self.windowSize]
            middlecount=correspondingcount[leftIndex]
            foldchange=(middlecount+1)/(sum(correspondingcount)/len(correspondingcount)+1)
            storels.append([onewindow[leftIndex],correspondingcount[leftIndex],foldchange])

            
        foldchangels=[i[2] for i in storels]
        sortindex=sorted(range(len(foldchangels)), key=lambda k: foldchangels[k],reverse=True)
        allsortls=[storels[i] for i in sortindex]

        return allsortls




    def _get_CTSS(self):
        alloneclusterdf=self.alloneclusterdf
        alloneclusterdf['gene_id']=alloneclusterdf['Unnamed: 0'].str.split('*',expand=True)[0]
        alloneclusterdf['TSS_start']=alloneclusterdf['Unnamed: 0'].str.split('*',expand=True)[1].str.split('_',expand=True)[0].astype('float')
        alloneclusterdf['TSS_end']=alloneclusterdf['Unnamed: 0'].str.split('_',expand=True)[1].astype('float')

        self.generefdf.reset_index(inplace=True)
        

        #print(self.generefdf)
        stranddf=self.generefdf[['Strand','gene_id']]
        alloneclusterdf=alloneclusterdf.merge(stranddf,on='gene_id')


        start_time=time.time()

        #pool = multiprocessing.Pool(processes=self.n_proc)
        allsortfddict={}

        for i in range(0,len(alloneclusterdf)):
            geneID=alloneclusterdf['gene_id'][i]
            #print(geneID)
            genereads=self.fetchadata[geneID]
            clusterID=alloneclusterdf['Unnamed: 0'][i]
            TSS_start=alloneclusterdf['TSS_start'][i]
            TSS_end=alloneclusterdf['TSS_end'][i]
            strand=alloneclusterdf['Strand'][i]
            windowreturn=self.window_sliding(genereads,TSS_start,TSS_end,strand)
            allsortfddict[clusterID]=windowreturn


        print('window sliding Time elapsed',int(time.time()-start_time),'seconds.')


        # allsortfddict={}
        # for geneid,resls in zip(alloneclusterdf['Unnamed: 0'],results):
        #     allsortfddict[geneid]=resls  

        
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

 
        #print(self.out_dir)
        self.ctss_out_dir=os.path.join(self.out_dir+'/CTSS/')
        if not os.path.exists(self.ctss_out_dir):
            os.mkdir(self.ctss_out_dir)




        allsortfddict=self._get_CTSS()
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
                for ele in self.fetchadata[geneid]:
                    if j[0]==ele[0]:
                        cellIDls.append(ele[1])
                cellIDdict[newid]=cellIDls

        #print(len(cellIDdict))


        #create a big matrix including cell ID
        cellidls=list(cellIDdict.values())
        cellidset = set([item for sublist in cellidls for item in sublist])
        ctssfinaldf=pd.DataFrame(index=list(cellidset))



        
        for clusterID in cellIDdict.keys():
            cellID,count=np.unique(cellIDdict[clusterID],return_counts=True)
            CTSSdf=pd.DataFrame({'cell_id':cellID,clusterID:count})
            CTSSdf.set_index('cell_id',inplace=True)
            ctssfinaldf[clusterID]=ctssfinaldf.index.map(CTSSdf[clusterID])


        ctssfinaldf.fillna(0,inplace=True)
        #print(ctssfinaldf)
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
    
