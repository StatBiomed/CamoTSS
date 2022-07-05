import brie
import numpy as np
import pandas as pd
from brie.utils import load_samfile, fetch_reads
import os
import pickle
from functools import reduce
import anndata as ad



class get_TSS_count():
    def __init__(self,generefPath,tssrefPath,bamfilePath,outdir,minCount=30,isoformNumber=2):
        self.generefdf=pd.read_csv(generefPath,delimiter='\t')
        self.tssrefdf=pd.read_csv(tssrefPath,delimiter='\t')
        self.bamFile=load_samfile(bamfilePath)
        self.count_out_dir=str(outdir)+'/count/'
        if not os.path.exists(self.count_out_dir):
            os.mkdir(self.count_out_dir)
        self.minCount=minCount
        self.isoformNumber=isoformNumber
        
        

    def _getreads(self,geneid):
        #fetch reads1 in gene 



        reads = fetch_reads(self.bamFile, 'chr'+str(self.generefdf.loc[geneid]['Chromosome']),  self.generefdf.loc[geneid]['Start'] , self.generefdf.loc[geneid]['End'],  trimLen_max=100)
        reads1_umi = reads["reads1"]

        #do filtering，selecting reads1 with cell barcode and umi+barcode
        reads1_umi = [r for r in reads1_umi if r.has_tag('CB')]  
        reads1_umi = [r for r in reads1_umi if r.has_tag('UR')]

        #filtering, unique umi and cell barcode
        _umi = [r.get_tag('CB') + r.get_tag('UR') for r in reads1_umi]
        _umi_uniq, _umi_idx = np.unique(_umi, return_index=True) 
        reads1_umi = [reads1_umi[x] for x in _umi_idx]

        #store start of reads,CB and UR
        reads_info=[(r.positions[0],r.get_tag('CB')) for r in reads1_umi]
        return reads_info
    
    
    def _get_gene_reads(self):
        self.generefdf.set_index('gene_id',inplace=True)
        readinfodict={}
        for i in self.generefdf.index:
            readinfodict[i]=self._getreads(i)
        return readinfodict
    
    
    #need to improve
    def _reads_distribution(self):
        readinfodict=self._get_gene_reads()
        selectdict={}
        for i in readinfodict.keys():
            tempgene=self.tssrefdf[self.tssrefdf['gene_id']==i]
            tempgene.set_index('transcript_id',inplace=True,drop=True)

            for j in tempgene.index:
                position=np.array([t[0] for t in readinfodict[i] ])[([t[0] for t in readinfodict[i] ]>tempgene['down'][j])&([t[0] for t in readinfodict[i] ]<tempgene['up'][j])]
                CBfile=np.array([t[1] for t in readinfodict[i] ])[([t[0] for t in readinfodict[i] ]>tempgene['down'][j])&([t[0] for t in readinfodict[i] ]<tempgene['up'][j])]
                #URfile=np.array([t[2] for t in readinfodict[i] ])[([t[0] for t in readinfodict[i] ]>tempgene['down'][j])&([t[0] for t in readinfodict[i] ]<tempgene['up'][j])]
                #selectdict[j]=(position,CBfile,URfile)
                selectdict[j]=(position,CBfile)
        return selectdict
    
    
    def _get_count(self):
        selectdict=self._reads_distribution()
        #print(selectdict)
        countdict={}
        for i in selectdict.keys():
            countdict[i]=selectdict[i][0].size
        countdf=pd.DataFrame(countdict.items(),columns=['transcript_id','count'])


        countdf=self.tssrefdf.merge(countdf,on='transcript_id')

        #do filtering,select gene which has two isoform, select count highest isoform
        countdf=countdf[countdf['count']>self.minCount]
        countdf=countdf[countdf.groupby('gene_id')['gene_id'].transform('size') > 1]
        countdf=countdf.sort_values(["gene_id", "count"]).groupby("gene_id").tail(self.isoformNumber)
        countdf.drop(labels='cluster',axis=1,inplace=True)

        #select corresponding UR/CB/position
        newdict={new_key:selectdict[new_key] for new_key in countdf.transcript_id}

        count_outputdf=self.count_out_dir+'Alt_TSS_count.tsv' 
        countdf.to_csv(count_outputdf,sep='\t',index=None)

        return countdf,newdict

    
    def produce_sclevel(self):
        countdf,newdict=self._get_count()
        transcriptdfls=[]
        for i in newdict.keys():
            #print(twoInfodict[i][1])
            cellID,count=np.unique(newdict[i][1],return_counts=True)
            transcriptdf=pd.DataFrame({'cell_id':cellID,i:count})
            transcriptdf.set_index('cell_id',inplace=True)
            transcriptdfls.append(transcriptdf)
        
        finaldf=reduce(lambda x,y: x.join(y,how='outer'), transcriptdfls )
        finaldf.fillna(0,inplace=True)
        adata=ad.AnnData(finaldf)  
        vardf=pd.DataFrame(adata.var.copy())
        vardf.reset_index(inplace=True)
        vardf.columns=['transcript_id']
        vardf=vardf.merge(countdf,on='transcript_id')
        adata.var=vardf
        #hope sanpy can fixed index problem and allow users equal directly by using dataframe
        # adata.var['transcript_id']=vardf['transcript_id']
        # adata.var['gene_id']=vardf['gene_id']
        # adata.var['gene_name']=vardf['gene_name']
        # adata.var['Chromosome']=vardf['Chromosome']
        # adata.var['Strand']=vardf['Strand']
        # adata.var['TSS']=vardf['Strand']
        # adata.var['down']=vardf['down']
        # adata.var['up']=vardf['up']  
        # adata.var['count']=vardf['count']   


        sc_output_h5ad=self.count_out_dir+'sc_TSS_count.h5ad'
        adata.write(sc_output_h5ad)
        return adata








    
    
    
    