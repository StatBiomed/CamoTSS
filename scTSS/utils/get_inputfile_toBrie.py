import pandas as pd
import scanpy as sc
import numpy as np
import anndata as ad
import os



class get_brie_input():
    def __init__(self,rawExpFilePath,splicingFilePath,cellInfoPath,quant_out_dir,cellnumThreshold):
        self.originadata=sc.read_10x_mtx(rawExpFilePath,var_names='gene_symbols')
        self.cellinfodf=pd.read_csv(cellInfoPath,delimiter='\t')
        countpath=os.path.join(splicingFilePath,'count','sc_TSS_count.h5ad')
        refpath=os.path.join(splicingFilePath,'ref_file','ref_TSS.tsv')
        self.adata=sc.read(countpath)
        self.refdf=pd.read_csv(refpath,delimiter='\t')
        self.quant_out_dir=quant_out_dir
        self.cellnumThreshold=cellnumThreshold


    def get_h5adFile(self):
        adata=self.adata

        adata=adata[:,adata.var['gene_id'].isin(self.originadata.var['gene_ids'])]
        originadata=self.originadata[self.originadata.obs.index.isin(adata.obs.index),self.originadata.var['gene_ids'].isin(adata.var['gene_id'])]
        splicedf=pd.DataFrame(adata.X,columns=adata.var.index)
        arr=np.arange(len(splicedf.columns))%2
        isoform1df=splicedf.iloc[:,arr==0]
        isoform2df=splicedf.iloc[:,arr==1]
   
   
        originadata.layers['isoform1']=isoform1df.values
        originadata.layers['isoform2']=isoform2df.values
        originadata.var['isoform1_name']=isoform1df.columns
        originadata.var['isoform2_name']=isoform2df.columns
        originadata.layers['ambiguous']=np.zeros((len(originadata.obs),len(originadata.var)))


        quant_input=self.quant_out_dir+'brie_input.h5ad'
        originadata.write(quant_input)
        return quant_input,originadata



    def get_cluster_cdrFile(self,mode,originadata):

        #pcdf=self.cellinfodf[['cell_id','PC1','PC2','PC3','PC4','PC5']]
        #pcdf=pd.DataFrame(originadata.obs.index,columns=['cell_id'])
        pcdf=self.cellinfodf[['cell_id']]

        #quant_input,addrawadata=self.get_h5adFile()
        cdr=np.array((originadata.X > 0).mean(1))[:,0]
        pcdf['detect_rate']=cdr
        cellinfodf=self.cellinfodf
        

        if mode=='cluster':
            clusterdf=cellinfodf[['cluster']]
            optiondf = pd.get_dummies(clusterdf.cluster, prefix='cluster')
        elif mode=="disease":
            diseasedf=cellinfodf[['disease']]
            optiondf=pd.get_dummies(diseasedf.disease, prefix='disease')
            optiondf=optiondf.iloc[:,0]

        cdrdf=pd.concat([pcdf,optiondf],axis=1)
        cdrdf.set_index('cell_id',inplace=True)
        cdrdf.loc['sum']=cdrdf.sum(axis=0)
        cdrdf=cdrdf.loc[:,cdrdf.loc['sum']>self.cellnumThreshold]



        print(len(cdrdf.columns))

        numls=[i for i in range(1,len(cdrdf.columns)-1)]
        

        cdr_input=self.quant_out_dir+'cdr.tsv'
        cdrdf.to_csv(cdr_input,sep='\t')

        return cdr_input,numls


