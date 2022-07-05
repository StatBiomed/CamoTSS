import pandas as pd
import scanpy as sc
import numpy as np
import anndata as ad



class get_brie_input():
    def __init__(self,rawExpFilePath,splicingFilePath,cellInfoPath,quant_out_dir):
        self.originadata=sc.read_10x_mtx(rawExpFilePath,var_names='gene_symbols')
        self.cellinfodf=pd.read_csv(cellInfoPath,delimiter='\t')
        self.splicingadata=sc.read(splicingFilePath)
        self.quant_out_dir=quant_out_dir

    def get_h5adFile(self):
        originadata=self.originadata[self.originadata.obs.index.isin(self.cellinfodf['cell_id']),:]
        # originadatacelldf=pd.DataFrame(originadata.obs.index,columns=['cell_id'])
        # originadata_adddf=originadatacelldf.merge(cellinfodf,on='cell_id')
        # originadata.obs['cluster']=originadata_adddf['cluster'].values
        # originadata.obs['disease']=originadata_adddf['disease'].values
        # originadata.obs['PC1']=originadata_adddf['PC1'].values
        # originadata.obs['PC2']=originadata_adddf['PC2'].values
        # originadata.obs['PC3']=originadata_adddf['PC3'].values
        # originadata.obs['PC4']=originadata_adddf['PC4'].values
        # originadata.obs['PC5']=originadata_adddf['PC5'].values

        # splicingadata=sc.read(splicingFilePath)
        # print(splicingadata)
        splicingadata=self.splicingadata[self.splicingadata.obs.index.isin(self.cellinfodf['cell_id']),:]
        splicinggenedf=pd.DataFrame(splicingadata.var['gene_id'])
        splicinggenedf.drop_duplicates('gene_id',inplace=True)
        #print(splicinggenedf)
        #print(originadata.var['gene_ids'])
        originadata=originadata[:,originadata.var['gene_ids'].isin(splicinggenedf['gene_id'])]
        #print(originadata)

        splicingadata=splicingadata[:,splicingadata.var['gene_id'].isin(originadata.var['gene_ids'])]
        splicedf=pd.DataFrame(splicingadata.X,columns=splicingadata.var['transcript_id'])
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



    def get_cluster_cdrFile(self,mode):

        pcdf=self.cellinfodf[['cell_id','PC1','PC2','PC3','PC4','PC5']]

        quant_input,addrawadata=self.get_h5adFile()
        cdr=np.array((addrawadata.X > 0).mean(1))[:,0]
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

        cdr_input=self.quant_out_dir+'cdr.tsv'
        cdrdf.to_csv(cdr_input,sep='\t')

        return cdr_input,cdrdf



