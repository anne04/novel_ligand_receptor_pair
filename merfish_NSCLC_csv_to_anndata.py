# 3:25
import pandas as pd
import anndata 
from collections import defaultdict
from scipy import sparse
import numpy as np

sample_name = 'Pt73'
for sample_name in ['Pt43', 'Pt67', 'Pt38']:
    print(sample_name)
    file = pd.read_csv('tx_merfish_'+ sample_name +'.csv')
    cell_id_dict = defaultdict(list)
    for i in range(0, len(file)):
        cell_id_dict[file['cellpose_segmentation'][i]].append([ file['global_x'][i], file['global_y'][i], file['gene'][i]])
    
    cell_metadata = defaultdict(list)
    gene_present = dict()
    coordinate = []
    for cell_id in cell_id_dict:
        gene_dict = defaultdict(list)
        x_list = []
        y_list = []
        for tupple in cell_id_dict[cell_id]:
            gene_dict[tupple[2]].append(1)
            x_list.append(int(tupple[0]))
            y_list.append(int(tupple[1]))
    
        for gene in gene_dict:
            gene_dict[gene] = len(gene_dict[gene])
            gene_present[gene] = 1
    
        center_x = (max(x_list)-min(x_list))/2 + min(x_list)
    
        center_y = (max(y_list)-min(y_list))/2 + min(y_list)
    
        cell_metadata[cell_id] = [gene_dict, center_x, center_y]
        coordinate.append([center_x, center_y])
    
    
    
    cell_id_dict = 0
    gene_ids = list(gene_present.keys())
    cell_id_list = sorted(list(cell_metadata.keys()))
    count_matrix = np.zeros((len(cell_id_list), len(gene_ids)))
    for i in range (0, len(cell_id_list)):
        for j in range(0, len(gene_ids)):
            if gene_ids[j] not in cell_metadata[cell_id_list[i]][0]: # or len(cell_metadata[cell_id_list[i]][0][gene_ids[j]])==0:
                count_matrix[i][j] = 0
            else:
                count_matrix[i][j] = cell_metadata[cell_id_list[i]][0][gene_ids[j]]
    
    
    for i in range(len(cell_id_list)):
        cell_id_list[i] = str(cell_id_list[i])
    
    #convert the count matrix to csr??
    count_matrix = sparse.csr_matrix(count_matrix)
    # make anndata using this count matrix, cell ids as obs and gene ids as var. 
    adata = anndata.AnnData(count_matrix)
    adata.obs_names = cell_id_list 
    adata.var_names = gene_ids
    adata.obsm['spatial'] = np.array(coordinate)
    # save it. 
    adata.write('NSCLC_'+ sample_name + '_cell_vs_gene.h5ad', compression="gzip")
