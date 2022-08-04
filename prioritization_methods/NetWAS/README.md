# NetWAS - network-wide association study

The following scripts are used to process the results produced by NetWAS:
1. parse_netwas_results.py
2. convert_gene_id_ensembl_id.R

The results of NetWAS need to be processed before they can be used for analysis. First of all NetWAS adds a header to the output files, this needs to removed for it to be able to be used for analysis. This can be done using the `parse_netwas_results.py` script, optionally it can also be used to filter out the gene symbols. 

To check how the python script works you can look at the help function:

```bash
python3 parse_netwas_results.py -h
```

#########################################################  
&#35; HumanBase NetWAS Analysis Results  
&#35;  
&#35; Job id:      d7732f19-916d-4458-97b5-936b8d6345cb  
&#35; Job title:  
&#35; Email:  
&#35; Created:     2017-08-21 17:07:33 EDT  
&#35; GWAS file:   bmi-2012.out.txt  
&#35; GWAS format: vegas  
&#35; Tissue:      adipose_tissue  
&#35; P-value:     0.01  
&#35;  
&#35; Result file format:  
&#35;  
&#35; Column 1) Gene symbol  
&#35; Column 2) Training label: 1 (+, nominally significant p-value)  
&#35; -1 (-, not nominally significant p-value)  
&#35;                           0 (not used in training)  
&#35; Column 3) NetWAS Score: Distance from the SVM separating hyperplane. 
&#35; Positive scores are in the positive direction (more like nominally  
&#35; significant), negative scores  are in the negative direction (more like non-  
&#35; significant)  
#########################################################  
&#35; NetWAS citation:  
&#35; Greene CS*, Krishnan A*, Wong AK*, Ricciotti E, Zelaya RA, Himmelstein   
&#35; DS, Zhang  
&#35; R, Hartmann BM, Zaslavsky E, Sealfon SC, Chasman DI, FitzGerald GA,   
&#35; Dolinski K,  
&#35; Grosser T, Troyanskaya OG. (2015). Understanding multicellular function   
&#35; and disease with human tissue-specific networks. Nature Genetics.   
&#35; 10.1038/ng.3259w.
######################################################### 

| --------- | ------- | --------- | 
| KRT6B     | -1      | 0.561327  |
| EMP1      | -1      | 0.541169  |
| ZBTB41    | -1      | 0.503238  |
| PNPLA8    | -1      | 0.454396  |
| ITGB4     | -1      | 0.440985  |


The next step is to convert the gene symbols