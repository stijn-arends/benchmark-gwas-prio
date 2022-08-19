## Description 
* * *
To be able to access whether or not a gene prioritization method has succeeded in prioritizing relevant genes is to compare them to genes that are known to be associated with a relevent human phenotype ontology (HPO) term. HPO is a standardize vocabulary of phenotypic abnormalities in human disease[[1]](#references). There are over 13000 term each describing a phenotypical abnormality. Each of these terms have genes that are known to be involved with that abnormality.

To get an overview of how well a prioritization method prioritized genes is to take a relevant HPO term and check how many of the prioritzed genes overlap with the genes documented for that HPO term. The **Fisher’s exact test** can be used to assess whether this enrichment is statistically significant. In addition, the odd ratios can be used to determine the direction of this association.

The scripts located in this directory can be used to peform fisher exact tests for a list of HPO terms([example](#example-hpo-list)) on the results of a prioritization method. Currently there are five different gene prioritzation methods supported:

1. NetWAS
2. PoPs
3. DEPICT
4. MAGMA
5. Downstreamer


## Getting Started
* * *
To be able to perform the fisher exact tests on the results of a gene prioritization method we need to run the [fisher_exact_test_prio_methods.py](fisher_exact_test_prio_methods.py) script. This script can process multiple results from different traits at the same time, which should be set in the configuration file.

Example:
```bash
python fisher_exact_test_prio_methods.py -c config.yaml -m NetWAS -o results/
```

### Requirements

* configuration file
    * Results of the prioritization method
    * HPO data
    * List of HPO terms
* Name of the prioritization method

### Config file

```yaml
traits:
  Height: "/path/to/height_results.txt"
  IBD: "/path/to/IBD_results.txt"
  PrC: "/path/to/prstcan_results.txt"

hpo_data: "/path/to/hpo_database.txt.gz"
hpo_info: "/path/to/hpo_list.csv
```

### Example HPO list

```csv
GWAS trait,Related HPO term,HPO ID
Body mass index,Abnormality of body mass index,HP:0045081
Height,Abnormality of body height,HP:0000002
Inflammatory bowel disease,increased inflamatory response,HP:0012649
Coeliac disease,increased inflamatory response,HP:0012649
Diastolic blood pressure,Abnormal systemic blood pressure,HP:0030972
Systolic blood pressure,Abnormal systemic blood pressure,HP:0030972
Pulse pressure,Abnormal systemic blood pressure,HP:0030972
Coronary artery disease,Abnormality of the cardiovascular system,HP:0001626
```



## References
* * *
**[1]** Köhler S, Carmody L, Vasilevsky N, et al. Expansion of the Human Phenotype Ontology (HPO) knowledge base and resources. Nucleic Acids Res. 2019;47(D1):D1018-D1027. doi:10.1093/nar/gky1105