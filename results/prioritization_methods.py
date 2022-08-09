"""
Module that contains functionality to process data from 
different prioritization methods to be able to perform fisher exact tests.
"""

from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
import scipy.stats as stats

class PrioritizationMethod(ABC):

    @abstractmethod
    def read_data(self):
        pass

    @abstractmethod
    def filter_data(self):
        pass
    
    @abstractmethod
    def get_overlap_genes(self):
        pass

    def get_overlap(self, hpo_data, genes):
        """
        Get the genes overlapping with the HPO database.

        :parameters
        -----------
        hpo_data - pd.DataFrame
            HPO data inside a pandas dataframe
        genes - pd.Series
            Series of gene IDs

        :returns
        --------
        overlap_hpo - pd.DataFrame
            HPO data overlapping with the supplied genes
        overlap_genes - pd.Series
            Genes overlapping with the HPO data
        total_overlap - int
            Total number of overlapping genes
        """
        overlapping_genes_data = hpo_data[hpo_data.index.isin(genes)]
        overlapping_genes = overlapping_genes_data.index
        total_overlap = overlapping_genes.shape[0]

        # Only keep the genes that overlap with HPO
        overlap_genes = genes[genes.isin(overlapping_genes)]
        
        # Only keep releveant HPO data
        overlap_hpo = hpo_data[hpo_data.index.isin(overlapping_genes)]
        return overlap_hpo, overlap_genes, total_overlap