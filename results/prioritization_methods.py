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


class NetWAS(PrioritizationMethod):

    def __init__(self, hpo, fisher):
        self.hpo = hpo
        self.fisher = fisher

    def read_data(self, data):
        """
        Read in data from a CSV file.

        :parameters
        -----------
        data - Path
            File containing the data

        :returns
        --------
        netwas_data - pd.DataFrame
            Data in a data frame
        genes - pd.Series
            Gene IDs
        """
        netwas_data = pd.read_csv(data, sep=",")
        genes = netwas_data.ensemble_id
        return netwas_data, genes

    def get_overlap_genes(self, data, genes):
        """
        Get the data that overlaps with a list of specified genes.

        :parameters
        -----------
        data - pd.DataFrame
            Data
        genes - pd.Series
            List of gene IDs

        :returns
        --------
        overlap_netwas - pd.DataFrame
            Data overlapping with specified genes
        """
        overlap_netwas = data[data.ensemble_id.isin(genes)]
        return overlap_netwas

    def filter_data(self, data, threshold=0.5):
        """
        Get the data that overlaps with a list of specified genes.

        :parameters
        -----------
        data - pd.DataFrame
            Data
        threshold - float
            Threshold to determine what significant is.

        :returns
        --------
        significant_netwas - pd.DataFrame
            Data containing only the significant genes
        significant_genes - pd.Series
            Gene IDs of the significant genes
        """
        significant_netwas = data[data.netwas_score > threshold]
        significant_genes = significant_netwas.ensemble_id
        return significant_netwas, significant_genes