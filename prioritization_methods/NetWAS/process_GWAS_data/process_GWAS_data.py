"""
This module is designed to prepare GWAS summary statistics data for running the
Versatile Gene-based Association Study (VEGAS2) software. 

VEGAS requires as input a GWAS summary statistics file containing rs SNP IDs and p values.

link: https://vegas2.qimrberghofer.edu.au/
"""


import pandas as pd
import numpy as np
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from pathlib import Path

__author__ = "Stijn Arends"
__version__ = "v.01"


class ExtractVEGASColumns:
    
    snp_columns = ["SNP", "rsid", "rs"]
    p_val_columns = ["P", "pvalue", "p-value", "p_value"]
    
    def __init__(self):
        self.snp = {"score":0, "name":None}
        self.pval = {"score":0, "name":None}
    
    def check_p_col(self, val) -> None: 
        """
        Tries to find the column containing the p-values by using the Levenshtein Distance by matching the column name 
        to a set of pre-defined possible p-value column names.
        
        :parameter
        ----------
        val - str
            A column name
        """
        score = process.extractOne(val, self.p_val_columns, scorer=fuzz.WRatio)[1]
        if score > self.pval["score"]:
            self.pval["score"] = score
            self.pval["name"] = val
            
    def check_SNP_col(self, val) -> None: 
        """
        Tries to find the column containing the SNP ID by using the Levenshtein Distance by matching the column name 
        to a set of pre-defined possible SNP ID column names.
        
        :parameter
        ----------
        val - str
            A column name
        """
        score = process.extractOne(val, self.snp_columns, scorer=fuzz.WRatio)[1]
        if score > self.snp["score"]:
            self.snp["score"] = score
            self.snp["name"] = val
    
    def find_column_names(self, df) -> None: 
        """
        Try to find the column names which are needed to run VEGAS.
        """
        for column in df.columns:
            self.check_p_col(column)
            self.check_SNP_col(column)
        
    def get_col_names(self) -> list:
        """
        Get the column names 
        
        :returns
        --------
        snp,pval - list
            Column names of the SNP and pvalue columns
        """
        if self.snp["name"] == None or self.pval["name"] == None:
            raise Exception("Run find column names")
        return [self.snp["name"], self.pval["name"]]



class PrepGWASData:
    
    def __init__(self, file, vegas):
        self.vegas = vegas
        self.df = self.prepare_data(Path(file))

    def prepare_data(self, file):
        """
        Prepare the GWAS summstats data by activating a chain of commands.

        :parameter
        ----------
        file - Path
            GWAS summstats file
        """
        df = self.read_data(file)
        
        df = self.select_columns(df)
        
        # Drop NaN
        df = self.drop_nan(df)
        
        df = self.filter_rs_id(df)
        
        df = self.set_index(df)
        print(df)
        
        
    @staticmethod
    def read_data(file):
        """
        Read a tab seperated file into a pandas data frame

        :parameter
        ----------
        file - Path
            GWAS summstats file
        """
        return pd.read_csv(file, sep="\t", low_memory=False)
    
    def select_columns(self, df):
        """
        Select the SNP ID and p value columns from the file.

        :parameter
        ----------
        df - pd.DataFrame
            GWAS summstats file

        :returns
        --------
        df - pd.DataFrame
            Dataframe only containing SNP and pvalue columns
        """
        self.vegas.find_column_names(df)
    
        cols = self.vegas.get_col_names()
        
        # Select SNP ID and p-val column
        df = df.loc[:, cols]
        
        return df
    
    @staticmethod
    def drop_nan(df):
        """
        Drop all the NaN values in a data frame.

        :parameter
        ----------
        df - pd.DataFrame
            data frame

        :returns
        --------
        df - pd.DataFrame
            data frame without NaN values
        """
        return df.dropna()
    
    @staticmethod
    def filter_rs_id(df):
        """
        Filter out all the rows that do not contain rs SNP identifiers.

        :parameter
        ----------
        df - pd.DataFrame
            data frame

        :returns
        --------
        df - pd.DataFrame
            Dataframe only containing rs SNP identifiers
        """
        return df[df.iloc[:, 0].str.startswith("rs")]
    
    @staticmethod
    def set_index(df):
        """
        Set the SNP column as the index and sort it.

        :parameter
        ----------
        df - pd.DataFrame
            data frame

        :returns
        --------
        df - pd.DataFrame
            Dataframe with the SNP column as the index
        """
        df.set_index(df.iloc[:,0], drop=True, inplace=True)
        df.drop(columns=df.columns[0], axis=1, inplace=True)
        df.sort_index(inplace=True)
        return df