'''
TODO: Expand search terms on clinicaltrials.gov
'''

import numpy as np
import pandas as pd
from functions import *
from bs4 import BeautifulSoup
import requests
import json
import re
from Bio import Entrez
from Bio import Medline

# Clean the autism trials .csv file, argument = location/name of file. 
df = clean_csv('autism_trials.csv')
print(len(df))
df.to_csv('autism_drug_trials.csv')

'''Overwrite DF, in a new column, if papers exist. List papers separated by '|'
'''

# articleList = []

# Entrez.email = 'nsmccarty3@gmail.com'

# handle = Entrez.esearch(db='pubmed',
#                         retmax='5000',
#                         term='("autism") AND (("2020/09/10"[Date - Publication] : "3000"[Date - Publication])')
# record = Entrez.read(handle)
# handle.close()
# idlist = record['IdList']

# print("ID List completed. Moving to Medline...")

# # Now that we have a list of PubMed IDs for all required articles, we will use Medline to extract information from them. 
# handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline",
#                        retmode="text")
# records = Medline.parse(handle)

# recordList = []
# for record in records:
#     affiliations = []
#     for i in record.get("AD", "NaN"):
#         if i not in affiliations:
#             affiliations.append(i)
#         else:
#             continue

#     for j in record.get("AID", "NaN"):
#         try:
#             start = j.find("10.") + len("10.")
#             end = j.find(" [doi]")
#             doi_string = "10." + j[start:end]
#             # Remove the added [pii in some strings, from the AID call. 
#             new_doi_string = doi_string.replace(" [pii", "")
#         except AttributeError:
#             new_doi_string = "NaN"

#     recordList.append({
#         "title:": str(record.get("TI", "NaN")),
#         "authors": str(record.get("FAU", "NaN")),
#         "doi": "https://doi.org/" + new_doi_string,
#         "affiliations": str(affiliations),
#         "journal_abbrev": str(record.get("TA", "NaN")),
#         "pub_date": str(record.get("DP", "NaN"))
#      })

# print("Medline information extracted. Exporting .csv file with data...")

# # # Generate Pandas DataFrame from list of dictionaries
# articlesPD = pd.DataFrame.from_dict(recordList)
# export_csv = articlesPD.to_csv('pubmed_papers_doi.csv', index = None, header=True)
# print("Task complete. Closing.")