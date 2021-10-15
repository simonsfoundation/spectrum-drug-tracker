

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

from pymed import PubMed
import numpy as np
import pandas as pd
from functions import *
from bs4 import BeautifulSoup
import requests
import time

print("Executing script...")
df = pd.read_csv('datasets/20211015_drug_trials_master.csv')
pubmed = PubMed(tool="PubMedSearcher", email="nmccarty@simonsfoundation.org")

NCTid_list = df['NCTId'].tolist()
print("The list of NCTids includes " + str(len(NCTid_list)) + " records.")

articleDois = []

def retrieve_DOIs(NCTid):
    result = pubmed.query(NCTid, max_results=20)

    articleList = []
    articleDoi_temp = []

    for article in result:
        articleDict = article.toDict()
        articleList.append(articleDict)

    for article in articleList:
        if len(articleList) == 0:
            articleDois.append({u'doi':'None Found'})
        elif len(articleList) == 1:
            articleDois.append({u'doi':article['doi']})
        elif len(articleList) > 1:
            articleDoi_temp.append(article['doi'] + '|')
            articleDois.append({u'doi':articleDoi_temp})
        else:
            continue
    time.sleep(1)
    print("Loop complete.")

for item in NCTid_list:
    retrieve_DOIs(item)

# Generate Pandas DataFrame from list of dictionaries
articlesPD = pd.DataFrame.from_dict(articleDois)
articlesPD['NCTids'] = pd.Series([NCTid_list])
export_csv = articlesPD.to_csv('list_of_pubmed_papers.csv', index = None, header=True) 