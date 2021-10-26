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
articleNCT = []

def retrieve_DOIs(NCTid):
    result = pubmed.query(NCTid, max_results=20)

    articleList = []
    articleDoi_temp = []

    for article in result:
        articleDict = article.toDict()
        articleList.append(articleDict)

    for article in articleList:
        if len(articleList) == 0 or type(articleList) == 'None':
            articleDois.append({u'doi':'None Found'})
        elif len(articleList) == 1:
            articleDois.append({u'doi':article['doi']})
        elif len(articleList) > 1:
            articleDoi_temp.append(article['doi'])
            articleDois.append({u'doi':articleDoi_temp})
        else:
            continue

    time.sleep(1)
    print("Loop for " + str(NCTid) + " complete.")

for item in NCTid_list:
    retrieve_DOIs(item)

print(len(articleDois))

# Generate Pandas DataFrame from list of dictionaries
articlesPD = pd.DataFrame.from_dict(articleDois)
articlesPD['NCTids'] = pd.Series([NCTid_list])
export_csv = articlesPD.to_csv('list_of_pubmed_papers.csv', index = None, header=True) 