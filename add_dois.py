from pymed import PubMed
import numpy as np
import pandas as pd
from functions import *
from bs4 import BeautifulSoup
import requests
import time
import json

print("Executing script...")
df = pd.read_csv('datasets/20211015_drug_trials_master.csv')
pubmed = PubMed(tool="PubMedSearcher", email="nmccarty@simonsfoundation.org")

NCTid_list = df['NCTId'].tolist()
print("The list of NCTids includes " + str(len(NCTid_list)) + " records.")

articleDois = []
articleNCTids = []

def retrieve_DOIs(NCTid):
    results = pubmed.query(NCTid, max_results=20)

    if pubmed.getTotalResultsCount(NCTid) == 0:
        articleDois.append('None found.')
        articleNCTids.append(NCTid)
    elif pubmed.getTotalResultsCount(NCTid) > 0:
        endstring = ''
        for article in results:
            endstring += article
            endstring += '|'
        articleDois.append(endstring)
    else:
        articleDois.append("Error.")

    articleNCTids.append(NCTid)

    print(articleDois)
    print(articleNCTids)

    # articleList = []
    # articleDoi_temp = []

    # if result:
    #     print(result)
        # for article in result:
        #     articleDict = article.toDict()
        #     articleList.append(articleDict) #now articleList contains the articles, in JSON format, possibly with multiple entries

        #     if len(articleList) == 0:
        #         articleDois.append('None Found')
        #         print("The length is zero")
        #     elif len(articleList) == 1:
        #         articleDois.append({article['doi']})
        #         print("The length is one.")
        #     elif len(articleList) > 1:
        #         articleDoi_temp.append(article['doi'])
        #         del articleDoi_temp[1:]
        #         articleDois.append({article[NCTid]: articleDoi_temp})
        #         print("The length is greater than one.")
        #         pass
        #     else:
        #         continue

    # else:  
    #     articleDois.append('None Found')
    #     print("The length is zero")

    # author_metrics.append({
    #         "headline": str(article_headline),
    #         "hyperlink": str(link),
    #         "article_date": str(article_date),
    #         "word_count": str(total_words),
    #         "number_authors": str(number_authors),
    #         "author_adjusted_words": str(adjusted_words)
    #     })

    # print(list(articleList))
    # print(articleDois)
    # time.sleep(1)
    # print("Loop complete.")

retrieve_DOIs('NCT00844753')

# Generate Pandas DataFrame from list of dictionaries
# articlesPD = pd.DataFrame.from_dict(articleDois)
# articlesPD['NCTids'] = pd.Series([NCTid_list])
# export_csv = articlesPD.to_csv('list_of_pubmed_papers.csv', index = None, header=True) 