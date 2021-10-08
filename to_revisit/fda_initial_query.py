'''
Notes to Add:
Use urllib.parse to input search terms using key:value dictionary pairs. 

Game plan:
Download all trials related to the search term 'autism'. 
If a clinical trial 'Has Results,' I should visit the URL, navigate to "Study Results" and download more detailed information about age, demographics, and so forth. 
Search PubMed for each trial NCT paper. IF a paper includes that term _exactly_, then grab the DOI of that paper and append it to a new column. 
Every time I run the script, I want to identify all NEW information and add it through the pipeline. 
'''
from bs4 import BeautifulSoup
import requests
import json
import re
from functions import *

print("Running script...")
print("Setting URL...")

# Base URL to search
base_url = 'https://clinicaltrials.gov/ct2/results?cond=autism'

headers = {
    'Connection': 'keep-alive',
    'Pragma': 'no-cache',
    'Cache-Control': 'no-cache',
    'DNT': '1',
    'Upgrade-Insecure-Requests': '1',
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/78.0.3904.97 Safari/537.36',
    'Sec-Fetch-User': '?1',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3',
    'Sec-Fetch-Site': 'none',
    'Sec-Fetch-Mode': 'navigate',
    'Accept-Encoding': 'gzip, deflate, br',
    'Accept-Language': 'ru,en-US;q=0.9,en;q=0.8,tr;q=0.7',
}

params = (
    ('term', 'autism'),
)

with requests.session() as s:
    response = s.get(base_url, headers=headers)
    #, params=params
    result_url = re.search(r'"url": "(.*)",', response.text).group(1)
    columns = re.search(r'"columns":\s+(.*),\s+"columnDefs"', response.text, re.DOTALL)
    #.group(1)
    columns = json.loads(re.sub(r"//\s+\d+", "", columns))
    # get required columns index from columns json.

    data = {
        'draw': '1',
        'columns[0][data]': '0',
        'columns[0][name]': '',
        'columns[0][searchable]': 'true',
        'columns[0][orderable]': 'false',
        'columns[0][search][value]': '',
        'columns[0][search][regex]': 'false',
        'columns[1][data]': '1',
        'columns[1][name]': '',
        'columns[1][searchable]': 'false',
        'columns[1][orderable]': 'false',
        'columns[1][search][value]': '',
        'columns[1][search][regex]': 'false',
        'columns[2][data]': '2',
        'columns[2][name]': '',
        'columns[2][searchable]': 'true',
        'columns[2][orderable]': 'false',
        'columns[2][search][value]': '',
        'columns[2][search][regex]': 'false',
        'columns[3][data]': '3',
        'columns[3][name]': '',
        'columns[3][searchable]': 'true',
        'columns[3][orderable]': 'false',
        'columns[3][search][value]': '',
        'columns[3][search][regex]': 'false',
        'columns[4][data]': '4',
        'columns[4][name]': '',
        'columns[4][searchable]': 'true',
        'columns[4][orderable]': 'false',
        'columns[4][search][value]': '',
        'columns[4][search][regex]': 'false',
        'columns[5][data]': '5',
        'columns[5][name]': '',
        'columns[5][searchable]': 'true',
        'columns[5][orderable]': 'false',
        'columns[5][search][value]': '',
        'columns[5][search][regex]': 'false',
        'columns[6][data]': '6',
        'columns[6][name]': '',
        'columns[6][searchable]': 'true',
        'columns[6][orderable]': 'false',
        'columns[6][search][value]': '',
        'columns[6][search][regex]': 'false',
        'columns[7][data]': '7',
        'columns[7][name]': '',
        'columns[7][searchable]': 'true',
        'columns[7][orderable]': 'false',
        'columns[7][search][value]': '',
        'columns[7][search][regex]': 'false',
        'columns[8][data]': '8',
        'columns[8][name]': '',
        'columns[8][searchable]': 'true',
        'columns[8][orderable]': 'false',
        'columns[8][search][value]': '',
        'columns[8][search][regex]': 'false',
        'columns[9][data]': '9',
        'columns[9][name]': '',
        'columns[9][searchable]': 'true',
        'columns[9][orderable]': 'false',
        'columns[9][search][value]': '',
        'columns[9][search][regex]': 'false',
        'columns[10][data]': '10',
        'columns[10][name]': '',
        'columns[10][searchable]': 'true',
        'columns[10][orderable]': 'false',
        'columns[10][search][value]': '',
        'columns[10][search][regex]': 'false',
        'columns[11][data]': '11',
        'columns[11][name]': '',
        'columns[11][searchable]': 'true',
        'columns[11][orderable]': 'false',
        'columns[11][search][value]': '',
        'columns[11][search][regex]': 'false',
        'columns[12][data]': '12',
        'columns[12][name]': '',
        'columns[12][searchable]': 'true',
        'columns[12][orderable]': 'false',
        'columns[12][search][value]': '',
        'columns[12][search][regex]': 'false',
        'columns[13][data]': '13',
        'columns[13][name]': '',
        'columns[13][searchable]': 'true',
        'columns[13][orderable]': 'false',
        'columns[13][search][value]': '',
        'columns[13][search][regex]': 'false',
        'columns[14][data]': '14',
        'columns[14][name]': '',
        'columns[14][searchable]': 'true',
        'columns[14][orderable]': 'false',
        'columns[14][search][value]': '',
        'columns[14][search][regex]': 'false',
        'columns[15][data]': '15',
        'columns[15][name]': '',
        'columns[15][searchable]': 'true',
        'columns[15][orderable]': 'false',
        'columns[15][search][value]': '',
        'columns[15][search][regex]': 'false',
        'columns[16][data]': '16',
        'columns[16][name]': '',
        'columns[16][searchable]': 'true',
        'columns[16][orderable]': 'false',
        'columns[16][search][value]': '',
        'columns[16][search][regex]': 'false',
        'columns[17][data]': '17',
        'columns[17][name]': '',
        'columns[17][searchable]': 'true',
        'columns[17][orderable]': 'false',
        'columns[17][search][value]': '',
        'columns[17][search][regex]': 'false',
        'columns[18][data]': '18',
        'columns[18][name]': '',
        'columns[18][searchable]': 'true',
        'columns[18][orderable]': 'false',
        'columns[18][search][value]': '',
        'columns[18][search][regex]': 'false',
        'columns[19][data]': '19',
        'columns[19][name]': '',
        'columns[19][searchable]': 'true',
        'columns[19][orderable]': 'false',
        'columns[19][search][value]': '',
        'columns[19][search][regex]': 'false',
        'columns[20][data]': '20',
        'columns[20][name]': '',
        'columns[20][searchable]': 'true',
        'columns[20][orderable]': 'false',
        'columns[20][search][value]': '',
        'columns[20][search][regex]': 'false',
        'columns[21][data]': '21',
        'columns[21][name]': '',
        'columns[21][searchable]': 'true',
        'columns[21][orderable]': 'false',
        'columns[21][search][value]': '',
        'columns[21][search][regex]': 'false',
        'columns[22][data]': '22',
        'columns[22][name]': '',
        'columns[22][searchable]': 'true',
        'columns[22][orderable]': 'false',
        'columns[22][search][value]': '',
        'columns[22][search][regex]': 'false',
        'columns[23][data]': '23',
        'columns[23][name]': '',
        'columns[23][searchable]': 'true',
        'columns[23][orderable]': 'false',
        'columns[23][search][value]': '',
        'columns[23][search][regex]': 'false',
        'columns[24][data]': '24',
        'columns[24][name]': '',
        'columns[24][searchable]': 'true',
        'columns[24][orderable]': 'false',
        'columns[24][search][value]': '',
        'columns[24][search][regex]': 'false',
        'columns[25][data]': '25',
        'columns[25][name]': '',
        'columns[25][searchable]': 'true',
        'columns[25][orderable]': 'false',
        'columns[25][search][value]': '',
        'columns[25][search][regex]': 'false',
        'start': '0',
        'length': '10',
        'search[value]': '',
        'search[regex]': 'false'
    }
    response = s.post(f'https://clinicaltrials.gov/{result_url}', headers=headers, data=data)
    data = response.json()
    print(data)
    # get total records count from recordsTotal field
    # get required data by index got from columns json