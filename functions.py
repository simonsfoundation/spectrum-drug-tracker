from bs4 import BeautifulSoup
import requests
import json
import re

def clinicalTrialsGov(nctid):
    data = defaultdict(list)
    soup = BeautifulSoup(requests.get("https://clinicaltrials.gov/ct2/show/" + nctid + "?displayxml=true").text, "xml")
    subset = ['intervention_type', 'study_type', 'allocation', 'intervention_model', 'primary_purpose', 'masking', 'enrollment', 'official_title', 'condition', 'minimum_age', 'maximum_age', 'gender', 'healthy_volunteers', 'phase', 'primary_outcome', 'secondary_outcome', 'number_of_arms']

    for tag in soup.find_all(subset):
        data['ct{}'.format(tag.name.capitalize())].append(tag.get_text(strip=True))

    for key in data:
        print('{}: {}'.format(key, ', '.join(data[key])))


'''
TODO:
Divide "Study Designs" column into four separate columns, including Allocation, Intervention Model, Masking, Primary Purpose
Clean up the Outcome Measures column
"Return" the filtered and cleaned DataFrame so that I can use it in other .py files. 
'''
import numpy as np
import pandas as pd

def clean_csv(file_location):
    # Load the data
    df = pd.read_csv(file_location, parse_dates=['Start Date', 'Primary Completion Date', 'Completion Date', 'First Posted', 'Results First Posted', 'Last Update Posted'])

    # Filter by drug or biological
    df_drugs = df.loc[(df['Interventions'].str.contains('Drug', na=False)) | (df['Interventions'].str.contains('Biolog', na=False))]

    return df_drugs