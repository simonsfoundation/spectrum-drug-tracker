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

def diff_pd(df1, df2):
    """Identify differences between two pandas DataFrames"""
    assert (df1.columns == df2.columns).all(), \
        "DataFrame column names are different"
    if any(df1.dtypes != df2.dtypes):
        "Data Types are different, trying to convert"
        df2 = df2.astype(df1.dtypes)
    if df1.equals(df2):
        return None
    else:
        # need to account for np.nan != np.nan returning True
        diff_mask = (df1 != df2) & ~(df1.isnull() & df2.isnull())
        ne_stacked = diff_mask.stack()
        changed = ne_stacked[ne_stacked]
        changed.index.names = ['id', 'col']
        difference_locations = np.where(diff_mask)
        changed_from = df1.values[difference_locations]
        changed_to = df2.values[difference_locations]
        return pd.DataFrame({'from': changed_from, 'to': changed_to},
                            index=changed.index)