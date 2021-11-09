# Imports
import numpy as np
import pandas as pd
from functions import *
import requests
import re
import io
import time
import csv
from datetime import datetime as dt

# ENTER DATES IN YYYYMMDD FORMAT.
current_date = dt.today().strftime('%-d-%b-%y')
current_date_formatted = dt.strptime(current_date, '%d-%b-%y').strftime('%Y%m%d')
print(f"The date today is: {str(current_date)}")
prior_date = str(20211108)

df_x, df_y, df_z = build_dataframes()
df = pd.concat([df_x, df_y, df_z], sort=False)
df.reset_index(inplace=True)

# Create an unfiltered .csv file 
df.to_csv(f'datasets/{str(current_date_formatted)}_drug_trials_unfiltered.csv')

'''
DATA FILTERING
'''

df_new = clean_dataframes(df)
df_new.to_csv(f'datasets/{str(current_date_formatted)}_drug_trials_filtered.csv')

# Read in the previous DataFrame. 
df_old = pd.read_csv(f'./datasets/{prior_date}_drug_trials_master.csv')

# Identify NEW studies, based on NCTId. 
df_new_NCTIds = df_new[~df_new['NCTId'].isin(df_old['NCTId'])]

# Find studies that have been updated, using the df_new database and date ranges. 
df_new['LastUpdatePostDate'] = df_new['LastUpdatePostDate'].astype('datetime64[ns]')
old_date = dt.strptime(prior_date, '%Y%m%d').strftime('%-d-%b-%y')

date_mask = df_new['LastUpdatePostDate'] > old_date
df_updated_trials = df_new.loc[date_mask]
print(f"The length of df_updated_trials is: {len(df_updated_trials)}")
print(f"The length of df_new_NCTIds is: {len(df_new_NCTIds)}")

df_updated_trials.to_csv(f'./datasets/{current_date_formatted}_updated_trials.csv', encoding='utf-8')
df_new_NCTIds.to_csv(f'./datasets/{current_date_formatted}_new_trials.csv', encoding='utf-8')