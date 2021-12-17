'''
IMPORTS
'''
import numpy as np
import pandas as pd
from functions import *
import requests
import re
import io
import time
import csv
from datetime import datetime as dt

'''
COLLECTING DATES
'''
# Enter the 'prior script execution' date in YYYYMMDD format. 
prior_date = str(20211109)

# Collect today's date. 
current_date = dt.today().strftime('%-d-%b-%y')
current_date_formatted = dt.strptime(current_date, '%d-%b-%y').strftime('%Y%m%d')
print(f"The date today is: {str(current_date)}")

'''
SCRAPE DATA
'''
df_x, df_y, df_z = build_dataframes()
df = pd.concat([df_x, df_y, df_z], sort=False)
df.reset_index(inplace=True)

# Create an unfiltered .csv file. Export it to the datasets folder.
df.to_csv(f'datasets/{str(current_date_formatted)}_drug_trials_unfiltered.csv')

'''
DATA FILTERING
'''
df_new = clean_dataframes(df)

# Export all filtered trials; useful if you're not looking solely for 'updated' or 'new' trials since last script execution.
df_new.to_csv(f'datasets/{str(current_date_formatted)}_drug_trials_filtered.csv')

'''
Read in the previous DataFrame and identify only studies between 'prior_date' and 'current_date'.
MODIFY SCRIPT / DATASET LOCATIONS PRIOR TO EXECUTING
'''
df_old = pd.read_csv(f'./datasets/November_2021/{prior_date}_drug_trials_filtered.csv')
df_new_NCTIds = df_new[~df_new['NCTId'].isin(df_old['NCTId'])]
df_new['LastUpdatePostDate'] = df_new['LastUpdatePostDate'].astype('datetime64[ns]')
old_date = dt.strptime(prior_date, '%Y%m%d').strftime('%-d-%b-%y')
date_mask = df_new['LastUpdatePostDate'] > old_date
df_updated_trials = df_new.loc[date_mask]
print(f"The length of df_updated_trials is: {len(df_updated_trials)}")
print(f"The length of df_new_NCTIds is: {len(df_new_NCTIds)}")

'''
EXPORT FINAL CSV (Updated trials only)
'''
df_updated_trials.to_csv(f'./datasets/{current_date_formatted}_updated_trials.csv', encoding='utf-8')
df_new_NCTIds.to_csv(f'./datasets/{current_date_formatted}_new_trials.csv', encoding='utf-8')