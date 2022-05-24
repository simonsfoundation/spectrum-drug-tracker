'''
IMPORTS
'''
import numpy as np
import pandas as pd
from functions import *
from datetime import datetime as dt

'''
ENTER DATA MANUALLY BEFORE EXECUTING SCRIPT
'''
# Enter the 'prior script execution' date in YYYYMMDD format. 
prior_date = str(20220413)
# Enter the local path to the prior dataset.
current_data = "./datasets/April_2022/20220413_drug_trials_filtered.csv"

'''
COLLECT DATES
'''
# Collect today's date. 
current_date = dt.today().strftime('%d-%b-%y')
current_date_formatted = dt.strptime(current_date, '%d-%b-%y').strftime('%Y%m%d')
print(f"The date today is: {str(current_date)}")

'''
SCRAPE DATA
'''
df = build_dataframes()

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
'''
# UPDATE THESE DATA. Read in the old dataset to compare the new data changes. We cannot use the Google Doc, because those data are selectively updated. Use the old, filtered data.
df_current = pd.read_csv(current_data)

# df_old = pd.read_csv(f'./datasets/December_2021/{prior_date}_drug_trials_filtered.csv')
df_new_NCTIds = df_new[~df_new['NCTId'].isin(df_current['NCTId'])]
df_new['LastUpdatePostDate'] = df_new['LastUpdatePostDate'].astype('datetime64[ns]')
old_date = dt.strptime(prior_date, '%Y%m%d').strftime('%d-%b-%y')
date_mask = df_new['LastUpdatePostDate'] > old_date
df_updated_trials = df_new.loc[date_mask]

print(f"There are {len(df_updated_trials)} updated trials.")
print(f"There are {len(df_new_NCTIds)} new trials.")

'''
Export final .csv files (New and updated trials only)
'''
df_updated_trials.to_csv(f'./datasets/{current_date_formatted}_updated_trials.csv', encoding='utf-8')
df_new_NCTIds.to_csv(f'./datasets/{current_date_formatted}_new_trials.csv', encoding='utf-8')

print("Script completed successfully. Closing...")