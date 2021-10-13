'''
TODO:
Add for loop to make additional requests, and append all lists as needed. Store for loop in a dictionary to run more easily. 
Compartmentalize code. 
Update DataFrame column with date that the script was run
Do sanity checks on my own scripts by manually doing these searches on clinicaltrials.gov, etc. Why are phases so infrequently listed, for example?
Transition to a Jupyter notebook?
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
import io
import time

'''
SECTION 1: RETRIEVE ALL CLINICAL TRIALS FOR AUTISM AND AUTISM-RELATED CONDITIONS
-- Autism + other terms
-- Fragile X
-- Rett syndrome
-- Tuberous sclerosis
-- Williams syndrome
-- Praeder Willi syndrome
-- Phelan McDermid syndrome
-- Dup15q syndrome 
-- Angelman syndrome
-- Timothy syndrome
-- 16p deletion syndrome
-- 16p duplication syndrome
'''

# Set a search expression, spaces are denoted as '+', but Boolean logic can still be used. 
search_terms = """autism+OR+autism+spectrum+disorder+OR+Fragile+X+OR+Rett+syndrome+OR+tuberous+sclerosis+OR+Williams+syndrome+OR+
                Praeder+Willi+syndrome+OR+Phelan+McDermid+syndrome+OR+Dup15q+OR+Angelman+OR+Timothy+syndrome+OR+16p+deletion+OR+16p+duplication"""

search_fields_1 = ("""Acronym,ArmGroupDescription,ArmGroupInterventionName,ArmGroupLabel,ArmGroupType,BriefSummary,BriefTitle,CentralContactName,CompletionDate,
                    CompletionDateType,Condition,ConditionMeshTerm,DesignAllocation,DesignInterventionModel,
                    DesignInterventionModelDescription,DesignMasking,DesignMaskingDescription,
                    DesignObservationalModel,DesignPrimaryPurpose,DesignWhoMasked""").replace(',','%2C')

search_fields_2 = ("""DetailedDescription,
                DispFirstSubmitDate,EligibilityCriteria,EnrollmentCount,EnrollmentType,Gender,
                InterventionArmGroupLabel,InterventionName,IsFDARegulatedDevice,LastUpdateSubmitDate,
                LeadSponsorClass,LeadSponsorName,LocationCity,LocationContactEMail,LocationContactName,
                MaximumAge,MinimumAge,StdAge,NCTId,OfficialTitle""").replace(',','%2C')

search_fields_3 = ("""OrgClass,OrgFullName,OrgStudyId,OutcomeMeasureDescription,
                OutcomeMeasureTitle,OutcomeMeasureType,Phase,ReferenceCitation,ReferencePMID,ResponsiblePartyInvestigatorAffiliation,
                ResponsiblePartyInvestigatorFullName,ResponsiblePartyInvestigatorTitle,ResponsiblePartyOldOrganization,
                ResponsiblePartyOldNameTitle,ResultsFirstSubmitDate,StudyPopulation,StudyType""").replace(',','%2C')

headers = {"User-Agent": "Mozilla/5.0 (X11; CrOS x86_64 12871.102.0) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/81.0.4044.141 Safari/537.36"}
format_type = 'csv'

min_rank = 1
max_rank = 1000

def compile_df():
    base_url_1 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_1}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_2 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_2}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_3 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_3}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    
    response_1 = requests.get(base_url_1, headers=headers).content
    time.sleep(1)
    response_2 = requests.get(base_url_2, headers=headers).content
    time.sleep(1)
    response_3 = requests.get(base_url_3, headers=headers).content

    df1 = pd.read_csv(io.StringIO(response_1.decode('utf-8')), error_bad_lines=False, skiprows=10)
    df2 = pd.read_csv(io.StringIO(response_2.decode('utf-8')), error_bad_lines=False, skiprows=10)
    df3 = pd.read_csv(io.StringIO(response_3.decode('utf-8')), error_bad_lines=False, skiprows=10)

    return df1, df2, df3

# df1, df2, df3 = compile_df()
# df_x = pd.concat((df1, df2, df3), axis=1)

# min_rank += 1000
# max_rank += 1000
# df1, df2, df3 = compile_df()
# df_y = pd.concat((df1, df2, df3), axis=1)

# min_rank += 1000
# max_rank += 1000
# df1, df2, df3 = compile_df()
# df_z = pd.concat((df1, df2, df3), axis=1)

# min_rank += 1000
# max_rank += 1000
# df1, df2, df3 = compile_df()
# df_xy = pd.concat((df1, df2, df3), axis=1)

# df = df_x.append(df_y, df_z, df_xy)
# print(len(df))
# print(df.head(5))
# print(df.tail(20))

# ADD FOR LOOP TO GO THROUGH MORE RANKS, BASED ON NUMBER OF AVAILABLE STUDIES


'''
TODO:
Descriptive summary of each drug
Clinical indication for the trial; e.g. core autism trait, co-occurring symptom, etc.
Mechanism of action (e.g. ‘blocks the dopamine D-2 receptor’)
When the trial ended
Trial design (e.g. Double blind, randomized, placebo, etc.)
Broad class of drug (e.g. antipsychotic, antidepressant)
Link to PubMed papers that include NCT number, where possible
(Optional: Number of intervention doses given in the trial)

COMPLETED:
Name of drug
Type of drug (e.g. biological)
Index (via Rank)
NCT number
Condition designed to treat (e.g. autism, Phelan McDermid syndrome…)
Phase of trial
Name of trial
When the trial began
Last trial update
Design intervention model
Design Allocation
Design Masking
Name of organization
Brief summary of trial
Detailed description of trial
Intervention description
Overall Status (e.g. Recruiting)
How many people were enrolled
Outcome measures
Minimum age
Maximum age
Indication / what it’s supposed to treat
Location of trial

'''

'''
SECTION 2: CLEAN AND FILTER THE DATAFRAME, BASED ON OUR CRITERIA. ADD 'DATE OF SCRIPT' COLUMN. Only add row if NCT not in original DataFrame column. 
Must be a DRUG, and not OTHER or Behavioral, etc. 
'''

'''
FILTER CRITERIA:
-- Phase 2, 3, 4
-- Only include single-arm trial if there’s a suitable outcome measure (e.g. fMRI)
-- Include combined modality trials (e.g. NeuroNext, oxytocin + behavioral intervention)
-- Sort DataFrame by drugs (after I tease out that column / make a Placebo column) and compare to existing drug database
'''


# Clean the autism trials .csv file, argument = location/name of file. 
# df = clean_csv('datasets/autism_trials.csv')
# print(len(df))
# df.to_csv('datasets/autism_drug_trials.csv')

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