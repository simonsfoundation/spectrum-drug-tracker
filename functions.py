"""
IMPORTS
"""
import numpy as np
import pandas as pd
import requests
import re
import io
import time
# import csv
from datetime import datetime as dt


'''
STUDY FIELDS TO QUERY
'''
# Set a search expression. Spaces are denoted as '+'. Boolean logic applies. 
search_terms = """autism+OR+autism+spectrum+disorder+OR+Fragile+X+OR+Rett+syndrome+OR+tuberous+sclerosis+OR+Williams+syndrome+OR+
                Praeder+Willi+syndrome+OR+Phelan+McDermid+syndrome+OR+Dup15q+OR+Angelman+OR+Timothy+syndrome+OR+16p+deletion+OR+16p+duplication"""

# Set study fields to query in the FDA's API. Each of these fields is a column in the final DataFrame.
search_fields_1 = ("""NCTId,Acronym,ArmGroupDescription,ArmGroupInterventionName,ArmGroupLabel,ArmGroupType,BriefSummary,BriefTitle,CentralContactEMail,CentralContactName,
CompletionDate,CompletionDateType,Condition,ConditionBrowseLeafAsFound,DesignAllocation,DesignInterventionModel,DesignMasking,DesignPrimaryPurpose,DesignWhoMasked""").replace(',','%2C')

search_fields_2 = ("""DetailedDescription,EligibilityCriteria,EnrollmentCount,EnrollmentType,Gender,HealthyVolunteers,
IPDSharing,InterventionArmGroupLabel,InterventionDescription,InterventionName,InterventionType,LastUpdatePostDate,LastUpdatePostDateType,LastUpdateSubmitDate,
LeadSponsorClass,LeadSponsorName,LocationCity,LocationCountry,LocationFacility""").replace(',','%2C')

search_fields_3 = ("""LocationState,MaximumAge,MinimumAge,OfficialTitle,OrgFullName,OverallOfficialAffiliation,OverallStatus,
OversightHasDMC,Phase,PrimaryCompletionDate,PrimaryCompletionDateType,PrimaryOutcomeDescription,PrimaryOutcomeMeasure,PrimaryOutcomeTimeFrame,
ReferenceCitation,ResponsiblePartyType,
ResultsFirstPostDate""").replace(',','%2C')

search_fields_4 = ("""ResultsFirstPostDateType,ResultsFirstSubmitDate,ResultsFirstSubmitQCDate,SecondaryOutcomeDescription,SecondaryOutcomeMeasure,SecondaryOutcomeTimeFrame,
StartDate,StartDateType,StatusVerifiedDate,StdAge,StudyFirstPostDate,StudyFirstPostDateType,StudyFirstSubmitDate,StudyFirstSubmitQCDate,StudyType,VersionHolder,WhyStopped""").replace(',','%2C')

# Set headers and format_type. CSV is used here, but could be changed to json. 
headers = {"User-Agent": "Mozilla/5.0 (X11; CrOS x86_64 12871.102.0) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/81.0.4044.141 Safari/537.36"}
format_type = 'csv'


"""
FUNCTIONS
"""
def compile_df(min_rank, max_rank):
    """
    This function executes the requests to the clinical trials API. It assembles the 'URLs' to query, saves requests to variables, 
    and then builds dataframes to store the data. The FDA's API is limited in several respects; one can only query 20 study fields at a time. 
    Additionally, one can only pull data on 1000 trials per request. Thus, DataFrames are concatenated in both dimensions to assemble the final
    unfiltered dataset. 

    Returns 4 DataFrames; one for each 'search_fields' string. 
    """
    base_url_1 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_1}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_2 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_2}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_3 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_3}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_4 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_4}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"

    response_1 = requests.get(base_url_1, headers=headers).content
    time.sleep(3)
    response_2 = requests.get(base_url_2, headers=headers).content
    time.sleep(3)
    response_3 = requests.get(base_url_3, headers=headers).content
    time.sleep(3)
    response_4 = requests.get(base_url_4, headers=headers).content
    time.sleep(3)

    df1 = pd.read_csv(io.StringIO(response_1.decode('utf-8')), skiprows=10)
    df2 = pd.read_csv(io.StringIO(response_2.decode('utf-8')), skiprows=10)
    df3 = pd.read_csv(io.StringIO(response_3.decode('utf-8')), skiprows=10)
    df4 = pd.read_csv(io.StringIO(response_4.decode('utf-8')), skiprows=10)

    return df1, df2, df3, df4

def build_dataframes():
    """
    This function determines the number of trials, based on a 'test request' and then executes compile_dataframes iteratively. 
    
    Returns three DataFrames. 
    """

    test_url = """https://clinicaltrials.gov/api/query/study_fields?expr=autism+OR+autism+spectrum+disorder+OR+Fragile+X+OR+Rett+syndrome+OR+tuberous+sclerosis+OR+Williams+syndrome+OR+
                Praeder+Willi+syndrome+OR+Phelan+McDermid+syndrome+OR+Dup15q+OR+Angelman+OR+Timothy+syndrome+OR+16p+deletion+OR+16p+duplication+OR+ADNP+OR+&fields=NCTId&min_rnk=1&max_rnk=3&fmt=csv"""

    response_test = requests.get(test_url, headers=headers).content.decode('utf-8')
    number_of_studies = int(re.search("(?<=NStudiesFound: )\d\d\d\d", response_test).group(0))
    print("Number of studies found: " + str(number_of_studies))

    min_rank = 1
    if 1000 <= number_of_studies:
        max_rank = 1000
    else:
        max_rank = number_of_studies

    print("Max rank set to: " + str(max_rank))

    df1x, df2x, df3x, df4x = compile_df(min_rank, max_rank)
    df_x = pd.concat((df1x, df2x, df3x, df4x), axis=1, ignore_index=False) # First 1000 studies

    print("Length of df_x: " + str(len(df_x)))

    if 2000 <= number_of_studies:
        max_rank += 1000
        min_rank += 1000
    else:
        min_rank += 1000
        max_rank = number_of_studies

    df1y, df2y, df3y, df4y = compile_df(min_rank, max_rank)
    df_y = pd.concat((df1y, df2y, df3y, df4y), axis=1, ignore_index=False) # Up to study 2000

    print("Max rank set to: " + str(max_rank))
    print("Length of df_y: " + str(len(df_y)))


    if 3000 <= number_of_studies:
        max_rank += 1000
        min_rank += 1000
    else:
        min_rank += 1000
        max_rank = number_of_studies

    df1z, df2z, df3z, df4z = compile_df(min_rank, max_rank)
    df_z = pd.concat((df1z, df2z, df3z, df4z), axis=1, ignore_index=False) # Up to study 3000

    print("Max rank set to: " + str(max_rank))
    print("Length of df_z: " + str(len(df_z)))

    return df_x, df_y, df_z

def clean_dataframes(df):
    """
    This function takes a DataFrame as input and cleans the data. Specifically, it returns only trials that are Phase II+, 
    that include Drugs as the intervention (not behavioral interventions), and that contain details on one of the queried conditions. 
    
    It also creates a 'Placebo' column, based on text in the 'ArmGroupInterventionName' column. 
    
    Returns 1 DataFrame. 
    """
    # Remove all trials that do not have a 'Phase' explicitly listed 
    df_phase = df[df['Phase'].notna()]

    # Remove trials that say 'Not Applicable' for Phase data. 
    df_phase = df_phase[~df_phase['Phase'].str.contains("Applicab")]

    # Remove Phase 1 clinical trials 
    df_phase = df_phase[~df_phase['Phase'].str.contains("1")]

    # Remove trials that do not include 'drug' as an intervention type 
    df_phase_drugs = df_phase[df_phase['InterventionType'].str.contains("Drug", na=False)]

    # Remove open label / non-placebo trials
    df_phase_drugs_placebo = df_phase_drugs[~df_phase_drugs['DesignMasking'].str.contains("None", na=False)]

    # Perform extra filtering, as some trials appear for certain lymphomas, etc. 
    targets = ['Autis', 'autis', 'Rett', 'Angelman', 'uberous', 'Phelan', '15q', 'Fragile X', 'Asperger', 'Williams', 'Pervasive Developmental Disorder', 'ADNP']
    extra_filter_mask = df_phase_drugs_placebo['ConditionBrowseLeafAsFound'].str.contains('Autis|autis|Rett|Angelman|uberous|Phelan|15q|Fragile X|Asperger|Williams|Pervasive Developmental Disorder|16p|ADNP', regex=True, na=False)
    s = pd.Series(extra_filter_mask)
    df_phase_drugs_extra_filter = df_phase_drugs_placebo[s.values]

    # Populate a placebo column automatically, based solely on ArmGroupInterventionName column (maybe not flawless)
    df_phase_drugs_extra_filter['Placebo'] = df_phase_drugs_extra_filter['ArmGroupInterventionName'].apply(lambda x: 'Yes' if 'Placebo' in str(x) else 'No')

    # Export to a .csv file, name it as 'master' with the current date.
    return df_phase_drugs_extra_filter