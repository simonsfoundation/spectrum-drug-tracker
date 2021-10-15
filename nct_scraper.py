'''
TODO:
Decide on 'final categories' that we want to pull -- I look at A through L
Compartmentalize code. 
Update DataFrame column with date that the script was run
Do sanity checks on my own scripts by manually doing these searches on clinicaltrials.gov, etc. Why are phases so infrequently listed, for example?
Transition to a Jupyter notebook?
Clean up code, obviously. Terrible approach. 
Add code to find "new" entries, based on whether NCTId exists, and when LastUpdate was posted, etc. Find way to highlight those (e.g. new updated column, export to new csv with JUST what's new...)
Search PubMed for NCTId and append DOIs to new column if papers exist?
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

TODO:
Descriptive summary of each drug
Clinical indication for the trial; e.g. core autism trait, co-occurring symptom, etc.
Mechanism of action (e.g. ‘blocks the dopamine D-2 receptor’)
When the trial ended
Whether drug has already been approved for something else and, if so, what. 
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
STUDY FIELDS:
Acronym,ArmGroupDescription,ArmGroupInterventionName,ArmGroupLabel,ArmGroupType,BriefSummary,BriefTitle,CentralContactEMail,CentralContactName,CentralContactRole,
CompletionDate,CompletionDateType,Condition,ConditionBrowseLeafAsFound,ConditionMeshId,ConditionMeshTerm,DesignAllocation,DesignInterventionModel,DesignMasking,
DesignObservationalModel,DesignPrimaryPurpose,DesignTimePerspective,DesignWhoMasked,DetailedDescription,DispFirstPostDate,DispFirstPostDateType,DispFirstSubmitDate,DispFirstSubmitQCDate,
EligibilityCriteria,EnrollmentCount,EnrollmentType,EventGroupDeathsNumAffected,EventGroupDeathsNumAtRisk,EventsFrequencyThreshold,EventsTimeFrame,FDAAA801Violation,Gender,HealthyVolunteers,
IPDSharing,InterventionArmGroupLabel,InterventionDescription,InterventionName,InterventionOtherName,InterventionType,IsFDARegulatedDrug,LastKnownStatus,LastUpdatePostDate,LastUpdatePostDateType,LastUpdateSubmitDate,
LeadSponsorClass,LeadSponsorName,LocationCity,LocationContactEMail,LocationContactName,LocationCountry,LocationFacility,LocationState,LocationStatus,MaximumAge,MinimumAge,
NCTId,OfficialTitle,OrgClass,OrgFullName,OtherEventStatsNumAffected,OtherEventStatsNumAtRisk,OtherEventStatsNumEvents,OtherEventTerm,OtherOutcomeDescription,OtherOutcomeMeasure,OtherOutcomeTimeFrame,
OutcomeAnalysisCILowerLimit,OutcomeAnalysisCINumSides,OutcomeAnalysisCIPctValue,OutcomeAnalysisCIUpperLimit,OutcomeAnalysisDispersionType,OutcomeAnalysisDispersionValue,
OutcomeAnalysisPValue,OutcomeAnalysisParamType,OutcomeAnalysisParamValue,OutcomeAnalysisStatisticalMethod,OutcomeClassDenomCountValue,OutcomeDenomCountValue,
OutcomeDenomUnits,OutcomeGroupDescription,OutcomeGroupTitle,OutcomeMeasureDescription,OutcomeMeasureDispersionType,OutcomeMeasureParamType,OutcomeMeasurePopulationDescription,OutcomeMeasureReportingStatus,OutcomeMeasureTimeFrame,OutcomeMeasureTitle,OutcomeMeasureType,
OutcomeMeasurementValue,OutcomeMeasureUnitOfMeasure,OutcomeMeasurementLowerLimit,OutcomeMeasurementSpread,OutcomeMeasurementUpperLimit,OverallOfficialAffiliation,OverallOfficialName,OverallOfficialRole,OverallStatus,
OversightHasDMC,Phase,PointOfContactEMail,PointOfContactOrganization,PointOfContactTitle,PrimaryCompletionDate,PrimaryCompletionDateType,PrimaryOutcomeDescription,PrimaryOutcomeMeasure,PrimaryOutcomeTimeFrame,
ReferenceCitation,ReferencePMID,ResponsiblePartyInvestigatorAffiliation,ResponsiblePartyInvestigatorFullName,ResponsiblePartyInvestigatorTitle,ResponsiblePartyType,
ResultsFirstPostDate,ResultsFirstPostDateType,ResultsFirstSubmitDate,ResultsFirstSubmitQCDate,SamplingMethod,SecondaryId,SecondaryIdType,SecondaryIdDomain,SecondaryIdLink,
SecondaryOutcomeDescription,SecondaryOutcomeMeasure,SecondaryOutcomeTimeFrame,SeriousEventAssessmentType,SeriousEventNotes,SeriousEventStatsNumAffected,SeriousEventStatsNumAtRisk,SeriousEventStatsNumEvents,SeriousEventTerm,
StartDate,StartDateType,StatusVerifiedDate,StdAge,StudyFirstPostDate,StudyFirstPostDateType,StudyFirstSubmitDate,StudyFirstSubmitQCDate,StudyPopulation,StudyType,VersionHolder,WhyStopped
'''

# DATE OF LAST SCRIPT EXECUTION
prior_date = pd.to_datetime("October 15, 2021")

# Set a search expression, spaces are denoted as '+', but Boolean logic can still be used. 
search_terms = """autism+OR+autism+spectrum+disorder+OR+Fragile+X+OR+Rett+syndrome+OR+tuberous+sclerosis+OR+Williams+syndrome+OR+
                Praeder+Willi+syndrome+OR+Phelan+McDermid+syndrome+OR+Dup15q+OR+Angelman+OR+Timothy+syndrome+OR+16p+deletion+OR+16p+duplication"""

search_fields_1 = ("""Acronym,ArmGroupDescription,ArmGroupInterventionName,ArmGroupLabel,ArmGroupType,BriefSummary,BriefTitle,CentralContactEMail,CentralContactName,CentralContactRole,
CompletionDate,CompletionDateType,Condition,ConditionBrowseLeafAsFound,ConditionMeshId,ConditionMeshTerm,DesignAllocation,DesignInterventionModel,DesignMasking,DesignObservationalModel""").replace(',','%2C')

search_fields_2 = ("""DesignPrimaryPurpose,DesignTimePerspective,DesignWhoMasked,DetailedDescription,DispFirstPostDate,DispFirstPostDateType,DispFirstSubmitDate,DispFirstSubmitQCDate,
EligibilityCriteria,EnrollmentCount,EnrollmentType,EventGroupDeathsNumAffected,EventGroupDeathsNumAtRisk,EventsFrequencyThreshold,EventsTimeFrame,FDAAA801Violation,Gender,HealthyVolunteers,
IPDSharing,InterventionArmGroupLabel""").replace(',','%2C')

search_fields_3 = ("""InterventionDescription,InterventionName,InterventionOtherName,InterventionType,IsFDARegulatedDrug,LastKnownStatus,LastUpdatePostDate,LastUpdatePostDateType,LastUpdateSubmitDate,
LeadSponsorClass,LeadSponsorName,LocationCity,LocationContactEMail,LocationContactName,LocationCountry,LocationFacility,LocationState,LocationStatus,MaximumAge,MinimumAge""").replace(',','%2C')

search_fields_4 = ("""NCTId,OfficialTitle,OrgClass,OrgFullName,OtherEventStatsNumAffected,OtherEventStatsNumAtRisk,OtherEventStatsNumEvents,OtherEventTerm,OtherOutcomeDescription,OtherOutcomeMeasure,OtherOutcomeTimeFrame,OutcomeAnalysisCILowerLimit,OutcomeAnalysisCINumSides,OutcomeAnalysisCIPctValue,OutcomeAnalysisCIUpperLimit,OutcomeAnalysisDispersionType,OutcomeAnalysisDispersionValue,
OutcomeAnalysisPValue,OutcomeAnalysisParamType,OutcomeAnalysisParamValue""").replace(',','%2C')

search_fields_5 = ("""OutcomeAnalysisStatisticalMethod,OutcomeClassDenomCountValue,OutcomeDenomCountValue,
OutcomeDenomUnits,OutcomeGroupDescription,OutcomeGroupTitle,OutcomeMeasureDescription,OutcomeMeasureDispersionType,OutcomeMeasureParamType,OutcomeMeasurePopulationDescription,
OutcomeMeasureReportingStatus,OutcomeMeasureTimeFrame,OutcomeMeasureTitle,OutcomeMeasureType,OutcomeMeasurementValue,OutcomeMeasureUnitOfMeasure,OutcomeMeasurementLowerLimit,OutcomeMeasurementSpread,OutcomeMeasurementUpperLimit,OverallOfficialAffiliation""").replace(',','%2C')

search_fields_6 = ("""OverallOfficialName,OverallOfficialRole,OverallStatus,
OversightHasDMC,Phase,PointOfContactEMail,PointOfContactOrganization,PointOfContactTitle,PrimaryCompletionDate,PrimaryCompletionDateType,PrimaryOutcomeDescription,PrimaryOutcomeMeasure,PrimaryOutcomeTimeFrame,
ReferenceCitation,ReferencePMID,ResponsiblePartyInvestigatorAffiliation,ResponsiblePartyInvestigatorFullName,ResponsiblePartyInvestigatorTitle,ResponsiblePartyType,
ResultsFirstPostDate""").replace(',','%2C')

search_fields_7 = ("""ResultsFirstPostDateType,ResultsFirstSubmitDate,ResultsFirstSubmitQCDate,SamplingMethod,SecondaryId,
SecondaryIdType,SecondaryIdDomain,SecondaryIdLink,SecondaryOutcomeDescription,SecondaryOutcomeMeasure,SecondaryOutcomeTimeFrame,SeriousEventAssessmentType,SeriousEventNotes,SeriousEventStatsNumAffected,SeriousEventStatsNumAtRisk,SeriousEventStatsNumEvents,SeriousEventTerm,
StartDate,StartDateType,StatusVerifiedDate""").replace(',','%2C')

search_fields_8 = ("""StdAge,StudyFirstPostDate,StudyFirstPostDateType,StudyFirstSubmitDate,StudyFirstSubmitQCDate,StudyPopulation,StudyType,VersionHolder,WhyStopped""").replace(',','%2C')


headers = {"User-Agent": "Mozilla/5.0 (X11; CrOS x86_64 12871.102.0) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/81.0.4044.141 Safari/537.36"}
format_type = 'csv'

test_url = """https://clinicaltrials.gov/api/query/study_fields?expr=autism+OR+autism+spectrum+disorder+OR+Fragile+X+OR+Rett+syndrome+OR+tuberous+sclerosis+OR+Williams+syndrome+OR+
                Praeder+Willi+syndrome+OR+Phelan+McDermid+syndrome+OR+Dup15q+OR+Angelman+OR+Timothy+syndrome+OR+16p+deletion+OR+16p+duplication&fields=NCTId&min_rnk=1&max_rnk=3&fmt=csv"""

response_test = requests.get(test_url, headers=headers).content.decode('utf-8')
number_of_studies = int(re.search("(?<=NStudiesFound: )\d\d\d\d", response_test).group(0))
print("Number of studies found: " + str(number_of_studies))

min_rank = 1
if 1000 <= number_of_studies:
    max_rank = 1000
else:
    max_rank = number_of_studies

print("Max rank set to: " + str(max_rank))

def compile_df(min_rank, max_rank):
    base_url_1 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_1}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_2 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_2}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_3 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_3}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_4 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_4}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_5 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_5}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_6 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_6}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_7 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_7}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"
    base_url_8 = f"https://clinicaltrials.gov/api/query/study_fields?expr={search_terms}&fields={search_fields_8}&min_rnk={min_rank}&max_rnk={max_rank}&fmt={format_type}"

    response_1 = requests.get(base_url_1, headers=headers).content
    time.sleep(3)
    response_2 = requests.get(base_url_2, headers=headers).content
    time.sleep(3)
    response_3 = requests.get(base_url_3, headers=headers).content
    time.sleep(3)
    response_4 = requests.get(base_url_4, headers=headers).content
    time.sleep(3)
    response_5 = requests.get(base_url_5, headers=headers).content
    time.sleep(3)
    response_6 = requests.get(base_url_6, headers=headers).content
    time.sleep(3)
    response_7 = requests.get(base_url_7, headers=headers).content
    time.sleep(3)
    response_8 = requests.get(base_url_8, headers=headers).content
    time.sleep(3)

    df1 = pd.read_csv(io.StringIO(response_1.decode('utf-8')), skiprows=10)
    df2 = pd.read_csv(io.StringIO(response_2.decode('utf-8')), skiprows=10)
    df3 = pd.read_csv(io.StringIO(response_3.decode('utf-8')), skiprows=10)
    df4 = pd.read_csv(io.StringIO(response_4.decode('utf-8')), skiprows=10)
    df5 = pd.read_csv(io.StringIO(response_5.decode('utf-8')), skiprows=10)
    df6 = pd.read_csv(io.StringIO(response_6.decode('utf-8')), skiprows=10)
    df7 = pd.read_csv(io.StringIO(response_7.decode('utf-8')), skiprows=10)
    df8 = pd.read_csv(io.StringIO(response_8.decode('utf-8')), skiprows=10)

    return df1, df2, df3, df4, df5, df6, df7, df8

df1x, df2x, df3x, df4x, df5x, df6x, df7x, df8x = compile_df(min_rank, max_rank)
df_x = pd.concat((df1x, df2x, df3x, df4x, df5x, df6x, df7x, df8x), axis=1) # First 1000 studies

print("Length of df_x: " + str(len(df_x)))

if 2000 <= number_of_studies:
    max_rank += 1000
    min_rank += 1000
else:
    min_rank += 1000
    max_rank = number_of_studies

df1y, df2y, df3y, df4y, df5y, df6y, df7y, df8y = compile_df(min_rank, max_rank)
df_y = pd.concat((df1y, df2y, df3y, df4y, df5y, df6y, df7y, df8y), axis=1) # Up to study 2000

print("Max rank set to: " + str(max_rank))
print("Length of df_y: " + str(len(df_y)))


if 3000 <= number_of_studies:
    max_rank += 1000
    min_rank += 1000
else:
    min_rank += 1000
    max_rank = number_of_studies

df1z, df2z, df3z, df4z, df5z, df6z, df7z, df8z = compile_df(min_rank, max_rank)
df_z = pd.concat((df1z, df2z, df3z, df4z, df5z, df6z, df7z, df8z), axis=1) # Up to study 3000

print("Max rank set to: " + str(max_rank))
print("Length of df_z: " + str(len(df_z)))

df = pd.concat([df_x, df_y, df_z], sort=False)
df = df.reset_index()
print(len(df))
print(df.head(5))
print(df.tail(5))

#MOVE THIS DOWN, AFTER APPENDING DOI?
df.to_csv('20211015_drug_trials_unfiltered.csv')

# df.to_csv('unfiltered_trials_all_conditions.csv')



''''''
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