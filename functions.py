import numpy as np
import pandas as pd
import requests
import re
import io
import time
import csv
from datetime import datetime as dt

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


def build_dataframes():
    response_test = requests.get(test_url, headers=headers).content.decode('utf-8')
    number_of_studies = int(re.search("(?<=NStudiesFound: )\d\d\d\d", response_test).group(0))
    print("Number of studies found: " + str(number_of_studies))

    min_rank = 1
    if 1000 <= number_of_studies:
        max_rank = 1000
    else:
        max_rank = number_of_studies

    print("Max rank set to: " + str(max_rank))

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

    return df_x, df_y, df_z

def clean_dataframes(df):
    # First, remove all trials that do not have a 'Phase' explicitly listed 
    df_phase = df[df['Phase'].notna()]

    # Remove trials that say 'Not Applicable' for Phase data. 
    df_phase = df_phase[~df_phase['Phase'].str.contains("Applicab")]

    # Remove Phase 1 clinical trials 
    df_phase = df_phase[~df_phase['Phase'].str.contains("1")]

    # Remove trials that do not include 'drug' as an intervention type 
    df_phase_drugs = df_phase[df_phase['InterventionType'].str.contains("Drug", na=False)]

    # Perform extra filtering, as some trials appear for certain lymphomas, etc. 
    targets = ['Autis', 'autis', 'Rett', 'Angelman', 'uberous', 'Phelan', '15q', 'Fragile X', 'Asperger', 'Williams', 'Pervasive Developmental Disorder']
    extra_filter_mask = df_phase_drugs['ConditionBrowseLeafAsFound'].str.contains('Autis|autis|Rett|Angelman|uberous|Phelan|15q|Fragile X|Asperger|Williams|Pervasive Developmental Disorder', regex=True, na=False)
    s = pd.Series(extra_filter_mask)
    df_phase_drugs_extra_filter = df_phase_drugs[s.values]

    # Populate a placebo column automatically, based solely on ArmGroupInterventionName column (maybe not flawless)
    df_phase_drugs_extra_filter['Placebo'] = df_phase_drugs_extra_filter['ArmGroupInterventionName'].apply(lambda x: 'Yes' if 'Placebo' in str(x) else 'No')

    # Export to a .csv file, name it as 'master' with the current date.
    return df_phase_drugs_extra_filter