# Spectrum Drug Tracker
A collection of Python scripts to query the FDA's clinical trials API and filter the data specifically for autism and autism-related conditions.

View the [Spectrum Drug Tracker]().

View the [ObservableHQ notebook] that powers the Drug Tracker. 

## How to use this code

This repository contains two Python files and two folders. The files `functions.py` and `nct_scraper.py` are all that is needed to scrape the FDA API and filter data for autism-related conditions. Executing that script, on a local machine running Python3, will result in an exported .csv file that contains your data. The script takes less than one minute to run on a modern MacBook machine.

Additional details are provided for each file below.

## How we filtered the data & what is included

The final data included in the _Spectrum_ Drug Tracker only include Phase 2, 3, and 4 trials. They only include trials which were placebo-controlled, too, and which were specifically conducted for autism or a related condition, such as Rett syndrome or Fragile X syndrome. A full list of included conditions:

- Autism (sometimes also referred to as pervasive developmental disorder or autism spectrum disorder) 
- Fragile X syndrome 
- Rett syndrome
- Tuberous sclerosis complex
- Williams syndrome 
- Praeder-Willi syndrome
- Phelan McDermid syndrome
- Dup15q
- Angelman syndrome
- Timothy syndrome
- 16p deletion
- 16p duplication

The FDA's clinical trials database is not exhaustive; the first trials appeared in the mid-2000s, and there are obviously trials that precede those dates. 

All data were filtered using the Pandas library (see function `clean_dataframes` for explicit code). Briefly, trials were removed if they did not have a 'Phase' listed, or if the listed 'Phase' was denoted as 'Not Applicable.' Phase 1 trials were also removed, as were non-drug interventions. This means that trials based on behavioral interventions are not included in the Drug Tracker or underlying datasets.

## Data columns, demystified

Not all data columns are included in the final Drug Tracker. The following data columns are those scraped using the Python code, or manually curated by the _Spectrum_ newsroom. Some of these columns may appear with spaces (rather than as camelCase) on the Drug Tracker, due to how the data is imported into the web application.

index
NCTId
TrialAcronym
ArmGroupDescription
ArmGroupInterventionName
DrugsTested
OtherDrugNames
DrugMechanism
PreviouslyApproved
ApprovedConditions
CombinedModality
PubMedPapers
ArmGroupLabel
ArmGroupType
BriefSummary
BriefTitle
CentralContactEMail
CentralContactName
CompletionDate
CompletionDateType
Condition
ConditionBrowseLeafAsFound
ConditionMeshId
ConditionMeshTerm
DesignAllocation
DesignInterventionModel
DesignMasking
DesignObservationalModel
DesignPrimaryPurpose
DesignWhoMasked
DetailedDescription
DispFirstPostDate
DispFirstPostDateType
DispFirstSubmitDate
DispFirstSubmitQCDate
EligibilityCriteria
EnrollmentCount
EnrollmentType
EventGroupDeathsNumAffected
EventGroupDeathsNumAtRisk
EventsFrequencyThreshold
EventsTimeFrame
FDAAA801Violation
Gender
HealthyVolunteers
IPDSharing
InterventionArmGroupLabel
InterventionDescription
InterventionName
InterventionOtherName
InterventionType
IsFDARegulatedDrug
LastKnownStatus
LastUpdatePostDate
LastUpdatePostDateType
LastUpdateSubmitDate
LeadSponsorClass
LeadSponsorName
LocationCity
LocationContactEMail
LocationContactName
LocationCountry
LocationFacility
LocationState
LocationStatus
MaximumAge
MinimumAge
OfficialTitle
OrgClass
OrgFullName
OtherEventStatsNumAffected
OtherEventStatsNumAtRisk
OtherEventStatsNumEvents
OtherEventTerm
OtherOutcomeDescription
OtherOutcomeMeasure
OtherOutcomeTimeFrame
OutcomeAnalysisCILowerLimit
OutcomeAnalysisCINumSides
OutcomeAnalysisCIPctValue
OutcomeAnalysisCIUpperLimit
OutcomeAnalysisDispersionType
OutcomeAnalysisDispersionValue
OutcomeAnalysisPValue
OutcomeAnalysisParamType
OutcomeAnalysisParamValue
OutcomeAnalysisStatisticalMethod
OutcomeClassDenomCountValue
OutcomeDenomCountValue
OutcomeDenomUnits
OutcomeGroupDescription
OutcomeGroupTitle
OutcomeMeasureDescription
OutcomeMeasureDispersionType
OutcomeMeasureParamType
OutcomeMeasurePopulationDescription
OutcomeMeasureReportingStatus
OutcomeMeasureTimeFrame
OutcomeMeasureTitle
OutcomeMeasureType
OutcomeMeasurementValue
OutcomeMeasureUnitOfMeasure
OutcomeMeasurementLowerLimit
OutcomeMeasurementSpread
OutcomeMeasurementUpperLimit
OverallOfficialAffiliation
OverallStatus
OversightHasDMC
Phase
PrimaryCompletionDate
PrimaryCompletionDateType
PrimaryOutcomeDescription
PrimaryOutcomeMeasure
PrimaryOutcomeTimeFrame
ReferenceCitation
ReferencePMID
ResponsiblePartyInvestigatorAffiliation
ResponsiblePartyInvestigatorFullName
ResponsiblePartyInvestigatorTitle
ResponsiblePartyType
ResultsFirstPostDate
ResultsFirstPostDateType
ResultsFirstSubmitDate
ResultsFirstSubmitQCDate
SecondaryOutcomeDescription
SecondaryOutcomeMeasure
SecondaryOutcomeTimeFrame
SeriousEventAssessmentType
SeriousEventNotes
SeriousEventStatsNumAffected
SeriousEventStatsNumAtRisk
SeriousEventStatsNumEvents
SeriousEventTerm
StartDate
StartDateType
StatusVerifiedDate
StdAge
StudyFirstPostDate
StudyFirstPostDateType
StudyFirstSubmitDate
StudyFirstSubmitQCDate
StudyPopulation
StudyType
VersionHolder
WhyStopped
Placebo - programmatically determined, using Pandas, by 

 (automated)
 (manually added)

## Repository contents
`nct_scraper.py`

`functions.py`
Limitations of API; e.g. limited to 1000 studies, 20 columns at once.

`data-analysis/data_cleaning.ipynb`

`data-analysis/spectrum_story_analysis.ipynb`

`datasets`

## Attributions and License
Code and analysis by Niko McCarty

Contributions
