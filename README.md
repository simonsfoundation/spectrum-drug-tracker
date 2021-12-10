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

- `index`: A unique value for each trial, based on number of trials returned in the search (automated).
- `NCTId`: A unique value assigned to each clinical trial when it is added to the FDA database. Append an NCTId to the end of `https://clinicaltrials.gov/ct2/show/` to view its data (automated).
- `TrialAcronym`: An acronym for the clinical trial, if specified (automated).
- `ArmGroupDescription`: A verbose explanation of each arm in the trial; describes which drugs were given and at which doses (automated).
- `ArmGroupInterventionName`: Lists the drugs used in the trial, and placebo, if indicated (automated).
- `DrugsTested`: Data derived from `ArmGroupInterventionName`; lists solely the experimental drugs tested for each trial, separated by a '|' character (manual).
- `OtherDrugNames`: Other names for the drugs listed in `DrugsTested`. Typically derived from a simple Google search term. Values for each drug are separated by a '|' character (manual).
- `DrugMechanism`: Curated description describing the biological mechanism by which the drug is thought to act on the human body (manual).
- `Sources`: A list of prior Spectrum articles regarding the clinical trial (manual).
- `PreviouslyApproved`: Whether each drug has been approved by the FDA, for any condition. Typically derived from a simple Google search. Marked as 'Yes' or 'No' for each drug, separated by a '|' character (manual).
- `ApprovedConditions`: If marked as 'Yes' in `PreviouslyApproved`, lists conditions for which the drug was approved. Values may not be exhaustive (manual).
- `CombinedModality`: Whether the trial employed multiple modalities (e.g. combined a drug and a behavioral intervention for an intervention arm, or mixed drugs). Possible values are 'Yes', 'No' and 'Combined drugs' (manual).
- `PubMedPapers`: DOI links to papers, listed on PubMed, based on an input `NCTId` (manual).
- `ArmGroupLabel`: Brief descriptions of each arm in the trial (automated).
- `ArmGroupType`: Brief descriptors for each arm in the trial, such as 'Placebo' or 'Experimental' interventions (automated).
- `BriefSummary`: A brief summary of the clinical trial, provided by the trial authors (automated).
- `BriefTitle`: A brief title of the clinical trial, provided by the trial authors (automated).
- `CentralContactEMail`: An email address for the clinical trial's leader, if provided (automated).
- `CentralContactName`: The name of the clinical trial's lead author, if provided (automated).
- `CompletionDate`: The actual or estimated completion date for the clinical trial. See also `CompletionDateType` (automated).
- `CompletionDateType`: Whether the `CompletionDate` provided is the actual or estimated end date (automated).
- `Condition`: The conditions that the clinical trial is intended to treat. Derived from `ConditionBrowseLeafAsFound`, consolidated to a single value (manual).
- `ConditionBrowseLeafAsFound`: The conditions that the clinical trial is intended to treat, provided as a list of values separated by a '|' character (automated).
- `ConditionMeshId`: A unique ID for each value listed in `ConditionBrowseLeafAsFound`; largely deprecated and replaced with ICD-10 values (automated).
- `ConditionMeshTerm`: A standardized 'Mesh' term for each value listed in `ConditionBrowseLeafAsFound` (automated).
- `DesignAllocation`: Whether the study was 'Randomized' or 'Non-randomized' (automated).
- `DesignInterventionModel`: The general design of the strategy for assigning interventions to participants in a clinical study. Types of intervention models include: single group assignment, parallel assignment, cross-over assignment, and factorial assignment. [See link for more details](https://clinicaltrials.gov/ct2/about-studies/glossary) (automated).
- `DesignMasking`: A clinical trial design strategy in which one or more parties involved in the trial, such as the investigator or participants, do not know which participants have been assigned which interventions. Types of masking include: open label, single blind masking, and double-blind masking. [See link for more details](https://clinicaltrials.gov/ct2/about-studies/glossary) (automated).
- `DesignObservationalModel`: N/A
- `DesignPrimaryPurpose`: The intended goal of the clinical trial; this is 'Treatment' or 'Prevention', largely (automated).
- `DesignWhoMasked`: The parties who were masked to who received the intervention or placebo; values are separated by a '| character (automated).
- `DetailedDescription`: An extended description of the clinical trial; see also `BriefSummary` (automated).
- DispFirstPostDate
- DispFirstPostDateType
- DispFirstSubmitDate
- DispFirstSubmitQCDate
- EligibilityCriteria
- EnrollmentCount
- EnrollmentType
- EventGroupDeathsNumAffected
- EventGroupDeathsNumAtRisk
- EventsFrequencyThreshold
- EventsTimeFrame
- FDAAA801Violation
- Gender
- HealthyVolunteers
- IPDSharing
- InterventionArmGroupLabel
- InterventionDescription
- InterventionName
- InterventionOtherName
- InterventionType
- IsFDARegulatedDrug
- LastKnownStatus
- LastUpdatePostDate
- LastUpdatePostDateType
- LastUpdateSubmitDate
- LeadSponsorClass
- LeadSponsorName
- LocationCity
- LocationContactEMail
- LocationContactName
- LocationCountry
- LocationFacility
- LocationState
- LocationStatus
- MaximumAge
- MinimumAge
- OfficialTitle
- OrgClass
- OrgFullName
- OtherEventStatsNumAffected
- OtherEventStatsNumAtRisk
- OtherEventStatsNumEvents
- OtherEventTerm
- OtherOutcomeDescription
- OtherOutcomeMeasure
- OtherOutcomeTimeFrame
- OutcomeAnalysisCILowerLimit
- OutcomeAnalysisCINumSides
- OutcomeAnalysisCIPctValue
- OutcomeAnalysisCIUpperLimit
- OutcomeAnalysisDispersionType
- OutcomeAnalysisDispersionValue
- OutcomeAnalysisPValue
- OutcomeAnalysisParamType
- OutcomeAnalysisParamValue
- OutcomeAnalysisStatisticalMethod
- OutcomeClassDenomCountValue
- OutcomeDenomCountValue
- OutcomeDenomUnits
- OutcomeGroupDescription
- OutcomeGroupTitle
- OutcomeMeasureDescription
- OutcomeMeasureDispersionType
- OutcomeMeasureParamType
- OutcomeMeasurePopulationDescription
- OutcomeMeasureReportingStatus
- OutcomeMeasureTimeFrame
- OutcomeMeasureTitle
- OutcomeMeasureType
- OutcomeMeasurementValue
- OutcomeMeasureUnitOfMeasure
- OutcomeMeasurementLowerLimit
- OutcomeMeasurementSpread
- OutcomeMeasurementUpperLimit
- OverallOfficialAffiliation
- OverallStatus
- OversightHasDMC
- Phase
- PrimaryCompletionDate
- PrimaryCompletionDateType
- PrimaryOutcomeDescription
- PrimaryOutcomeMeasure
- PrimaryOutcomeTimeFrame
- ReferenceCitation
- ReferencePMID
- ResponsiblePartyInvestigatorAffiliation
- ResponsiblePartyInvestigatorFullName
- ResponsiblePartyInvestigatorTitle
- ResponsiblePartyType
- ResultsFirstPostDate
- ResultsFirstPostDateType
- ResultsFirstSubmitDate
- ResultsFirstSubmitQCDate
- SecondaryOutcomeDescription
- SecondaryOutcomeMeasure
- SecondaryOutcomeTimeFrame
- SeriousEventAssessmentType
- SeriousEventNotes
- SeriousEventStatsNumAffected
- SeriousEventStatsNumAtRisk
- SeriousEventStatsNumEvents
- SeriousEventTerm
- StartDate
- StartDateType
- StatusVerifiedDate
- StdAge
- StudyFirstPostDate
- StudyFirstPostDateType
- StudyFirstSubmitDate
- StudyFirstSubmitQCDate
- StudyPopulation
- StudyType
- VersionHolder
- WhyStopped
- Placebo: programmatically determined, using Pandas, by 

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
