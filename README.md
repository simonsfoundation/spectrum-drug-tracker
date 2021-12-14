# Spectrum Drug Tracker
A collection of Python scripts to query the FDA's clinical trials API and filter the data specifically for autism and autism-related conditions.

View _Spectrum's_ interactive [Autism Drug Tracker]().

View the [ObservableHQ notebook](https://observablehq.com/@spectrumnews/drugtracker) that powers the Drug Tracker.

## How to use this code

This repository contains Python files and folders. The files `functions.py` and `nct_scraper.py` are all that is needed to scrape the FDA API and filter data for autism-related conditions. Executing that script, on a local machine running Python3, will result in an exported .csv file that contains your data. The script takes less than one minute to run on a modern MacBook machine.

Additional details are provided for each file below.

## How we filtered the data & what is included

The final data included in the _Spectrum_ Drug Tracker only includes Phase 2, 3, and 4 trials. These data only include trials which were placebo-controlled, too, and which were specifically conducted for autism or a related condition, such as Rett syndrome or Fragile X syndrome. A full list of included conditions:

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

The FDA's clinical trials database is not exhaustive; the first trials were not posted until September 2008.

All data were filtered using the Pandas library (see function `clean_dataframes` for explicit code). Briefly, trials were removed if they did not have a 'Phase' listed, or if the listed 'Phase' was denoted as 'Not Applicable.' Phase 1 trials and joint Phase 1/2 trials were also removed, as were non-drug interventions. This means that trials based on behavioral interventions are not included in the Drug Tracker or underlying datasets.

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
- `SpectrumCoverage`: A list of prior Spectrum articles regarding the clinical trial, where available (manual).
- `PreviouslyApproved`: Whether each drug has been approved by the FDA, for any condition. Typically derived from a simple Google search. Marked as 'Yes' or 'No' for each drug, separated by a '|' character (manual).
- `ApprovedConditions`: If marked as 'Yes' in `PreviouslyApproved`, lists conditions for which the drug was approved. Values may not be exhaustive (manual).
- `CombinedModality`: Whether the trial employed multiple modalities (e.g. combined a drug and a behavioral intervention for an intervention arm, or mixed drugs). Possible values are 'Yes', 'No' and 'Combined drugs' (manual).
- `Placebo`: Whether the study was placebo-controlled, 'Yes' or 'No'. (manual).
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
- `DesignAllocation`: Whether the study was 'Randomized' or 'Non-randomized' (automated).
- `DesignInterventionModel`: The general design of the strategy for assigning interventions to participants in a clinical study. Types of intervention models include: single group assignment, parallel assignment, cross-over assignment, and factorial assignment. [See link for more details](https://clinicaltrials.gov/ct2/about-studies/glossary) (automated).
- `DesignMasking`: A clinical trial design strategy in which one or more parties involved in the trial, such as the investigator or participants, do not know which participants have been assigned which interventions. Types of masking include: open label, single blind masking, and double-blind masking. [See link for more details](https://clinicaltrials.gov/ct2/about-studies/glossary) (automated).
- `DesignPrimaryPurpose`: The intended goal of the clinical trial; this is 'Treatment' or 'Prevention', largely (automated).
- `DesignWhoMasked`: The parties who were masked to who received the intervention or placebo; values are separated by a '| character (automated).
- `DetailedDescription`: An extended description of the clinical trial; see also `BriefSummary` (automated).
- `EligibilityCriteria`: The criteria, used by the trial sponsor, to determine which participants can enroll (automated).
- `EnrollmentCount`: The number of participants enrolled in the trial (automated).
- `EnrollmentType`: Whether the value given in EnrollmentCount is the actual value, or just the expected value (automated).
- `Gender`: The gender of trial participants; 'Male', 'Female', or 'All' (automated).
- `HealthyVolunteers`: Whether the trial accepted only healthy volunteers; typically 'No' (automated).
- `IPDSharing`: Whether there is a plan to make individual participant data (IPD) collected in this study, including data dictionaries, available to other researchers (automated).
- `InterventionArmGroupLabel`: The 'arms' tested in the trial; each arm is separated by a '|' character (automated).
- `InterventionDescription`: A description of how each 'arm' in the trial was actually conducted (automated).
- `InterventionName`: A brief name for each arm; typically this is similar to `InterventionArmGroupLabel` (automated).
- `LastUpdatePostDate`: The date at which the last update for the trial was posted to ClinicalTrials.gov; DD-Mon-YY (automated).
- `LastUpdatePostDateType`: Whether the date provided in `LastUpdatePostDate` was the actual or estimated value (automated).
- `LastUpdateSubmitDate`: The data at which the authors submitted data for their last update; precedes `LastUpdatePostDate` (automated).
- `LeadSponsorClass`: Whether the lead sponsor of the trial was affiliated with a company, academia, or other entity (manual).
- `LeadSponsorName`: The name of the lead sponsor; typically a company or individual associated with an academic lab (automated).
- `LocationCity`: The cities in which the clinical trial was conducted, with values separated by '|' (automated).
- `LocationCountry`: The countries in which the clinical trial was conducted, with values separated by '|' (automated).
- `LeadSponsorCountry`: The country in which the lead sponsor is headquartered, determined by Google searches (manual).
- `LocationFacility`: The facility name(s) in which the clinical trial was performed, with values separated by '|' (automated).
- `LocationState`: The states or, often, cities in which the clinical trial was performed, with values separated by '|' (automated).
- `MaximumAge`: The maximum age of participants in the trial, in years (automated).
- `MinimumAge`: The minimum age of participants in the trial, in years (automated).
- `OfficialTitle`: The official title of the clinical trial, as described by the authors (automated).
- `OrgFullName`: The full name of the lead sponsor's organization or affiliation (automated).
- `OverallOfficialAffiliation`: Full name of the official's, or lead sponsor's representative, affiliation (automated).
- `OverallStatus`: The current status of the clinical trial; 'Completed' or 'Recruiting' are common (automated).
- `OversightHasDMC`: Whether the clinical trial is using a data monitoring committee to track its progress (automated).
- `Phase`: The phase of the trial, as a string (automated).
- `PrimaryCompletionDate`: The date at which the primary outcomes of the clinical trial were or are expected to be completed (automated).
- `PrimaryCompletionDateType`: Whether the `PrimaryCompletionDate` is the actual or expected value (automated).
- `PrimaryOutcomeDescription`: A text-based description of the clinical trial's primary outcome measures; e.g. the goal of the trial (automated).
- `PrimaryOutcomeMeasure`: How the trial is measuring whether the primary outcomes have been achieved (automated).
- `PrimaryOutcomeTimeFrame`: The time frame of the trial for evaluating primary measures (automated).
- `ReferenceCitation`: Relevant citations for the clinical trial, provided by the authors (automated).
- `ResultsFirstPostDate`: The date at which results were first posted, if applicable (automated).
- `ResultsFirstPostDateType`: Whether the `ResultsFirstPostDate` is the actual or estimated value (automated).
- `ResultsFirstSubmitDate`: The date on which the results were submitted by the authors (automated).
- `ResultsFirstSubmitQCDate`: The date on which the results were subjected to quality control measures at the FDA (automated).
- `SecondaryOutcomeDescription`: A text-based description of the clinical trial's secondary outcome measures; e.g. the secondary goal of the trial (automated).
- `SecondaryOutcomeMeasure`: How the trial is measuring whether the secondary outcomes have been achieved (automated).
- `SecondaryOutcomeTimeFrame`: The time frame of the trial for evaluating secondary measures (automated).
- `StartDate`: When the clinical trial began (automated).
- `StartDateType`: Whether the `StartDate` is an actual or expected value (automated).
- `StatusVerifiedDate`: The date on which the responsible party last verified the clinical study information in the entire ClinicalTrials.gov record for the clinical study, even if no additional or updated information is being submitted (automated).
- `StdAge`: The 'type' of participants enrolled in the trial; e.g. Adult, Child, or Older Adult.
- `StudyFirstPostDate`: The date on which the trial was first posted on the ClinicalTrials.gov database (automated).
- `StudyFirstPostDateType`: Whether `StudyFirstPostDate` is an actual or expected value (automated).
- `StudyFirstSubmitDate`: The date on which the trial was first submitted to ClinicalTrials.gov (automated).
- `StudyFirstSubmitQCDate`: The date on which the trial was first submitted for quality control measures (automated).
- `StudyType`: The type of study; almost always 'Interventional' (automated).
- `VersionHolder`: The date on which _Spectrum_ last updated the clinical trial data for our drug tracker (manual).
- `WhyStopped`: Why the study was stopped, if applicable (automated).

## Repository contents
`nct_scraper.py`: 

`functions.py`:
Limitations of API; e.g. limited to 1000 studies, 20 columns at once.

`data-analysis/data_cleaning.ipynb`:

`data-analysis/spectrum_story_analysis.ipynb`:

`datasets`: This folder contains all of the datasets, generated by the _Spectrum_ newsroom, while creating the drug tracker. Each file is a .csv filetype.

## Attributions and License
Code and analysis by Niko McCarty for _Spectrum_.
