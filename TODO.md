2020-04-26
- add: fix biomet gapfill plots, do not appear to be writing to PDF properly
- add: consider drastically simplifying biomet QC (e.g. replace auto flags with
       simple Hampel filter)

2020-04-24
- add: [DONE] separate output types (i.e. data, plots, documentation) into 
       different folders where applicable?
- add: [DONE] correct wind direction only if not already done in eddypro
       1. import eddypro settings file in 02_correct_eddypro
       2. pull setting for geographic or magnetic north
       3. make wd correction dependent on this setting