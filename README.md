# Code documentation

## Project Structure

Below is a template for documenting the structure of your codebase. Each entry should include the folder/file name and a brief description of its contents or purpose.

```
kd_looms
├── code/         
│   └── archives.R  
│   └── functions.R
│   └── init.R
│   └── weight_db.R
├── data/              
│   ├── kd_level1.nex   
│   ├── kd_level2.nex
│   ├── kd_level3.nex
│   ├── kd_level4.nex
│   └── Kra-DaiLooms28master-2.xlsx
└── outputs    
```

## File Descriptions

### Inputs

1. `data/`: This directory contains the five database we use. You can reproduce all the  `kd_level`'s files with  `kd_by_level` function.