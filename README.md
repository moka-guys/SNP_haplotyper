# BASHer (BiAllellic SNP Haplotyper)

This tool can be used to evaluate the risk haplotypes associated with genetic conditions and assess risk embryo. It is applicable for different MOI (AD, AR, AR consanguineous and XL) 

## Inputs
- SNP array file (Recommended to generate this by using https://github.com/moka-guys/pre_basher_filter )
- Sample sheet

## Outputs
- html 
- pdf 

## Background logic of BASHer
https://biallelic-snp-haplotyper-basher.readthedocs.io/en/latest/internal_logic.html

## How to run Basher

### Run web-app locally with docker-compose

Clone the github repo and run docker-compose file 
```
docker-compose up -d
```
This should load the web app in `http://127.0.0.1:5000/basher`. 

### Run web-app on GSTT server

Go to Genapp server and load Basher web app

