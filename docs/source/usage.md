# BASHer Usage

There are two ways to run BASHer:

1. Directly via the command line.

2. Running the helper script, excel_parser.py, which parses a user supplied Excel sheet containing the relevant metadata (sample names, inheritance mode, genomic co-ordinates of ROI etc) and will auto generate & run the relevant commands.

## Command line usage

```bash
# Example command for autosomal dominant case
python3 snp_haplotype.py --input_file F4_BRCA2_AD.txt --output_folder output/ --output_prefix F4_BRCA2_AD --mode_of_inheritance autosomal_dominant --male_partner 22.F4.MP.rhchp --male_partner_status affected --female_partner 23.F4.FP.rhchp --female_partner_status unaffected --reference 19.F4.PGF.rhchp --reference_status affected --reference_relationship grandparent --embryo_ids 24.F4.EMB11.rhchp 25.F4.EMB12.rhchp 26.F4.EMB13.rhchp --embryo_sex unknown unknown unknown --gene_symbol BRCA2 --gene_start 32315086 --gene_end 32400268 --chr 13

# Example command for autosomal recessive case
python3 snp_haplotype.py --input_file F3_HBB_AR.txt --output_folder output/ --output_prefix F3_HBB_AR --mode_of_inheritance autosomal_recessive --male_partner 13.F3.MP.rhchp --male_partner_status carrier --female_partner 12.F3.FP.rhchp --female_partner_status carrier --reference 15.F3.EMB7.rhchp --reference_status affected --reference_relationship child --embryo_ids 16.F3.EMB8.rhchp 17.F3.EMB9.rhchp 18.F3.EMB10.rhchp --embryo_sex unknown unknown unknown --gene_symbol HBB --gene_start 5225464 --gene_end 5227197 --chr 11

# Example command for X-linked case
python3 snp_haplotype.py --input_file F10_FMR1_XL.txt --output_folder output/ --output_prefix F10_FMR1_XL --mode_of_inheritance x_linked --male_partner 61.F10.MP.rhchp --male_partner_status unaffected --female_partner 62.F10.FP.rhchp --female_partner_status carrier --reference 63.F10.CCM.rhchp --reference_status affected --reference_relationship child --embryo_ids 64.F10.EMB33.rhchp 65.F10.EMB34.rhchp 66.F10.EMB35.rhchp --embryo_sex male female female --gene_symbol FMR1 --gene_start 147911919 --gene_end 147951125 --chr x
```

## Running Basher via excel_parser.py helper script

### Running excel_parser.py from the command line

Rather than manually generating the command with all the correct parameters users can fill in an excel template with the analysis details and use the `excel_parser.py` script to parse the details and automatically generate the command to run BASHer.

`excel_parser.py` can either be run from the command line:

```powershell
python3 "\path\to\excel_parser.py" --input_file "\path\to\excel_template.xlms
```

### Running excel_parser.py via a GUI

The simplest way to launch BASHer is for users to right-click on the powershell script `launch_SNP_analysis.ps1`, select "Open with powershell", and use the file browser which opens to select the relevant excel template which will then be passed to `excel_parser.py` causing BASHer to be run with the analysis details in the supplied spreadsheet.
