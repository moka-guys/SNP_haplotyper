# Deployment

## Get BASHer from GitHub Repo

BASHer can be deployed in two separate ways. Check what release version of BASHer you are downloading to ensure that it is the one you want.  Previously released versions are available from the GitHub repo if required.

* If you have `git` installed then, follow these steps:

1. Navigate to the directory you want to install BASHer in.
2. Use the `git clone` command to clone the repository from GitHub.

```bash
clone git@github.com:moka-guys/SNP_haplotyper.git
```

* To deploy BASHer using the second method, you will need to have a tool installed on your system that can unpack ZIP files. Then, follow these steps:

1. Go to the GitHub page for the [BASHer project](https://github.com/moka-guys/SNP_haplotyper) and download the ZIP file for the repository by clicking on the green `Code` dropdown box and selcting the zip file.

2. Unpack the ZIP file in the desired directory using your preferred tool. This will create a new directory containing the contents of the repository.

## Set Configuration via Config.py

There is a config.py file in the main directory with several settings which need configuring before deployment:

Ensure the correct genome build is entered - this will be recorded on the HTML reports produced.

```python
# The genome build used by the SNP array
genome_build = "GRCh38"
```

When debugging the output can be streamed to stdout in a JSON format.  Before deploying to production ensure that stream_results = False, otherwise HTML reports will not be produced.

```python
# If this flag is set to True, the program will stream the results of the analysis in JSON
# format to the standard output intead of creating a HTML report. This is useful for
# debugging issues with the test suite as pytest captures this output during testing and
# compares it to the expected output

stream_results = False
```

Several flags can be set to prevent features which have not been properly validated from being used in production.  Set the flags accordingly before a release.

```python
# These flags can be used to prevent the program from running certain parts of the analysis
# for example, if you have validated specific modes of inheritance, you can set the flags
# to skip the analysis of other modes of inheritance
allow_autosomal_dominant_cases = True
allow_autosomal_recessive_cases = True
allow_x_linked_cases = True
allow_consanguineous_cases = False
```

## Configure helper scripts, if used

### Background

There are three helper scripts which can be implemented to help uses run BASHer.  They require the user to fill in an excel template with the analysis details.

* The helper script `excel_parser.py` can parse this excel template and generate the arguments to pass to `snp_haplotyper.py`.  This script can be run from the command line, by passing the excel file as an argument, or via the helper scripts below which both open a file browser and allow the excel spreadsheet to be selected and passed to `excel_parser.py`:
    1. `launch_SNP_analysis.ps1`, a powershell script which opens a file browser. This can be used with the standalone distribution of python which does not support the TKinter module.
    2. `select_input_file.py`, a python script which opens a file browser.  This requires that TKinter module is installed.

### Configuring Helper Scripts

In `launch_SNP_analysis.ps1` there are two lines of code that will need manually updating:

* Set the file path in the line below to where you want the filebrowser to initially open:

    ```powershell
    $FileBrowser.InitialDirectory = "C:\filepath\to\test_data"
    ```

* Set the file paths in the line below, so that:

   1. The first filepath points to the `python.exe` for the version of python3 you want to use.
   2. The second filepath points to the location of `excel_parser.py`
  
```powershell
    # Call excel_parser.py with the provided excel template
    & C:\filepath\to\Software\python-3.10.0-embed-amd64\python.exe C:filepath\to\snp_haplotyper\excel_parser.py --input_file $FilePath
```

Note that the command has to start  with `&` so that powershell calls the following command.

Configuration of `select_input_file.py` is not currently supported.

### Test Deployment


