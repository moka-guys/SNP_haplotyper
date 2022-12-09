# Developer

The project has the following structure:

```bash
.SNP_haplotype
    ├── docs
    ├── htmlcov
    ├── output
    ├── snp_haplotyper
    ├── test_data
    └── tests
```

The source code for BASHer is contained in the snp_haplotyper directory:

```bash
SNP_haplotype
    └── snp_haplotyper
            ├── autosomal_dominant_logic.py
            ├── autosomal_recessive_logic.py
            ├── debug_data.ipynb
            ├── excel_parser.py
            ├── exceptions.py
            ├── json2excel.py
            ├── sample_sheet_reader.py
            ├── snp_haplotype.py
            ├── snp_plot.py
            ├── stream_output.py
            ├── templates
            │   └── report_template.html
            └── x_linked_logic.py
```

**snp_haplotype.py** contains the majority of the code for processing samples

**autosomal_dominant_logic.py, autosomal_recessive_logic.py, x_linked_logic.py** contain the logic for assigning informative SNPs for each mode of inheritance.

**sample_sheet_reader.py** parses the command line arguments from the supplied excel sheet template containing the sample's metadata.

**report_template.html** is the HTML template for the report.

**stream_output.py** streams the relevant output as a JSON when in testing mode so that the output is easily readable by the test suite.

## Automated Testing Using Pytest

BASHer has two types of tests

```bash
SNP_haplotype
    └── tests
        ├── functional
        │   └── test_snp_filtering.py
        └── unit
            └── test_functions.py
```

### Unit Tests

Unit tests are used to test the functionality of individual functions within the script to ensure they return expected results.

### Functional Tests

Functional tests run test data through the software and compare the output to benchmark results manually curated by the PGD team.
This benchmark covers all possible scenarios (reference type, mode of inheritance, embryo sex, and consanguinity).  These were converted to machine readable JSON format and integrated into the test suite.  Before any push to production these tests should be run to ensure that results match the expected output.

### Coverage

The amount of code covered by the testing is calculated using coverage.py, the results of which can be found as an HTML report in:

```bash
SNP_haplotype
    └── htmlcov
        └── index.html
```

### Run Tests and Coverage

```bash
# Run the tests use the following command within the local repo:
python -m pytest

# Calculate test coverage and produce HTML report (SNP_haplotype/htmlcov/index.html)
pytest --cov=snp_haplotyper
```

### Debugging Tests

There is a known issue where the VS Code debugger doesn't stop at breakpoints set in pytest test modules if pytest-cov is used to calculate coverage.  This issue also affects PyCharm.  When debugging tests you can manually set the "--no-cov" flag in the VS Code's settings.json as shown below.

```json
{
    "python.testing.pytestArgs": [
        ".",
        "--no-cov",
    ]
}
```

You can then debug the tests using breakpoints.  Remember to change the settings.json back again once the debugging is completed.

## Sphinx Documentation

The documentation for Basher is generated using Sphinx.
The docstrings in the python modules are parsed using sphinx-apidoc and are automatically added to the documentation.

### Editing Sphinx Documentation

The source files for the docs are located in the source folder as shown below.  The majority of these source files are in markdown, with some being in rst or HTML.  These files can be edited in a text editor and compiled following the instructions in the next section.  This will build the corresponding HTML files from the source files and store the output in the build folder.

```bash
SNP_haplotype
    └── docs
        ├── apidocs
        ├── build
        └── source
```

### Compiling Sphinx Documentation

```bash
# From within the local repo run the following commands.

# Auto generate documentation from docstrings
sphinx-apidoc -f -o docs/apidocs snp_haplotyper/

# Build documentation from source files
(cd docs && make clean) # Optional - removes any cached html files
sphinx-build -b html docs/source docs/build
```

To view the output open index.html in the build folder.

## Release Process

TODO