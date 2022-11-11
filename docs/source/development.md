# Developer

The project has the following structure:

```
.SNP_haplotype
    ├── docs
    ├── htmlcov
    ├── output
    ├── snp_haplotyper
    ├── test_data
    └── tests
```

The source code for BASHer is contained in the snp_haplotyper directory:

```
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

## Running pytest

`BASHer` has two types of tests

```
SNP_haplotype
    └── tests
        ├── functional
        │   └── test_snp_filtering.py
        └── unit
            └── test_functions.py
```

```bash
# Run the tests use the following command within the local repo:
python -m pytest

# Calculate test coverage and produce HTML report (SNP_haplotype/htmlcov/index.html)
pytest --cov=snp_haplotyper
```

## Compiling Sphinx Documentation

```bash
# From within the local repo run the following commands.

# Auto generate documentation from docstrings
sphinx-apidoc -f -o docs/apidocs snp_haplotyper/

# Build documentation from source files
make clean -C docs # Optional - removes any caches html files
sphinx-build -b html docs/source /docs/build
```

## Editing Sphinx Documentation


```
SNP_haplotype
    └── docs
        ├── apidocs
        ├── build
        └── source
```

## Release Process

