# Deployment

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
allow_cosanguineous_cases = True
```
