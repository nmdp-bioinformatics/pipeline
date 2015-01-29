ngs-validation-report
========
This tool takes in the output of the NGS pipeline and produces a 

[![Build Status](https://travis-ci.org/nmdp-bioinformatics/pipeline.svg?branch=master)](https://travis-ci.org/nmdp-bioinformatics/pipeline)


Using ngs-validation-report
========
1) Running with no parameters
```bash
/path/to/tool/ngs-validation-report
```
If you pass it no parameters then it looks in the current directory for the expected, validated, and observed file. The output directory is defaulted to your current directory and the report title is defaulted to report.

2) Using only an experiment name
```bash
/path/to/ngs-validation-report -x ex014
```
If you pass the tool only an experiment name, then it looks in the current directory for an expected, validated, and observed file with that experiment name.

3) Using an experiment name and an input directory
```bash
/path/to/ngs-validation-report -x ex014 -d /path/to/directory/of/data
```
4) Providing the experiment, validated, and observed files.
```bash
/path/to/ngs-validation-report -l ex00_validated.txt -e ex00_expected.hml -o ex00_observed.txt
```

Output
========

```
report/
├── index.html
├── results.html
├── errors.html
├── experiments.html
├── help.html
├── log.html
├── subjects/
│   ├── subject1.html
│   └── ...
├── css/
│   ├── dashboard.css
│   ├── bootstrap.min.css
│   └── default.css
├── js/
│   ├── Chart.js
│   ├── docs.min.js
│   ├── ie-emulation-modes-warning.js
│   ├── ie10-viewport-bug-workaround.js
│   ├── ie8-responsive-file-warning.js
│   ├── init.js
│   ├── jquery.js
│   ├── raphael.js
│   └── bootstrap.min.js
└── img/
    └── bethematch.jpeg

```

Caveats
========

The email option will only work if you're on a linux machine with mailx properly installed. This option DOES NOT work with nmdp email addresses from the private cloud. In order to do the allele code expansion you must have LWP::UserAgent installed on your machine. 

Debugging
========
Elapsed times and STDERR messages are stored in the log.html file if your run with the verbose flag.


Copyright and License
=====================
Code released under the MIT license.
Copyright (c) 2014 National Marrow Donor Program (NMDP)


