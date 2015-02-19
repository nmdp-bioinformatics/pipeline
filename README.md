pipeline
========
This is a prototypical pipeline.  Greater detail on how this works, why it's supposed to work, and much, much more can be found in the wiki at the github page.

[![Build Status](https://travis-ci.org/nmdp-bioinformatics/pipeline.svg?branch=master)](https://travis-ci.org/nmdp-bioinformatics/pipeline)


DESCRIPTION
========
splitter is the master script which figures out what files to work on, and divides the workload up -- into one file per core on the machine.  It then starts process_fastq on each of these files.  Note that it isn't particularly clever about dividing up the workloads.  If there are 100 files, and 4 cores, you get 4 files of 25 lines each -- regardless of sample size.  Still, some parallelization is better than none, right?

USING pipeline
========

 - clone the repository to your working directory
 - edit splitter.bash to include the correct path to your copy of the pipeline
    note -- splitter.bash will not run unless you set this value.
 - ensure your data is in the format that splitter and process_fastq expect:
     - all data files to be processed must be in fastq format, with a filename ending in .fastq or .fq
     - your directory of data must have all fastq files in a subdirectory named raw/ .  Optionally, you can organize your fastq files in subdirectories of raw/, but this may add complexity.
     - avoid having underscores in subdirectory names under raw/
     - any non-fastq files in the raw/ subdirectory will be ignored
     - you must have two writable directories under the same parent as final, named intermediate/ and final/
     - that is, if your data is in /mnt/foo/mydata, there must be /mnt/foo/mydata/raw, /mnt/foo/mydata/intermediate and /mnt/foo/mydata/final

 - run splitter, thusly:
```bash
/path/to/splitter.bash /path/to/my/compliant/directory/of/data  [ DEBUG ]
```
    - note that DEBUG is an optional flag and will GREATLY increase the verbosity of the tools _any_ string value will activate the debugging.

  - check the result files in final/
  - intermediate files are stored in intermediate/ and can safely be deleted after processing has completed

CAVEATS
========
This tool requires a fair number of other tools to be installed.  The shell script checks for them, and will try to give you breadcrumbs as to how to solve any of the issues it finds (e.g. where to find bwa).

If your machine/instance has less than 4GB of RAM, you'll want to modify the script so it uses chr6.fa as the reference genomic, install of all_chr.fa.   Make sure those files are indexed.  This script does not do the indexing for you.

This is a work in progress, but we welcome improvements and questions.  Please feel free to file issues in github.

This has been primarily tested on Debian-based Linux systems (Ubuntu 14.04.01 LTS, to be exact).  However, we do try to support other Unix-like operating systems, such as Apple OSX.  If you run into issues on your distribution of choice, please file an issue on github.  OSX, mind you, has some VERY old shell utilities, so it's still a bit prototypical.

DEBUGGING
========
Elapsed times, STDOUT and STDERR messages are stored in files named timedata.[0-9]+ -- that is, timedata, followed by a numeric string derived from the PID of the splitter.bash script. 

All timedata files may safely be removed after  run validation.


ngs-validation-report
========
This tool takes in the output of the NGS pipeline and produces a html report of the results. This allows for an easily traversible report for viewing the validation results. 

Using ngs-validation-report
========
1) Running with no parameters
```bash
/path/to/tool/ngs-validation-report
```

2) Running with only an experiment name
```bash
/path/to/ngs-validation-report -x ex014
```

3) Running with only an experiment name and an input directory
```bash
/path/to/ngs-validation-report -x ex014 -d /path/to/directory/of/data
```

4) Running with the experiment, validated, and observed files.
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
│   ├── subject2.html
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