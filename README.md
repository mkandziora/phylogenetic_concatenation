# Concat     

Concat is a command-line tool written in python3 to automatically concatenate alignments and update phylogenies.

As input it needs  single locus alignments and files with the information about the sequence names and the corresponding taxon names. 
If a phylogeny is provided, this can either be used as starting tree or as a backbone constraint.

## Short tutorial:

### Before you can start

#### 1. install the dependencies:

* [RAxML-NG](https://github.com/amkozlov/raxml-ng/archive/master.zip) - tree estimation program
* [gappa](https://github.com/lczech/gappa/archive/master.zip)  - transforms the output from EPA-NG into a readable output
* [ModelTest-NG](https://github.com/ddarriba/modeltest/archive/master.zip) - to determine the substitution model


make sure the programs are accessible from everywhere, thus add them to your PATH using the command line:
* UNIX: `export PATH=$PATH:/path/to/my/program`
* Windows: `set PATH=%PATH%;C:\path\to\my\program`
* MAC: `export PATH=$PATH:~/path/to/program`

(! set PATH=%PATH%:  it takes the current path and sets PATH to it.)

#### 2.a) download Concat using the command line:
* as a normal package: `git clone https://github.com/blubbundbla/Concat.git`
* as a git repository: `git clone 'git@github.com:blubbundbla/Concat.git'`

#### 2.b) install a virtual environment
  This is very useful if you want to run it on a cluster and/or do not want to change already installed python packages on your computer.
  A virtual environment will locally install the packages needed.

  `pip install virtualenv` 
  `virtualenv -p python3 NameOfYourENV`  # you may need to just say `python` instead of `python3`, depending on your system

  To use the virtual machine you need to activate it before doing anything else. 
  This needs to be done before you start installing software in your virtual maschine or before running PhylUp.

  `source NameOfYourENV/bin/activate`

  and to deactivate it: `deactivate`

#### 3. install python requirements and dependencies:

run from within the Concat main folder:

* `python setup.py install`
* `pip install -r requirements.txt`

#### 4. install a local instance of the BLAST taxonomy database: 

Should be done automatically. If not follow description below.
   
   *  `wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'`
   *  `gunzip  -cd taxdump.tar.gz | (tar xvf - names.dmp nodes.dmp)`  
   *  move files into `tests/data/`
   *  `rm taxdump.tar.gz`

   
### Set up a run

There are several example files in the main folder. Edit them for your purpose. 
The file called `example_concat_phylup.py` shows how to use the package in combination with the PhylUp package.
   
The program randomly decides which sequences of the same otu to concatenate if there are more sequences available for one loci.
 In the future, the user can also specify a file, which sequences shall be concatenated.

