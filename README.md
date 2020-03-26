# Concat     

Concat is a command-line tool written in python3 to automatically concatenate alignments and calculate the corresponding phylogeny.

As input it needs  single locus alignments and files with the information about the sequence names and the corresponding taxon names. 
If a phylogeny is provided, this can either be used as starting tree or as a backbone constraint.
There is also an option available to run it in combination with the package [PhylUp](https://github.com/blubbundbla/PhylUp_remote.git).


### Before you can start

#### 1. install the dependencies:

* [RAxML-NG](https://github.com/amkozlov/raxml-ng/archive/master.zip) - tree estimation program
* [ModelTest-NG](https://github.com/ddarriba/modeltest/archive/master.zip) - to determine the substitution model


make sure the programs are accessible from everywhere, thus add them to your PATH using the command line:
* UNIX: `export PATH=$PATH:/path/to/my/program`
* Windows: `set PATH=%PATH%;C:\path\to\my\program`
* MAC: `export PATH=$PATH:~/path/to/program`

(! set PATH=%PATH%:  it takes the current path and sets PATH to it.)

#### 2. download Concat using the command line:
* as a normal package: `wget https://github.com/blubbundbla/Concat.git`
* as a git repository: `git clone 'git@github.com:blubbundbla/Concat.git'`


#### 3. install python requirements and dependencies:

run from within the Concat main folder:

* `python setup.py install`
* `pip install -r requirements.txt`

   
### Set up a run

There are several example files in the main folder. Edit them for your purpose. 
The file called `example_concat_phylup.py` shows how to use the package in combination with the 
[PhylUp](https://github.com/blubbundbla/PhylUp_remote.git) package.
   
The program randomly decides which sequences of the same otu to concatenate if there are more sequences available for one loci.
 In the future, the user can also specify a file, which sequences shall be concatenated.

