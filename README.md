# Concat     

Concat is a command-line tool written in python3 to automatically concatenate alignments and calculate the 
corresponding phylogeny. Process can be split into two main steps:
1. The program concatenates sequences from the same ncbi taxon id. If more sequences of the same taxon id for the same locus are available those can be either randomly selected, or the user can define it via file. For details see below. 
2. If enabled a concatenated phylogeny is being calculated. If a phylogeny is provided, this can either be used as starting tree or as a backbone constraint using RAxML-NG.

As input the program needs multiple single locus alignments and files with the information about the sequence names and 
the corresponding taxon names. The file has to be provided as a comma-delimited file with sequence tip names 
and a taxon name matching the ncbi taxonomic database. 



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

* as a normal package: `wget https://github.com/mkandziora/Concat.git`
* as a git repository: `git clone 'git@github.com:mkandziora/Concat.git'`


#### 3. install python requirements and dependencies:

run from within the Concat main folder:

* `python setup.py install`
* `pip install -r requirements.txt`

   
### Set up a run

There are several example files in the main folder.

To concatenate sequences from multiple loci, there are three options evailable.
1. The program randomly decides which sequences of the same otu/ncbi taxon id to concatenate if there are 
more sequences available for one locus. 
2. The program first suggests which taxa could be concatenated and in a second steps concatenated them.
 During the first step (see `example_concat_table.py`), suggested taxa to concatenate are written to the file 'concattable_selfselect.txt'. The association by `concat_id` can then be modified by the user. In the second step (see `example_concat_calc.py`), the concatenated alignment file is being written.
3. The user specifies which sequences shall be concatenated in `concattable_selfselect.txt` and needs to enable the `suggest` and `select` option. For examples please refer to `example_concat_table.py`. 

There is also a method available to run it in combination with the package [PhylUp](https://github.com/mkandziora/PhylUp.git): `example_concat_phylup.py`.

