# *tombRaider* - an algorithm to identify artefacts from metabarcoding data

## 1. Introduction

*tombRaider* is an algorithm capable of identifying and removing artefact sequences from metabarcoding data.

## 2. Installation

*tombRaider* is a command-line only toolkit running on typical Unix/Linux environments and is exclusively written in Python3. *tombRaider* can be installed manually from GitHub. Below are details for the manual installation process (section 2.1 Manual installation). Due to the native implementation of the incorporated functions, only a limited number of dependencies restricted to non-built-in Python3 modules are required for *tombRaider* to successfully execute. A complete list of the dependencies can be found below (section 2.2 Dependencies). Remember to install the dependencies separately, since the GitHub installation of *tombRaider* will not install the dependencies for you!

### 2.1 Manual installation

To manually install *tombRaider*, download the files from GitHub using the `git clone` command in the Terminal.

```{code-block} bash
git clone https://github.com/gjeunen/tombRaider.git
```

Dependent on your OS settings, the downloaded files might need to be made executable, which can be achieved by running the `chmod +x tombRaider` command. You can make *tombRaider* globally accessible in your OS using the following command:

```{code-block} bash
export PATH="/path/to/tombRaider/folder:$PATH"
```

Substitute `/path/to/tombRaider/folder` with the actual path to the repo folder on your system. Adding this line to the `.bash_profile`, `.bashrc`, or `.zshrc` file in your home directory will keep *tombRaider* globally accessible, even when closing the Terminal application.

### 2.2 Dependencies

Once the GitHub files are downloaded, ensure the following Python3 modules with the correct or compatible versions are accessible to *tombRaider*:

1. [Python3](https://www.python.org/downloads/) (*v 3.11.5*)
2. [rich](https://rich.readthedocs.io/en/stable/introduction.html) (*v 13.3.5*)
3. [rich-click](https://github.com/ewels/rich-click) (*v 1.6.1*)
4. [numpy](https://numpy.org) (*v 1.25.2*)
5. [pandas](https://pandas.pydata.org/docs/getting_started/install.html) (*v 2.1.4*)

### 2.3 Check the installation

To check if the installation was successful, type the following command into the Terminal to prompt the help information:

```{code-block} bash
tombRaider -h
```

![tombRaider help prompt](figures/help-prompt.png)

Additionally, you can run *tombRaider* on example files in the `exampleFiles` subdirectory by executing:

```{code-block} bash
tombRaider --example-run
```

## 3. Code execution

With *tombRaider* only a single line of code is required to identify and remove artefacts from your metabarcoding data!

To run *tombRaider* with default values, use the following line of code:

```{code-block} bash
tombRaider --method 'taxon-dependent co-occurrence' --frequency-input countTable.txt --sequence-input sequences.fasta --blast-input blastnResults.txt --frequency-output countTableNew.txt --sequence-output sequencesNew.fasta --blast-output blastnResultsNew.txt --condensed-log condensedLog.txt --detailed-log detailedLog.txt --count 0 --sort 'total read count'
```

Detailed information about all input file structures (section 4. Input files) and parameters (section 5. Parameters) can be found below.

## 4. Input files

An example of all input files can be found in the `exampleFiles` subdirectory on GitHub or where *tombRaider* is installed on your OS. This file list can help determine the file structures required by *tombRaider* when formatting errors are preventing a successful execution on your local files.

### 4.1 Count table

### 4.2 Sequence list

### 4.3 Taxonomy assignment

#### 4.3.1 BLAST

#### 4.3.2 BOLD

#### 4.3.3 SINTAX

#### 4.3.4 IDTAXA

## 5. Parameters

### 5.1 Main function

#### 5.1.1 --method

### 5.2 Parameters

#### 5.2.1 --occurrence-type

#### 5.2.2 --detection-threshold

#### 5.2.3 --similarity

#### 5.2.4 --negative

#### 5.2.5 --global-ratio

#### 5.2.6 --local-ratio

#### 5.2.7 --count

### 5.3 Input files

#### 5.3.1 --frequency-input

#### 5.3.2 --sequence-input

#### 5.3.3 --blast-input

#### 5.3.4 --bold-input

#### 5.3.5 --sintax-input

#### 5.3.6 --idtaxa-input

### 5.4 Output files

#### 5.4.1 --frequency-output

#### 5.4.2 --sequence-output

#### 5.4.3 --blast-output

#### 5.4.4 --bold-output

#### 5.4.5 --sintax-output

#### 5.4.6 --idtaxa-output

#### 5.4.7 --condensed-log

#### 5.4.8 --detailed-log

### 5.5 Frequency table details

#### 5.5.1 --taxa-are-rows

#### 5.5.2 --omit-rows

#### 5.5.3 --omit-columns

#### 5.5.4 --sort

### 5.6 Taxonomy file details

#### 5.6.1 --blast-format

#### 5.6.2 --use-accession-id

#### 5.6.3 --bold-format

#### 5.6.4 --sintax-threshold

### 5.7 Options

#### 5.7.1 --example-run

#### 5.7.2 --help
