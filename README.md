# *tombRaider* - an algorithm to identify artefacts from metabarcoding data

## 1. Introduction

*tombRaider* is an algorithm capable of identifying and removing artefact sequences from metabarcoding data.

## 2. Installation

*tombRaider* is a command-line only toolkit running on typical Unix/Linux environments and is exclusively written in Python3. Due to the native implementation of the incorporated functions, only a limited number of dependencies restricted to non-built-in Python3 modules are required for successful code execution. A complete list of the dependencies can be found below. *tombRaider* can be installed manually from GitHub. Remember to install the dependencies separately! Below are details for the manual installation process and a list of current dependencies with compatible versions.

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

## 4. Input files

An example of all input files can be found in the `exampleFiles` subdirectory on GitHub or where *tombRaider* is installed on your OS. This file list can help determine the file structures required by *tombRaider* when formatting errors are preventing a successful execution of *tombRaider* on your local files.
