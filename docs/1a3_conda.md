---
layout: default
title: Conda environment
parent: 1. Introduction
nav_order: 4
description: A comprehensive guide to understanding epigenetics.
published: true
---


FINAL VERSION
{: .label .label-green }

# Conda

From the website, `conda` provides ["Package, dependency and environment management for any language"](https://docs.conda.io/en/latest/).

Conda is a package manager allows specific versions of programs to be installed, alongside their dependencies.
Different sets of programs can be installed to different [virtual environments](https://www.anaconda.com/moving-conda-environments/).
A virtual environment is basically a set of programs.

## Installing `conda`

Conda is part of [Anaconda](https://www.anaconda.com/distribution/), which is available for free.
Conda is also available through [Miniconda](https://docs.conda.io/en/latest/miniconda.html), a free minimal installer for conda.

Conda can be installed on a 64-bit Linux system with the following commands...

```bash
# Downloading miniconda
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# Installing miniconda
$ bash miniconda.sh -b -p $HOME/miniconda
# Updating conda
$ conda update -q conda
```

## Cloning and activating a `conda` environment

Conda virtual environments can be shared, either as a `.yml` file or a `.txt` file.
A `.yml` copy of a conda environment can be used to recreate that environment on another machine, regardless of the operating system platform used.
A `.txt` copy of a conda environment is more explicit: it can be used to create an identical copy of a conda environment using the same operating system platform as the original machine.


A conda environment can be activated using `$ conda activate name_of_environment`.
Once activated, the programs installed in this environment are available.
Conda can be deactivated using `$ conda deactivate`.



## Further reading
1. Downloading conda: <https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html>
2. Conda packages: <https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/packages.html>
3. Conda environments: <https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>


## Configure the prompt 

In order to configure the prompt and be able to use conda, you need to run the following command

```bash
eval "$(/course/tadb/miniconda3/bin/conda shell.bash hook)"
```

Now you will be able to use teh conda command.

In order to use the conda environments, you need to run the following command:

```bash
conda init bash
```

![alt text](conda-init.png)