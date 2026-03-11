Requirements and Set-Up
=======================

*SU2 DataMiner* is a python-based software which was originally developed for python 3.11 and currently only supports Linux distributions. 
The software can be used through terminal commands for easy use, but can be controlled more effectively through the python API. 

Storage Location and Configuration Information
----------------------------------------------

During the various processes in *SU2 DataMiner*, data are generated, processed, and analyzed. All information regarding these
processes is stored on the current hardware in a user-defined location. *SU2 DataMiner* configurations can be saved locally under 
different names in order to keep track of various data sets and manifolds at once. 
The following functions can be used to manipulate and access the storage location for fluid data and manifolds of the *SU2 DataMiner* configuration 
and save and load configurations.


Requirements
------------
The following python modules are required by the various modules in *SU2 DataMiner*:

* argparse
* cantera
* CoolProp
* gmsh
* joblib
* keras
* matplotlib
* numpy
* pickle
* pyfiglet
* random
* scipy
* scikit-learn
* tensorflow
* time
* tkinter
* tqdm

A suitable conda environment can also be created using the [environment recipe](environment.yml) through the following command:

:: 

    conda env create -f environment.yml

Set-Up 
------

After cloning `the repository`_, and checking out the desired branch, add the following lines to your ```~/.bashrc``` to update your pythonpath and access the *SU2 DataMiner* modules through your terminal and python API:

::

    export PINNTRAINING_HOME=<PATH_TO_SOURCE>
    export PYTHONPATH=$PYTHONPATH:$PINNTRAINING_HOME
    export PATH=$PATH:$PINNTRAINING_HOME/bin


where ```<PATH_TO_SOURCE>``` is the path to where you cloned the repository.

Check whether your setup was successful running the following command in your terminal:

::

    which GenerateConfig.py 

which should return ```<PATH_TO_SOURCE>/bin/GenerateConfig.py```.

Tutorials can be found under ```TestCases```, proper documentation will follow soon.

.. _the repository : https://github.com/EvertBunschoten/SU2_DataMiner.git