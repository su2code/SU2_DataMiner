.. _configurations:

***************************
SU2 DataMiner Configuration
***************************
*SU2 DataMiner* uses a configuration class in order to store important information regarding the data generation, data mining, and 
manifold generation processes. This page lists some of the important functions of the *Config* class which acts as the :ref:`base<base>` 
class for configurations specific to the application such as :ref:`NICFD<NICFD>` and FGM :ref:`FGM<FGM>`, for which additional settings can be specified.

.. contents:: :depth: 2

.. _base:

Base Config Class 
=================

*SU2 DataMiner* uses a configuration class in order to store important information regarding the data generation, data mining, and 
manifold generation processes. This page lists some of the important functions of the *Config* class which acts as the base 
class for configurations specific to the application such as :ref:`NICFD<NICFD>` and FGM :ref:`FGM<FGM>`, for which additional settings can be specified.


Storage Location and Configuration Information
----------------------------------------------

During the various processes in *SU2 DataMiner*, data are generated, processed, and analyzed. All information regarding these
processes is stored on the current hardware in a user-defined location. *SU2 DataMiner* configurations can be saved locally under 
different names in order to keep track of various data sets and manifolds at once. 
The following functions can be used to manipulate and access the storage location for fluid data and manifolds of the *SU2 DataMiner* configuration 
and save and load configurations.

.. autofunction:: Common.Config_base.Config.SetConfigName 

.. autofunction:: Common.Config_base.Config.GetConfigName 

.. autofunction:: Common.Config_base.Config.SaveConfig 


.. code-block::

       from su2dataminer.config import Config 

       c = Config()
       c.SetConfigName("test")
       c.SaveConfig()

       
.. autofunction:: Common.Config_base.Config.SetOutputDir 

.. autofunction:: Common.Config_base.Config.GetOutputDir 

.. autofunction:: Common.Config_base.Config.PrintBanner
    
.. autofunction:: Common.Config_base.Config.SetConcatenationFileHeader 

.. autofunction:: Common.Config_base.Config.GetConcatenationFileHeader 


Settings for Machine Learning Applications 
------------------------------------------

The data-driven fluid modeling applications of *SU2 DataMiner* involve the use of multi-layer perceptrons to calculate the thermodynamic state of the fluid during flow calculations in SU2. 
The values of the weights and biases, the hidden layer architecture(s) and other parameters needed to train the network can be stored in and retrieved from the *SU2 DataMiner* configuration. 

.. autofunction:: Common.Config_base.Config.SetControllingVariables

.. autofunction:: Common.Config_base.Config.GetControllingVariables

.. autofunction:: Common.Config_base.Config.SetTrainFraction

.. autofunction:: Common.Config_base.Config.GetTrainFraction

.. autofunction:: Common.Config_base.Config.SetTestFraction

.. autofunction:: Common.Config_base.Config.GetTestFraction   

The following functions are used to specify the parameters used for training the networks on the fluid data. Currently, *SU2 DataMiner* uses supervised, gradient-based methods for training where the learning 
rate follows the exponential decay method. The value of the initial learning rate and decay parameter can be accessed through 

.. autofunction:: Common.Config_base.Config.SetAlphaExpo

.. autofunction:: Common.Config_base.Config.GetAlphaExpo

.. autofunction:: Common.Config_base.Config.SetLRDecay

.. autofunction:: Common.Config_base.Config.GetLRDecay
   
The networks are trained using batches of training data. The size of the training batches and the maximum number of epochs for which the networks are trained can be accessed through

.. autofunction:: Common.Config_base.Config.SetBatchExpo 

.. autofunction:: Common.Config_base.Config.GetBatchExpo 

.. autofunction:: Common.Config_base.Config.SetNEpochs

.. autofunction:: Common.Config_base.Config.GetNEpochs

Additionally, the number of nodes in the hidden layers of the networks and activation function can be accessed. *SU2 DataMiner* currently supports the use of a single activation function to all hidden layers which can be selected from a list of options and the input and output layers use a linear function.

.. autofunction:: Common.Config_base.Config.SetHiddenLayerArchitecture

.. autofunction:: Common.Config_base.Config.GetHiddenLayerArchitecture

.. autofunction:: Common.Config_base.Config.SetActivationFunction

.. autofunction:: Common.Config_base.Config.GetActivationFunction

Finally, the weights and biases data can be accessed in the configuration using 

.. autofunction:: Common.Config_base.Config.SetWeights 

.. autofunction:: Common.Config_base.Config.SetBiases

.. autofunction:: Common.Config_base.Config.GetWeightsBiases

All the settings regarding the training parameters, architecture, weights, and biases of the trained network can be stored automatically after training a network from the :ref:`Trainer` class 

.. autofunction:: Common.Config_base.Config.UpdateMLPHyperParams


To use the network stored in the configuration in SU2 simulations, the network information needs to be written to an ASCII file such that it can be loaded in SU2 through the `MLPCpp`_ module. 
All relevant information about the network is automatically written to a properly formatted ASCII file using 

.. autofunction:: Common.Config_base.Config.WriteSU2MLP

.. _NICFD:

Configuration for real-gas applications 
=======================================

The following section describes the most important functions of the Config_NICFD class.

.. autofunction:: Common.DataDrivenConfig.Config_NICFD.__init__

.. autofunction:: Common.DataDrivenConfig.Config_NICFD.SetFluid

.. autofunction:: Common.DataDrivenConfig.Config_NICFD.SetEquationOfState 

Example
-------
The code snippet below demonstrates the 
:: 

    from su2dataminer.config import Config_NICFD 

    config = Config_NICFD()
    config.SetFluid("MM")
    config.SetEquationOfState("HEOS")
    config.SetName("siloxane_MM_heos")
    config.UseAutoRange(True)
    config.PrintBanner()


.. _FGM:

Configuration for combustion applications 
=========================================

The following section describes the most important functions of the SU2 DataMiner configuration class for FGM applications. 

SU2 DataMiner supports the generation of premixed flamelets using `Cantera <CanteraLink_>`_ for FGM applications. 


Reaction Mechanism and Species Transport Model 
----------------------------------------------


The reaction mechanism used for flamelet calculations is specified through the following function:


.. autofunction:: Common.DataDrivenConfig.Config_FGM.SetReactionMechanism

.. autofunction:: Common.DataDrivenConfig.Config_FGM.GetReactionMechanism

.. autofunction:: Common.DataDrivenConfig.Config_FGM.SetTransportModel

.. autofunction:: Common.DataDrivenConfig.Config_FGM.GetTransportModel


Reactants and Flamelet Types
----------------------------


The fuel and oxidizer used for flamelet calculations can be defined through the following functions: 


.. autofunction:: Common.DataDrivenConfig.Config_FGM.SetFuelDefinition

.. autofunction:: Common.DataDrivenConfig.Config_FGM.SetOxidizerDefinition


The types of flamelet data that are currently supported by SU2 DataMiner are adiabatic flamelets, burner-stabilized flamelets, and chemical equilibrium data. 
Each of these flamelet types can be included or excluded from the manifold through the following functions: 

.. autofunction:: Common.DataDrivenConfig.Config_FGM.RunFreeFlames 

.. autofunction:: Common.DataDrivenConfig.Config_FGM.RunBurnerFlames

.. autofunction:: Common.DataDrivenConfig.Config_FGM.RunEquilibrium


FGM Controlling Variables 
-------------------------

.. autofunction:: Common.DataDrivenConfig.Config_FGM.SetControllingVariables

.. autofunction:: Common.DataDrivenConfig.Config_FGM.GetControllingVariables

.. autofunction:: Common.DataDrivenConfig.Config_FGM.SetProgressVariableDefinition 

.. autofunction:: Common.DataDrivenConfig.Config_FGM.ResetProgressVariableDefinition 

.. autofunction:: Common.DataDrivenConfig.Config_FGM.GetUnburntScalars 


Preferential Diffusion 
----------------------

.. autofunction:: Common.DataDrivenConfig.Config_FGM.EnablePreferentialDiffusion 

.. autofunction:: Common.DataDrivenConfig.Config_FGM.SetAverageLewisNumbers 

.. autofunction:: Common.DataDrivenConfig.Config_FGM.ComputeBetaTerms 


Thermochemical State Variables 
------------------------------

.. autofunction:: Common.DataDrivenConfig.Config_FGM.AddOutputGroup

.. autofunction:: Common.DataDrivenConfig.Config_FGM.RemoveOutputGroup 

    
.. autofunction:: Common.DataDrivenConfig.Config_FGM.EditOutputGroup



.. _MLPCpp : https://github.com/EvertBunschoten/MLPCpp.git 

.. _CanteraLink: https://cantera.org/