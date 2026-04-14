.. _NICFD_LUT: 

.. sectionauthor:: Evert Bunschoten

Table Generation for NICFD Applications 
=======================================

SU2 DataMiner supports the creation of look-up table methods for thermophyscial state evaluations in NICFD simulations in SU2. 
This tutorial demonstrates how to generate thermophyisical tables for NICFD applications in SU2. 

To get started, you will need to have installed SU2 DataMiner according to the :ref:`installation instructions <label_setup>`. 


.. contents:: :depth: 2


1. Config Generation
--------------------

As for any process within the SU2 DataMiner workflow, all settings regarding the setup of the fluid data generation and tabulation are stored in an SU2 DataMiner :ref:`configuration object <NICFD>`.
The tutorial for setting up a generic SU2 DataMiner configuration can be found :ref:`here <tutorialconfigs>`. 

In this example, a look-up table will be created for the application of modeling fluid properties of carbondioxide in supercritical conditions.
The following Python code snippet shows the initial set-up. 


.. code-block::

    #!/usr/bin/env python3
    from su2dataminer.config import Config_NICFD 

    config = Config_NICFD()
    config.SetFluid("CarbonDioxide")
    config.SetEquationOfState("HEOS")

    




Step 2:
-------



.. image:: /Tutorials/SU2DataMiner_logo.png
   :height: 200 px
   :width: 200 px
   :scale: 50 %
   :alt: this is a detailed caption of the image
   :align: left
   :loading: embed


and some text here

Step 3:
-------


.. _literature_link_1:


.. _literature_link_2:
