.. _doc_nicfd_tabulation:

.. sectionauthor:: Evert Bunschoten 

Tabulation methods for NICFD applications 
=========================================


This page documents the tabulation methods for NICFD applications in *SU2 DataMiner* 

.. contents:: :depth: 2


.. autofunction:: Manifold_Generation.LUT.LUTGenerators.SU2TableGenerator_NICFD.__init__ 


Refinement settings 
-------------------

.. autofunction:: Manifold_Generation.LUT.LUTGenerators.SU2TableGenerator_NICFD.SetTableDiscretization 
.. autofunction:: Manifold_Generation.LUT.LUTGenerators.SU2TableGenerator_NICFD.SetCellSize_Coarse
.. autofunction:: Manifold_Generation.LUT.LUTGenerators.SU2TableGenerator_NICFD.SetCellSize_Refined 
.. autofunction:: Manifold_Generation.LUT.LUTGenerators.SU2TableGenerator_NICFD.SetRefinement_Radius 
.. autofunction:: Manifold_Generation.LUT.LUTGenerators.SU2TableGenerator_NICFD.AddRefinementCriterion

Table Generation 
----------------

.. autofunction:: Manifold_Generation.LUT.LUTGenerators.SU2TableGenerator_NICFD.SetTableVars 
.. autofunction:: Manifold_Generation.LUT.LUTGenerators.SU2TableGenerator_NICFD.GenerateTable 

Output of tabulation files 
--------------------------

.. autofunction:: Manifold_Generation.LUT.LUTGenerators.SU2TableGenerator_NICFD.WriteTableFile

.. autofunction:: Manifold_Generation.LUT.LUTGenerators.SU2TableGenerator_NICFD.WriteOutParaview

