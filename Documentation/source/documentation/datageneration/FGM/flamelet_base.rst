.. _flamelet_solver_base_class:

Cantera Flamelet Solver Base Class 
==================================

This page describes the general methods of the flamelet solver base class from which all other flamelet types are derived. Go to :ref:`this page <flamelettypes>` to read about the types of flamelets supported by *SU2 DataMiner*.

.. contents:: :depth: 2

.. _flamelet_base_init:

Initialization
--------------

The flamelet solver class is initialized from the :ref:`SU2 DataMiner configuration <FGM>`. When initializing the solver, the reaction mechanism and storage folder are read from the configuration. In addition, the species transport mechanism is retrieved.

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.__init__


.. _flamelet_inflow_conditions_base: 

Inflow Conditions
-----------------

The inflow conditions and initial condition largely depend on the flamelet type, but the base class allows for the temperature, pressure, and mixture of the premixed reactants to be specified.
 
.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.setReactantTemperature

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.setPressure

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.setMixtureStatus


.. _flamelet_grid_parameters:

1D Grid Parameters
------------------

The stability of the flamelet simulation and the resolution of the solution depend on the refinement of the one-dimensional grid.
The one-dimensional grid is automatically refined using the methods documented on the `Cantera webpage <https://cantera.org/dev/reference/onedim/grid-refinement.html>`_.
The parameters used to automatically refine the grid can be specified with the following method.

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.setGridRefinementCriteria

The flamelet simulation is initialized from an initial grid. The length and resolution of the initial grid are specified with the following method. 

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.setInitialGrid


.. _output_info:

Output Information
------------------

The following functions can be used to retrieve data from the flamelet solution and the simulation process. 
The verbosity level of the Cantera solver can be specified with 

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.setCanteraVerbose 

If the Cantera verbosity level is higher than zero, information on the convergence process is printed to the terminal.

The verbosity level of the general solution process can be specified with 

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.setSolverVerbose 

If the verbosity level of the solver is set to one, messages are displayed to the terminal on whether the flamelet solution converged or not.

The solution data of the converged flamelet simulation can be retrieved in the python API with the following method. 
The solution data frame contains information on the thermochemical state, as well as the grid (under the label "Distance"), and the inflow conditions. 

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.getSolution

The solution data can be saved to a csv file using the following method. The storage location is retrieved from the configuration class used to initialize the flamelet solver, with the sub-directory depending on the specific type of flamelet.

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.saveFlameletSolution


.. _running_simulation:

Running the Solver 
------------------

The methods described in this section regard running simulations of individual flamelets and in batches. 

Flamelet simulations can be initialized by manually specifying the inflow conditions described in the previous sections on this page, or they can be initialized by reading data from a converged flamelet solution file using the following method.
By calling this function, the initial guess of the flamelet solution is set according to the thermochemical state information extracted from the file.  

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.loadSolution

After specifying the flamelet-specific inflow conditions and solver settings, the simulation of an individual flamelet is initialized with the following function. 
Information on the solution process is displayed in the terminal, depending on the verbosity level.
When running a single simulation, the solution is initialized from the initial grid specified by :ref:`these settings <flamelet_grid_parameters>`. Subsequent calls of :code:`startSolver()` will initialize the solution from the **previously converged solution** to reduce computation time.

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.startSolver

Inflow conditions can also be specified through arguments. Specific information regarding the function arguments can be found on the documentation pages of the respective flamelet types.

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.solveFor

To automatically store the solution after convergence, the following function can be used.

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.solveAndSaveFor

For premixed flamelets, simulations can be run in batches with the following function. Flamelet simulations are run and saved over a range of inflow conditions, depending on the specific flamelet type.

.. autofunction:: Data_Generation.FlameletSolvers.FlameletSolver_Cantera.solveForMixtureStatus



