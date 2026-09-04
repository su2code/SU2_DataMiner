.. _flamelet_solver_adiabatic:

Adiabatic Flamelet Solver
=========================

Adiabatic or 'free' flamelets are one-dimensional, unstretched, premixed flames for which a zero heat flux boundary condition is imposed at the inflow and outflow boundaries. 

.. autofunction:: Data_Generation.FlameletSolvers.FreeFlameSolver.__init__ 

Adiabatic flamelets are included in the manifold **by default**, but can be enabled or disabled manually by inserting the label "FREEFLAME" in the :code:`includeFlameletType` or :code:`excludeFlameletType` in the configuration.

The adiabatic flamelet solver is a wrapper for the Cantera `FreeFlame <https://cantera.org/3.0/sphinx/html/cython/onedim.html#freeflame>`_ module and is derived from the :ref:`.FlameletSolver_Cantera base class <flamelet_solver_base_class>`.

When running batches of adiabatic flamelet simulations, the reactant temperature is linearly varied between the maximum and minimum reactant temperature specified in the :ref:`SU2 DataMiner configuration <FGM>`. 

Several tutorials on running adiabatic flamelet simulations can be found on :ref:`this page <freeflame_tutorial>`.


The value of the adiabatic mass flux must be retrieved when :ref:`burner-stabilized <flamelet_solver_burnerstabilized>` flamelets are included in the manifold.
The adiabatic mass flux of the converged adiabatic flamelet solution can be retrieved with the following function:

.. _adiabatic_mass_flux:

.. autofunction:: Data_Generation.FlameletSolvers.FreeFlameSolver.getMassFlowRate
