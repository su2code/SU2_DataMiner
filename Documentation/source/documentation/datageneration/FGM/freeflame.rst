.. _flamelet_solver_adiabatic:

Adiabatic Flamelet Solver
=========================

Adiabatic or 'free' flamelets are one-dimensional, unstretched, premixed flames for which a zero heat flux boundary condition is imposed at the inflow and outflow boundaries. 

.. autofunction:: Data_Generation.FlameletSolvers.FreeFlameSolver.__init__ 

Adiabatic flamelets are included in the manifold **by default**, but can be enabled or disabled manually by inserting the label "FREEFLAME" in the :code:`includeFlameletType` or :code:`excludeFlameletType` in the configuration.

The adiabatic flamelet solver is a wrapper for the Cantera `FreeFlame <https://cantera.org/3.0/sphinx/html/cython/onedim.html#freeflame>`_ module and is derived from the :ref:`.FlameletSolver_Cantera base class <flamelet_solver_base_class>`.

When running batches of adiabatic flamelet simulations, the reactant temperature is linearly varied between the maximum and minimum reactant temperature specified in the :ref:`SU2 DataMiner configuration <FGM>`. 

