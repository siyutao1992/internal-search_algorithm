# Internal-Search Algorithm for Constrained Optimization

* **Purpose**: This project aims to develop an "internal-search" algorithm for solving constrained optimization problems. A test case problem of reliabitility-based design optimization (RBDO) will be solved by the developed algorithm for demonstration.
* **Programming language**: Initially the algorithm is developed in MATLAB, but I may develop the Python version of it in the near future.
* **File descriptions**:
  * "InterSearchFunc.m": main function of the developed algorithm
  * "InterSearchFunc_rect.m": entrance function which rectifies the inputs (normalize the objective functions and variable bounds, etc.) before applying the internal search algorithm
  * "main_RBDO_case_test.m": main script for testing the optimization algorithm on an RBDO case problem.