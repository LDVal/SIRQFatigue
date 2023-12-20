## FFiles Description

### Fortran Files (.f90)

1. **Strategy_A_Scatter_vs_f.f90**
   - Description: Performs stochastic simulations of Strategy A on ER networks with cliques.
   - Output:  Generates a table where the first column represents the detection probability f, and the second column represents the corresponding fraction of recovered people R at the final stage. Each row corresponds to an individual realization.
   
2. **Strategy_A_Scatter_vs_I0.f90**
   - Description: Performs stochastic simulations of Strategy A on ER networks with cliques.
   - Output: Generates a table where the first column represents the initial fraction of infected individuals (I0), and the second column represents the corresponding fraction of recovered people R at the final stage. Each row corresponds to an individual realization.

3. **Strategy_B_Scatter_vs_f.f90**
   - Description: Performs stochastic simulations of Strategy B on ER networks with cliques.
   - Output:  Generates a table where the first column represents the detection probability f, and the second column represents the corresponding fraction of recovered people R at the final stage. Each row corresponds to an individual realization

4. **Strategy_B_Scatter_vs_I0.f90**
   - Description: Performs stochastic simulations of Strategy B on ER networks with cliques.
   - Output: Generates a table where the first column represents the initial fraction of infected individuals (I0), and the second column represents the corresponding fraction of recovered people R at the final stage. Each row corresponds to an individual realization.

   *Note: To model a random regular network (RR), set kminC equal to kmaxC and kminI equal to kmaxI in the Fortran programs.*

### Mathematica Notebooks (.nb)

5. **R0_StrategyA.nb**
   - Output: Generates a table where the first column corresponds to the values of beta, and the second column contains the corresponding values of f that satisfy the equation R0(beta,f)=1 for Strategy A.

6. **R0_StrategyB.nb**
   - Output: Generates a table where the first column corresponds to the values of beta, and the second column contains the corresponding values of f that satisfy the equation R0(beta,f)=1 for Strategy B.
