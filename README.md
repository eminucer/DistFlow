# DistFlow Solver
This function solves the distribution power flow equations for a given, generic, radial, single-phase distribution network based on the DistFlow equations introduced by Baran&Wu (1989).

-Example.m file showcases the solution for a 33-bus distribution grid network. (Run `Example.m`) <br />
-The following are the resulting plots for the network graph and node-voltage profiles <br />
-The solver function `SolveDistFlow(Grid,PlotGraph,MaxIter,tol)` takes 4 inputs <br />
-Grid is a structure that contains grid-related parameters, i.e., number of nodes `(N)`, grid graph `(L)`, line impedences `(rk, xk)`, node powers `(pN, qN)` and substation voltage `(Vo)` <br />
-`PlotGraph (1,0)` enables plotting the results after execution <br/>
-`MaxIter` sets the maximum number of iterations for the Newton-Raphson method <br/>
-`tol` is the error tolerance <br/>

-There are also different test feeders under `\Grid Data`. You can uncomment the corresponding lines in `Example.m` and run it to see the results for yourself<br/>
-You can use the first code section in `Example.m` to create your own custom distribution grid and solve it<br/>

<img src="https://user-images.githubusercontent.com/85322612/143732190-018a04e8-1ef0-457b-9036-591229d37bc6.png" width="700"> 
<img src="https://user-images.githubusercontent.com/85322612/143732191-92687466-974d-4fce-8ccf-7b0a7c1b7a9b.png" width="700">

```
Iteration: 9
Terminal node error: 0.000002< tolerance = 0.000100 -> Y

Total Num of Iterations: 9

Node voltages and phase angles:
V[0] = 1.0000 /_ 0.00
V[1] = 0.9970 /_ 0.01
V[2] = 0.9829 /_ 0.10
V[3] = 0.9755 /_ 0.16
V[4] = 0.9681 /_ 0.23
V[5] = 0.9497 /_ 0.13
V[6] = 0.9462 /_ -0.10
V[7] = 0.9413 /_ -0.06
V[8] = 0.9351 /_ -0.13
V[9] = 0.9292 /_ -0.20
V[10] = 0.9284 /_ -0.19
V[11] = 0.9269 /_ -0.18
V[12] = 0.9208 /_ -0.27
V[13] = 0.9185 /_ -0.35
V[14] = 0.9171 /_ -0.38
V[15] = 0.9157 /_ -0.41
V[16] = 0.9137 /_ -0.49
V[17] = 0.9131 /_ -0.50
V[18] = 0.9965 /_ -0.51
V[19] = 0.9929 /_ -0.57
V[20] = 0.9922 /_ -0.59
V[21] = 0.9916 /_ -0.61
V[22] = 0.9794 /_ -0.65
V[23] = 0.9727 /_ -0.73
V[24] = 0.9694 /_ -0.78
V[25] = 0.9477 /_ -0.74
V[26] = 0.9452 /_ -0.68
V[27] = 0.9337 /_ -0.60
V[28] = 0.9255 /_ -0.52
V[29] = 0.9220 /_ -0.42
V[30] = 0.9178 /_ -0.50
V[31] = 0.9169 /_ -0.52
V[32] = 0.9166 /_ -0.53
```

