# DistFlow
This function solves the distribution power flow equations for a given, generic, radial, single-phase distribution network based on the Baran(1989) paper in Matlab.

-The following 10-node grid is solved for the node voltages as an example. (Run Example.m) <br />
-The only imputs are grid parameters, i.e., number of nodes (N), grid graph (L), line impedences (rk, xk), node powers (pN, qN) and substation voltage (Vo). <br />
-For now, there are only two test feeders, i.e, KekatosGrid.mat and TwoNodeGrid.mat  <br/>
-The script that constructs the graph (L) for any generic distribution grid will be included later.  <br/>

<img src="https://user-images.githubusercontent.com/85322612/143179629-2f3a13f9-c82e-4697-ac8d-dddb0a558148.png" width="500"> 
<img src="https://user-images.githubusercontent.com/85322612/143179631-a30210c0-07bc-454b-a3a3-f3c76a9d260f.png" width="500">

```
V[0] = 1.0000 /_ 0.00
V[1] = 0.9858 /_ -0.18
V[2] = 0.9786 /_ -0.32
V[3] = 0.9804 /_ -0.39
V[4] = 0.9765 /_ -0.42
V[5] = 0.9729 /_ -0.44
V[6] = 0.9781 /_ -0.46
V[7] = 0.9785 /_ -0.49
V[8] = 0.9777 /_ -0.50
V[9] = 0.9644 /_ -0.66
Error: 0.000065 < tolerance = 0.000100 -> Y
```

