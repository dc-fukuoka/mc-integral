mc-integral - Monte Carlo integration test.  
intel MKL is requited.
------

<img src="https://latex.codecogs.com/gif.latex?\int_{a}^{b}dxp(x)f(x)&space;\approx&space;\frac{\left&space;(&space;b&space;-&space;a&space;\right&space;)}{N}\sum_{n=1}^{N}f(x_n)">

```
$ module load intel
$ make
$ ./a.out
 npts:   536870912
integrate gaussian function from  -1.00000E+05 to    1.00000E+05
 Monte Carlo integration:   1.77602888859703
 Exact solution:            1.77245385090552
```
