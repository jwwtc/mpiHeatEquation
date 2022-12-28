# mpiHeatEquation

## Abstract
Overall, MPI parallelization is a powerful tool for improving the efficiency and accuracy of heat conduction simulations, and it has numerous applications in various fields of engineering and science. Here, we are presenting the parallelized solution of the two-dimensional heat conduction equation on a square plate. We leverage collective communications and certain MPI functions to achieve convergence to the analytical solution. Metrics of scalability are widely used to illustrate the ability of both hardware and software to deliver greater compute with more hardware. Thus, we will also perform a weak and a strong scaling analysis of our program.

### Theory
For the numerical solution, we time-match:

$\frac{\partial \theta}{\partial t}=\kappa\left(\frac{\partial^2 \theta}{\partial x^2}+\frac{\partial^2 \theta}{\partial y^2}\right)$

to a steady-state solution on a uniform mesh, with the following finite difference approximation:


## Results

#### Numerical and Analytical Solutions Comparison

Surface Plot            |  Contour Plot
:-------------------------:|:-------------------------:
<img src="figures/numerical-temp-surf.png" style="width:600px;"/>  |  <img src="figures/comparison_temp.png" style="width:400px;"/>

#### Strong Scaling
<img src="figures/strong-scaling_time.png" style="width:400px;"/>

