# mpiHeatEquation

## Abstract
Overall, MPI parallelization is a powerful tool for improving the efficiency and accuracy of heat conduction simulations, and it has numerous applications in various fields of engineering and science. Here, we are presenting the parallelized solution of the two-dimensional heat conduction equation on a square plate. We leverage collective communications and certain MPI functions to achieve convergence to the analytical solution. Metrics of scalability are widely used to illustrate the ability of both hardware and software to deliver greater compute with more hardware. Thus, we will also perform a weak and a strong scaling analysis of our program.

### Theory
For the numerical solution, we time-match:
![equation]([img]http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7B%5Cpartial%20%5Ctheta%7D%7B%5Cpartial%20t%7D%3D%5Ckappa%5Cleft%28%5Cfrac%7B%5Cpartial%5E2%20%5Ctheta%7D%7B%5Cpartial%20x%5E2%7D%2B%5Cfrac%7B%5Cpartial%5E2%20%5Ctheta%7D%7B%5Cpartial%20y%5E2%7D%5Cright%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0[/img])
to a steady-state solution on a uniform mesh, with the following finite difference approximation:


## Results

#### Numerical and Analytical Solutions Comparison

Surface Plot            |  Contour Plot
:-------------------------:|:-------------------------:
<img src="figures/numerical-temp-surf.png" style="width:600px;"/>  |  <img src="figures/comparison_temp.png" style="width:400px;"/>

#### Strong Scaling
<img src="figures/strong-scaling_time.png" style="width:400px;"/>

