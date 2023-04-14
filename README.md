# Boundary_Value_Problem-TMA_Method

In the present work, a second-degree equation in a Boundary Value Problem will be solved. 
The following situation is posed: a thick-walled container subjected to internal and external pressure with inner radius a and outer radius b. 
Additionally, the equation that models the radial displacement μ of a point along its thickness is given:

d**2μ/dr**2 + 1/r * dμ/dr - μ/r**2 = 0

The corresponding data for "a" and "b" as well as the values of μ at these points are provided. 
The problem will be solved using different numerical methods, and then compared against an analytical solution to draw conclusions.

In the first step, the values of μ are requested to be obtained. The data will be obtained using the Thomas method or Tridiagonal Matrix Algorithm (TMA).
At the same time, the analytical solution is provided by the following equation:

μ = C1 * r + C2/r

By performing the calculations and plotting the analytical μ against the one obtained by the Thomas method, 
it can be observed that the values obtained by the numerical method are quite consistent.

Subsequently, the maximum normal stress (σ_max) is to be obtained. 
The derivative is calculated using the numerical method of indeterminate coefficients, for which a centered derivative of order 1 and precision order 2 is used.
Note: During the procedure, as we have a vector of 10000 values of μ but we are using a centered derivative that requires a previous and a subsequent point 
and do not operate with ghost nodes, we will obtain a vector of 9998 numerical derivatives. Therefore, the same number of analytical derivatives is calculated.
For this solution, it is noted that a reasonable fit is obtained, and the proportional error is also plotted.

On the other hand, it is proposed to calculate the values of μ again using the shooting method.
This method involves finding the solution to a BVP by obtaining the result of two IVPs associated with it. 
For solving the initial value problems, both Laplace for obtaining the analytical solution and a numerical method for obtaining numerical approximations can be used; 
in the present resolution, a 4th order Runge-Kutta method was used as it is recommended and developed by the Burden-Faires book.
Note: The steps detailed in the previously mentioned bibliographic source were followed for the development of the calculation algorithm.
Once the numerical solutions are obtained, the relationship between them is known.
The solutions are solved and compared with the analytical solution, and the same indeterminate coefficients algorithm is developed to obtain the maximum stresses 
(again, the proportional errors are also plotted).

Finally, the maximum error values for each case are printed.


Discussion of results and conclusions:

First, as an interesting fact, it can be observed that the Thomas method offers a good fit in comparison to the short execution time required by its code. 
Solving the same matrix system using numpy functions requires times more than 10 times longer than using this method. 
On the other hand, obtaining results by the shooting method is much more accurate but requires a slightly longer execution time than Thomas' method.

As requested by the assignment, it can be seen that the modeling by the shooting method offers a modeling that meets the condition of relative error being less than 5%, 
as for Poisson = 0.3 and Poisson = 0.5, the maximum error values are 0.15% and 0.26% respectively.
In contrast, the Thomas method has a fast but less precise modeling, as it does not meet the condition of the assignment for Poisson = 0.5, 
but does so for Poisson = 0.3 with values of 7.73% and 4.12% respectively.

Finally, it is observed that obtaining the numerical derivatives using the method of indeterminate coefficients through a simple approach (second-order accuracy, requiring only 2 points) offers a very good fit when combined with the shooting method, yielding relative error values below unity. 
It should be noted that in almost all cases, at some point, the relative error becomes zero as the numerical curves "intersect" with the analytical curves. 
The only case where this does not happen is in the curve obtained by the shooting method for the values of μ corresponding to Poisson = 0.5.

Finally, as always, it is concluded that if the appropriate resources are available, it is convenient to use a more memory-demanding method such as the shooting method to obtain more accurate results. 
On the other hand, if resources are limited, the Thomas algorithm is a good tool to achieve relatively accurate results (according to the requirements of each case) without consuming a large amount of memory or requiring a long execution time.
