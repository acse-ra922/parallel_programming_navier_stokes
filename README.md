## Parallel Solver for Naver Stokes Equation
Using Method Programming Interface (MPI) to enhance the speedup of the Navier Stokes Equation

**Introduction and Code Overview**


The aim of this project is to develop a parallel solver for the Navier Stokes equation. Navier-Stokes equation, in fluid mechanics, a partial differential equation that describes the flow of incompressible fluids.

At each time step there will be an inner loop in which the pressure is solved using a projection method followed by an updating of the velocity.

<img width="545" alt="eqn1" src="https://github.com/acse-ra922/parallel_programming_navier_stokes/assets/110524155/a1c5c5bd-428a-419d-959a-7f1a6f1a52f2">

<img width="484" alt="eqn2" src="https://github.com/acse-ra922/parallel_programming_navier_stokes/assets/110524155/8040e325-04fb-489d-93b4-53d2f1815685">

<img width="467" alt="eqn3" src="https://github.com/acse-ra922/parallel_programming_navier_stokes/assets/110524155/574533b2-1c28-4637-917e-0addd6e5a3a6">

The code implements a parallelization strategy using the MPI library.
The parallelization is achieved by distributing the computational domain across multiple processes. Each 
process is responsible for a portion of the domain and performs computations on its local data. The domain 
is divided into multiple rows along the x-direction, and each process is assigned a subset of these rows.
The main steps involved in the parallelization of the code are as follows:
1. Domain decomposition:
 - The computational domain is divided into equal-sized subdomains along the x-direction based on the 
number of processes.
 - The number of rows and the starting row index are determined for each process.
2. Data structures:
 - Two-dimensional arrays are used to store the pressure (P) and velocity (vel) fields, respectively.
 - Each process creates its local arrays corresponding to its subdomain.
3. Pressure calculation:
 - The pressure values at the boundaries between subdomains are exchanged using non-blocking 
communication (MPI_Isend and MPI_Irecv).
4. Velocity calculation:
 - Each process updates the intermediate velocity field by performing advection and diffusion calculations 
on its local subdomain.
 
5. Output and visualization:
 - The computed results are written to output files for further analysis and visualization.
 - Each process writes its local portion of the velocity and pressure fields to separate files


**Performance Analysis**


The parallelized code was executed on the Imperial College HPC to evaluate its performance with varying numbers of processors and analyze the corresponding changes in execution time. The graph demonstrates a gradual decrease in execution time as the number of processors increases. This reduction in execution time can be attributed to the availability of more computing resources when using additional processors, thereby enhancing the computational power.
By distributing the workload among multiple processors, the computational tasks can be executed concurrently, further contributing to the reduction in execution time. This approach allows each processor to handle a smaller portion of the problem, preventing any single processor from becoming a bottleneck and maximizing resource utilization.
Additionally, the speedup factor was calculated to quantify the performance improvement achieved by parallelization. The speedup is determined by dividing the time taken by the serial code by the execution time of the parallelized code. In this case, a speedup of 35 was obtained for 128 processors, indicating that the parallelized code performs over thirty times better than the serial code.
Additionally, the efficiency of the parallelized program is calculated to be 0.27 for 128 processors. The graphs below display the speedup and efficiency for different number of processors.

<img width="540" alt="results1" src="https://github.com/acse-ra922/parallel_programming_navier_stokes/assets/110524155/142da7dd-fc52-478f-be6f-0b6334690852">
<img width="557" alt="results2" src="https://github.com/acse-ra922/parallel_programming_navier_stokes/assets/110524155/0af0d0bb-931b-4e18-9de0-5407622668be">


Additionally, the code's performance was analyzed using a fixed number of processors (8) across different grid sizes. The findings reveal that as the grid size increases, the computational workload intensifies, leading to longer execution times. Conversely, smaller grid sizes require less memory for data storage, resulting in faster computations and reduced communication overhead. This indicates that the code's efficiency is influenced by the interplay between grid size, computational workload, memory usage, and communication requirements.

<img width="500" alt="results3" src="https://github.com/acse-ra922/parallel_programming_navier_stokes/assets/110524155/b7c9a3c9-049e-49c7-8572-19d4dbf036e9">

Furthermore, the individual execution time for each Jacobi iteration was recorded and depicted in the accompanying graph. It demonstrates that as the number of processors increases, the time required to execute a single Jacobi iteration decreases noticeably. This behavior is anticipated since the function responsible for determining or executing the Jacobi iterations (pressure_jacobi) has been parallelized across the striped domain. By harnessing additional processors, the computational process is accelerated, resulting in faster computations and reduced iteration times.

<img width="514" alt="results4" src="https://github.com/acse-ra922/parallel_programming_navier_stokes/assets/110524155/c8e9eaf0-29e3-4859-adc2-17343ff657ce">

Finally, the performance analysis focused on the time taken per time step for varying numbers of processors. The results clearly demonstrate that as the number of processors increases, the time taken per time step decreases.This observation can be attributed to the implementation of parallelization, which enables the concurrent 
execution of computational tasks associated with each time step. By distributing the workload among multiple processors, the calculations can be performed simultaneously, resulting in a significant reduction in the overall time required to complete each time step.The graph visually depicts this relationship, highlighting the benefits of parallelization. As more processors are utilized, the computational tasks are efficiently divided and executed in parallel, leading to faster completion of each time step. This increased concurrency and workload distribution among processors effectively optimize the computational resources, resulting in improved performance.

<img width="466" alt="results5" src="https://github.com/acse-ra922/parallel_programming_navier_stokes/assets/110524155/1e318786-ea43-4176-bfb7-f4371976b46e">


**Animation**

Video/animations of the simulation output

https://github.com/acse-ra922/parallel_programming_navier_stokes/assets/110524155/f2d7c260-577a-4ed5-bbf5-6065acc929be



**Conclusion**


The present study employed the parallelized code on the Imperial College HPC to assess its performance under varying processor configurations and investigate the resultant changes in execution time. The results exhibited a gradual decrease in execution time as the number of processors increased. This decrease can be attributed to the availability of additional computing resources when employing extra processors, which enhances computational power.
Furthermore, the speedup factor was computed to quantify the performance enhancement achieved through parallelization. The obtained speedup of 35 indicated that the parallelized code performs over thirty times better than its serial counterpart. The analysis encompassed the individual execution time for each Jacobi iteration, visually depicted in the accompanying graph. The results indicated a noticeable decrease in the time required to execute a single Jacobi iteration as the number of processors increased. 
Ultimately, the study focused on the time taken per time step for varying numbers of processors. The findings convincingly demonstrated that as the number of processors increased, the time taken per time step decreased. Consequently, the utilization of more processors led to enhanced concurrency, workload distribution, and 
optimization of computational resources, consequently yielding improved performance



