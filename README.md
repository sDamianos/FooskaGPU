# FooskaGPU

FooskaGPU is a mulit-GPU accelerated CFD software developed forom scratch be me. The code is written in CUDA using cuda aware MPI foro communication and thus works only in NVIDIA GPUs. The code provided here is just a first version of the one phase solver, while to two phase solver is under development. 

The solver is density-based for 3D unstructure meshes, using an HLLC approximate Riemann solver for flux evaluation, with the WALE model for turbulence. High-order accuracy is achieved through a MUSCL reconstruction and a high-order Runge–Kutta time integration scheme.

A stream-based, GPU-aware, non-blocking MPI communication layer ensures excellent scalability, performing slightly better than P2P NCCL in my tests.

For validation cases and updates check my linkedIn account https://www.linkedin.com/in/sotiris-damianos-0272b3225/
