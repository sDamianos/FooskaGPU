# FooskaGPU

FooskaGPU is a mulit-GPU accelerated CFD software developed from scratch by me. The code is written in CUDA using cuda aware MPI for communication and thus works only in NVIDIA GPUs. You can find preliminary versions of both single phase and mulitphase solver here. In the mulitphase solver you can find two versions one based in double precision and one in mixed floating point precission for higher performance.

The solver is density-based for 3D unstructure meshes, using an HLLC approximate Riemann solver for flux evaluation, with the WALE model for turbulence. Mass transfer is computed through gibbs free neergy equilibrium. High-order accuracy is achieved through a MUSCL reconstruction and a high-order Runge–Kutta time integration scheme.

A stream-based, GPU-aware, non-blocking MPI communication layer ensures excellent scalability, performing slightly better than P2P NCCL in my tests.

For validation cases and updates check my linkedIn account https://www.linkedin.com/in/sotiris-damianos-0272b3225/
