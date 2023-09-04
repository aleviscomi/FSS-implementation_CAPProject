# "Fish School Search" in x86-32+SSE assembly, x86-64+AVX assembly, and OpenMP

*Project related to the '**Computer Architectures and Programming**' exam @ **University of Calabria.***   
***TEAM MEMBERS:** Alessandro Viscomi, Alessia Ramogida, Gianluca Massara*

## Project Description

Fish School Search (FSS) is an optimization algorithm inspired by the behavior of fish schools. The basic idea of the algorithm arises from the feeding and coordinated movement mechanisms typical of this context. In particular, fish swim towards the 'positive gradient' to feed and gain weight. Since heavier fish have a greater influence on the school in the collective search for food, the center of mass of the school tends to move towards the best locations in the search space.

FSS consists of three main steps:
1. Individual movement: during this phase, each fish makes a random movement. If the new position is better, it keeps it; otherwise, it returns to the previous one.

2. Instinctive movement: during this phase, fish are attracted to the most promising areas.

3. Volitive movement: this phase is responsible for regulating the balance between exploration and exploitation in the search process.

## Design Implementation

The implementation of the algorithm has been divided into three main phases:

1. **High-Level Implementation**: This phase involves implementing the algorithm in the C language.

2. **Performance Evaluation**: An intermediate phase during which the execution times of various functions implemented in C were analyzed to determine which ones were most suitable for optimization.

3. **Optimization Implementation**: This phase involves implementing the slower functions in Assembly language with the support of x86 instruction set extensions such as SSE and AVX.

Regarding the third phase, the optimization techniques used include:

- **Code Vectorization**: This is the primary SIMD-based technique. SIMD-based techniques allow you to work on multiple data with a single instruction.

- **Loop Unrolling**: This is an ILP technique that involves "unrolling of loops", meaning executing multiple instructions in the same iteration, avoiding multiple conditional operations. ILP (Instruction Level Parallelism) techniques involve executing sequential instructions, that have no data dependencies, in parallel.

- **OpenMP Paradigm**: This is a MIMD-based technique that allows, through certain directives, specifying that certain portions of code can be processed independently, enabling the compiler to translate the code to run on multiple CPUs simultaneously. MIMD-based techniques specifically aim to leverage the presence of multi-cores in CPUs.
