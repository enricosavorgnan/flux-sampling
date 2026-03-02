# Notes on Flux Sampling Papers


The following represent my personal notes on papers shared into `./papers` folder. \
Most of the papers are related to mathematical optimisation for flux sampling, i.e. ways for uniformly sampling values  
from a bounded space, \
This topic then finds applications in the field of metabolic flux analysis.


## R.L. Smith, 1984,  Efficient Monte Carlo Procedures for Generating Points Uniformly Distributed over Bounded Regions
Mathematical Optimisation. \

Problem: given a bounded $k$-dimensional surface $S \subset \mathbb{R}^n$, with $k < n$, we want to generate pseudorandomly a vector $(X_1, ..., X_n) \in S \sim \text{Unif}(S)$, i.e. s.t. $\mathbb{P}[X \in A \subset S] = V(A) / V(S) $, with $V$ content of Jordan measurable sets in $S$.

3 main techniques for uniform sampling:

1. Transformation Techniques \
    Idea: if you sample form $[0, 1]$ and then map into $S$ you have done.; this approach can be seen as composition of mappings. The problem is that to build the map is not trivial. 
2. Rejection Techniques \
    Idea: You can uniformly sample from a (well-behaved) region $D \supset S$, then keep only $X \in S$. The problem is that this approach suffers the curse of dimensionality, and the number of rejections grows with $k$.
3. Random Mixing Algorithms \
    Idea: choose a point $x_o \in S$ and then slightly modify it at random, so that $\forall i, x_i \in S$. The algorithm should show **strong reversibility**, meaning that the probability of moving from $y$ to $z$ in the region $S$ is equal to the one of moving from $z$ to $y$. The difficult part here is how to choose the set of the possible directions $D$:
    - Symmetry Mixing Algorithm: $D = {d in \mathbb{R}^n s.t. ||d|| = 1}$
    - Coordinate Directions Algorithm: $D = {d \in \mathbb{R}^n s.t. d=\pm e^i}$. This approach is clearly much convenient.
   Both the approaches are proved to work.

The latter approach will be called in literature *Hit-and-run*.


## L. Lovász, 1999, Hit-and-run mixes fast
Proves that **Hit-and-Run** approach *mixes fast*, i.e. the steps required for the method to approximate a uniform distribution are just a few. \
In particular, the walk mixes in $\mathcal{O}(n^2R^2) = \mathcal{O}(n^3)$, where $n$ is the dimension of the space defined by $S$, and $R$ is the radius of a ball containing $S$; the equation derives from the fact that it is possible to prove that $R = \mathcal{O}(\sqrt{n})$. \
The fact that the convex space needs to be **isotropic** (well-rounded), represents the main caveat to the result.

Future works established that an ideal approach consists into finding ways to make the system isotropic. Cousins, 2014 and Vempala, 2016 proposed ways for *approximate* isotropic transformations without solving the Semidefinite Program used to find the optimal affine transformation from a polytope $S$ to a rounded polytope $\mathcal{S}$.



## H. S. Haraldsdóttir et al., 2017, CHRR: Coordinate Hit-and-Run with rounding for uninform sampling of constraint-based models
An idea for having well-rounded surface $S$ is to ensure that $S$ is **isotropic**, i.e. is s.t. $\mathbb{E}_{X \sim S}[X] = 0$ and $\mathbb{E}_{X \sim S}[X^T X] = 1$. The transformation from a polytope to an isotropic polytope exists always. The problem is that the operation is expensive.

An alternative is to place the polytope $S$ in *John's position* (i.e. the maximal inscribed ellipsoid is the unit ball) by calculating the maximum volume inscribed ellipsoid and transforming it to the unit ball. Thus, the space is *rounded* from $S$ to $\mathcal{S}$ by applying a transformation $\mathcal{T}$. \
The right transformation $\mathcal{T}$ can be found using simple Linear Programs in $\mathcal{O}(n)$.
Then take $q$ steps of this Coordinate Hit-and-Run:
1. Pick a random coordinate direction $e_i$
2. Move from current point $v_k \in \mathcal{S}$ to a random point $v_{k+1} \in \mathcal{S}$ along $v_k + \alpha e_i \cap \mathcal{S}$.
3. Map samples back from $\mathcal{S}$ to $S$ by applying the inverse transformation $\mathcal{T}^{-1}$.


## W. Megchelenbrink et al., 2013, optGpSampler: An Improved Tool for Uniformly Sampling the Solution-Space of Genome-Scale Metabolic Networks
**OptGpSampler** merges the efficiency of GpSampler initialization with the ACHR sampling operations. Nothing more to notice here.






