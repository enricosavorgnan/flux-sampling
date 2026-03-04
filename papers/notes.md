# Notes on Flux Sampling Papers


The following represent my personal notes on papers shared into `./papers` folder. \
Most of the papers are related to mathematical optimisation for flux sampling, i.e. ways for uniformly sampling values  
from a bounded space, \
This topic then finds applications in the field of metabolic flux analysis.



### Overall Idea
Flux Sampling Analysis is a tool that aims at understand better the metabolism of a given organism or cell. \

First of all, to better describe the topic, it is needed to introduce the Stoichiometric Matrix $S \in \mathbb{R}^{m \times n}$. \
Each row of the matrix represents a different metabolite, while each column represents a different reaction. \
If a metabolite $i$ is consumed in reaction $j$, then $S_{ij} < 0$; if a metabolite $i$ is produced by reaction $j$, then $S_{ij} > 0$. If metabolite is not involved, then $S_{ij} = 0$. \
Since different reactions interest different metabolites, $S$ is typically a highly sparse matrix.

Usually, the interesting element of a reaction $j$ regards the concentration $x_i$ of the metabolite $i$ over time. The rate at which reaction $j$ converts substrates into products is called the flux $v_j$ , expressed as amount of metabolite processed per unit time. The two, together, are linked by the system of ODE $d\mathbf{X}/dt = S \mathbf{v}(t)$. \
However, usually $\mathbf{v}$ is a highly non-linear function, thus exact solution of the system is practically impossible. \
An interesting step, however, is given by the steady-state space, i.e. values for which $S \mathbf{v} = 0$. \
In such scenario, a solution is theoretically possible; however, given the high dimension of the elements and the fact that the system is intrinsically undetermined, it is still difficult to get relevant insights.

What's more is that usually $\mathbf{v}$ values are bounded by some physiological/thermodynamical factors, so that $V_{i,inf} \leq v_i \leq V_{i,sup}$.
The bound, together with the steady-state assumption, defines a convex polytope of feasible solutions, which is the main object of interest: indeed, it to write a linear optimization problem; the analysis of the steady-state space is called Flux Balance Analysis (FBA). \
The limitation of this approach is that it gives only a single solution, while the space of feasible solutions is typically huge. \

A solution for higher dimension problems is to base the analysis on the volume of the bounded space; the volume and its properties can be estimated by sampling point belonging to it via Monte Carlo algorithms.
The idea of a naive one is to pick points in a given space, calculate how many of them are in the volume of the polytope and take the ratio. This approach suffers of the curse of dimensionality, and is still practically not feasible.
Another idea is instead to consider that ratio between two different bounded regions: the original one has the boundaries defined by the previous introduced V_inf, V_sup arrays, the tighter ones are instead bounded by ulterior restrictions on V_inf and V_sup.

Another intelligent approach is to approximate uniform sampling via Artificial Centering Hit and Run: starting from a point inside the feasible region, move it sampling on a direction in the space and moving randomly on it. \
The crucial note is that the direction is chosen towards the "center" of the space, i.e. the average of the previously sampled points, which actually breaks the Markovian property of the sampling.
By sampling uniformly, it is possible to reconstruct the probability density function (PDF) for every single reaction, leading to uncertainty measures  quantification, like variances and covariances of the flux values, thus learning hidden biological modules between different reactions.

Other improvements are proposed by the GpSampler and OptGpSampler approach, consisting in several Monte Carlo Markov Chains run in parallel with different starting points and starting point conditions, and then points are taken after a given number of step, when their position is approximately uniform - the mixing problem.
Further improvement is provided by the Coordinate Hit and Run with Rounding, which applies the HR method but in a rounded (isotropic!) space. 
Indeed, CHRR finds the Maximal Volume Inscribed Ellipsoid of the original space, transforms it via an affine transformation and applies the HR method in the new space. Final points are then mapped back to the original space. CHRR is proved to efficiently lead to approximation of a uniform distribution, while previous approaches were not.


### R.L. Smith, 1984,  Efficient Monte Carlo Procedures for Generating Points Uniformly Distributed over Bounded Regions
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


### L. Lovász, 1999, Hit-and-run mixes fast
Proves that **Hit-and-Run** approach *mixes fast*, i.e. the steps required for the method to approximate a uniform distribution are just a few. \
In particular, the walk mixes in $\mathcal{O}(n^2R^2) = \mathcal{O}(n^3)$, where $n$ is the dimension of the space defined by $S$, and $R$ is the radius of a ball containing $S$; the equation derives from the fact that it is possible to prove that $R = \mathcal{O}(\sqrt{n})$. \
The fact that the convex space needs to be **isotropic** (well-rounded), represents the main caveat to the result.

Future works established that an ideal approach consists into finding ways to make the system isotropic. Cousins, 2014 and Vempala, 2016 proposed ways for *approximate* isotropic transformations without solving the Semidefinite Program used to find the optimal affine transformation from a polytope $S$ to a rounded polytope $\mathcal{S}$.



### H. S. Haraldsdóttir et al., 2017, CHRR: Coordinate Hit-and-Run with rounding for uninform sampling of constraint-based models
An idea for having well-rounded surface $S$ is to ensure that $S$ is **isotropic**, i.e. is s.t. $\mathbb{E}_{X \sim S}[X] = 0$ and $\mathbb{E}_{X \sim S}[X^T X] = 1$. The transformation from a polytope to an isotropic polytope exists always. The problem is that the operation is expensive.

An alternative is to place the polytope $S$ in *John's position* (i.e. the maximal inscribed ellipsoid is the unit ball) by calculating the maximum volume inscribed ellipsoid and transforming it to the unit ball. Thus, the space is *rounded* from $S$ to $\mathcal{S}$ by applying a transformation $\mathcal{T}$. \
The right transformation $\mathcal{T}$ can be found using simple Linear Programs in $\mathcal{O}(n)$.
Then take $q$ steps of this Coordinate Hit-and-Run:
1. Pick a random coordinate direction $e_i$
2. Move from current point $v_k \in \mathcal{S}$ to a random point $v_{k+1} \in \mathcal{S}$ along $v_k + \alpha e_i \cap \mathcal{S}$.
3. Map samples back from $\mathcal{S}$ to $S$ by applying the inverse transformation $\mathcal{T}^{-1}$.


### W. Megchelenbrink et al., 2013, optGpSampler: An Improved Tool for Uniformly Sampling the Solution-Space of Genome-Scale Metabolic Networks
**OptGpSampler** merges the efficiency of GpSampler initialization with the ACHR sampling operations. Nothing more to notice here.


### J. D. Orth, 2010, What is flux balance analysis?
A review on Flux Balance Analysis (FBA), a mathematical method for analysing the flow of metabolites through a metabolic network. FBA is based on the assumption that the system is in a steady state, meaning that the concentration of metabolites does not change over time. FBA uses linear programming to find the optimal flux distribution that maximizes or minimizes a given objective function, such as biomass production or energy generation.

FBA approaches the problem in 2 steps:
1. **Model Construction**: The first step is to construct a metabolic model, which consists of a set of reactions and their associated stoichiometry. The model can be represented as a matrix, where rows correspond to metabolites and columns correspond to reactions. The entries in the matrix represent the stoichiometric coefficients of the reactions.
2. **Optimization**: The second step is to optimize the flux distribution by solving a linear programming problem. The objective function is typically defined as a linear combination of the fluxes, and the constraints are defined by the stoichiometry of the reactions and any additional constraints, such as bounds on the fluxes or the requirement for certain reactions to be active or inactive.

FBA has several limitations, including the assumption of steady state, the reliance on accurate stoichiometric models, and the inability to capture dynamic changes in metabolism.


### S. J. Wilback, 2004, Monte Carlo  sampling can be used to determine the size and shape of the steady-state flux space
The paper is seminal on the specific topic of Flux Sampling Analysis.
The main reason for the introduction is related to the fact that calculating the size of the space of feasible solutions for the optimization problem of metabolic network modeling.

This constraint-based problem define a polytope $S$ by a system of (in)equalities. \
However, linear programming models rarely have the power of solving such complex problems: the problem is indeed #P-hard, so an exact solution is not possible. \
A possible solution is then presented by using Monte Carlo sampling into a (large) given space, and then counting the number of values in the feasible region.

The **steady-state flux space** is a problem defined as follows:
$$
    \begin{align}
        & S \mathbf{v} = 0 \\
        & 0 \leq v_i \leq V_{\text{max}, i}
    \end{align}
$$
A common assumption is that $S \in \mathbb{R}^{m \times n, +}$ has full row rank. 

The main drawback of a naïve MC approach is that it suffers the curse of dimensionality, meaning that in higher dimensions the numbers of points required for an exact approximation explodes. 
A more intelligent approach is to define an hypercube bounded by the "effective" $\mathbf{V}_{\max, i}$: the reference space is the original flux space, whether the volume is calculated after subsequent space reductions due to tighter new values of $\mathbf{V}_{\max, i}$.

**Steps** \
1. Determine the independent variables, i.e. project the flux polytope into a full-dimensional object.
2. Sample from the full-dimensional object; dependent-variables values are back-calculated based on the dependencies in the stoichiometric matrix. The result is a **random flux distribution** $\mathbf{v}$.
3. Determine whether $\mathbf{v}$ is valid or not, i.e. it is inside the space defined by the original $\mathbf{V}_{\text{max}}$ values. If valid, is stored.
4. Extra constraints are then applied in order to consider the fraction of reactions of interest (or including many other relevant physical properties). The resulting polytope $\mathcal{S}$ can be measured as
$$
    vol(\mathcal{S}) \approx vol(S) \frac{M}{N},
$$
where $M$ is the number of points inside $\mathcal{S}$, and $N$ the total valid ones.
5. The variance can also be calculated as $\mathbb{V} = \frac{p(1-p)}{N}$,
but the measure is minimal both for $p \approx 1$ (good) and for $p \approx 0$ (bad). To decouple, is okay to consider the ratio "variance to mean":
$$
    \mathbb{V}/\mathbb{M} = \frac{p(1-p)}{Np} = \frac{(1-p)}{N}  
$$

To find the right values of $\mathbb{V}_{\text{max}}$ is however a complex problem; usually physiologic values are too large and unuseful in practice. An approach is to use the $\alpha$-spectrum method.

The solution space is then characterized as a function of the extreme pathways by an SVD:
$$
    v_{\text{random}} = U \Sigma V^T,
$$
where the first column of $U$ captures the most information about the polytope and gives the principal direction, and $\Sigma$ is a diagonal matrix of singular values. 



