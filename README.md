# PIMC_demo
A demo for Path Integral Monte Carlo

This presentation will be dedicated to the primary application of the path integral. 
https://pitp.phas.ubc.ca/confs/sherbrooke2014/archives/pimc_notes_DelMaestro.pdf

https://arxiv.org/pdf/1712.08508

https://people.physics.illinois.edu/Ceperley/papers/095.pdf

https://www-zeuthen.desy.de/~kjansen/lattice/qcd/miscellaneous/CreutzFreedman.pdf (M. Creutz and B. Freedman, 1980)

# Preliminary
In PHYS 551 Quantum Theory (McGill University) we have discussed of the propagator

$$
\braket{ x_{f}, t_{f} |x_{i},t_{i}  }= \braket{ x_{f} |\exp(-i\hat{H}\Delta t/\hbar) |x_{i} } 
$$

which gives the probability amplitude of particle at position $x_{i}$ transiting to $x_{f}$ in the time $\Delta t=t_{f}-t_{i}$. 

Recall when we derived it we considered infinitesimal time $\Delta t=\delta t$ and calculated the short time propagator, by inserting completeness $\mathbf{1}=\int \frac{dp}{2\pi \hbar}\ket{p} \bra{p}$.

Which gives

$$
\braket{ x_{f}, t_{i}+\delta t |x_{i},t_{i}}=\sqrt{ \frac{m}{2\pi \hbar i\delta t} }\exp\left( \frac{iL\delta t}{\hbar}\right)  
$$

where,

$$
L=\frac{m}{2}\left( \frac{dx}{dt} \right)^2-V(x_{i})
$$

We do it by slicing time into many small $\delta t$ slices to write the propagator as a series of short time propagator, and insert position completeness to them:

$$
\braket{ x_{f},t_{f} | x_{i}, t_{i} }=\int \prod_{i}dx_{i}\braket{ x_{i+1},t_{i+1} |x_{i},t_{i}  }
$$

Expanding the short-time propagator with momentum completeness, taking a limit as$N\to \infty$ give,

$$
\braket{ x_{f} ,t_{f}| x_{i},t_{i} } = \int \mathcal Dx(t)\exp(-iS/\hbar)

$$
where $\mathcal Dx(t)=\lim_{ n \to \infty }\prod_{n=1}^{N-1}dx_{i}$ and

$$
S=\int_{t_{i}}^{t_{f}}L(x(t))dt=\int_{t_{i}}^{t_{f}}\left( \frac{m}{2}\left( \frac{dx}{dt} \right)^2-V(x(t)) \right)dt
$$

It is the sum of the Lagrangians of short-time propagators. Path integral is all nice, but notoriously difficult to solve for. And unfortunately, even using Monte Carlo, this is still not feasible due to the oscillatory convergence of the exponential. $\exp(-iS/\hbar)$.

# Imaginary Time Path Integral
So, the path integral for the purpose of multiple particles and how particles interact with their environment, so we should think about things in the stat mech way. The goal here is to reformulate stat mech into a path integral.

Recall the stat mech probability of a system in a state with energy $E_{n}$ called the Boltzmann factor,

$$
P(E_{n})=\frac{1}{Z}\exp(-E_{n}\beta )
$$

You probably have seen the partition function, (where $\beta=1/k_{B}T$), which is the normalization factor

$$
Z=\sum_{n}\exp(-\beta E_{n})
$$

And the observable,

$$
\braket{ \hat{O}  }= \frac{1}{Z}\sum_{n} \exp(-E_{n}\beta)\braket{ n | \hat{O}|n } 
$$

More generally,

$$
\braket{ \hat{O}  }= \frac{1}{Z}\sum_{n} P(E_{n})\braket{ n | \hat{O}|n } 
$$

We are working in the energy Eigenspace. Usually, we will start with some energy spectrum $E_{n}$ and do calculations with it. But usually solving for $E_n$ requires finding a solution to Schrodinger's equation with some crazy Hamiltonian:

$$
\hat{H} = -\frac{\hbar^2}{2m}\sum_{i}\nabla^2_{i}+\sum_{i}V_{ext}(\vec{r}_{i})+\sum_{i<j} \frac{q_{i}q_{j}}{|\vec{r}_{i}-\vec{r}_{j}|}
$$

(Urrh, ugly...) The first term is kinetic energy for n particles. The second term is the external potential applied to any particle. The third term is the Columb interaction between the charges.

$$
i\hbar  \frac{ \partial  }{ \partial t }\psi(\vec{x},t)=\hat{H}\psi(\vec{x},t)
$$

(You don't want to solve this)

The partition function is particularly useful because it tells us the free energy

$$
F=-\frac{1}{\beta}\ln Z
$$

Internal energy

$$
\braket{  E }_{\beta} =-\frac{ \partial }{ \partial \beta}\ln Z 
$$

BUT, the partition function is almost always impossible to find (unless we have some really good systems analytically). And no, of course, our monte Carlo algorithm that I will talk about later is impractical to find $Z$ as well. There is, however, a way to find expectation values bypassing the need for $Z$.

We can write $Z$ and $\braket{  \hat{O} }$ in canonical form where $\exp(-\beta E_{n})=\braket{ n | \exp(-\beta \hat{H})|n }$ hence

$$
Z=\sum_{n}\braket{ n | \exp(-\beta \hat{H}) |n} =\mathrm{Tr}(\exp(-\beta \hat{H}))
$$

and, to be more general, we can insert arbitrary basis completeness, in case $\hat{O}$ depends on other bases.

$$
\braket{ \hat{O}  }= \frac{1}{Z}\sum_{n,m} \braket{ n | \exp(-\beta \hat{H}) |m}\braket{ m| \hat{O}|n } =\frac{\mathrm{Tr}(\exp(-\beta \hat{H})\hat{O})}{\mathrm{Tr}(\exp(-\beta \hat{H}))}
$$

Now let's do it in a position basis, because the trace is invariant under a change of basis.

$$
Z=\int \braket{ x | \exp(-\beta\hat{H})|x }dx = \mathrm{Tr}(\exp(-\beta\hat{H}))
$$

and,

$$
\braket{ \hat{O}  }= \frac{1}{Z}\int dx' \int dx \braket{ x' | \exp(-\beta \hat{H}) |x}\braket{ x'| \hat{O}|x } =\frac{\mathrm{Tr}(\exp(-\beta \hat{H})\hat{O})}{\mathrm{Tr}(\exp(-\beta \hat{H}))}
$$

If $\hat{O}(\hat{x})$, i.e. $\hat{O}$ is a function of $\hat{x}$ then we can easily show $\braket{ x' | \hat{O}(\hat{x}) |x }=\braket{ x' | (O(x)|x })=O(x)\braket{ x' | x }=O(x)\delta_{xx'}$ 
is diagonal. So we can simplify to,

$$
\braket{ \hat{O}  }= \frac{1}{Z} \int dx \braket{ x | \exp(-\beta \hat{H}) |x}O(x)=\frac{\int dx \braket{ x | \exp(-\beta \hat{H}) |x}O(x) }{\int dx \braket{ x | \exp(-\beta\hat{H})|x }}
$$

Notice in these expression it contains $\braket{ x | \exp(-\beta\hat{H})|x }$ if a subsitution: we let $\tau=\hbar/k_{B}T=\hbar\beta$

$$
\braket{ x | \exp(-(\tau/\hbar)\hat{H})|x } 
$$

 $\braket{ x | \exp(-\hat{H}\tau/\hbar)|x }$ looks like a propagator with $\tau=-it$. i.e. the __density matrix__ is the imaginary time propagator. So the partition function is the integral over real space of all CLOSED ($x_{i}=x_{f}=x$) path integrals in the space. Although, it is still almost impossible to find $Z$ numerically, say we discretize $x$ into $X$ points and for each path we have $N_{\tau}$ points, the integration domain is $O(X^{N_{\tau}})$ and even a stochastic approach of Monte Carlo cant nicely estimate exponential path it.

Let's look at the integrand... Well, we can follow the same logic as the regular PI, with the substitutions with __Wick Rotation__:

$$
\begin{align}
t  & \to - i \tau\\ 
  it & \to\tau\\
dt  & \to - i d\tau
\end{align}
$$

Now

$$
\frac{dx}{dt}\to\frac{dx}{d\tau} \frac{d\tau}{dt}\to i \frac{dx}{d\tau}
$$

so 

$$
\left( \frac{dx}{dt} \right)^2 \to-\left( \frac{dx}{d\tau} \right)^2
$$

The short-time propagator becomes,

$$
\braket{ x_{f},\tau_{i}+\delta \tau | x_{i},\tau_{i} }=\sqrt{ \frac{m}{2\pi \hbar \tau} }\exp\left( \frac{L_{E}\delta \tau}{\hbar}\right)  
$$

where $L$ become $L_{E}$ defined as as,

$$
L_{E}=\frac{m}{2}\left( \frac{dx}{d\tau} \right)^2+V(x(\tau))
$$

is the Euclidean Lagrangian. Now the action becomes,

$$
S_{E}=\int_{\tau_{i}}^{\tau_{f}}L_{E}d\tau=\int_{\tau_{i}}^{\tau_{f}}\left( \frac{m}{2}\left( \frac{dx}{d\tau} \right)^2+V(x(\tau)) \right)d\tau
$$

Is the Euclidean action.

$$
\braket{ n | \exp(-\beta \hat{H})|n }=\braket{ x | \exp(-(\tau/\hbar)\hat{H})|x } = \int \mathcal Dx(\tau)\exp(-S_{E}/\hbar)
$$

Note that this path integral 

Together,

$$
Z=\int \braket{ x | \exp(-(\tau/\hbar)\hat{H})|x } dx = \int dx \int \mathcal Dx(\tau)\exp(-S_{E}/\hbar)
$$

similarly, 

$$
\braket{  \hat{O} } =\frac{\int dx \int \mathcal Dx(\tau)\exp(-S_{E}/\hbar)O(x)}{\int dx \int \mathcal Dx(\tau)\exp(-S_{E}/\hbar)}=\int dx \int\mathcal Dx(\tau)P[x(\tau)]O(x)
$$

We can merge the $dx$ into the measure. Previously, we considered a specific $x_{i}=x_{f}$. Merging the measure now means we are considering all $x$'s.

$$
\braket{  \hat{O} } =\frac{\int \mathcal Dx(\tau)\exp(-S_{E}/\hbar)O(x)}{\int \mathcal Dx(\tau)\exp(-S_{E}/\hbar)}=\int\mathcal Dx(\tau)P[x(\tau)]O(x)
$$

Now we effectively turned this energy problem into a path problem. The observation here is that the operator expected value went from the probability of an energy eigenstate to a probability of paths. 

$$
P(E_{n})= \frac{1}{Z} \exp(-E_{n}\beta )
$$

Becomes

$$
P[x(\tau)]=\frac{1}{Z} \exp(-S_{E}[x]/\hbar)
$$

But solving the path integral is hard. How can we do it?

# Motivation of Monte Carlo Integration
The idea of the Monte Carlo integration method is that many computations, like path integrals, are hard to carry out even numerically. Some of them just have too many terms to sum up. 

Let's begin with some Monte Carlo 
As a simple example:
This algorithm helps to integrate functions numerically. Recall the average value of a function is,

$$
\braket{  f } =\frac{1}{b-a}\int _{a}^b \, dx f(x) 
$$

So we have,

$$
\braket{  f } =\int _{a}^b \, dx \frac{f(x)}{(b-a)} 
$$

If $f(x_{1},\dots,x_{n})$ then its computationally expensive to do RHS directly, involving $\prod_{n} dx^n$ many differential. Use law of large number, sampling uniformly
$$\braket{  f }\approx\frac{1}{N}\sum_{i}f(x_{i})$$
then

$$
 \frac{1}{N}\sum_{i} \frac{f(x_{i})}{P_{\text{uniform}}(x_{i})} \approx \int_{a}^b dxf(x) \;\;\; \text{where} \, P_{\text{uniform}}(x_{i})=\frac{1}{b-a}
$$

There is no reason to sample uniformly for an arbitrary probability.

$$
 \frac{1}{N}\sum_{i} \frac{f(x_{i})}{P(x_{i})} \approx \int_{a}^b dx\frac{f(x)}{P(x)}P(x)
$$

This form can be matched with $\braket{  \hat{O} }=\int\mathcal Dx(\tau)P[x(\tau)]O(x)$. The idea is that we sample as many path as possible.
Hence,

$$
\braket{ \hat{O}  } _{\beta}\approx \frac{1}{N_{\text{samples}}}\sum_{i=1}^{N_{\text{samples}}}\braket{  \hat{O}(x_{i}(\tau)) } =\braket{ \hat{O}  } _{\beta}\approx \frac{1}{N_{\text{samples}}}\sum_{i=1}^{N_{\text{samples}}}\frac{1}{N_{\tau}} \sum_{\tau=1}^{N_{\tau}}\braket{  \hat{O} }_{\text{time slice}} 
$$

$$
\braket{  \hat{O}(x_{i}(\tau)) }=
$$

# Metropolis-Hastening Algorithm
Our goal is to find $Z$; in doing so, we need to generate many paths via the Metropolis algorithm. These path should have probability of $P=\exp\left( -\frac{S_{E}}{\hbar} \right)$.
1. Pick a random path $x^{(0)}=(x_{0},\dots,x_{N_{\tau}})$ with periodic condition $x_{0}=x_{N_\tau}$
2. Propose a perturbed path, $x'_{j}=x_{j}+\delta$ so $x^{(k)}=(x_{0},\dots,x_{j-1},x_{j}',x_{j+1},\dots,x_{N_{\tau}})$.
3. Calculate the acceptance probability $A=\min\left( 1, \exp\left( -\frac{\Delta S_{E}}{\hbar} \right) \right)$
If we proposed a path with a lower action than exp > 1, then the path is a more favourable configuration.
If the path has a higher action than exp < 1. The path is less favourable, so we accept it with lower probability
$$
x^{(0)}\to x^{(1)}\to x^{(2)}\to x^{(3)}\to x^{(4)}\to x^{(5)}\to x^{(6)}\to x^{(7)}\to x^{(8)}\to x^{(9)}\to x^{(10)}\to\dots
$$
Implementation:
- Thermalization (remove the first few paths) because the  original path is random, we almost never start with the correct distribution.
- Reduce correlation between paths, so they become more independent.
- Dynamic change rate.

# Verify correctness of Metropolis
For a sanity check, we try to test the result on a quantum H.O.

## Quantum H.O.
The Hamiltonian is,

$$
\hat{H}=\frac{\hat{p}^2}{2m}+\frac{m\omega^2\hat{x}^2}{2}
$$

The energy values are $E_{n}=\hbar\omega\left( n+\frac{1}{2} \right)$ for $n=0,1,2,\dots$

The wave function came from solving Schrödinger's equation.

$$
\psi_{n}=\left( \frac{m\omega}{\pi \hbar} \right)^{1/4} \frac{1}{\sqrt{ 2^nn! }}H_{n}\left( \sqrt{ \frac{m\omega}{\hbar} }x \right) \exp\left( -\frac{m\omega x^2}{2\hbar} \right)
$$

The expectation position is

$$
\braket{ x }_{n}=\braket{ n | x|n } =\int _{-\infty}^\infty \psi_{n}^*(x)x\psi_{n}(x)dx=0
$$

because the integral is odd. (Or we can argue the same using ladder operator with $\hat{x}=\sqrt{ \frac{\hbar}{2m\omega} }(\hat{a}+\hat{a}^{\dagger})$)

Using ladder operator $\hat{a}^{\dagger}, \hat{a},\hat{x}^2$ we can also derive the expectation variance,

$$
\braket{ x^2 }_{n}=\frac{\hbar}{2m\omega}(2n+1) 
$$

The easiest way to derive this is not through the integral, but instead through the ladder operator.
## Thermal H.O.
Put one harmonic oscillator in a thermal ensemble of particles with temperature $T$, $\beta=\frac{1}{k_{B}T}$. You have probably seen the derivation average position, $\braket{  x }_{\beta}$, energy $\braket{  E }_{\beta}$ and fluctuation $\braket{ x^2 }_{\beta}$ in stat mech before. Luckily, for the harmonic oscillator partition function can be solved explicitly.

$$
Z=\sum_{n}\exp(-\beta E_{n})=\sum_{n}\exp\left( -\beta \hbar\omega\left( n+\frac{1}{2} \right) \right)=\exp\left( -\frac{\beta \hbar\omega}{2} \right)\sum_{n}\exp\left( -\beta \hbar\omega n \right)=\frac{\exp\left( -\frac{\beta \hbar\omega}{2} \right)}{1-\exp(-\beta \hbar\omega)}
$$

The second term is a geometric serie,s so

$$
Z=\frac{\exp\left( -\frac{\beta \hbar\omega}{2} \right)}{1-\exp(-\beta \hbar\omega)}
$$

### Average Thermal Position
The average position is 

$$
\braket{  x }_{\beta}=\frac{1}{Z}\sum_{n}\braket{ n | x| n }\exp(-\beta E_{n})=0
$$

Not interesting... potential is symmetric.
### Average Thermal Energy
The average energy is $\braket{  E }=\frac{1}{Z}\sum_{n=0}^\infty E_{n}\exp(-\beta E_{n})$.
The numerator is just a geometric series, together

$$
\braket{  E }_{\beta} =\frac{1}{Z}\sum_{n=0}^\infty E_{n}\exp(-\beta E_{n})=\frac{\hbar\omega}{2}\coth\left( \frac{\beta \hbar\omega}{2} \right)
$$

### Average Thermal Fluxation
Another measurable we are interested in is $\braket{  x^2 }$. Where 

$$
\braket{  x^2 }_{\beta}=\frac{1}{Z}\sum_{n=0}^\infty \exp(-\beta E_{n})\underbrace{ \int_{-\infty}^\infty dx\psi_{n}(x)^2x^2 }_{ \braket{  x^2 }_{n}=\frac{\hbar}{2m\omega}(2n+1 ) }=\frac{\hbar}{2m\omega}\coth\left( \frac{\beta \hbar\omega}{2} \right)
$$

This is thermal variance. This is a classical result from simplifying a bunch of geometric series. 

In summary

| Quantity                  | Analytical Expression                                                   |
| ------------------------- | ----------------------------------------------------------------------- |
| $\braket{  x }_{\beta}$   | $0$                                                                     |
| $\braket{  E }_{\beta}$   | $\frac{\hbar\omega}{2}\coth\left( \frac{\beta \hbar\omega}{2} \right)$  |
| $\braket{  x^2 }_{\beta}$ | $\frac{\hbar}{2m\omega}\coth\left( \frac{\beta \hbar\omega}{2} \right)$ |

## PIMC on Thermal H.O.
Now lets do it with PIMC

$$
Z=\int\prod_{i}dx(\tau_{i}) \sqrt{ \frac{m}{2\pi \hbar\delta \tau} }\exp\left( -\frac{\delta \tau}{\hbar}\sum_{i} \left( \frac{m}{2}\left( \frac{x_{i+1}-x_{i}}{\delta \tau} \right)^2+\frac{m\omega^2x_{i}^2}{2} \right) \right)
$$

The energy estimator is the derivative of $Z$.

$$
\braket{  \hat{E} }_{\beta} =-\frac{ \partial }{ \partial \beta}\ln Z 
$$

$$
\braket{  \hat{T}  }_{\text{time slice}} = \frac{1}{2\delta \tau}-\frac{m}{2\delta \tau^2}\braket{  (x_{i+1}-x_{i})^2 } 
$$

$$
\braket{   x^2}_{\text{time slice}} = x_{i}^2 
$$

Put them all together,

$$
\braket{ x^2 } _{\beta} \approx \frac{1}{N_{\text{sample}}}\sum_{s=1}^{N_{\text{sample}}}\left( \frac{1}{N_{\tau}}\sum_{i=0}^{N_{\tau}-1}(x_{i}^{(s)})^2 \right)
$$

$$
\braket{  H } =\frac{1}{N_{\text{sample}}}\sum_{s=1}^{N_{\text{sample}}}\left( \frac{1}{N_{\tau}}\sum_{i=0}^{N_{\tau}-1}\left( \frac{1}{2\delta\tau} -\frac{m}{2\delta \tau^2}(x_{i+1}^{(s)} - x_{i}^{(s)})^2+\frac{m\omega^2}{2}(x_{i}^{(s)})^2\right) \right)
$$

## Misc.
The correlation is calculated with the Wiener-Khinchin Theorem, which is fast with large samples.

$$
\rho(m) = 
\frac{
\frac{1}{N} \sum_{k=0}^{N-1}
\left|
\sum_{n=0}^{N-1} (x_n - \bar{x}) e^{- 2\pi i k n / N}
\right|^{2}
\, e^{2\pi i k m / N}}{
\left[
\frac{1}{N}
\sum_{k=0}^{N-1}
\left|
\sum_{n=0}^{N-1} (x_n - \bar{x}) e^{- 2\pi i k n / N}
\right|^{2}
\right]}
$$

Bosons?

$$
P=\frac{1}{Z}\exp(-n_{s}E_{s}\beta+n_{s}\mu\beta)
$$

where $n_{s}=0, 1, 2, \dots$

$$
Z_{s}=\frac{1}{1-\exp(-\beta(E_{s}-\mu))}
$$

Fermions? $n_{s}=0\text{ or }1$

$$
Z_{s}=1+\exp(-E\beta+\mu\beta)
$$

More thorough analysis needed...
