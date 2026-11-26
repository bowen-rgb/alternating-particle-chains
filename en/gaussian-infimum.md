# The energy infimum for alternating Gaussian chains and the sharp cost of a minimum spacing

**GPT-6.2 Astra**

[English PDF](gaussian-infimum.pdf) · [中文 PDF](../zh/gaussian-infimum.pdf) · [中文全文](../zh/gaussian-infimum.md)

**Abstract.** We consider a one-dimensional periodic chain with an even number of particles of alternating sign. For a fixed period and particle count, we determine the energy infimum in an explicit range of Gaussian widths and characterize all configurations attaining it when collisions are allowed. If neighboring particles must remain at least $`h`$ apart, the minimum energy exceeds this infimum by a term quadratic in $`h`$ as $`h\downarrow0`$. We compute its sharp coefficient. Every minimizing configuration then concentrates near two points half a period apart. Each cluster contains an odd number of particles, and the two odd counts are as close as possible.

**Keywords.** Alternating particle chains, periodized Gaussian potential, rearrangement on the circle, energy minimization

## 1. The periodic model and its energy infimum

Bétermin, Knüpfer and Nolte [1, Conjecture 2.7] conjectured that if $`f(r)=F(r^2)`$ with $`F`$ completely monotone, then for the interaction triple $`\mathcal F=(-f,-f,f)`$ the equidistant chain $`e^\rho`$ uniquely minimizes the energy in the periodic alternating class $`\mathcal A_N^\rho`$ at every density $`\rho`$, up to translation. The Gaussian corresponds to $`F(t)=e^{-\alpha t}`$ and satisfies this hypothesis. Yet the four-particle example in [5] has lower energy than the equidistant chain. We fix $`N`$ and $`L=N/\rho`$ and ask how low the energy can fall over the entire class, and what it costs to exclude collisions.

For related work on equidistant ground states, high-density limits and competition between equidistant and bipartite structures, see [6–8]. Here the quadratic coefficient depends on the internal gaps of alternating clusters. A comparison restricted to equidistant or bipartite configurations does not determine it over the full admissible class.

Fix $`L,\alpha>0`$ and an even particle count $`N=2m`$. Write $`0\le x_0<\cdots<x_{N-1}<L`$, extend by $`x_{j+N}=x_j+L`$, and set $`\sigma_j=(-1)^j`$. Like signs interact through the negative Gaussian and unlike signs through the positive Gaussian. Counting ordered pairs as in [1, Definition 2.1], define

```math
G(x)=\sum_{k\in\mathbb Z}e^{-\alpha(x+kL)^2},\qquad Q=\sum_{i,j=0}^{N-1}\sigma_i\sigma_jG(x_i-x_j),\qquad E=1-\frac{Q}{N}.\qquad\text{(1)}
```

The constant $`1`$ removes each particle’s zero-distance self-interaction in its own cell. Images in other periods are included. Equation (1) extends continuously to coincident particles. Put

```math
s=\alpha L^2,\qquad q=e^{-\pi^2/s},\qquad c=\frac{\sqrt{\pi/\alpha}}{L},\qquad g(\theta)=G\!\left(\frac{L\theta}{2\pi}\right).\qquad\text{(2)}
```

**Theorem 1 (Infimum and equality).** Suppose $`q\le1/3`$. The energy per particle has infimum

```math
I_N=1-\frac{2}{N}\,[G(0)-G(L/2)]=1-\frac{8c}{N}\sum_{j=0}^{\infty}q^{(2j+1)^2}.\qquad\text{(3)}
```

over strictly ordered configurations. In the collision closure, equality holds exactly when cancellation of opposite signs at each point leaves a positive and a negative particle half a period apart. Neutral coincident pairs may lie elsewhere in the cell. Among strictly ordered configurations the antipodal pair attains the minimum for $`N=2`$; the infimum is unattained for $`N\ge4`$.

## 2. Rearrangement on the circle and quantitative stability

Use angular coordinates $`\theta=2\pi x/L`$. Let $`A`$ be the union of the positively oriented arcs from each positive particle to the next negative particle. Write $`|A|`$ for angular length. The distributional derivative of $`\mathbf1_A`$ on the circle is the signed particle measure $`\mu`$. For a smooth, even periodic function $`v`$, integration by parts in each variable gives

```math
Q_v(A):=\iint v(u-w)\,d\mu(u)\,d\mu(w)=\int_A\!\int_A[-v''(u-w)]\,du\,dw.\qquad\text{(4)}
```

If $`A`$ is a single arc of length $`a`$, its endpoints carry opposite unit masses, so $`Q_v(A)=2[v(0)-v(a)]`$. At collisions, zero-length arcs disappear or adjacent arcs merge. The same identity (4) continues to hold.

**Lemma 1 (Rearrangement of the signed interaction).** Let $`v`$ be smooth, even and $`2\pi`$-periodic. Suppose both $`v`$ and $`-v''`$ are nonincreasing on $`[0,\pi]`$. For any finite union of arcs $`A`$ and any semicircle $`H`$, one has $`Q_v(A)\le Q_v(H)`$.

**Proof.** Choose $`C=\max\{0,\max v''\}`$. Then $`K=C-v''`$ is nonnegative, even and nonincreasing on the half-period. The circle rearrangement inequality [2, Theorem 4.1] applies to this kernel and two nonnegative functions. Take both functions to be $`\mathbf1_A`$. Their symmetric decreasing rearrangement is the indicator of an arc $`J`$ with $`|J|=|A|`$. Hence

```math
\int_A\!\int_A K(u-w)\,du\,dw\ \le\ \int_J\!\int_J K(u-w)\,du\,dw.
```

Subtract $`C|A|^2`$ from both sides. This yields $`Q_v(A)\le Q_v(J)`$. By periodic symmetry and monotonicity, $`Q_v(J)=2[v(0)-v(|A|)]\le2[v(0)-v(\pi)]=Q_v(H)`$. $`\square`$

### The first Fourier mode and monotonicity of the residual kernel

Poisson summation gives $`g(\theta)=c[1+2\sum_{n\ge1}q^{n^2}\cos(n\theta)]`$, with absolute convergence after differentiation to any order. Set $`\beta=q/3`$ to retain an explicit positive margin, without optimizing this choice, and put $`v(\theta)=g(\theta)-2c\beta\cos\theta`$. For $`0<\theta<\pi`$, the bound $`|\sin(n\theta)|\le n\sin\theta`$ controls the higher modes in both $`v'`$ and $`(-v'')'`$ by $`2cqR(q)\sin\theta`$, where

```math
R(q)=\sum_{n=2}^{\infty}n^4q^{n^2-1}\le\frac{16}{27}\cdot\frac{48}{47}<\frac23.\qquad\text{(5)}
```

Indeed, the first term is at most $`16/27`$ and successive term ratios are at most $`(3/2)^4(1/3)^5=1/48`$. The first mode in both derivatives is $`-2c(q-\beta)\sin\theta`$, so Lemma 1 applies to $`v`$. Put $`z=\int_A e^{i\theta}\,d\theta`$. The mode we set aside contributes $`2c\beta|z|^2`$ to $`Q_g(A)`$, and a semicircle has $`|z|=2`$. Thus

```math
Q_g(H)-Q_g(A)\ge 2c\beta(4-|z|^2).\qquad\text{(6)}
```

**Theorem 2 (Stability relative to a semicircle).** Let $`d(A)=\min_H|A\triangle H|`$, where $`H`$ ranges over semicircles and $`\triangle`$ denotes symmetric difference. Under the assumptions of Theorem 1,

```math
E-I_N\ge\frac{16c\beta}{N}\left[1-\cos\!\left(\frac{d(A)}4\right)\right].\qquad\text{(7)}
```

**Proof of Theorems 1 and 2.** Set $`\phi=\arg z`$, with any choice if $`z=0`$, and write $`H_\phi=\{\cos(\theta-\phi)>0\}`$. Taking the real part of $`e^{-i\phi}z`$ gives

```math
2-|z|=\int_{A\triangle H_\phi}|\cos(\theta-\phi)|\,d\theta.\qquad\text{(8)}
```

Among sets of length $`a`$, the integral on the right is smallest on the four intervals of width $`a/4`$ flanking the two zeros. Its value there is $`4[1-\cos(a/4)]`$. Take $`a=|A\triangle H_\phi|\ge d(A)`$. Since $`4-|z|^2\ge2(2-|z|)`$, (6) implies (7) with semicircle energy $`1-2[g(0)-g(\pi)]/N`$.

The Fourier series of the semicircle energy gives (3). If this energy is attained, then $`Q_g(H)-Q_g(A)=0`$. The residual-kernel deficit and the first-mode deficit $`2c\beta(4-|z|^2)`$ are both nonnegative, so both vanish. Since $`\beta>0`$, we obtain $`|z|=2`$. The integral in (8) is then zero, and its weight vanishes only at the two endpoints. Thus $`|A\triangle H_\phi|=0`$. The two indicators have the same distributional derivative, an antipodal pair of opposite unit masses. Conversely, if the derivatives agree, the indicators differ by a constant; since both the semicircle and its complement have positive measure, this constant is zero. Strictly ordered particles give $`N`$ jumps. For any even $`N`$, fix a positive particle at $`0`$ and let the other $`N-1`$ alternating particles approach $`\pi`$ while remaining distinct. Their energies approach (3). $`\square`$

## 3. Minimum energy under a spacing constraint

Fix $`N=2m\ge4`$, $`L`$ and $`\alpha`$, with $`q\le1/3`$. Impose a lower bound $`0<h\le L/N`$ on every neighboring gap, including the gap across the period boundary. Denote the minimum energy by $`M_N(h)`$ and put

```math
\eta=\frac{2\pi h}{L},\qquad D=g''(\pi)-g''(0)=4c\sum_{\substack{n\ge1\\n\ \mathrm{odd}}}n^2q^{n^2}>0.\qquad\text{(9)}
```

**Theorem 3 (Small-spacing asymptotics and cluster sizes).** As $`h\downarrow0`$,

```math
M_N(h)=I_N+\frac{\lfloor m^2/2\rfloor}{N}D\,\eta^2+O_{N,L,\alpha}(\eta^4).\qquad\text{(10)}
```

For sufficiently small $`h`$, every minimizer can be rotated so that its particles lie in two clusters within angular distance $`O(\eta)`$ of a pair of antipodal points. The clusters contain $`2k+1`$ and $`2l+1`$ particles, with $`k+l=m-1`$ and $`|k-l|\le1`$.

### Existence and a trial upper bound

Fix one particle by translation. The admissible cyclic gaps form a nonempty compact simplex, so the continuous energy attains its minimum. Place one positive particle at $`0`$ and the remaining $`N-1`$ particles symmetrically about $`\pi`$ with spacing $`\eta`$. The two gaps between the clusters are $`\pi-(m-1)\eta\ge\eta`$. Each angle has the form $`\theta_i=b_i+t_i\eta`$, with $`b_i\in\{0,\pi\}`$ and $`|t_i|\le m-1`$. Since $`g'(0)=g'(\pi)=0`$, a second-order Taylor bound on each ordered pair gives

```math
0\le M_N(h)-I_N\le C_{\mathrm{tr}}\eta^2,\qquad C_{\mathrm{tr}}=2N(m-1)^2\|g''\|_\infty .\qquad\text{(11)}
```

### 3.1. Uniform localization of every particle

**Lemma 2 (A localization constant).** Set $`B=NC_{\mathrm{tr}}/(4c\beta)`$ and $`C_0=\max\{1,\tfrac12+\pi B\}`$. Every minimizer can be rotated so that all its particles lie within angular distance $`C_0\eta`$ of $`\{0,\pi\}`$.

**Proof.** Rotate the semicircle selected in (8) to $`H=(0,\pi)`$. Equations (6), (8) and (11), together with $`4-|z|^2\ge2(2-|z|)`$, give

```math
\int_{A\triangle H}|\sin\theta|\,d\theta\ \le\ B\eta^2.\qquad\text{(12)}
```

Consider any particle and let $`r`$ be its angular distance to the nearer endpoint. If $`r<\eta`$, the claim follows. If $`\eta\le r\le\pi/2`$, take an interval of length $`\eta/2`$ on each side of the particle. Neither interval contains another particle in its interior or crosses an endpoint of $`H`$. At the particle, $`\mathbf1_A`$ jumps while $`\mathbf1_H`$ is constant. One of the two intervals therefore belongs to the symmetric difference, up to endpoints.

Every point in that interval is at least $`r-\eta/2`$ from the nearer endpoint. Since $`\sin u\ge2u/\pi`$ for $`0\le u\le\pi/2`$, its contribution to (12) is at least $`(\eta/\pi)(r-\eta/2)`$. Comparison with the upper bound gives $`r\le(\tfrac12+\pi B)\eta`$. This constant is uniform over all particles and minimizers. $`\square`$

Choose $`\eta`$ so that $`C_0\eta<\pi/4`$ and $`B\eta^2<\sqrt2`$. The indicator $`\mathbf1_A`$ is constant on each of the two long arcs outside the clusters. If either value differed from $`\mathbf1_H`$, the whole of $`[\pi/4,3\pi/4]`$ or its antipodal arc would lie in the symmetric difference. The integral of the weight on either arc is $`\sqrt2`$, contradicting (12). Thus the cluster near $`0`$ has net sign $`+1`$, and that near $`\pi`$ has net sign $`-1`$. Alternation forces both counts to be odd.

### 3.2. The signed second moment of an alternating cluster

Normalize the signs in each odd cluster to $`\epsilon_i=(-1)^i`$. Write $`u_0,\ldots,u_{2k}`$ for the ordered offsets from the corresponding endpoint and define

```math
P=\sum_{i=0}^{2k}\epsilon_i u_i,\qquad V=\sum_{i=0}^{2k}\epsilon_i u_i^2-P^2.\qquad\text{(13)}
```

Since $`\sum_i\epsilon_i=1`$, translating the entire cluster leaves $`V`$ invariant. Set $`u_0=0`$ and write its successive gaps as $`a_1,b_1,\ldots,a_k,b_k`$. Then $`P=\sum_j b_j`$. Subtracting adjacent squares gives $`\sum_i\epsilon_i u_i^2=\sum_j(2u_{2j-1}b_j+b_j^2)`$. Substitution into (13) and cancellation of the terms from $`P^2`$ yield

```math
V=2\sum_{j=1}^k b_j\sum_{i=1}^j a_i\ \ge\ k(k+1)\eta^2.\qquad\text{(14)}
```

Equality holds when all internal gaps equal $`\eta`$. For a singleton both sides vanish.

### 3.3. A uniform fourth-order remainder and selection of cluster sizes

**Proof of Theorem 3, continued.** Use subscripts $`0`$ and $`\pi`$ for the two clusters. Normalize their internal signs as in (13), restoring the opposite net signs in the cross interaction. Set $`C_*=\max\{C_0,m-1\}`$. By Lemma 2, every offset of every minimizer satisfies $`|u_i|\le C_*\eta`$. The symmetric trials used below satisfy the same bound.

Each ordered pair contributes a signed copy of $`g(b+t)`$ at $`b=0`$ or $`\pi`$ modulo the period, with $`|t|\le2C_*\eta`$. Odd derivatives vanish at both centers. With $`K_4=\|g^{(4)}\|_\infty`$, Taylor’s theorem gives

```math
\left|g(b+t)-g(b)-\tfrac12 g''(b)t^2\right|\le \frac{K_4}{24}(2C_*\eta)^4.\qquad\text{(15)}
```

Within a cluster, $`\sum_{i,j}\epsilon_i\epsilon_j(u_i-u_j)^2=2V`$. The corresponding sum between normalized clusters is $`V_0+V_\pi+(P_0-P_\pi)^2`$. Counting both directions and using $`NE=N-Q`$, we obtain

```math
\begin{aligned}N(E-I_N)&=D(V_0+V_\pi)+g''(\pi)(P_0-P_\pi)^2+\mathcal R,\\ |\mathcal R|&\le C_R\eta^4,\qquad C_R=\frac{N^2K_4}{24}(2C_*)^4.\end{aligned}\qquad\text{(16)}
```

There are at most $`N^2`$ ordered pairs, so $`C_R`$ is independent of the minimizer and of its cluster counts. The Fourier series and (5) give $`g''(\pi)\ge2cq[1-R(q)]>0`$. Equations (14) and (16) now imply

```math
N(E-I_N)\ge D[k(k+1)+l(l+1)]\eta^2-C_R\eta^4 .
```

For nonnegative integers $`k+l=m-1`$ with $`k\ge l+2`$, replacing $`(k,l)`$ by $`(k-1,l+1)`$ reduces the bracket by $`2(k-l-1)`$. Its minimum is therefore $`\lfloor m^2/2\rfloor`$, attained precisely when $`|k-l|\le1`$.

Choose these counts and center the clusters symmetrically at $`0`$ and $`\pi`$, with all internal gaps equal to $`\eta`$. Both first moments vanish, (14) is an equality in each cluster, and the two outer gaps $`\pi-(m-1)\eta`$ are admissible. Applying (16) gives the matching upper bound and proves (10). An unbalanced split increases the quadratic coefficient of $`N(E-I_N)`$ by at least $`2D`$, so it cannot minimize when $`2D\eta^2>2C_R\eta^4`$. Combining this condition with feasibility and the two localization conditions after Lemma 2 gives the sufficient threshold

```math
h_0=\frac{L}{2\pi}\min\!\left\{\frac{2\pi}{N},\,\frac{\pi}{4C_0},\,\sqrt{\frac{\sqrt 2}{B}},\,\sqrt{\frac{D}{C_R}}\right\}.
```

The stated cluster counts hold for $`0<h<h_0`$. $`\square`$

For $`N=L=4`$ and $`\alpha=1/4`$, this gives $`h_0\approx1.19263\times10^{-5}`$. This sufficient bound comes from uniform remainder estimates and need not be the spacing at which the cluster pattern changes.

### Four particles versus six

With four particles the only optimal split is $`1+3`$. With six particles both $`1+5`$ and $`3+3`$ can collapse to the same antipodal dipole, but their tightly spaced symmetric trials have quadratic coefficients in the ratio $`6:4`$. The collision-limit equality classification alone does not select between them.

## 4. Numerical checks and the range of validity

At $`L=4`$, $`\alpha=1/4`$, the equidistant four-particle cell $`(0,1,2,3)`$ and the example $`(0,1.5,2,2.5)`$ in [5] have energies $`0.999633292157728`$ and $`0.924751872822186`$, respectively. Equation (3) gives $`I_4=0.849687099565508`$. Table 1 compares the coefficient in physical distance $`h`$ with the energy quotient of the symmetric trial.

**Table 1. Quadratic energy cost and finite-spacing check**

| $`N`$ | Cluster sizes | Theoretical coefficient | Finite-spacing quotient |
|--:|:--:|--:|--:|
| 4 | 1 + 3 | 0.37088222 | 0.37088192 |
| 6 | 3 + 3 | 0.49450963 | 0.49450892 |
| 8 | 3 + 5 | 0.74176445 | 0.74176200 |
| 10 | 5 + 5 | 0.89011734 | 0.89011313 |

The coefficient is $`\lfloor m^2/2\rfloor D(2\pi/L)^2/N`$. The last column gives $`(E-I_N)/h^2`$ at $`h=0.001`$ for the symmetric trial.

Exhaustive six-particle searches on grids of $`16,24,32`$ equispaced nodes give minima $`0.927980536870`$, $`0.912984112405`$ and $`0.907345330218`$, respectively. Each agrees with the symmetric $`3+3`$ trial on that grid. The repository script also compares real-space and Fourier sums and checks (7) on 675 fixed-seed random configurations.

### Competing configurations outside the width range

The condition $`q\le1/3`$ is equivalent to $`s=\alpha L^2\le\pi^2/\log3\approx8.984`$. Write $`E_{\rm dip}^{(N)}=1-2[G(0)-G(L/2)]/N`$ for the antipodal dipole energy, which remains defined outside this range. At $`L=1`$ and $`s=22`$, consider the six-particle cell in the collision closure

```math
X=(0,\tfrac18,\tfrac18,\tfrac14,\tfrac12,\tfrac34).\qquad\text{(17)}
```

The opposite signs at $`1/8`$ cancel. The remaining signed measure is that of four equidistant particles. Denoting their energy by $`E_{\rm eq}^{(4)}`$, the normalization in (1) gives

```math
E(X)=1-\tfrac46(1-E_{\rm eq}^{(4)})\approx0.665009396 < 0.669391181\approx E_{\rm dip}^{(6)}.\qquad\text{(18)}
```

Both values follow directly from (1). Truncation to $`|k|\le3`$, with $`|x_i-x_j|\le1`$, incurs an error per particle of at most $`12e^{-198}/(1-e^{-154})`$. By continuity, separating the coincident pair slightly preserves the strict inequality. Thus the dipole has lost global optimality at $`s=22`$. The competitor does not establish the six-particle global minimum.

For $`N=4`$ we fixed a positive particle at zero on a unit-period grid of 120 nodes and enumerated the other three positions in nondecreasing order, allowing collisions. At $`s=18`$ the grid minimum and the dipole energy both equal $`0.5111089813`$ to the digits shown. At $`s=22`$ the equidistant energy $`0.4975140939`$ is below the dipole energy $`0.5040867712`$. The energy difference of these two trials has opposite signs at $`s=21.7707`$ and $`21.7709`$. Bisection gives a crossing at $`s\approx21.77085`$; another configuration may win earlier. At fixed $`N\ge4,L`$, as $`\alpha\to\infty`$, the equidistant energy tends to zero while the dipole energy tends to $`1-2/N>0`$.

Inserting a coincident pair of opposite signs in an even cell preserves its signed measure and $`Q`$. If $`\widetilde E`$ denotes the energy after insertion, then

```math
\widetilde E-E_{\rm dip}^{(N+2)}=\frac{N}{N+2}\big(E-E_{\rm dip}^{(N)}\big).\qquad\text{(19)}
```

Let $`S_N`$ be the set of parameters $`s`$ at which some configuration has energy strictly below the dipole energy, and put $`s_*(N)=\inf S_N`$. By (19), $`S_N\subseteq S_{N+2}`$, hence $`s_*(N+2)\le s_*(N)`$. This does not require the sets to be intervals or imply strict monotonicity.

## 5. Conclusions

For $`q\le1/3`$, unconstrained minimizers in the collision closure can have neutral pairs at different points of the cell. Once a small positive minimum gap is imposed, Theorem 3 localizes every particle near two antipodal points and fixes the two cluster sizes. Equation (16) explains the selection. Each cluster has a lower bound on its signed second moment that depends on its size, and the relative displacement of the clusters contributes a nonnegative term. A single family of trial configurations attains both moment bounds and makes the displacement term vanish.

The period, particle count and Gaussian parameter are fixed throughout, and signs alternate strictly. The constants in (15)–(16) depend on these parameters, so these estimates do not give a thermodynamic limit at fixed density. Higher-order corrections to particle positions at positive spacing, the critical widths for different counts and the classification of minimizers outside the proven range remain open.

Quantitative rearrangement inequalities on groups are studied more generally in [3]. Here retaining part of the first Fourier mode yields an explicit constant even as the occupied length $`|A|`$ varies. For one-dimensional block minimization with competing interactions, see [4]. The small-spacing coefficient here comes from the signed second moment in (14).

### Numerical reproducibility

The standard-library script verify_gaussian_infimum.py accompanying [5], PR #1, reproduces Table 1 and the three six-particle grid minima. To check the four-particle enumeration independently, tabulate the periodic kernel (1) at grid differences, set $`x_0=0`$ and enumerate $`x_j=i_j/120`$ with $`0\le i_1\le i_2\le i_3<120`$. To reproduce (18), substitute the six positions in (17) directly into (1). These calculations check the values; the proofs of Theorems 1–3 are analytical.

### Research note

AI tools were used in the derivations, preparation of the manuscript and numerical checks. Theorems 1–3 have not undergone independent peer review or formal verification.

## References

[1] Bétermin L, Knüpfer H, Nolte F. Note on crystallization for alternating particle chains. Journal of Statistical Physics, 2020, 181: 803–815. DOI: [10.1007/s10955-020-02603-2](https://doi.org/10.1007/s10955-020-02603-2).

[2] Cabré X, Csató G, Mas A. Periodic solutions to integro-differential equations: variational formulation, symmetry, and regularity. Communications in Partial Differential Equations, 2025, 50(1–2). DOI: [10.1080/03605302.2024.2441851](https://doi.org/10.1080/03605302.2024.2441851).

[3] Christ M, Iliopoulou M. Inequalities of Riesz–Sobolev type for compact connected Abelian groups. arXiv:1808.08368, 2018.

[4] Buttà P, Esposito R, Giuliani A, Marra R. Froth-like minimizers of a nonlocal free energy functional with competing interactions. Communications in Mathematical Physics. DOI: [10.1007/s00220-013-1740-z](https://doi.org/10.1007/s00220-013-1740-z).

[5] A Gaussian counterexample and a Riesz bound for alternating particle chains. Repository manuscript, revision a48fe4c. The proof and numerical script accompanying this paper are in [PR #1](https://github.com/xiujiang4/alternating-particle-chains/pull/1).

[6] Ventevogel W J. On the configuration of a one-dimensional system of interacting particles with minimum potential energy per particle. Physica A, 1978, 92: 343–361.

[7] Sütő A. Ground state at high density. Communications in Mathematical Physics, 2011, 305: 657–710.

[8] Bétermin L, Šamaj L, Travěnec I. Equidistant versus bipartite ground states for 1D classical fluids at fixed particle density. Analysis and Mathematical Physics, 2025, 15: 88. DOI: [10.1007/s13324-025-01076-4](https://doi.org/10.1007/s13324-025-01076-4).
