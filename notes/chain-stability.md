# Equal spacing for exponential-power potentials and the Gaussian stability threshold

This contributed note uses the signed chain and pair-counting convention of the [main manuscript](../en/paper.tex). It extends the four-site example by allowing arbitrary even periods. The derivations and verification program were prepared with OpenAI Codex. [Full Chinese manuscript and figure](https://github.com/gamehhh7/alternating-chain-stability/blob/3a7012e285a7482a9caadc5130cc21bbe8057172/paper-zh.md).

## Which exponents give equal spacing at every density?

Let $h$ be the mean spacing. For $f_\alpha(r)=e^{-|r|^\alpha}$, $\alpha>0$, equal spacing uniquely minimizes energy up to translation **at every density and every even period if and only if $0<\alpha\le1$**. The positive direction also follows from the exponential-mixture theorem in the main manuscript. For $\alpha>1$, at each sufficiently high density there is a strictly ordered competitor of some finite even period.

To see why the latter holds, define $T_f(h)=\sum_{j\in\mathbb Z}(-1)^j f(|j|h)$, so $E_f(U_h)=f(0)-T_f(h)$. Assume $f$ is continuous, $f(0)>0$ and $|f(r)|=O((1+r)^{-1-\epsilon})$ for some $\epsilon>0$. At fixed $0<h<t$, the infimum over all periodic chains is at most

$$
f(0)-\frac{h}{t}T_f(t).
$$

Indeed, if $M/N=h/t$ with $M,N$ even, place $M$ equally spaced sites in length $Nh$ and insert $(N-M)/2$ pairs of opposite signs. At coincidence each pair cancels from the signed quadratic form, including every periodic image. The self-term in the energy remains $f(0)$ and the denominator is $N$. Separating the pairs by positive distances gives admissible chains whose energies converge to the bound. For arbitrary ratios, take $M/N\to h/t$ while fixing $h$ and using $t_N=Nh/M$. The decay assumption gives continuity of the sums. Thus all-density optimality requires $T_f(h)/h$ to be nonincreasing.

If $f'$ is locally absolutely continuous up to zero, $f''\in L^1(0,\infty)$, and $f'(\infty)=0$, then

$$
\lim_{h\downarrow0}\frac{T_f(h)}h=-\frac12 f'(0+).
$$

Group the signed sum as $\sum_{j\ge0}[f(2jh)-2f((2j+1)h)+f((2j+2)h)]$. Each difference integrates $f''$ against a triangle of height $h$ and base $2h$. After division by $h$ the periodic triangular wave has mean $1/2$. Its integral against each interval indicator tends to half the interval length, and approximation by step functions gives the limit for $f''$. If $f'(0+)=0$, choose $t$ large enough that $T_f(t)>0$, since $T_f(t)\to f(0)$. The bound then improves on uniform energy by at least $ch$ for some $c>0$ and all sufficiently small $h$. This applies to every $\alpha>1$.

For the positive direction there is an explicit remainder. Write $q_i=e^{-\kappa d_i}$, $t_i=(1-q_i)/(1+q_i)$ for gap $d_i$, and let $p_i=1+2\sum_{m\ge1}(-1)^m\prod_{\ell=0}^{m-1}q_{i+\ell}$. The exponential energy is $1-\frac1N\sum_i p_i$. With $\Delta_i=p_i-p_{i+1}$, the recurrence $p_i+q_i p_{i+1}=1-q_i$ gives $p_i+p_{i+1}=t_i(2-\Delta_i)$. Summing, and separately multiplying by $\Delta_i$ and telescoping, yields

$$
\sum_i(t_i-p_i)=\tfrac14\sum_i t_i\Delta_i^2 .
$$

Strict concavity of $\tanh(\kappa d_i/2)=t_i$ gives uniqueness. For $0<\alpha\le1$, the positive Laplace mixture of these exponential kernels yields the claim. The mixture representation is a standard consequence of the composition theorem for Bernstein functions. [Schilling, Song and Vondraček, Chapter 3](https://doi.org/10.1515/9783110269338).

## When does the Gaussian Hessian become indefinite?

For $f(r)=e^{-ar^2}$ set $b=ah^2$. With $k=2\pi l/N$, the Fourier eigenvalues of the Hessian have the sign of

$$
\Lambda_b(k)=\sum_{m\ge1}A_m(b)(1-\cos(mk)),\qquad
A_m(b)=(-1)^{m+1}(2bm^2-1)e^{-bm^2}.
$$

This follows by differentiating the pair energy twice and using discrete Fourier orthogonality. Put $D(b)=\sum_{m\ge1}m^2A_m(b)$. At $k=0$ the reduced symbol $\Lambda_b(k)/(1-\cos k)$ extends to $D(b)$. A negative $D(b)$ therefore gives a negative mode on a sufficiently large finite cell.

For $b\ge1$ the reduced symbol is strictly increasing with $k$ on $[0,\pi]$. To verify this, put $x=\cos k$ and expand $R_m=(1-\cos(mk))/(1-\cos k)=m+2\sum_{j=1}^{m-1}(m-j)\cos(jk)$. The derivative in $x$ obeys $|R'_m|\le m^2(m^2-1)/6$, since $|\sin(jk)/\sin k|\le j$. The $m=2$ term in the derivative is $-2(8b-1)e^{-4b}\le-14be^{-4b}$. The absolute sum of derivatives for $m\ge3$ is at most $(be^{-4b}/3)\sum_{m\ge3}m^6e^{-(m^2-4)}<2be^{-4b}$. The summands in the last series have consecutive ratios bounded by $(4/3)^6e^{-7}<0.006$.

To locate the unique zero of $D$, define $V(s)=\sum_{m\in\mathbb Z}(-1)^m e^{-s^2m^2}$. Then $V''(s)=-4D(s^2)$.

- For $s^2\le1$, Poisson summation expresses $V$ as $\sqrt\pi\,s^{-1}\sum_{j\in\mathbb Z}\exp[-\pi^2(j+1/2)^2/s^2]$. Each summand including the prefactor has positive second derivative, since $4z^2-10z+2>0$ for $z=\pi^2(j+1/2)^2/s^2\ge\pi^2/4$. Thus $D<0$.
- On $1\le b\le1.1$, the first term in $D'(b)$ exceeds $0.26$, the second is positive, and the absolute tail is at most $2.2\sum_{m\ge3}m^6 e^{-m^2}<0.21$.
- For $b\ge1.1$, the alternating tail beginning at $m=3$ is nonnegative, since $z(2z-1)e^{-z}$ decreases for $z\ge4$. The first two terms give $D(b)\ge (2b-1)e^{-b}-4(8b-1)e^{-4b}>0$. The last inequality holds at $1.1$, and $e^{3b}(2b-1)/(8b-1)$ increases thereafter.

The zero $b_s$ lies in $(1.08182608,1.08182609)$. For $b\ge b_s$ every fixed even cell has positive Hessian modulo translation. At $b_s$ the reduced symbol has strictly negative derivative in $x$ even at $x=1$, giving positive quartic stiffness at long wavelengths.

## An interval with lower-energy finite competitors

The derivative of $V(s)/s$ has the sign of $J(s)=sV'(s)-V(s)$ and $J'(s)=-4sD(s^2)$. By the preceding sign change $J$ first strictly increases and then strictly decreases. The Poisson formula gives $J(0+)=0$, and the direct series gives $J(+\infty)=-1$. Its unique positive zero is $\sqrt{b_d}$, with $2.47416700<b_d<2.47416702$. Since $T_f(h)=V(\sqrt a\,h)$, choosing $t=\sqrt{b_d/a}$ in the insertion bound yields a lower-energy chain for every $b<b_d$. Thus $b_s\le b<b_d$ is an interval of strict local minimality in every fixed cell and failure of global optimality across cells. The value $b_d$ describes this family of competitors. The global ground-state boundary is open.

At $a=3/2$ and $h=1$, repeat either six-site cell with period six:

| Positions | Energy per particle |
| --- | ---: |
| 0, 1, 2, 3, 4, 5 | 0.441305557786 |
| 0, 0.74, 0.76, 1.5, 3, 4.5 | 0.391612365072 |

Run the following from the repository root with Python 3.9 or later.

~~~sh
python3 code/verify_stability.py
~~~

The verifier uses rational bounds for the root endpoint signs, the tail inequalities and the energy improvement (strictly greater than 0.04969). It also checks the gap identity with fractions and compares three Fourier Hessian values with real-space finite differences. The arguments above give the continuum conclusions.
