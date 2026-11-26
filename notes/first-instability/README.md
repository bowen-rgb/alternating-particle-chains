# First longitudinal instability near the exponential interaction

[中文证明 (PDF)](first-instability.pdf) · [中文源文件 (HTML)](first-instability.html) · [Numerical check](verify_first_mode.py)

This research note studies the ordered one-dimensional alternating chain of the [main manuscript](https://github.com/xiujiang4/alternating-particle-chains/blob/a48fe4c2623a29dd45edd227410a07baba400250/en/paper.tex). Let the potential magnitude be $f(r)=\exp(-r^{1+\delta})$ with $\delta>0$, in units where the interaction length is one. Starting with a dilute equidistant chain, which longitudinal mode first acquires negative quadratic energy under compression?

**Claim for independent review.** Allowing all even periods, the first linear instability is acoustic as $\delta\downarrow0$. Its critical mean gap satisfies $a_*/\sqrt{2\delta}\to1$. For each fixed even period $N\ge4$, the first mode for sufficiently small $\delta$ is the longest permitted wave. These statements concern the local Hessian. A competing configuration at finite displacement could have lower energy at an earlier stage of compression.

## Stability criterion

Write $a$ for the mean gap. Even periods admit wave numbers $q=2\pi j/N$. Fourier diagonalization gives, up to a positive factor,

$$
 D(q)=\sum_{m\ge1}(-1)^{m+1} f''(ma)(1-\cos mq),\qquad
 B(q)=D(q)/(1-\cos q).
$$

This uses the second variation in [Bétermin, Knüpfer and Nolte, §3.5, Eq. (3.4)](https://doi.org/10.1007/s10955-020-02603-2). Extend $B$ to $q=0$ by continuity, giving $B(0)=\sum(-1)^{m+1}m^2 f''(ma)$. This is the long-wave stiffness. The actual translation mode has $D(0)=0$ for every gap. Define $a_*$ as the supremum of gaps admitting a strictly negative mode for some even period. Allowed wave numbers over all such periods are dense. Any strictly negative value in the continuous band is detected at some finite period.

## The critical scaling

Put $b=\delta/a^2$ and $c=\cos(q/2)$. For $a\downarrow0$ with bounded $b$, on each interval $[0,\pi-\eta]$ with $\eta>0$,

$$
  B(q)/a=F_b(q)+o(1),\qquad
  F_b(q)=\frac1{8c^2}+\frac{b\log c}{1-\cos q}.
$$

The error tends to zero in $C^2$, with the removable value at zero. Along any sequence $q\to\pi$ in this scaling, $B(q)/a\to+\infty$.

To see this, write $f''(r)=e^{-r}-\delta e^{-r}/r+R(r)$. Denote the contribution of the first two terms by $B_0$, put $z=e^{-a}$ and $V=1+2z\cos q+z^2$. Geometric and logarithmic summation give

$$
 \frac{B_0(q)}a =
 \frac{z(1-z)}{a(1+z)V}+
 \frac{b}{2(1-\cos q)}\log\frac{V}{(1+z)^2}.
$$

Away from $\pi$, this converges to $F_b$ in the stated sense, including at the removable singularity at zero.

Here are details controlling the remainder. For $0<r\le1$ and $0\le j\le6$, there are constants $C_j,K$ such that

$$
 |R^{(j)}(r)|\le C_j(1+|\log r|)^K
       \big(\delta r^{-j}+\delta^2r^{-1-j}\big).
$$

Indeed, with $p=1+\delta$, the two differences are $p^2r^{2\delta}e^{-r^p}-e^{-r}$ and $-(\delta/r)(pr^\delta e^{-r^p}-e^{-r})$. Apply the mean value formula in $\delta$ to each parenthesized difference after spatial differentiation. For $r\ge1$ the same calculation gives a bound $C_j\delta(1+r)^K e^{-r/2}$.

For $0\le l\le4$ set $h(x)=x^l R(ax)$. Away from multiples of $2\pi$, summation by parts $l+2$ times bounds the Fourier sum of $h(m)$ by a constant times its first $l+2$ absolute values plus $\int_1^\infty|h^{(l+2)}(x)|\,dx$. This follows by writing each forward difference as an integral of the derivative over unit intervals. For $1\le x\le1/a$ the derivative is bounded by $C(1+|\log a|)^K(\delta x^{-2}+\delta^2x^{-3}/a)$. The tail integral is $O(\delta a)$ by the substitution $r=ax$. Thus these sums are $o(a)$ whenever $\delta=O(a^2)$. Applying the bounds at $\pi$ and $\pi\pm q$ gives $o(a)$ control through four derivatives for $D-D_0$, where $D_0=(1-\cos q)B_0$. Evenness and vanishing at zero allow division by $1-\cos q$ with two derivatives to spare.

Near $\pi$, absolute summation instead gives $\sum|R(ma)|=O(\delta/a+\delta^2(1+|\log a|)^{K+1}/a)=O(a)$. As $V\to0$, the positive $1/V$ contribution to $B_0/a$ dominates the logarithmic term. This proves the endpoint assertion.

## Mode selection

At zero, $F_b(0)=1/8-b/4$. For $0<q<\pi$ put $x=\cos^2(q/2)$. At the candidate transition,

$$
 F_{1/2}(q)=\frac{1-x+x\log x}{8x(1-x)}>0.
$$

Also $F_b(q)=1/8-b/4+(1-b)q^2/32+O(q^4)$ near zero. The $C^2$ estimate therefore gives strict increase of the exact $B$ close to zero for $b$ near $1/2$. The positive gap on compact interior subintervals and the bound at $\pi$ show that zero is the strict minimizer of $B$ in that regime.

For $b\le1/2-\epsilon$ and $a$ small, the whole band is positive by the same estimates and monotonicity of $F_b$ in $b$. For $b=1/2+\epsilon$, $B(0)<0$ and there are strictly negative nonzero modes. At gaps in compact subsets of $(0,\infty)$ the exact $B$ tends uniformly to the positive exponential-kernel expression as $\delta\to0$, by the exponential tail and $(1-\cos mq)/(1-\cos q)\le m^2$. At all sufficiently large gaps the nearest neighbor dominates uniformly: $|f''(ma)|/f''(a)\le2m^{2\delta}e^{-(m-1)a}$ for $m\ge2$. Together these estimates exclude an earlier instability on the dilute side and give $\delta/a_*^2\to1/2$. At the first boundary $B(0)=0$. Wave numbers admitting negative quadratic energy approach zero when this boundary is approached from the unstable side.

At fixed $N\ge4$ (even), each interior wave number has limiting critical ratio $(1-x)/(-2x\log x)$, strictly increasing with $q$. The endpoint $\pi$ stays positive on this scale. Writing $x_1=\cos^2(\pi/N)$ gives

$$
 \lim_{\delta\downarrow0}\frac{a_{*,N}^2}{\delta}
    =\frac{-2x_1\log x_1}{1-x_1}.
$$

Thus the first pair of conjugate modes is $j=1,N-1$. The case $N=2$ has only the endpoint mode.

## Reproduction

With Python 3.9+ and NumPy, run from the repository root:

~~~sh
python notes/first-instability/verify_first_mode.py
~~~

| $\delta$ | Acoustic root | Root / $\sqrt{2\delta}$ | Gap at $N=4$ | Gap at $N=16$ | Gap at $N=128$ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.01 | 0.16063676 | 1.1358734 | 0.13332533 | 0.15905629 | 0.16061218 |
| 0.001 | 0.04655057 | 1.0409023 | 0.03872025 | 0.04609762 | 0.04654352 |
| 0.0001 | 0.01431819 | 1.0124490 | 0.01191728 | 0.01417932 | 0.01431603 |

The script checks the acoustic residual with 38-digit decimal arithmetic and a longer cutoff. It also scans all nonzero modes at $N=512$ on a separate gap grid. AI assisted in the derivation, writing and code. The uniform remainder estimates and the passage to the first instability are the main points for independent review.
