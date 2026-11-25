# A Gaussian counterexample and a Riesz bound for alternating particle chains

**Codex (OpenAI AI model)** · Corresponding author / 通讯作者 **guowei jiang**  
Research draft / 研究手稿

| Manuscript / 手稿 | PDF | LaTeX |
| --- | --- | --- |
| English | [Read the paper](en/paper.pdf) | [Source](en/paper.tex) |
| 中文 | [阅读中文稿](zh/paper.pdf) | [源文件](zh/paper.tex) |

## Results

Must a periodic chain of two alternating species have equally spaced particles to minimize its energy? Bétermin, Knüpfer and Nolte proposed this for pair potentials completely monotone in the **squared distance**, with attraction between like species and repulsion between unlike species. The paper gives a four-particle counterexample to Conjecture 2.7 as stated in the published article. [Original conjecture](https://doi.org/10.1007/s10955-020-02603-2)

For the Gaussian potential $f(r)=\exp(-r^2/4)$, moving the two particles of one species inward lowers the energy while keeping the other two fixed. Both cells repeat with period four. The table uses the original paper's energy per particle, counting each pair in both directions.

| Positions in one cell / 周期内位置 | Energy per particle / 每粒子能量 |
| --- | ---: |
| $(0,1,2,3)$ | 0.999633292158 |
| $(0,1.5,2,2.5)$ | 0.924751872822 |

The Fourier proof gives a strict inequality and shows that arbitrarily small displacements in this direction also lower the energy. For positive mixtures of exponentials in the **distance**, the paper gives a direct proof that equal spacing uniquely minimizes energy. This includes every inverse power $r^{-p}$ with $p>0$, the range proposed in Remark 2.6. The inverse-power result also follows from earlier bounds for one-dimensional spin domains. [Earlier bound, Eqs. (2.4)–(2.6)](https://doi.org/10.1007/s00220-013-1740-z)

## 中文摘要

一维周期链上有两类交替排列的粒子，同类相吸，异类相斥。Bétermin、Knüpfer 和 Nolte 在期刊版猜想 2.7 中提出，当势函数关于**距离的平方**完全单调时，等间距排列在任意密度下都唯一地使能量最低。上表中的高斯势反例保持周期、密度和粒子顺序，仅移动同一类的两个粒子，每粒子能量就下降约 0.0748814。手稿用傅里叶级数严格证明能量下降，并证明沿同一方向任意小的位移也能降低能量。[原猜想](https://doi.org/10.1007/s10955-020-02603-2)

对于关于**距离本身**的衰减指数函数的正混合，手稿直接从粒子间距出发证明等间距排列是唯一最优解。结论覆盖原文评注 2.6 所问的全部反幂势 $r^{-p}$（$p>0$），也可由此前对一维自旋畴的估计推出。[已有估计，式 (2.4)–(2.6)](https://doi.org/10.1007/s00220-013-1740-z)

## Lean proof / 形式化证明

Lean 4.28.0 verifies the four-site Gaussian counterexample in two ways. The Fourier argument proves the strict inequality for the specified sites and the Gaussian satisfies the conjecture’s squared-distance hypothesis. A direct real-space proof also gives $E(U)-E(X)>1/36$. Both proofs include every periodic image. 两条机器验证的证明覆盖完整无穷链的高斯反例。手稿中的局部不稳定性与里斯势结论尚未形式化。[Proofs and build instructions / 证明与复现说明](lean/README.md) · [Main theorems / 主定理](lean/Counterexample.lean)

## Numerical verification / 数值复核

Run from the repository root with Python 3. The script uses only the standard library. 在仓库根目录运行下列命令即可，只需 Python 3 标准库。

```sh
python3 code/verify_chains.py
```

The script compares direct and Fourier energy sums for both cells, checks four sample displacements, and evaluates the exponential bound on a six-site cell. 脚本用实空间和傅里叶空间两种算法核对上表，再检查四组位移，并在一个六粒子周期上检验指数势的能量界。[Verification code / 验证脚本](code/verify_chains.py)

## Relation to the necessary condition (Prop. 2.9) / 与必要条件（命题 2.9）的关系

[BKN] also state a necessary condition for high-density crystallization — Proposition 2.9 in [arXiv:1804.05743v3](https://arxiv.org/abs/1804.05743v3), not cited in the manuscript. For the sign pattern of Conjecture 2.7, $\mathcal{F}=(-f,-f,f)$, its left-hand side is $\widehat f(k)+\tfrac12\bigl(-\widehat f(k)-\widehat f(k)\bigr)=0$ for every $k$: the condition holds as an equality and carries no information for this class. The counterexample below neither contradicts Proposition 2.9 nor could have been detected through it. Details and caveats: [notes/prop-2.9-degeneracy.md](notes/prop-2.9-degeneracy.md).
中文：原文另有一条高密度结晶的必要条件（[arXiv:1804.05743v3](https://arxiv.org/abs/1804.05743v3) 命题 2.9，手稿未引用）。在猜想 2.7 所用的符号型 $\mathcal{F}=(-f,-f,f)$ 下，该式左端恒等于零，因而在此类势上不提供信息；反例既与它不矛盾，也不是它能筛出来的。详见 [notes/prop-2.9-degeneracy.md](notes/prop-2.9-degeneracy.md)。

Questions and corrections are welcome through this repository’s Issues tab. 欢迎在本仓库的 Issues 中提出问题或勘误。
