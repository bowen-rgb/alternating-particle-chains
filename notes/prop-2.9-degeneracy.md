# Note: BKN's high-density necessary condition (Prop. 2.9) is identically degenerate for ℱ = (−f, −f, f)

**Scope.** This is a positioning comment on top of the counterexample in this repository. It does not
question any numerical or Lean-checked result. It records one check that, as far as I can tell, is not
made in `en/paper.tex` / `zh/paper.tex`: the necessary condition for high-density crystallization
stated in [BKN] as Proposition 2.9 is satisfied *with equality* by every potential in the class of
Conjecture 2.7, and therefore gives no information there.

## What BKN state

In arXiv:1804.05743v3, §2 (journal version: *J. Stat. Phys.* **181**), Proposition 2.9 reads:

> **Proposition 2.9** (Necessary condition for high density crystallization). Suppose that
> f_αβ ∈ C²(ℝ) ∩ L¹(ℝ) are mirror symmetric and strongly tempered, i.e. |f_αβ(x)| < C|x|^(−1−η) for
> |x| > r₀. If the equidistant configuration is a minimizer of E_ℱ at high density, then
>
>   f̂₁₂(k) + ½ ( f̂₁₁(k) + f̂₂₂(k) ) ≥ 0   for all k ∈ ℝ,
>
> where f̂(k) = (1/√(2π)) ∫_ℝ f(x) cos(kx) dx.

## The check

Conjecture 2.7 fixes the sign pattern ℱ = (−f, −f, f), i.e. (f₁₁, f₂₂, f₁₂) = (−f, −f, f). Substituting:

  f̂₁₂(k) + ½ ( f̂₁₁(k) + f̂₂₂(k) ) = f̂(k) + ½ ( −f̂(k) − f̂(k) ) = 0.

The quantity is **identically zero — for every f, not only for the completely monotone ones** — so the
inequality holds as an equality for all k. Consequently Proposition 2.9 cannot exclude any
counterexample with this sign pattern: in the terminology used in this repository, the leading-order
(LO) test is degenerate, and the N = L = 4 configuration could not have been detected through it.

## Why this may deserve one sentence in the paper

The manuscript cites Conjecture 2.7 and Theorem 2.4 of [BKN] but not Proposition 2.9. A reader who
knows Proposition 2.9 might reasonably either (a) assume the necessary condition rules out the four-site
configuration, or (b) think the counterexample contradicts the proposition. Neither is the case: the
proposition is a one-sided *high-density* implication whose functional vanishes identically for
ℱ = (−f, −f, f), so the implication imposes no constraint on this class. Stating this explicitly would
make the "why was this missed" part of the story complete.

## Checks performed

- **Reference.** Text above quoted from arXiv:1804.05743v3, §2 (HTML rendering). The numbering in the
  *J. Stat. Phys.* version was **not** verified — if it differs, the citation should follow the
  journal version.
- **Numbers.** `code/verify_chains.py` re-derived independently (no Python on the checking machine, so
  a small Node script was used): equidistant energy per particle 0.999633292158, four-site
  0.924751872822, difference 0.074881419336, and q = e^(−π²/(16a)) ≈ 0.084805 with 32q⁴/(1−q⁴) < 4q.
  All consistent with the values stated in this repository; four additional displacements
  (s = 0.01, −0.125, 0.5, 0.9) also fall below the equidistant value.
- **Formalization.** The Lean development covers the two proofs over the full infinite chain; the
  local-instability and Riesz parts are, per the README, not formalized. That is consistent with the
  files present and is outside the scope of this note.

## Authorship

Drafted by 照澜 (an AI agent) on a local clone, at the suggestion of Lan Zhao. Please treat it as a
suggestion rather than an authoritative claim; corrections are welcome.
