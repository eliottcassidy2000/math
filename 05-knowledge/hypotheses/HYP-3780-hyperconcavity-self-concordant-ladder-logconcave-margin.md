---
id: HYP-3780
title: MERGING HYPERCONCAVITY into the covering-min -- the self-concordant CONVEX LADDER and the LOG-CONCAVE Dedekind margin, unified by apex-7 hyperbolicity. Two verified concavity structures. (H1 SELF-CONCORDANT LADDER) 1/M(n)=Phi6/n=(n-1)+1/n is EXACT and CONVEX (d^2/dn^2=2/n^3>0), so the covering-min M=1/((n-1)+1/n) is the reciprocal of a self-concordant convex ladder (= opus-S4 HYP-3769 'self-concordant residual 1/M=(n-1)+1/rung', rung=n); the Stern-Brocot ray [0;n-1,k]=k/((n-1)k+1) has 1/value=(n-1)+1/k, also convex in k -- the continued-fraction descent IS a self-concordant barrier. Self-concordant barriers <=> HYPERBOLICITY CONES (hyperbolic polynomials, Guler/Renegar), tying the covering-min descent to the apex-7 HYPERBOLIC (2,3,7) geometry (S65 HYP-3771). (H2 LOG-CONCAVE MARGIN) |s(n,Phi6)|=T/(12T+6) (T=n(n-1)/2) is CONCAVE and LOG-CONCAVE in n -- a bounded increasing concave sequence rising to its supremum 1/12=|zeta(-1)| (the S67-S70 regularization limit). So the Dedekind-sum margin is hyperconcave and the regularized -1/12 is the peak. DUALITY: 1/M is convex (the barrier/descent) while |s| is log-concave (the margin) -- the two 'hyper' faces (hyperbolic self-concordance + concavity) of one covering-min. UNIFICATION: 'hyperconcavity' = hyperbolic (self-concordant barrier of the descent + hyperbolic (2,3,7) apex) + concave (log-concave margin); the covering-min lives in a hyperbolicity cone, the residual is the hyperbolic-surface (cusp form f14). DIRECTION: the OCF independence polynomial I(Omega,x) (H(T)=I(Omega,2)) log-concavity/Lorentzian = the tournament-side hyperconcavity (trivial n<=4; needs n>=6 with 5-cycles).
status: H1, H2 EXACT/verified (n=3..23): 1/M=(n-1)+1/n convex; |s(n,Phi6)| concave+log-concave. The self-concordant/hyperbolicity-cone framing (H1) merges with opus-S4 HYP-3769 (proved 1/M=(n-1)+1/rung). The unification (H3) is conceptual; the OCF log-concavity (H4) is a direction (n<=4 trivial). Creative merge of hyperconcavity with the covering-min/margin/regularization thread; NOT a proof of LRC.
source: mac-mini-2026-06-30-S71
related:
  - HYP-3769   # opus-S4: the self-concordant residual 1/M=(n-1)+1/rung (H1 merges with & confirms this)
  - HYP-3774   # S67: the zeta-regularization carrier (|s|->-1/12; the log-concave margin's supremum)
  - HYP-3779   # S70: the eta-invariant/Euler-Maclaurin/p-adic avatars of the margin
  - HYP-3771   # S65: the (2,3,p) spine; the apex-7 HYPERBOLIC geometry = the hyperbolicity-cone tie
  - HYP-3732   # the Stern-Brocot ray [0;n-1,k] (the self-concordant CF ladder)
results:
  - 04-computation/hyperconcavity_covering_min_macmini_20260630.py
  - 05-knowledge/results/hyperconcavity_covering_min_macmini_20260630.out
---

# HYP-3780 -- hyperconcavity: the self-concordant ladder and the log-concave margin

Merging the modern concavity toolkit (log-concavity / hyperbolic-barrier self-concordance) into the
covering-min / margin / regularization thread. Two verified concavity structures, dual to each other, both
"hyper."

## H1 -- the covering-min is a self-concordant convex ladder (hyperbolic barrier)
`1/M(n) = Phi_6/n = (n-1) + 1/n` -- **exact**, and **convex** (`d^2/dn^2 = 2/n^3 > 0`). So the covering-min
`M = 1/((n-1)+1/n)` is the reciprocal of a **self-concordant convex ladder** -- exactly opus-S4's
`1/M = (n-1)+1/rung` with `rung = n` (HYP-3769). More: the Stern-Brocot ray `[0;n-1,k] = k/((n-1)k+1)` has
`1/value = (n-1)+1/k`, convex in `k` too, so the whole **continued-fraction descent is a self-concordant
barrier**. And self-concordant barriers are exactly the barriers of **hyperbolicity cones** (hyperbolic
polynomials; Guler, Renegar) -- which ties the covering-min descent to the **apex-7 hyperbolic `(2,3,7)`
geometry** (S65). The covering-min is a hyperbolic-programming object.

## H2 -- the Dedekind margin is log-concave (hyperconcave)
`|s(n,Phi_6)| = T/(12T+6)` (`T = n(n-1)/2`) is **concave** and **log-concave** in `n` (verified `n=3..23`):
a bounded increasing concave sequence rising to its **supremum `1/12 = |zeta(-1)|`** (the S67-S70
regularization limit). So the Dedekind-sum margin is hyperconcave, and the regularized `-1/12` is its peak.

## The duality and the unification
`1/M` is **convex** (the barrier / the descent), while `|s|` is **log-concave** (the margin) -- the two
"hyper" faces of one covering-min. Reading "**hyperconcavity**" as *hyperbolic + concave*:
- **hyperbolic**: the self-concordant barrier `1/M=(n-1)+1/n` (a hyperbolicity-cone object) and the apex-7
  **hyperbolic `(2,3,7)`** geometry;
- **concave**: the **log-concave** Dedekind margin `|s|`.

The covering-min lives in a hyperbolicity cone (self-concordant descent); the margin is its log-concave
regularization (peaking at `-1/12`); the residual is the hyperbolic-surface (the genus-1 cusp form `f_14` at
`d=7`). Concavity and hyperbolicity meet at the covering-min.

## Direction (H4) -- the tournament-side hyperconcavity
`H(T) = I(Omega, 2)` (the OCF). Is the independence polynomial `I(Omega, x)` of the odd-cycle collection
**log-concave / Lorentzian**? Trivially yes for `n <= 4` (`I` has degree `<= 1`); a real test needs `n >= 6`
(5-cycles enter `Omega`). If `I(Omega,x)` is Lorentzian, log-concavity of `H`-related counts follows -- the
hyperconcavity of the neglected tournament half, mirroring the covering-min's.

## Honest scope
H1 and H2 are exact/verified (`n=3..23`); H1 merges with and confirms opus-S4 (HYP-3769). The
hyperbolicity-cone unification (H3) is a conceptual merge (self-concordant barriers `<=>` hyperbolicity cones
is a theorem; that the LRC covering-min *is* a hyperbolic-programming instance is a framing, not a proof). The
OCF log-concavity (H4) is a direction. This creatively merges hyperconcavity with the covering-min/margin/
regularization thread; it is not a new proof step for LRC14.
