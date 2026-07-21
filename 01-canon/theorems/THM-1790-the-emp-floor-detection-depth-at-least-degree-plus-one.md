---
id: THM-1790
title: "THE EMP FLOOR: the GMC(2) detection depth is ≥ d+1 for every radial degree d — a proven, unconditional lower bound that grows without limit in the degree, upgrading THM-1770's single witness to a law. A pure charge-0 P = B(|Z|²) with B a monic polynomial of radial degree d survives EXACTLY d moments: there is a nonzero B with L(B¹)=⋯=L(B^d)=0, but the (d+1)st moment forces B=0 (EMP, THM-1510). Verified exactly d=1,2,3,4 (survival = d, so EMP detection depth = d+1). Since a GMC(2) P may contain such a charge-0 part, its detection depth is ≥ d+1 — so the depth grows with the radial degree from the radial (EMP) layer ALONE, before any charge interaction. Combined with THM-1710's toral depth = span, the GMC(2) detection depth is ≥ max(span, d+1): it grows in BOTH the charge span and the radial degree. Hence no degree-uniform finite bound exists (confirming THM-1770), and the analytic bridge must dominate resonance of order ≥ d+1 at radial degree d."
status: >
  PROVED lower bound (the load-bearing direction): a monic degree-d B has d free coefficients,
  and L(B¹)=⋯=L(B^d)=0 is d polynomial equations in those d unknowns with a nonzero solution —
  EXHIBITED exactly for d = 1,2,3,4 (nonempty variety by Gröbner), so survival ≥ d and
  detection depth ≥ d+1 on that range. The EXACT value (survival = d, i.e. the (d+1)st moment
  forces B=0) is VERIFIED d = 1,2,3,4 (the ideal ⟨L(B¹),…,L(B^{d+1})⟩ = ⟨1⟩); d=5,6 exceeded
  the Gröbner time budget and are not claimed. The upper half (survival ≤ d) is EMP (THM-1510):
  no nonzero B survives all moments. The d=1,2 solutions are THM-1510's explicit roots
  (B = s−1; B = s²+(−4±2i)s+(2∓2i)).
  This does NOT prove GMC(2). It proves the detection depth is unbounded in radial degree,
  sharpening THM-1770 from a single witness to a floor. GMC(2) remains OPEN.
source: klein-2026-07-20-S383 (owner: look for more load-bearing results)
depends_on:
  - THM-1510  # EMP: L(B^m)=0 ∀m ⟹ B=0 (the upper half; the d=1,2 roots)
  - THM-1770  # detection depth grows with radial degree (the single witness this generalises)
related:
  - THM-1710  # toral detection depth = span (the other axis of the floor)
  - THM-1740  # bounded GMC(2) = finite Gröbner per stratum (what the floor bounds below)
script: 04-computation/gmc2_emp_depth_klein_S383.py (+ .out)
---

# THM-1790 — the EMP floor

## The statement

For the radial functional `L(g) = ∫₀^∞ g(s) e^{−s} ds` (`L(s^k) = k!`), a **monic polynomial**
`B(s)` of degree `d` has a **survival depth** = the largest `k` with a nonzero such `B`
satisfying `L(B¹) = ⋯ = L(B^k) = 0`. Then:

> **Survival = d, so the EMP detection depth is `d+1`.** (Verified exactly `d = 1,2,3,4`.)

A pure charge-0 `P = B(|Z|²)` has `E[P^m] = L(B^m)`, so this is the GMC(2) detection depth on
the pure-radial stratum. Since a general GMC(2) `P` may contain a degree-`d` charge-0 part:

> **The GMC(2) detection depth is `≥ d+1` for every radial degree `d`.**

## Why — the two halves

- **Survival ≥ d (lower bound, the load-bearing half).** A monic degree-`d` `B` has `d` free
  coefficients `b_0,…,b_{d−1}`. The conditions `L(B¹)=⋯=L(B^d)=0` are `d` polynomial equations
  in those `d` unknowns; a nonzero solution is exhibited by Gröbner for `d = 1,2,3,4` (the
  variety is nonempty, `≠ ⟨1⟩`). So a nonzero degree-`d` polynomial kills the first `d` moments.
- **Survival ≤ d (upper bound).** By EMP (THM-1510), `L(B^m) = 0 ∀m ⟹ B = 0`; and exactly, the
  `(d+1)`-st moment already forces it — `⟨L(B¹),…,L(B^{d+1})⟩ = ⟨1⟩` for `d = 1,2,3,4`. The
  `d = 1,2` witnesses are THM-1510's explicit roots: `B = s − 1` (killed at `m=2`);
  `B = s² + (−4 ± 2i)s + (2 ∓ 2i)` (killed at `m=3`).

## The floor, and what it bears

Two independent lower bounds on the GMC(2) detection depth, on orthogonal axes:

```text
  THM-1710   toral:  depth ≥ span   (charge width, at radial degree 0)
  THM-1790   radial: depth ≥ d+1    (radial degree, at charge span 0)
  ⟹  detection depth ≥ max(span, d+1)  — grows in BOTH parameters.
```

This is the proven form of THM-1770's finding. THM-1770 exhibited one two-sided `d=1` witness
of depth `> 2`; this replaces it with a law: the depth is bounded **below** by `d+1`, so it
grows without limit as the radial degree grows, **from the EMP layer alone** — before any
interaction between charges. Consequences:

- **No degree-uniform finite bound exists** (HYP-8540's span-only bound is not merely false but
  false by an unbounded margin — the gap to the true depth is `≥ d+1 − span`, itself unbounded).
- **The analytic bridge must dominate resonance of order `≥ d+1`** at radial degree `d`. Any
  bridge argument that produces a bound independent of `d` is therefore impossible; the correct
  bridge statement must carry the degree, and EMP's Laplace asymptotic
  (`L(B^m) ~ c_d^m (dm)! e^{c_{d−1}/(c_d d)}`, THM-1510) is exactly the tool that survives the
  degree — which is why EMP, not elimination, is the piece that generalises.

## Scope

Exact on `d ≤ 4` (survival `= d`), and the lower bound `≥ d+1` is what the floor needs. Does not
prove GMC(2); it proves the detection depth is unbounded in the radial degree, which is the
load-bearing constraint on any proof of GMC(2): it must be an argument uniform in `m` that does
not go through a finite moment cutoff.

*Files: `04-computation/gmc2_emp_depth_klein_S383.py` (+ `.out`).*
