---
id: THM-2836
title: "Three-slot SFC(1) restriction through support 12 and window 8 (certified census)"
status: >
  FINITE-EXACT with per-cell rigorous certificates + independent-engine
  control.  For every 3-term f = a z^p + b z^q + c z^r with
  0 <= p < q < r <= 12 and every k <= 8: the three consecutive factorial
  moments L(f^{k+1}), L(f^{k+2}), L(f^{k+3}) (L(z^n)=n!) have no common
  nonzero zero — 2574/2574 cells certified empty by full Macaulay graded
  rank (one-way, escalation never needed), extending the independently
  controlled 840-cell three-prime core.  Extends THM-2812
  (consecutive supports, first window) in both open directions inside the
  box.  This is a three-slot restriction of ambient SFC(1), not SFC(3);
  MISTAKE-350 corrects the old slot-index notation.  The unbounded
  three-slot SFC(1) restriction remains OPEN.
source: mac-mini-2026-07-28-S171 (external open-problem raid; Epoch
  FrontierMath factorial-conjecture cluster = the repo GMC lane)
depends_on: []
related:
  - THM-2812-consecutive-three-slot-factorial-detection (slug approximate; see file)
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-2173-sparse-projective-factorial-moment-floor
script: 04-computation/sfc3_moment_window_census_macmini_S171.py
output: 05-knowledge/results/sfc3_moment_window_census_macmini_S171.out
script_sha256: bce6d58130a98c310c2825981736f5489ea93a950fae333547663f9974ab6549
output_sha256: e24ee09943713780ac294ce5e334cdc0db9066d85b9335451cd66946e8946f51
hash_basis: LF-normalized bytes
---

# THM-2836 — certified three-slot SFC(1) census: arbitrary supports, shifted windows

## Statement

Let `L : C[z] -> C`, `L(z^n) = n!`.  For all integers `0 <= p < q < r <= 12`
and `0 <= k <= 8`, the only `(a,b,c) in C^3` with

    L(f^{k+1}) = L(f^{k+2}) = L(f^{k+3}) = 0,   f = a z^p + b z^q + c z^r,

is `a=b=c=0`.  Boundary cases with one or two nonzero coefficients are
excluded by the inherited two-slot univariate result and by nonvanishing of
pure-power moments.  Equivalently, the three-slot restriction of ambient
`SFC(1)` holds on this box.  It is not a result about ambient `SFC(3)`.

## Certificate

For fixed cell, `F_m(a,b,c) = sum_{i+j+l=m} m!/(i!j!l!) (ip+jq+lr)!
a^i b^j c^l` is homogeneous of degree m.  With `D = 3k+4`, the census
verifies that the multiples `{x^alpha F_m : |alpha| = D - m}` span the whole
degree-`D` graded piece modulo a 30-bit prime.  Full rank mod p lifts to Q,
so `(F_{k+1},F_{k+2},F_{k+3})` contains `m^{>= D}` and the projective
variety is empty over `C`.  This is a one-way certificate: a deficient cell
would only have been *flagged* (none were).  Since the two-slot univariate
subproblem is proved in the inherited lane and pure-power moments are
nonzero factorials, any
common projective zero would necessarily have `abc != 0`, i.e. be a genuine
three-slot `SFC(1)` counterexample; emptiness therefore decides the cell
completely.

## Controls

* Consecutive-support cells (p,p+1,p+2, k=0) reproduce THM-2812's theorem.
* Independent engine: sympy Groebner over Q with saturation `t*abc = 1` on
  the shifted non-consecutive cell `(0,2,5,1)` returns basis `[1]`.
* Smoke boxes (r<=4, k<=2) and (r<=5, k<=3): all empty.

## Boundary / loss ledger

* Finite box only: exponents `r>=13` and windows `k>=9` are open here
  (matrix sizes grow as C(3k+6,2); the method extends mechanically).
* One-way certificates: emptiness is proved; a hypothetical deficient cell
  would need exact elimination (never triggered in this box).
* No claim about ambient `SFC(m)` for `m>=2`, the multivariate factorial
  conjecture, or the Furter rigidity transfer -- those remain exactly as
  before.

## Extension (same session)

The identical certified pipeline was rerun on the larger box
`0 <= p < q < r <= 12`, `k <= 8`: **2574/2574 cells certified empty**
(no flags), recorded in
`05-knowledge/results/sfc3_moment_window_census_ext_macmini_S171.out`
(sha256 c048322f9900ff6fbbe467c1b1ef990966ce2fddb6d374d38729e0031c8046ae).  The theorem statement therefore holds with `9 -> 12` and
`6 -> 8`.
