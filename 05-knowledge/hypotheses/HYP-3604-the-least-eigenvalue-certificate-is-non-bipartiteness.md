---
id: HYP-3604
title: The LEAST-EIGENVALUE CERTIFICATE for the LRC apex floor -- the binding value 4cos^2(3pi/7) is lambda_min(2I+A(C_7)), the least eigenvalue of the doublet's autocorrelation Gram (= twice-identity plus the 7-cycle adjacency); and its POSITIVITY is exactly the NON-BIPARTITENESS of C_p (the apex cycle has odd length p, so A(C_p) has no eigenvalue -2, so 2I+A(C_p) is positive definite with lambda_min = 2-2cos(pi/p) = 4sin^2(pi/2p) > 0). Over all 127 non-full Z_p cores the least nonzero-mode eigenvalue is >= this doublet value (THM-590). The certificate is SET-INDEPENDENT (depends only on the apex prime p), is the Bochner/SOS minorant (its worst eigenvector is the middle Fourier mode k=(p+-1)/2), and certifies the sigma-ODD / discrete / EXISTENCE content (klein-S18), not the measure. So: the LRC floor is positive BECAUSE the apex cycle is odd
status: VERIFIED (exact: lambda_min(2I+A(C_p))=4sin^2(pi/2p) for odd p, =0 for even p; 127-core census). The certificate FORM + the non-bipartiteness mechanism are the contribution; the value/positivity is klein's THM-590. Finite & explicit (Lean-able: forall O subsetneq Z_p, lambda_min^nonzero(Gram O) >= 4sin^2(pi/2p)).
source: mac-mini-2026-06-30-S37
related:
  - THM-590   # klein: the apex gap 4cos^2(3pi/7)>0 (this is its least-eigenvalue/PSD certificate form)
  - HYP-3590  # mac-mini S31: 4cos^2(3pi/7)=2+lambda_min(A(C_7)); the doublet = the 7-cycle
  - HYP-3599  # klein-S18: the certificate is the sigma-odd/discrete side; measure-rho_j binds at the cusp
  - HYP-3602  # mac-mini S34: intransitivity; non-bipartite = has an odd cycle = the danger relation doesn't factor
  - HYP-3535  # klein: the Fejer-Bochner cyclotomic SOS (= this minorant)
results:
  - 04-computation/least_eigenvalue_certificate_macmini_20260630.py
  - 05-knowledge/results/least_eigenvalue_certificate_macmini_20260630.out
---

# HYP-3604 -- the least-eigenvalue certificate is non-bipartiteness

Working the least-eigenvalue certificate (the owner's ask). The apex floor is a positive-definiteness
statement about one explicit matrix, and its positivity has a one-line reason: the apex cycle is odd.

## The certificate
For an apex core `O subset Z_p` (p the apex prime), the autocorrelation **Gram** is the circulant with
eigenvalues `|sum_{x in O} w^{kx}|^2` (k=0..p-1). The mean mode (k=0) is `|O|^2`; the floor is the least
NONZERO-mode eigenvalue, `gap(O) = min_{k!=0}|sum w^{kx}|^2`. For the binding **doublet** `O={a,b}` this Gram
is exactly `2I + A(C_p)`, where `C_p = Cay(Z_p, {+-(b-a)})` is the **p-cycle**, and
> `lambda_min(2I + A(C_p)) = 2 + lambda_min(A(C_p)) = 2 - 2cos(pi/p) = 4 sin^2(pi/2p)`.
For `p=7`: `= 4sin^2(pi/14) = 0.19806 = 4cos^2(3pi/7)` (HYP-3590, THM-590).

## The positivity IS non-bipartiteness (the mechanism)
The eigenvalues of `A(C_p)` are `2cos(2pi k/p)`, all in `(-2, 2]`. The bottom value `-2` is ATTAINED iff
`k=p/2` exists iff **`p` is even iff `C_p` is bipartite**. So:
> `2I + A(C_p) succ 0  <=>  lambda_min(A(C_p)) > -2  <=>  C_p has no eigenvalue -2  <=>  C_p NON-BIPARTITE
> <=>  p ODD.`
- `p` ODD: `lambda_min(A) = -2cos(pi/p) > -2`, so the Gram is positive definite, floor `= 4sin^2(pi/2p)>0`.
- `p` EVEN: `lambda_min(A) = -2`, so the Gram is SINGULAR, floor `= 0` -- the disproof boundary.
**The least-eigenvalue certificate is the odd-cycle / non-bipartiteness certificate.** An odd cycle is the
obstruction to a 2-coloring; spectrally, that is `lambda_min(A) > -2`; analytically, that is `2I+A succ 0`;
in the LRC, that is the floor `> 0`. The apex prime `7` is odd, so `C_7` is non-bipartite, so the floor is
positive -- the same odd-cycle fact (HYP-3602: the danger relation is intransitive = has an odd cycle =
does not factor) read spectrally.

## The worst eigenvector (the Bochner/SOS minorant)
The least eigenvalue sits at the **middle Fourier mode** `k = (p+-1)/2` (the near-`pi/2` frequency), with
eigenvector `v_j = cos(2pi k j / p)`. This is the most-oscillatory resonance, the one that for odd `p` never
hits the `-1` phase exactly. The least-eigenvalue bound is precisely the **Fejer-Bochner SOS minorant**
(klein HYP-3535/3581): the certificate is `Gram - 4sin^2(pi/2p) I succeq 0`, an SOS/PSD witness, and the
worst mode is where it is tight.

## Uniform, set-independent, finite
- **Uniform** (THM-590): over all 127 non-full `Z_7` cores the gap values are `{0.198, 0.308, 1, 2}`; the
  least nonzero-mode eigenvalue is `0.198` at the doublets. So `inf_{O subsetneq Z_p} gap(O) = 4sin^2(pi/2p)`.
- **Set-INDEPENDENT**: the bound depends ONLY on `p` (the apex prime) -- the covering set never enters. This
  is the set-independent floor the gatekeeper (THM-579/OPEN-Q-108) needs, in certificate form, NOT a per-set
  estimate (which klein-S4 showed is unbounded).
- **FINITE / Lean-able**: `forall O subsetneq Z_p, lambda_min^{nonzero}(Gram(O)) >= 4sin^2(pi/2p)` is a
  finite (127-case at p=7) explicit cyclotomic check.

## Placement (klein-S18) and what it does NOT do
This certifies the **sigma-ODD / discrete / Bochner** content -- the apex Gram gap, binding at the doublet:
the odd cycle is present and non-degenerate. Per klein-S18 (HYP-3599) that is the EXISTENCE-supporting side
(the discrete skeleton), and it is what the proof needs once klein-S16 retired the measure (inf=0). It does
NOT bound the measure-`rho_j` (which binds at the opposite, cusp, end). So the least-eigenvalue certificate
discharges the DISCRETE floor (the odd cycle is there, spectrally non-degenerate); the remaining open piece
is the measure->existence passage at the top level (klein-S18's `rho_0 > 0`).

## What it buys
THM-590 gets a crisp certificate FORM (a PSD/Bochner statement about `2I+A(C_p)`) and its deepest REASON
(non-bipartiteness = odd length). The floor's positivity is identified with the apex cycle being odd -- one
line, set-independent, finite, Lean-able. It unifies the spectral (least eigenvalue), combinatorial (odd
cycle / non-bipartite), and analytic (PSD/SOS) faces of the same fact.
