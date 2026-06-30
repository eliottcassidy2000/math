---
id: HYP-3604
title: The LEAST-EIGENVALUE CERTIFICATE for the sigma-EVEN LRC floor -- THM-590's g(O) is EXACTLY the smallest eigenvalue of the core's autocorrelation circulant C(O) (Bochner-PSD = SOS = sigma-even; verified 0 mismatches over all 127 cores); the BINDING doublet's Gram is C=2I+A(C_7)= the SIGNLESS LAPLACIAN Q(C_7) of the apex 7-cycle, and lambda_min(Q(C_p)) = 2-2cos(pi/p) is POSITIVE iff p is ODD (C_p non-bipartite); so 4cos^2(3pi/7)=2-2cos(pi/7)=0.198062 > 0 PRECISELY because the apex prime 7 is odd (an even apex gives 0). 14=2*7 = (even 2-adic descent) x (odd apex cycle's positive signless-Laplacian gap)
status: VERIFIED (g(O)=lambda_min(C(O)) for all 127 cores, 0 mismatches; doublet Gram = Q(C_7); lambda_min(Q(C_p))=2-2cos(pi/p), =0 iff p even, verified p=3..14). Integrates mac-mini-S36 (HYP-3601: doublet autocorrelation = 2I+A(C_7)) + THM-590 + the sigma-even frame (mac-mini-S6).
source: klein-2026-06-29-S19
depends_on:
  - THM-590   # g(O) = the apex cyclotomic gap (here = lambda_min of the autocorr circulant)
related:
  - HYP-3601   # mac-mini S36: the binding doublet autocorrelation = 2I+A(C_7); Q1/Q2
  - HYP-3594   # mac-mini S34: the truth is the ODD cycle (here = the odd C_7 that makes lambda_min>0)
  - THM-588    # algebraic connectivity / graph-spectral least eigenvalue
  - HYP-3599   # the bridge (the floor is sigma-even; this is its certificate)
  - HYP-3606   # mac-mini-S37: CONVERGENT independent rediscovery (non-bipartiteness form); cites this
  - HYP-3609   # klein-S20: chips the remaining gap (deeper levels closed, open part = level 0)
results:
  - 04-computation/least_eigenvalue_certificate_klein.py
  - 05-knowledge/results/least_eigenvalue_certificate_klein.out
  - 05-knowledge/results/signless_laplacian_odd_cycle_klein.out
---

# HYP-3604 — the least-eigenvalue certificate is the signless Laplacian of the apex odd cycle

## The certificate is a least eigenvalue (verified)

The LRC floor is entirely sigma-EVEN (mac-mini-S6: "lonely is even, Redei is odd"; the sigma-odd witness
does NOT apply -- a lonely tournament is not self-converse). Its certificate is therefore a Bochner-positive
(SOS = sum of squares = sigma-even) object: the cyclotomic Gram. Concretely, for a core `O subset Z_7` form
the **autocorrelation circulant** `C(O)_{ij} = a((i-j) mod 7)`, `a(d) = #{x in O: (x+d) in O}`. Then:
> `C(O)` is real-symmetric PSD (eigenvalues `|Ohat(k)|^2 >= 0` -- Bochner), and THM-590's
> `g(O) = min_{k!=0}|Ohat(k)|^2` is EXACTLY the smallest eigenvalue of `C(O)`.
Verified: `lambda_min(C(O)) = g(O)` for all 127 proper cores (0 mismatches); the five values
`{0, 0.198062, 0.307979, 1, 2}` are THM-590's. The least eigenvector is the explicit SOS/dual certificate
(for the doublet, `Cv = lambda v` to 1e-16).

## The binding atom is the signless Laplacian of the apex ODD cycle

The binding core is the **doublet** `O={0,1}`, with autocorrelation `a = [2,1,0,0,0,0,1]`, i.e.
> `C({0,1}) = 2I + A(C_7) = Q(C_7)` -- the **SIGNLESS LAPLACIAN** `Q=D+A` of the 7-cycle
(integrating mac-mini-S36/HYP-3601: the doublet autocorrelation IS `2I+A(C_7)`; `C_7` is an even graph, the
apex of the even-graph dual). Its eigenvalues are `2 + 2cos(2pi k/7)`, so
> `lambda_min(Q(C_7)) = 2 + 2cos(6pi/7) = 2 - 2cos(pi/7) = 4cos^2(3pi/7) = 0.198062`
(and `sin(pi/14) = cos(3pi/7)`, so this is also `4sin^2(pi/14)`). The binding mode is `k=3` (the sigma-orbit
`{3,4}={3,-3}`): `|1 + zeta^3|^2 = 2 + 2cos(6pi/7)` -- which is why the angle is `3pi/7`.

## Positive PRECISELY because the apex prime is ODD

`lambda_min` of a signless Laplacian is `0` iff the graph is **bipartite**. `C_p` is bipartite iff `p` is
EVEN. So:
> `g(doublet on Z_p) = lambda_min(Q(C_p)) = 2 - 2cos(pi/p) > 0  <=>  p is ODD` (`C_p` non-bipartite).
Verified p=3..14: positive `1, 0.382, 0.198, 0.121, 0.081, 0.058` at odd `p=3,5,7,9,11,13`; exactly `0` at
even `p=4,6,8,10,12,14`. **The LRC(14) apex obstruction `4cos^2(3pi/7)` is positive precisely because the
apex prime `7` is ODD** -- the odd cycle `C_7` is non-bipartite. An even apex prime would give `0` (a
degenerate, bipartite cusp). This is exactly mac-mini-S34's "the truth is the odd cycle" (HYP-3594), made
quantitative: the odd cycle `C_7` is the graph whose signless-Laplacian least eigenvalue IS the positive
apex atom.

## The 2*7 split

`14 = 2 * 7 = (the even 2-adic descent) x (the odd apex cycle's positive signless-Laplacian gap)`:
- the **2-part** is the sigma = `t -> t+1/2` half-translation / 2-adic descent (THM-580, the even-category
  degree) -- it peels the even structure down to the apex prime;
- the **7-part** is the cyclotomic least eigenvalue `lambda_min(Q(C_7)) > 0` (this) -- positive because 7
  is odd.
The certificate lives at the APEX (mod 7), NOT the full grid: the full lonely-measure power spectrum has
zeros at most frequencies (verified `min |Lhat(k)|^2 ~ 0`), so the naive full least eigenvalue is `~0`. The
"danger does not factor" means the certificate is the apex cyclotomic gap, reached after the descent
removes the 2-part -- not a product over the whole circle.

## Net

The least-eigenvalue certificate for the sigma-even LRC floor is the smallest eigenvalue of the core
autocorrelation circulant (Bochner-PSD = SOS). Its binding value is `lambda_min` of the signless Laplacian
of the apex 7-cycle, `2-2cos(pi/7) = 4cos^2(3pi/7)`, positive iff the apex prime is odd. This unifies
THM-590 (the gap), mac-mini-S36 (the `C_7` atom), mac-mini-S34 (the odd cycle), and the sigma-even frame
(mac-mini-S6): the apex obstruction is the positive signless-Laplacian gap of an ODD cycle.
