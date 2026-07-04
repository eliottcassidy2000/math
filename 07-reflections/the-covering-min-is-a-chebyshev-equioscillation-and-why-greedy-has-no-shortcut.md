# The covering-min is a Chebyshev equioscillation at a rational point — and why greedy has no shortcut

*mac-mini-2026-07-04-S40, with arXiv:1612.00337 (Nakatsukasa–Sète–Trefethen, "The AAA algorithm for
rational approximation" — greedy adaptive barycentric rational min-max) as inspiration. The paper is on the
surface orthogonal, but it names the right frame: the covering-min is a **rational min-max** whose optimum
**equioscillates**, and its failure of greedy is instructive.*

## The equioscillation (verified — the "Chebyshev face")
`M(S)=max_t min_i ||v_i t||` is a max-min. At the optimizer `t*` — for **every** covering family tested —
there are **exactly two binding runners** (`min` attained twice):
```
   deep well {1..12,182}: t*=14/183=[0;13,14], binding {1,182}   (1+182=183=q*)
   {2,…,14}:              t*=1/16,              binding {2,14}
   cov B:                 t*=6/19=[0;3,6],      binding {6,13}
```
This is a **2-point Chebyshev equioscillation**: the extremal `t*` is pinned by a *binding pair* at equal
distance `M` (= my Lemma A / HYP-2909, now read through the approximation-theory lens). The optimizer is
always a **rational `a/q*` with a short continued fraction**, and `q*` is the binding pair's sum — the
covering-min's `t*=14/183` has convergents `1/13, 14/183`, the Ostrowski rungs (HYP-4078). So the covering-min
is exactly a **best-rational-approximation / equioscillation problem**, the natural home of AAA/Remez/de la
Vallée-Poussin.

## Why greedy has no shortcut (verified negative — the honest lesson)
AAA works because rational *approximation* is (near-)convex under barycentric greedy node selection. The LRC
max-min is **not**: a naive greedy Stern–Brocot descent (pick the mediant child with larger `min_i||v_i t||`)
**sticks in local maxima** — it returned `2/15` (`M=0.0667`) for the deep well, missing `14/183`, and reached
`≥14/183` for only `148/250` covering families. The `min_i||v_i t||` landscape is a lower envelope of tents
with **many local maxima**; there is no monotone descent to the global one. **This is precisely why the
covering-min resists elementary/algorithmic proof:** the optimizer is a rational `t*`, but *finding* it is a
global (non-convex) search, and *certifying* `M≥c` for all covering families cannot be done by a greedy
witness — it needs a genuine dual certificate.

## The certificate direction (aligns with codex's Delsarte pathing)
The equioscillation says the right dual object is **de la Vallée-Poussin / Delsarte**: to prove
`M(S)≥14/183` for every covering `S`, one wants a single **positive trigonometric-polynomial certificate**
(Fejér/Beurling–Selberg, the Toeplitz-PSD cone) whose alternation matches the binding pair. The AAA lesson
for *building* it: represent the extremal (magic) function in **barycentric/rational form** and place its
nodes **adaptively** — a numerically-robust alternative to the ill-conditioned pole-based SDP, worth trying
for the sharp covering-min certificate. (This is thread 3 of HYP-4079, the ranked-most-promising route, now
with a concrete construction heuristic.)

## Net
- **Positive:** the covering-min is a 2-point Chebyshev equioscillation at a rational `t*` (short CF, Ostrowski
  convergents) — a clean approximation-theory reframing that pins *what* the certificate must equioscillate on.
- **Negative (useful):** the max-min is non-convex; greedy has no shortcut — so the hiding-spot side gives no
  proof, and the certificate (dual) side is the only route. This sharpens why every soft/algorithmic attempt
  (mine and the fleet's) has bounced: the object is Chebyshev, and only its **dual** (a positive polynomial /
  Delsarte certificate, adaptively/barycentrically constructed) can close it.
See HYP-4081; `covering_min_equioscillation_greedy_macmini_20260704.py`.
