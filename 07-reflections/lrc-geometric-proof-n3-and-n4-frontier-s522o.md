---
source: oracle-2026-06-01-S522o
status: result (geometric proof of LRC n=3; n=4 frontier mapped)
tags: [lonely-runner, permutohedron, geodesic, n3-proof, center-grid, view-obstruction]
---

# Proving Small LRC Cases with the Permutohedron/Geodesic Methodology

Using ONLY the S521o framing (runners = a line on the permutohedral torus;
loneliness = the line meets the box B=[1/n,1-1/n]^{n-1}; per-cell witnesses):
**n=3 is proven outright**, and n=4 is reduced to a sharp residual that exposes
why n grows hard.

## n=3: a clean geometric proof (the center grid)

THEOREM. For distinct nonzero integers v1,v2 there is t with ||v1 t||,||v2 t||>=1/3.

PROOF. By scaling invariance (Lonely is invariant under v->cv, t->t/c) assume
gcd(v1,v2)=1; relabel so 0<v1<v2.
- If v2>=3: the v2 times t_k=(2k+1)/(2 v2), k=0..v2-1, have ||v2 t_k||=1/2>=1/3
  (runner 2 at its far-band CENTER). The positions {v1 t_k mod 1} step by
  v1 t_{k+1}-v1 t_k = v1/v2, so they are an arithmetic progression of common
  difference v1/v2; since gcd(v1,v2)=1 they form a full coset of the spacing-1/v2
  grid: v2 EQUALLY-SPACED points, spacing 1/v2. The far-band [1/3,2/3] has length
  1/3 >= 1/v2 (as v2>=3), so a closed interval of length >= the spacing contains a
  grid point. Hence some t_k has ||v1 t_k||>=1/3, and is lonely.
- If v2<=2: then {v1,v2}={1,2}, and t=1/3 gives ||1/3||=||2/3||=1/3. []

VERIFIED: 1101/1101 coprime pairs v1<v2<=60 (the single "t=1/3 fallback" is {1,2}).
Geometrically: the line's positions at the largest runner's far-centers form a
1-D equally-spaced grid that is FORCED to meet the 1-D far-band once spacing
(1/v_max) <= band length (1/n). This is the methodology proving the first
nontrivial case.

## n=4: the frontier (1-D grid vs 2-D box)

For 1336 primitive triples (speeds <= 22):
- single center-grid witness t_k=(2k+1)/(2 v_max): 1275/1336 (95%);
- center-grid OR n-gon t=a/4: 1299/1336 (97%);
- residual 37 (all lonely, 0 counterexamples) are exactly the family {1, k, k+1}
  (a consecutive pair plus the unit) -- the resonant near-extremal triples.

Why it does not close by one witness: fixing the largest runner far (the center
trick) leaves the OTHER TWO runners to be far simultaneously -- a 1-D grid of
center-positions must now meet a 2-D far-box [1/4,3/4]^2, and a 1-D set need not
hit a 2-D region. Exactly here the per-cell lonely-interval step (S521o) is
required, and it bites precisely on {1,k,k+1}.

## The pattern: codimension grows with n

The center trick reduces n runners to "n-2 others far on a 1-D grid." For n=3 that
is a 0-codimension hit (grid vs band, automatic). For n=4 it is a 1-D grid vs 2-D
box (codimension 1 -- not automatic). Each added runner raises the grid-vs-box
codimension by one. So:

> The methodology proves n=3 cleanly and uniformly; from n=4 on, a single
> equally-spaced witness grid undershoots the box dimension, and one must iterate
> (per-cell lonely intervals / a multi-scale grid). The resonant {1,...} families
> are the geometric carriers of the residual difficulty.

This both validates the geometry (it settles the first case from scratch) and
locates its limit: the open problem at n=14 is the high-codimension version of the
n=4 "1-D grid vs 2-D box" gap, on the resonant near-extremal families.

## Artifacts
04-computation/lrc_geometric_small_n_s522o.py
05-knowledge/results/lrc_geometric_small_n_s522o.out
