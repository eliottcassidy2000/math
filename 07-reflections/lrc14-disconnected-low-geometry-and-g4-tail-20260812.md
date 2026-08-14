# Disconnected-low geometry and reflected pair-floor reductions

Status: PROVED structural/analytic reductions + VERIFIED-EXACT finite
residual, subsequently completed by THM-3355. This file preserves the route
into that closure and the repaired affine sidecar.

## Structural reduction

Join two of six distinct levels when their reduced ratio `(P,Q)` satisfies
`P+Q<=7`.  If this low graph is disconnected, its components have primitive
connected shapes of sizes at most five.  The connected-shape counts from
THM-3350 are

`c_1,...,c_5 = 1,8,94,1295,19389`.

Consequently the number of unordered component-shape profiles is

`[x^6] product_(n=1)^5 (1-x^n)^(-c_n)=36520`.

This coefficient counts realizable profiles, not merely formal multisets.
For each selected normalized connected component, clear denominators and
divide its internal gcd.  Place the first component at its primitive scale;
if all coordinates already placed are at most `M`, multiply the next
component by `6M+1`.  Every new/old ratio is then strictly greater than six,
so no cross pair is low, while scaling preserves every internal low graph.
The companion reconstructs the low components of all `36,520` resulting
six-level sets with zero mismatch (semantic digest
`9ad5d29562c8aa99fd6b852c8f9786042495cc93cfc5c968df2a42947970e324`).

The cross-component graph is complete multipartite.  Choosing vertices
`x in V_1`, `y in V_2`, joining `x` to every vertex outside `V_1`, and
joining `y` to every other vertex of `V_1` gives a five-edge spanning tree.
Its homogeneous credit is at least `5/105=1/21`, exceeding the universal
debt by

`1121969414539/35319431272125`.

On a body-safe cell the first clause has exactly `p` full teeth.  The
mean-zero primitive of the `q`-clause indicator has oscillation
`6L/[49(qL-f)]`, so interval-by-interval integration gives

`I(p,q) >= pL/[49(pL-e)] - 6pL/[49(qL-f)]`.

Thus `q>=8p` implies

`I(p,q) >= 23/4655 = 1/294 + 43/27930`.

The companion `lrc14_disconnected_low_geometry_verify_20260812.py` checks
the connected-shape enumeration, all profile realizations, multipartite tree
counts, and all displayed exact gaps.

## The large-ruler floor

The cutoff `L>=4592` is now frozen rather than inherited implicitly.  There
are `1,514` feasible large-ruler contexts on twelve rulers and `180` ordered
endpoint lanes.  For `P>=4`, the linked projective floor and the same
one-variable midpoint envelope used below already close at `g=1`; the weakest
exact envelope margin is

`1179868965463/947350767068160 > 0`

at lane `(L,e,f)=(5096,14,13)`.  Every loss decreases with `P` and with the
common dilation `g`, so this one corner certifies all `P>=4,g>=1`.

There are only nineteen primitive high channels with `P<4` and `Q<8P`.
The contextwise midpoint bound closes all `28,766` channel-context rows
already at `g=1`; its weakest tail margin is
`23885167/27322254960`.  As an independent positive control, literal
THM-3352 mass evaluation on the same `28,766` rows has weakest margin
`349997/61905522`, at `(P,Q;L,j,e,f)=(3,5;5096,2772,14,4)`.
Thus every large-ruler primitive high pair has physical mass at least
`1/294` at every dilation.  Reproduce with
`lrc14_disconnected_large_ruler_floor_20260812.py`.

## Common dilation g>=4

Universe: all 2,530 ordered upper-median contexts `(L,j,e,f)` with `L<4592`; every coprime `P<Q<8P`, `P+Q>=8`, `(P,Q)!=(3,5)`; every integer dilation `g>=4`. Claim: exact reflected physical overlap is at least `1/294`.

The projective fibre satisfies both `Phi>=1/105` and `Phi>=1/49-12/(49PQ)`. The latter follows from the exact 14x14 linked tent-correction bank; its minimum is `-12/49` (attained on eight residue pairs). We use their maximum.

THM-3350 gives error

`B = gamma_P eP/(gLP-e) + gamma_Q fQ/(gLQ-f) + |C|(floor(|C|/L)+1)/(2g^2LPQ)`,

where `C=Qe-Pf` and `gamma_k=(k^2+1)/(2k^2)` is a lawful common upper bound (equal to the theorem constant for odd k and above it for even k).

For `P<=Q<=8P`, `gamma_Q<=gamma_P`, `Q/(gLQ-f)<=P/(gLP-f)`, and

`C^2/(PQ) <= max((e-f)^2,(8e-f)^2/8) <= max((e+f)^2,(8e+f)^2/8)`.

Also `|C|<=P(8e+f)` and `floor(|C|/L)+1 <= (|C|+L)/L`. The script therefore uses the conservative one-variable envelope

`U_P = gamma_P eP/(gLP-e)+gamma_P fP/(gLP-f)+max((e+f)^2,(8e+f)^2/8)/(2g^2L^2)+(e+f)/(2g^2LP)`.

The final `e+f` term is conservative because `|C|/(PQ)<=e/P+f/Q<=(e+f)/P`. Each endpoint term is strictly decreasing in real `P>0`, since

`d/dP [gamma_P eP/(gLP-e)] = -e(2gLP+eP^2-e)/(2P^2(gLP-e)^2)<0`.

The remaining losses `1/P` and `12/[49P(P+1)]` decrease. Hence the envelope margin is increasing. Exact base checks give cutoffs: g=4: P>=12; g=5: P>=7; g=6: P>=6; g=7..12: P>=5; g=13..15: P>=4. At g=16 the actual midpoint bound handles P<4 and the envelope handles P>=4; all errors decrease thereafter in g.

Below these cutoffs, exact contextwise midpoint evaluation leaves 118 `(g,P,Q)` rows. Literal arbitrary-ratio physical evaluation on all 2,530 contexts gives 298,540 checks and no failure. Weakest gap over `1/294` is `1354081/141260826`, at `(g,P,Q)=(10,2,11)` and `(L,j,e,f)=(4368,2236,8,1)`; overlap is `6240/480479`.

Reproduce from the repository root:

`python3 04-computation/lrc14_disconnected_non35_g4_tail_20260812.py`

`python3 -O 04-computation/lrc14_disconnected_non35_g4_tail_20260812.py`

The two outputs are byte-identical and equal the frozen result.  The script
pins the canonical THM-3350 tail and THM-3352 exact mass engine by SHA-256.

## The all-dilation 3:5 lane

The companion
`lrc14_disconnected_35_small_ruler_symbolic_20260812.py` proves on all 2,530
small-ruler contexts that

`I(3g,5g) >= 1/294` for every integer `g>=1`.

It splits by the exact period `|5e-3f|`, certifies affine branch
stabilization, checks every finite head, and minimizes the resulting exact
quadratic on every residue ray.  There are `56,191` residue classes and the
largest stabilization point is `g=134`; all cleared margins are strictly
positive.  The weakest context is `(L,j,e,f)=(336,174,12,3)`, where the
physical overlap at `g=1` is `158/46397`.

At this historical stage these were frontier theorems rather than a closure.
After the far-ratio, `3:5`, and `g>=4` certificates, the remaining physical
bank was exactly the 29 small rulers, moderate non-`3:5` ratios `P<Q<8P`, and
primitive common scales `g=1,2,3`. A repaired arbitrary-ratio Dirichlet block
argument reduced one route to `22,890` nonzero-resonance affine rays and a
finite raw head `p<264`; see
[`lrc14-generalized-dirichlet-resonance-reduction-20260812.md`](lrc14-generalized-dirichlet-resonance-reduction-20260812.md).
That reduction is proved.  The later
[`head-263 finite-exact certificate`](lrc14-disconnected-head263-finite-exact-certificate-20260812.md)
closes all 201,377 raw head channels by 509,483,810 exact physical
comparisons.

THM-3355 later bypasses the affine scan. Its centered primitive-grid lemma
and coefficient-positive tail leave four genuine scale-one horns, and a
weighted complete-multipartite forest plus a debt-sensitive `K_(1,5)` repair
closes them. Hence the disconnected-low branch is now closed. MISTAKE-377
also records that the old many-turn routing omitted `9|c|<=p`; after repair
the honest affine residual is `8,079,264` occurrences. The affine carrier
inequalities themselves therefore remain an optional OPEN sidecar, not a
proof obligation for the reflected closure.
