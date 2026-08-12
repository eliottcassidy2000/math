# Disconnected-low geometry and the non-3:5 small-ruler g>=4 tail

Status: PROVED structural/analytic reductions + VERIFIED-EXACT finite
residual.  The disconnected branch itself remains open.

## Structural reduction

Join two of six distinct levels when their reduced ratio `(P,Q)` satisfies
`P+Q<=7`.  If this low graph is disconnected, its components have primitive
connected shapes of sizes at most five.  The connected-shape counts from
THM-3350 are

`c_1,...,c_5 = 1,8,94,1295,19389`.

Consequently the number of unordered component-shape profiles is

`[x^6] product_(n=1)^5 (1-x^n)^(-c_n)=36520`.

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
the profile ledger, multipartite tree counts, and all displayed exact gaps.

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

This is a frontier theorem, not a closure of the disconnected branch.  The
remaining physical bank consists of the 29 small rulers, moderate ratios
`P<Q<8P`, and primitive common scales `g=1,2,3`.  A repaired arbitrary-ratio
Dirichlet block argument gives a further analytic cutoff at raw smaller level
`p>=264`; its finite affine-resonance compiler remains open.
