# Independent audit: formal cycle defects and the restricted signed-square obstruction

**Verdict: PASS.** The universal formal locality argument, exact degree-nine
expressions, homogeneous signed-square witness, and component-additivity
obstruction in the [producer report](overnight4_20260906_no3line_cycle_defect.md)
are sound. The third-moment geometric coefficients are inherited from the
previous independently audited census, not new counts in this round.

This audit uses a new BFS classifier of literal edge-subset graphs and an
Euler-derivative recurrence for formal logarithms. It imports neither the
producer nor its cyclic-run algorithm. No 45-million-event census is repeated.

## 1. Universal locality audit

The coefficient ring records every component of a non-induced edge subset:
odd paths have one variable, even paths retain the larger shore, and even
cycles retain their circumference. There are no isolated-vertex factors.
This makes disjoint skeleton components multiply exactly. For a simple
bipartite two-regular skeleton with `n` vertices on each shore,

```text
Z_G=product_L P_L^(c_L),
log Z_G=sum_L c_L log P_L.
```

The claimed universal path series `b` is justified coefficientwise. Use
nonempty connected edge subsets as polymers, incompatible when they share
a vertex and incompatible with themselves. Every compatible family is
exactly the component decomposition of a unique edge subset. The formal
logarithm retains only connected incompatibility clusters, with repeated
polymers allowed in its expansion.

I checked the cited primary source directly:
[Fernández--Procacci, *Cluster Expansion for Abstract Polymer Models*,
equations (2.2)--(2.4)](https://webspace.science.uu.nl/~ferna107/papers/ferpro07.pdf).
It includes self-incompatibility and states that the logarithmic coefficient
vanishes for a disconnected incompatibility graph. Only that formal identity
is used; no convergence bound or positivity theorem is imported.

A cluster of total edge weight `m` has a connected union with at most `m`
distinct edges. For `L>m`, it cannot fill the cycle. Its union is a proper
path arc, and lifting that arc to the infinite alternating path preserves
every polymer incidence, including repeated polymers. Fixing the initial
edge of the arc in an oriented cycle removes rotational ambiguity. Each
shore orientation then has `L/2` placements; the cluster weight and its
combinatorial coefficient are independent of `L`. This proves

```text
[log P_L]_m=L[b]_m  when L>m.
```

The completion by positive edge weight is sufficient: only finitely many
component variables and cluster patterns enter a fixed weight. Full-cycle
polymers cannot enter below their circumference. Therefore

```text
Delta_L=log P_L-Lb starts in weight L,
Z_G=exp(2n b)*exp(sum_L c_L Delta_L).
```

The equality uses the exact edge count `sum_L L*c_L=2n`. The locality
argument concerns the multiplicative component-type ring before geometric
averaging; it does not apply a logarithm through that averaging map.

## 2. Independent coefficient and truncation audit

The referee enumerates all edge combinations of weights at most nine in
cycles of lengths `4,6,8,10,12,16`. For each subset it constructs the
vertex adjacency graph, finds connected components by BFS, and counts
edges and both shores. This independently reconstructs every monomial.

It obtains the formal logarithm from `P*Euler(log P)=Euler(P)`:

```text
L_m=P_m-(1/m)sum_(i=1)^(m-1) (m-i)P_i L_(m-i),
where log P=sum_m L_m.
```

This differs from the producer's powers-of-`P-1` logarithm. The length-10,
12 and 16 logarithms divided by circumference agree through weight nine,
and the four/six/eight-cycle defects start at their stated circumferences.

Writing `D=[Delta_4]_4`, `p=p1`, `u=p2^L`, `v=p2^R`, and
`r=p4^L,s=p4^R`, the independent coefficients are exactly

```text
D=z4+p^4-2p^2(u+v)+4p*p3+u^2+v^2-2(r+s),
T=p^2*p3-p*u*v+p*z4-p(r+s)+(u+v)*p3-p5,
[Delta_4]_5=-8p*D+4T.
```

The four leading/next defect sizes are respectively `9,13,45,74` monomials
for `D4,D5,D8,D9`. Expanding the exact exponential gives

```text
[c4^2]Z_G = D^2/2+(n-8)p*D^2+4D*T through weight9,
[c8]Z_G   = D8+D9+2n*p*D8 through weight9.
```

For the first formula, `Delta_4^2/2` contributes `D^2/2+D*D5`, while
the bulk's weight-one term contributes `n*p*D^2`. This accounts for both
the factor one-half and the cancellation at `n=8`. Mixed `c4*c6` terms
start at weight ten and cannot enter this computation.

At `n=8` the full resulting polynomials, not merely geometrically surviving
terms, agree with the independent skeleton contrasts

```text
[c8] = (P8^2-P16)/2,
[c4^2] = (P4^4-4P4*P12+3P16)/12.
```

They have 108 and 103 monomials. Specializing every component variable to
the corresponding power of one edge variable annihilates each defect,
consistent with `Z_G -> (1+t)^(2n)`.

## 3. Geometric normalization and the meaningful square obstruction

Only the frozen weight-eight/nine geometric bank is used. For each of its
150 profiles the referee reconstructs the number of copies in `K_(8,8)`
from `(8)_left*(8)_right` divided by shore-preserving automorphisms:
an odd path contributes one, an even path two, a cycle its edge count,
and repeated identical components contribute their permutation factorial.
Every grid-copy denominator and every formal coefficient agrees with the
inherited certificate. Restoring the factor six for ordered triples gives

```text
E=172483/529200,
F=3W8(D^2)+24W8(D*T)
 =456371/2116800+42631/2116800
 =11881/50400.
```

Global nonpositivity on the entire unital ring is already immediate from
`W8(1)=0` while `W8` is nonzero. The substantive new witness is narrower:
`Q=u*(p1^2-u)` is homogeneous of edge weight four, and every monomial of
its square is a **realizable weight-eight forest type**. In that retained
sector the exact independently normalized values give

```text
W8(Q^2)=16747/1058400-2*(2749/235200)+1/147
        =-397/529200<0.
```

Thus positivity cannot justify the `D^2` term even when restricted to
the actual weight-four/weight-eight setting of the proposed argument.
There is no negative expected square of an actual random variable:
formal multiplication joins component types, while `W8` counts exactly
three whole geometric events under global injective shore labeling.

The map preserves each complete union-type count and its ordered-event
factor. It loses event incidence/count under multiplication and independent
allocation of the global labels. An event-number grading is a plausible
necessary sidecar for a repaired construction, but adding that grading alone
is not proved sufficient for positivity or component independence.

## 4. Third-cumulant and fixed-n identifiability audit

The ordinary/factorial conversion is correct:

```text
K3=M3-3mu*M2+2mu^3,
kappa3=M3+3(1-mu)M2+mu-3mu^2+2mu^3.
```

At fixed `n`, the inherited short-edge theorem makes `mu` cycle independent
and `M2` affine in `c4,c6`. Thus neither cumulant subtraction can change
the `c8` or `c4^2` coefficients. This uses the complete moment operators,
not an assumed independence of cycle components.

The five valid `n=8` skeletons `C16,C4+C12,C6+C10,2C8,4C4` give determinant
24 on the basis `1,c4,c6,c8,c4^2`, independently checked by the full signed
permutation determinant. There is no fixed-n basis ambiguity here.

More decisively, the contrast

```text
stat(2C8)-2stat(2C4+C8)+stat(4C4)=8F=11881/6300
```

annihilates every additive sum of component contributions: its combined
coefficient of each cycle length is zero. All three skeletons have the same
`n=8`, so allowing the per-component contribution to depend on `n` does not
repair additivity. The positive contrast applies to `M3` and both third
cumulants above. A connected expansion involving multiple skeleton
components remains possible and is not refuted.

## 5. Reproduction and audit boundary

[Independent source](../../04-computation/overnight4_20260906_no3line_cycle_defect_audit.py),
[normal output](overnight4_20260906_no3line_cycle_defect_audit.out), and
[optimized output](overnight4_20260906_no3line_cycle_defect_audit_optimized.out):

```text
python 04-computation/overnight4_20260906_no3line_cycle_defect_audit.py
python -O 04-computation/overnight4_20260906_no3line_cycle_defect_audit.py
```

Both pass **278,884 explicit gates**, which remain active under optimization.
The evaluator is deliberately restricted to retained geometric weights eight
and nine. It does not treat absent lower-weight data as zero geometry.

The inherited complete geometry certificate is frozen at SHA-256

```text
f2c566ac1b2bcedb530af72f3a290479841e7db87844b331558bed68c93ba727
```

The all-circumference conclusion is the formal connected-cluster proof;
the finite graph enumeration independently challenges its low-weight
consequences. No convergence, all-size sign, zero-defect probability,
or geometric event-independence conclusion is asserted.

Audit source SHA-256:

```text
3a0af6e0adb422a2db494d8fe60f9d23fcefc1aef663122de2191142ad2c74c0
```

Normal and optimized outputs are byte-identical after explicit LF
normalization, with SHA-256:

```text
8fe7462a695b3506f370823e8b4f2c64c5ac253706a618c28edcdce59c84074a
```
