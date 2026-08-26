---
id: THM-4202
title: "Vertex-transitive ordinal-remainder positivity"
status: >
  PROVED RELATIVE TO THM-4181/4184/4187 + VERIFIED-EXACT PRIMARY +
  ALGORITHMICALLY INDEPENDENT FINITE AUDIT. For all finite nonempty
  vertex-transitive tournaments A,B with |B|>=3, R_+(A,B)>0. No
  self-converse hypothesis is used. If only A is vertex-transitive, an exact
  all-order identity remains, and at every right order its complete defect
  from the symmetric baseline is a variance-covariance penalty. Regular score first
  fails to force those sidecars at order seven. General (OS+), the all-order
  no-sink gate law, LRC(14), JC(2), and DC(2) remain OPEN.
source: open-frontiers-incoming-20260826b / forgotten-tournament signal
depends_on:
  - THM-4181-ordinal-sum-capacity-transfer-and-parity-component-exchange
  - THM-4184-path-cover-parity-ordinal-cocycle-and-lollipop-positivity
  - THM-4187-universal-source-padding-supermodularity-and-transitive-left-ordinal-positivity
related:
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4193-cycle-first-transitive-tail-crossing-and-transitive-context-positivity
script: 04-computation/tournament_vertex_transitive_ordinal_remainder_positivity_thm4202.py
output: 05-knowledge/results/tournament_vertex_transitive_ordinal_remainder_positivity_thm4202.out
independent_audit_script: 04-computation/tournament_vertex_transitive_ordinal_remainder_positivity_independent_audit_thm4202.py
independent_audit_output: 05-knowledge/results/tournament_vertex_transitive_ordinal_remainder_positivity_independent_audit_thm4202.out
minimal_hostile_script: 04-computation/tournament_regular_order7_minimal_hostile_thm4202.py
minimal_hostile_output: 05-knowledge/results/tournament_regular_order7_minimal_hostile_thm4202.out
script_sha256: 954a2a0b85036d446a1bbf5ec1c0f5f46e9f5b2646a7d72a63c17775569dbafd
output_sha256: 96d7a18c1aaff89556ecae3456471d37d910aa6cd3edaa71555c765389384634
primary_dependency_lrc_sha256: 0a60b3133d567eea99da3d29e85c70352a918f90d919647724e3c97f12ba6ba7
primary_dependency_thm4184_sha256: 9ab09bf8b70ee5dcf3d86698180a456d67f012655b49a16dfea9903caefbb39c
independent_audit_script_sha256: 37e0555155ec6715a238df62cb6a9327156631ee84214de70a0c61edd5be1b18
independent_audit_output_sha256: 7962a5404785bc06399e84495783b585503dee6f3301c196797f3dae24f97d8a
minimal_hostile_script_sha256: 0f03b29fe3a15902a44ea63170d7f5d77d9ecd5c02748fa161b9805d949e9244
minimal_hostile_output_sha256: e56e2f8d0d80c354b38b7e25a8c00d4517b05b084437a76c0c6c4dcb67c3f08e
hash_basis: raw LF bytes
primary_audit: >
  PASS. The inherited exact engine checks 930 vertex-transitive pairs, 90
  transitive-prefix rows, and the one-sided identity on all 1,098 labelled
  right factors of orders two through five. It also freezes an exact
  regular order-nine hostile. Semantic digest
  faf77ea0f456efdfe0cdfdafac5fa880c401196c0d056796494fc1b1d8c3b22a.
independent_audit: >
  ACCEPT. A standalone literal-permutation path rebuilds Hamilton counts,
  exposure capacities, rooted states, and 26 direct ordinal children without
  importing tournament code. It independently reconstructs 74 even/odd
  one-sided defects and the order-nine hostile. Semantic digest
  9428843ca20a1b09d9dbc6c47d5f0a89d65a37c722117e637a0dbfc7f0c46056.
minimal_hostile_audit: >
  ACCEPT. A standalone literal enumeration proves every regular tournament
  below order seven vertex-transitive, then reconstructs the first regular
  non-vertex-transitive hostile, its order-ten child, and the exact 3,816
  variance-covariance penalty. Semantic digest
  87a5682f9b258a3ee766dbc476c585c0aea6ca48f16a4551b6e03e957ccb5df2.
---

# THM-4202 -- vertex-transitive ordinal-remainder positivity

**PROVED RELATIVE TO THM-4181/4184/4187 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** General `(OS+)`, the all-order no-sink gate law, LRC(14), JC(2),
and DC(2) remain open.

Closest proved mechanisms: THM-4177's root capacity/current and odd-path
decomposition, THM-4181's rank-two ordinal transfer and exact block-gate
formula, THM-4184's normalized remainder cocycle, and THM-4187's
universal-source supermodularity. The main guardrail is MISTAKE-013: a
vertex-transitive tournament need not be self-converse. No anti-automorphism
or self-converse hypothesis is used below.

## 1. Theorem and strict scope

Let `A,B` be finite nonempty vertex-transitive tournaments and suppose
`n=|B|>=3`. With THM-4181's notation

```text
R_+(A,B)
 =G_+(c(A triangleright B))
  -H(B)^2G_+(c(A))-H(A)^2G_+(c(B)),
```

the theorem is

```text
R_+(A,B)>0.                                             (1)
```

There is no equality case on this scope. A singleton `A` is allowed. The
condition `|B|>=3` is essential: `B=1` includes THM-4187's exact negative
family `R_+(C3,P_1)=-72`.

An immediate proved-relative corollary of `(1)` and THM-4187 is

```text
R_+(P_a triangleright A,B)>R_+(A,B)>0,                 (2)
```

for every `a>=1`. Indeed THM-4187 gives

```text
R_+(P_a triangleright A,B)
 =R_+(A,B)+Theta(P_a,A,B),
```

and its equality case for `Theta(P_a,A,B)` is incompatible with `|B|>=3`.
Thus `(2)` supplies nontransitive-left all-order families, for example
`P_a triangleright C3` against every nontrivial vertex-transitive right
factor.

## 2. Uniform rooted and capacity coordinates

Write

```text
m=|A|,       n=|B|,
h=H(A),      s=H(B),
w=W(A),      v=W(B).                                   (3)
```

A vertex-transitive tournament has constant outdegree `(order-1)/2`, so its
order is odd (apart from the already odd singleton). Automorphisms preserve
Hamilton paths, marked exposed words, directed simple paths, and their
complements. Therefore all coordinates in each of the following vertex
families are constant:

```text
d_a(c(A)), o_a(c(A)), r_a(c(A)), U_a^0(A), U_a^1(A),
d_b(c(B)), o_b(c(B)), r_b(c(B)), V_b^0(B), V_b^1(B).    (4)
```

Summing capacity degrees and directed capacity masses gives

```text
d_a=2w/m,        o_a=r_a=w/m,
d_b=2v/n,        o_b=r_b=v/n.                           (5)
```

The equalities `o=r` use only transitivity of the automorphism action:
outgoing and incoming capacity masses are separately constant, while each
sums globally to the same total edge mass. They do not use an
anti-automorphism.

THM-4181's path-cover identities at odd order give

```text
U^0+U^1=w+h,        U^1-U^0=h,
V^0+V^1=v+s,        V^1-V^0=s.                          (6)
```

Dividing the total coordinates uniformly over the vertex orbits yields

```text
U_a^0=w/(2m),             U_a^1=(w+2h)/(2m),
V_b^0=v/(2n),             V_b^1=(v+2s)/(2n).            (7)
```

These formulas include the singleton: `(U^0,U^1)=(0,1)` there.

## 3. Uniform cross capacity and exact remainder formula

THM-4181's rank-two transfer now makes every cross capacity in
`A triangleright B` equal:

```text
c_ab=2(U_a^0V_b^0+U_a^1V_b^1)
    =(wv+ws+hv+2hs)/(mn).                               (8)
```

Put

```text
K=wv+ws+hv+2hs.                                        (9)
```

In THM-4181's block-gate notation, the cross row/column/total square data are

```text
p_a=K/m,       q_b=K/n,       Z=K,
Z_2=K^2/(mn),  sum_a p_a^2=K^2/m,  sum_b q_b^2=K^2/n. (10)
```

The scaled internal tensors are `x=s c(A)` and `y=h c(B)`. Equations `(5)`
give their two linear brackets:

```text
W_x-d_a^x+4o_a^x = sw(1+2/m),
W_y-d_b^y-4r_b^y = hv(1-6/n).                           (11)
```

Substitution of `(10)--(11)` into THM-4181's exact generic block formula,
after removing the two internal gates as in the definition of `R_+`, gives

```text
R_+(A,B)
 =h s w v
  +K s w(1+2/m)
  +K h v(1-6/n)
  +(K^2/2)(1+1/(mn)+3/m-5/n).                          (12)
```

Normalize by

```text
x=w/h,       y=v/s,       k=xy+x+y+2=K/(hs).           (13)
```

Then the exact dimensionless form is

```text
R_+(A,B)/(h^2s^2)
 =xy+kx(1+2/m)+ky(1-6/n)
  +(k^2/2)(1+1/(mn)+3/m-5/n).                          (14)
```

### 3.1 One-sided symmetry and the exact defect address

There is a stronger exact identity before imposing transitivity on `B`.
Assume only that `A` is vertex-transitive and let `B` be arbitrary. Since the
left rooted states are uniform, the cross capacity is independent of
`a in A`; write

```text
t_b=c_ab(A triangleright B)
   =[w V_b^0(B)+(w+2h)V_b^1(B)]/m,
S=sum_b t_b,        Q=sum_b t_b^2.                     (14a)
```

The block formula gives the all-order one-sided identity

```text
R_+(A,B)
 =h s w v +(m+2)s w S
  +m h sum_b t_b(v-d_b-4r_b)
  +(m/2)[(m+3)S^2-(5m-1)Q].                           (14b)
```

THM-4184's all-order parity balance holds without any symmetry of `B`, so for
every `n=|B|`

```text
S=K/m.                                                  (14c)
```

Put

```text
tbar=S/n,       dbar=2v/n,       rbar=v/n,
Var_t=sum_b(t_b-tbar)^2,
Cov=sum_b(t_b-tbar)[(d_b-dbar)+4(r_b-rbar)].            (14d)
```

Then `Q=S^2/n+Var_t`, and exact centering in `(14b)` gives

```text
R_+(A,B)
 =R_sym(h,w;s,v;m,n)
  -m h Cov -(m/2)(5m-1)Var_t,                          (14e)
```

where `R_sym` is precisely the right side of `(12)`. Thus the unresolved
one-sided `(OS+)` problem at **every right order** has an exact address: bound
the rooted-state variance and its covariance with `d+4r` by the positive
symmetric baseline. This is more precise than asking for pointwise
cross-capacity dominance, which THM-4181 already refutes. THM-4177's
incoming-fan injection is the natural least-used sidecar for attacking the
covariance term.

## 4. Positivity in the three right-order regimes

First note the marked-Hamilton-path lower bound

```text
W(X)>=(|X|-1)H(X).                                     (15)
```

Indeed every Hamilton path of `X`, with any one of its `|X|-1` actual
adjacencies marked, is a distinct exposed marked word counted by `W(X)`.
Thus

```text
x>=m-1,          y>=n-1.                               (16)
```

Only the first inequality is needed below.

### 4.1 `n>=7`

Every displayed term in `(14)` is nonnegative, because

```text
1-6/n>0,
1+1/(mn)+3/m-5/n > 1-5/7>0.                           (17)
```

The sum is strict. If `x=0` (only possible here at `A=1`), the positive
`ky(1-6/n)` and quadratic terms remain.

### 4.2 `n=5`

The quadratic coefficient in `(14)` is `16/(5m)` before multiplication by
`1/2`. Combine it with the only negative term:

```text
R_+/(h^2s^2)
 =xy+kx(1+2/m)+ k(8k-my)/(5m).                         (18)
```

Using `k=xy+x+y+2` and `x>=m-1`,

```text
8k-my
 =y(8x+8-m)+8x+16
 >=7my+8x+16>0.                                        (19)
```

Therefore `(18)` is strict.

### 4.3 `n=3`

The only vertex-transitive tournament of order three is `C3`, for which

```text
s=3,       v=6,       y=2.                             (20)
```

Equation `(14)` becomes

```text
R_+/(h^2s^2)
 =[63x^2+(144-24m)x+80-40m]/(3m).                     (21)
```

Put `u=x-m+1>=0`. The numerator is

```text
63u^2+(102m+18)u+39m^2+2m-1.                          (22)
```

Every term except the final `-1` is nonnegative and
`39m^2+2m-1>=40` for `m>=1`. Hence `(22)` is strict. This completes the
proof of `(1)`.

## 5. What this changes, and what it does not

This supplies an all-order `(OS+)` family with a nontransitive, source-free
left factor: take any nontrivial vertex-transitive `A` and any nontrivial
vertex-transitive `B`. It complements THM-4187's transitive-left theorem and
THM-4193's cycle-first prefixes in transitive contexts rather than
duplicating either.

For a hypothetical counterexample to general `(OS+)`, THM-4187 already
permits left-order minimization to remove universal sources. The theorem
further says that a counterexample with vertex-transitive right
factor cannot also have vertex-transitive left factor. With `(2)`, no
transitive prefix placed before such a left core can create a counterexample.

This does **not** prove general `(OS+)`: an arbitrary no-sink right factor
does not have uniform rooted `V` states. It also does not prove THM-4177's
no-sink/no-source gate law, the order-eleven asymmetric Johnson bank, exact
Johnson cosets, or actual response maximizers.

## 6. Regular score is not enough: sharp order-seven boundary

### 6.1 First possible failure

The standalone hostile audit exhausts every labelled tournament below order
seven. The regular labelled/class counts at orders `1,...,6` are

```text
(1/1),(0/0),(2/1),(0/0),(24/1),(0/0),                 (23a)
```

and every regular row is vertex-transitive. At order seven, start with the
circulant tournament on `Z/7` with connection set `{1,2,3}` and reverse the
directed triangle `(0,2,4)`. The resulting pair-bit label is

```text
101100111001010111111,                                 (23b)
```

with score `(3,3,3,3,3,3,3)`, but its automorphism group has order three and
vertex orbits

```text
(0,3,5), (1,4,6), (2).                                (23c)
```

Literal enumeration gives

```text
H=171, W=1638,
V states=((118,142),(116,141),(117,141),(118,142),
          (116,141),(118,142),(116,141)),
capacity degrees=(468,468,468,468,468,468,468),
incoming masses=(236,232,234,236,232,236,232).         (23d)
```

Thus regular score and even uniform capacity degree do not force the rooted
and directed sidecars. With left factor `C3`, the cross row and exact
remainders are

```text
t=(804,796,798,804,796,804,796),
actual R_+=233,570,700,
false uniform value=233,574,516,
difference=3,816.                                      (23e)
```

The all-order defect identity accounts for the discrepancy exactly:

```text
Var(t)=696/7, Cov(t,d+4r)=192,
9 Cov+21 Var=3,816.                                    (23f)
```

The independent audit constructs the order-ten child both literally and by
a separate transfer assembly. Hence order seven is the first possible
regular-score failure of the uniform sidecars used in the proof.

### 6.2 Order-nine simultaneous sidecar failure

Start with the consecutive circulant tournament on `Z/9`, connection set
`{1,2,3,4}`, and reverse the directed triangle `(0,2,5)`. A directed-cycle
reversal preserves one incoming and one outgoing edge at each of its three
vertices, so the result still has regular score `(4,...,4)`. Its pair-bit
label is

```text
101110001111000110100111101111111111.                  (23)
```

The independent referee reconstructs

```text
H=3255,

V states:
(2402,2767),(2390,2758),(2387,2749),
(2391,2748),(2394,2753),(2378,2738),
(2390,2747),(2381,2745),(2376,2739),

capacity degrees:
9584,9584,9548,9534,9550,9536,9536,9550,9534,

incoming capacity masses:
4804,4780,4774,4782,4788,4756,4780,4762,4752.         (24)
```

Thus all three symmetry sidecars used in `(5)--(7)` are nonuniform. For the
left factor `C3`, the actual exact remainder and the value obtained by
illegally imposing the uniform formula are

```text
actual R_+(C3,B)       =169,021,426,548,
false uniform value    =169,022,173,596,
difference             =747,048.                       (25)
```

The exact defect ledger `(14e)` explains every unit of this discrepancy:

```text
(t_b)=(15872,15812,15770,15774,15800,15708,15768,15742,15708),
Var_t=21,496,              Cov=32,848,
9 Cov+21 Var_t=747,048.                                 (25a)
```

The first failed implication is

```text
regular outdegree => uniform rooted/capacity sidecars. (false)
```

The strongest survivor is exactly vertex transitivity, or more generally an
explicit hypothesis that all coordinates in `(4)` are uniform. A score
vector, regularity, `H`, or `W` alone is insufficient.

## 7. Exact audits

### 7.1 Primary inherited-engine audit

The primary certificate imports the current THM-4184 subset-DP and
ordinal-transfer implementation. It checks every labelled circulant
presentation in the banks

```text
n=1,3,5,7,9:       1,2,4,8,16 presentations.           (26)
```

Against every nontrivial right presentation this gives `930` exact formula
checks, all strict, with minimum

```text
R_+(1,C3)=120.                                          (27)
```

It also checks `90` prefixed-cycle rows from `(2)`, `(14b)--(14e)` on all
`1,098` labelled right factors of orders two through five with left `C3`,
the order-nine hostile, and the LRC projection hostile in Section 8. The
frozen semantic digest is

```text
faf77ea0f456efdfe0cdfdafac5fa880c401196c0d056796494fc1b1d8c3b22a.
```

### 7.2 Algorithmically independent referee

An independent standard-library script imports no tournament computation
module. It rebuilds:

- `H` and every exposure capacity by literal permutations with one marked
  gap;
- `U,V` by literal directed simple paths weighted by a separately enumerated
  complement Hamilton count;
- `26` ordinal children of total order at most eight by direct exposed-word
  enumeration, not by the rank-two transfer;
- the one-sided defect identity on all eight labelled order-three right
  factors;
- the order-nine hostile sidecars literally; and
- the order-nine hostile child gate with its own rank-two assembler only
  after those literal factor sidecars are fixed (direct order-twelve
  permutation enumeration would be wasteful).

All `26` direct child remainders equal `(12)` and are strict, again with
minimum `120`. Its semantic digest is

```text
9428843ca20a1b09d9dbc6c47d5f0a89d65a37c722117e637a0dbfc7f0c46056.
```

Thus an algorithmically independent audit is not merely feasible; the
bounded referee is implemented. The all-order result still rests on the
elementary derivation `(4)--(22)`, not on finite enumeration.

### 7.3 Minimal regular-score hostile audit

A third standalone program imports no repository module. It exhausts all
labelled tournaments through order six to establish the minimal boundary in
`(23a)`, computes the order-seven automorphism orbits, and reconstructs every
entry in `(23d)--(23f)`. Unlike the order-nine child audit, it can enumerate
the complete order-ten ordinal child literally; an independent transfer
assembly gives the same `233,570,700`. Its semantic digest is

```text
87a5682f9b258a3ee766dbc476c585c0aea6ca48f16a4551b6e03e957ccb5df2.
```

Artifact byte hashes are

```text
primary script       954a2a0b85036d446a1bbf5ec1c0f5f46e9f5b2646a7d72a63c17775569dbafd
primary output       96d7a18c1aaff89556ecae3456471d37d910aa6cd3edaa71555c765389384634
independent script   37e0555155ec6715a238df62cb6a9327156631ee84214de70a0c61edd5be1b18
independent output   7962a5404785bc06399e84495783b585503dee6f3301c196797f3dae24f97d8a
hostile script       0f03b29fe3a15902a44ea63170d7f5d77d9ecd5c02748fa161b9805d949e9244
hostile output       e56e2f8d0d80c354b38b7e25a8c00d4517b05b084437a76c0c6c4dcb67c3f08e
```

## 8. LRC deletion/resonance map: exact failure and survivor

There are two distinct precise failures.

### 8.1 Half-turn tournament does not preserve the LRC safe predicate

Take moving speeds `(1,3,4)` with the stationary observer, so `N=4`. At

```text
t=41/96:
positions=(0,41/96,9/32,17/24),
distances to 0=(41/96,9/32,7/24),       safe;

t=7/48:
positions=(0,7/48,7/16,7/12),
distances to 0=(7/48,7/16,5/12),        unsafe.         (28)
```

The two half-turn tournaments have the same unmarked canonical word
`0001100011000110` and the same observer-pointed canonical word
`0011100001000110`; both have `H=5`, score `(1,1,2,2)`, and observer score
two. Thus the source, target, and loss contract is

```text
source:       exact LRC phase cell with threshold 1/4,
map:          observer-pointed half-turn tournament class,
preserved:    circular pair orientations, H, score, observer score,
destroyed:    metric distances and threshold slack,
target bit:   all moving phases have distance at least 1/4,
hostile:      the two rational cells in (28),
needed sidecar: exact gap lengths/endpoint slacks.       (29)
```

So even the pointed tournament does not preserve the actual LRC predicate.

### 8.2 A static pool/deletion tournament does not preserve THM-4188 repair

For THM-4188's fixed pool, keep the same anchor and deletion set

```text
A={40,63,84,286},
R={16,20,60,80,95,170,240},       A intersect R=empty. (30)
```

The exact theorem gives

```text
R in E_7(50),             R notin E_7(48).              (31)
```

Any static tournament built only from pool order, anchor membership,
deletion membership, or their pairwise incidence is identical in the two
rows, while the actual repair predicate

```text
mu(G_((P union {q})\R))>=4/63                          (32)
```

flips. The destroyed coordinate is the newcomer joint-wall geometry. The
needed sidecar is the exact `lcm(D,14q)` atom refinement (or an equivalent
exact measure/resonance word). A tournament on repair routes may still be a
scheduler, but hyperedge membership, anchor disjointness, and exact Haar mass
must remain attached; it is not itself a repair certificate.

THM-4190 gives the positive controlled-forgetting model. Its map

```text
R |-> R intersect V_i
```

is legal because it proves the exact iff

```text
K hits every projected edge
iff
K hits every anchor-disjoint repair.                   (32a)
```

The faithful target is therefore a projected **hypergraph**, not a bare
tournament. The four depth-five blocker rows at `q=25,105,210,256` sharpen
the warning: a shallow certificate can fail while every blocker body has
positive direct Haar margin, and the exact depth-six hypergraph repairs it.
This is the same procedural pattern as `(14e)`: expose the defect coordinate
before deciding the target predicate.

Even the pairwise co-occurrence graph of a repair hypergraph is insufficient.
On labels `{1,2,3}`, the hypergraphs

```text
{{1,2},{1,3},{2,3}}       and       {{1,2,3}}
```

have the same complete two-section, but for anchor `{1}` the first has the
disjoint edge `{2,3}` and the second has none. Orienting that common
two-section cannot restore the lost higher-arity predicate. A useful
tournament must therefore remain a scheduler over the hypergraph, with the
full projected edge incidence attached.

## 9. Planar Jacobian root-splitting map: exact projection obstruction

The tempting bridge is to take critical points as vertices and orient them by
a root/projection coordinate. THM-4189 supplies an exact hostile. At the
algebraic coefficient specialization `(39)--(42)` of that theorem, two
distinct reduced Morse points

```text
(X,T)=(0,tau),(1,tau)                                  (33)
```

share the same normalized `T` coordinate, while the source gauge separates
them as

```text
(s,p)=(0,tau),(tau,tau+tau^2).                          (34)
```

Therefore:

```text
source:       full reduced critical scheme,
map:          normalized T projection / T-root list,
preserved:    projected root location and resultant support,
destroyed:    fibre cardinality, point ownership, Hessian/source value,
target data:  affine critical length L=23 and response supports,
hostile:      gcd_X(f,h)=X(X-1) over T=tau,
needed sidecar: source p (and ultimately full fibre, Hessian, value). (35)
```

If equal `T` values are merged, the affine critical length is wrong. If they
are kept as two vertices, their pairwise projection observable ties and does
not define an intrinsic tournament; an imposed `X` tie-break is a coordinate
gauge and does not preserve the geometric response predicate. THM-4189's
second hostile, where an extra point joins the universal `T=-1/6` fibre,
shows that this is structural rather than isolated.

The genuine overlap with THM-4177 is a controlled-forgetting principle:

```text
tournament: corrected restriction of the child tensor != card capacity;
Jacobian:   projected resultant root != full critical fibre.             (36)
```

Both require a root/fibre sidecar before splitting. No Johnson sign or
ordinal-remainder predicate is transported to the Jacobian wall.

THM-4192 sharpens this from a projection warning to a specialization
warning. On the `K=0` wall the generic quadratic carrier is not two limiting
labelled vertices: it becomes one primitive rational toric branch, and there
is no finite carrier response. Its three coefficient strata have residual
lengths `18,15,14` (affine lengths `22,19,18`) and distinct packet shapes.
Thus a tournament whose vertices were the generic carrier branches would not
even have a conserved vertex set under the wall specialization. The legal
order of operations is

```text
specialize -> saturate coordinate artefacts -> factor primitive branches
 -> attach fibre/Hessian/value data -> only then consider pairwise structure.
```

This is exactly the tournament corrected-restriction lesson at a different
type: split the **actual specialized object**, not the generic object's
projected shadow.

## 10. Replay

From the repository worktree:

```powershell
python -B 04-computation/tournament_vertex_transitive_ordinal_remainder_positivity_thm4202.py
python -B 04-computation/tournament_vertex_transitive_ordinal_remainder_positivity_independent_audit_thm4202.py
python -B 04-computation/tournament_regular_order7_minimal_hostile_thm4202.py
```

The streams line-match the three frozen `.out` files. The all-order conclusion
rests on the symbolic proof, while these finite paths provide positive and
hostile controls.
