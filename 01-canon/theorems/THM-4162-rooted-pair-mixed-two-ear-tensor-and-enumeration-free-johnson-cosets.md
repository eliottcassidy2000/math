---
id: THM-4162
title: "Rooted-pair mixed-two-ear tensor and enumeration-free Johnson cosets"
status: >
  PROVED ELEMENTARY ROOTED ENDPOINT TENSOR + ENUMERATION-FREE EXACT JOHNSON
  COSETS + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For every rooted
  tournament (Q,r), the parent endpoint table and a double-clone endpoint
  slice reconstruct the complete directed-exposure capacity tensor of the
  homogeneous-pair expansion P_r(Q). At quotient order ten the new recursive
  slice has 4,608 entries; the child cut field has 55 new integer coordinates
  plus the inherited internal capacity. For any capacity tensor, all
  cardinality-layer moments, Johnson lattices, and exact coset floors follow
  in O(N^3) arithmetic without enumerating cuts. This theorem alone makes no
  order-eleven centrality or actual-maximizer claim.
source: codex-frontier-synthesis-creative-20260825ar
depends_on:
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
  - THM-4131-strong-tournament-centrality-through-order-eight
  - THM-4145-rooted-homogeneous-pair-expansion-two-defect-formula
related:
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4137-strong-tournament-centrality-complete-order-ten
  - THM-4144-order-eleven-large-homogeneous-module-johnson-centrality
  - THM-4163-order-eleven-homogeneous-pair-johnson-centrality
script: 04-computation/tournament_rooted_pair_mixed_twoear_tensor_thm4162.py
output: 05-knowledge/results/tournament_rooted_pair_mixed_twoear_tensor_thm4162.out
script_sha256: 2ebe3e743e0479b420dc3443adee71a729efb756aa13e93d5f45509a59353d57
output_sha256: e10247ea988ffc01031691dd0f6c361655822ed41e2ed76fd4462f970f75d2f1
finite_engine: 04-computation/tournament_order_eleven_pair_module_centrality_engine_thm4163.cpp
finite_engine_sha256: 2cc8a6ef46a189db9b77c8cc929fd02659e571d0a63410c339f8949c17d8dec0
independent_literal_audit: 04-computation/tournament_order_eleven_pair_module_centrality_independent_literal_audit_thm4163.cpp
independent_literal_audit_sha256: f88d4b3414e2f5326b84ddea49c1f097f4409c7deb01048fac06e54d971b2d9b
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standard-library Python implementation separately builds the
  literal child subset-DP and the rooted double-clone tensor, then compares
  every capacity, moment, lattice, and clone-decomposition formula on all
  5,404 rooted labelled quotients of orders two through five. It also checks
  the directed-triangle one-defect hostile and a strong order-seven row with
  every interior lattice equal to six. Normal, optimized, and hash-seeded
  executions byte-match the frozen output.
independent_audit: >
  ACCEPT. A warning-clean C++ tensor engine implements the parent-once
  recurrences and O(N^3) lattice formulas, while a separate literal child-DP
  engine enumerates all 2^11 cuts. On the pinned 0/1024 strong order-ten
  shard, all 20,580 rooted rows have identical semantic sum/xor, extrema,
  zero failures, and fourteen non-parity lattice hostiles.
---

# THM-4162 -- rooted-pair mixed-two-ear tensor and enumeration-free Johnson cosets

**PROVED ELEMENTARY ROOTED ENDPOINT TENSOR + ENUMERATION-FREE EXACT JOHNSON
COSETS + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement, inheritance, and information budget

Let `(Q,r)` be a rooted tournament of order `q>=2`, put
`U=V(Q)-{r}`, and let `T=P_r(Q)` replace `r` by clones `a->b` with the same
external orientation as `r`. THM-4145 proves strongness equivalence,
`H(T)=F_Q(N_Q^+(r))`, and the internal capacity

```text
c_ab=2w_ab=2H(Q).                                      (1)
```

For an arbitrary tournament `T` of order `N`, orient each unordered edge as
its tournament arc and write `c_ij=2w_ij`. THM-4128's directed-cut identity is

```text
F_T(S)=H(T)+sum_(i in S, j notin S, i->j)c_ij.         (2)
```

> **Theorem.** The parent endpoint tables of `Q` and the double-clone
> recurrences in Section 2 reconstruct the complete child endpoint oracle and
> hence every capacity `c_ij` of `P_r(Q)`. At `q=10`, the genuinely new
> recursive state has `4,608` entries. The full child cut field is determined
> by `H(T)` and its `55` capacities; `(1)` is inherited, so a rooted row needs
> `55` new integers (`H(T)` plus the other `54` capacities), together with
> the parent value `H(Q)`. For every capacity tensor, Sections 4--5 compute
> all fixed-cardinality first and second moments, exact Johnson-layer
> lattices, and THM-4123 coset floors in `O(N^3)` arithmetic without
> enumerating the `2^N` cuts.

The representation is coefficient-minimal within `(2)`, because

```text
c_ij=F_T({i})+F_T({j})-F_T(empty)-F_T({i,j}),
F_T({i})-H(T)=sum_(i->j)c_ij.                           (3)
```

Thus the constant and every edge coefficient are uniquely recovered from the
field, while the singleton linear terms are forced by flow conservation.

The typed connection is

```text
source:       rooted quotient endpoint tables (Q,r)
target:       mixed two-ear cut tensor (H(T),c)
map:          double-clone endpoint recurrence and exposed-block convolution
preserved:    every F_T(S), packet, rational floor, layer lattice/coset
destroyed:    individual Hamilton paths and presentation multiplicity
sidecar:      root r and the ordered clone bit a->b
hostiles:     C3 one-defect repair; strong q=7 all-layer gcd six.
```

## 2. The 4,608-entry double-clone endpoint tensor

For `A subseteq U` and `u in A`, define

```text
D^-_r(A,u) = # Hamilton paths on A union {a,b} ending at u,
D^+_r(A,u) = # Hamilton paths on A union {a,b} starting at u.       (4)
```

Let `E_Q(M,v)` and `S_Q(M,v)` be the usual parent ending/starting subset
tables.  The double-clone endpoint values at the clones themselves are not
new recursive state.  They are

```text
E_a(A)=sum_(p in A, p->r) E_Q(A union {r},p),
E_b(A)=E_a(A)+E_Q(A union {r},r),                         (5)

S_b(A)=sum_(p in A, r->p) S_Q(A union {r},p),
S_a(A)=S_b(A)+S_Q(A union {r},r).                         (6)
```

The internal orientation explains the asymmetric extra terms: only `a->b`,
so a double-clone path can end at `b` after ending at `a`, and can start at
`a` before starting at `b`.

For `u in A`, deleting the endpoint gives the recurrences

```text
D^-_r(A,u)
 =sum_(p in A-{u}, p->u)D^-_r(A-{u},p)
  +1_(r->u)[E_a(A-{u})+E_b(A-{u})],                      (7)

D^+_r(A,u)
 =sum_(p in A-{u}, u->p)D^+_r(A-{u},p)
  +1_(u->r)[S_a(A-{u})+S_b(A-{u})].                      (8)
```

These are ordinary endpoint recurrences, with the two possible clone
predecessors/successors written explicitly.  They are initialized by the
empty-set values implicit in `(5)--(6)`:

```text
E_a(empty)=0, E_b(empty)=1,
S_a(empty)=1, S_b(empty)=0.                               (9)
```

Every child subset has zero, one, or two clones.

- Zero clones: use the parent table on `A`.
- One clone: identify that clone with `r` and use the parent table on
  `A union {r}`.
- Two clones: use `(5)--(8)`.

This reconstructs every child endpoint number.  For `|U|=q-1=L`, each
direction in `(4)` has

```text
sum_(A subseteq U)|A| = L*2^(L-1)
```

entries.  Hence the two directions need `2L2^(L-1)`, which is `4,608` at
`L=9`.  The parent ending/starting tables have `2q2^(q-1)=10,240` meaningful
entries and are computed once per quotient, not once per root.

## 3. Capacity reconstruction

For a child subset `M` and a vertex `x` outside it, define

```text
Before_x(M)=1_(M=empty)+sum_(p in M,p->x) E_T(M,p),
After_x(M) =1_(M=empty)+sum_(p in M,x->p) S_T(M,p).       (10)
```

For an arc `x->y`, put `R=V(T)-{x,y}`.  Contracting the good adjacent block
`x,y` and the reversed one-defect block `y,x` gives

```text
c_xy=sum_(L subseteq R) [
       Before_x(L) After_y(R-L)
      +Before_y(L) After_x(R-L)].                         (11)
```

This is exactly THM-4131's independent exposed-block formula.  Equations
`(5)--(11)` therefore reconstruct the whole 55-edge capacity matrix without
a fresh order-eleven endpoint DP.  As an internal invariant, `(11)` returns
`c_ab=2H(Q)` on every row.

## 4. Exact Johnson moments without cut enumeration

Let the child order be `N`, and let its oriented capacity edges be `e`.  Put

```text
C =sum_e c_e,                 Q2=sum_e c_e^2,
O =sum_(e<f, common tail) c_e c_f,
I =sum_(e<f, common head) c_e c_f,
D =sum_(e<f, disjoint) c_e c_f.                           (12)
```

For `1<=m<=N-1`, write

```text
A_m=sum_(|S|=m)F_T(S),       B_m=sum_(|S|=m)F_T(S)^2.
```

With out-of-range binomial coefficients zero, exact double counting gives

```text
A_m=binom(N,m)H+binom(N-2,m-1)C,                          (13)

B_m=binom(N,m)H^2
   +binom(N-2,m-1)(2HC+Q2)
   +2[binom(N-3,m-1)O
      +binom(N-3,m-2)I
      +binom(N-4,m-2)D].                                  (14)
```

Proof: one directed edge is selected by a size-`m` directed cut in
`binom(N-2,m-1)` subsets.  Two disjoint arcs require their two tails in and
two heads out.  Adjacent arcs can be selected simultaneously only when they
have a common tail or common head; a directed two-edge path imposes the
contradictory requirement that its shared vertex be both in and out.  This
gives exactly the four terms in `(14)`.

The rational support floor is then recovered directly as

```text
J_m=(B_m-H A_m)/(A_m-binom(N,m)H).                        (15)
```

The five aggregates in `(12)` are computed in `O(N^2)` from degrees and
squared degrees; in particular

```text
D=(C^2+Q2-sum_i d_i^2)/2.                                (16)
```

## 5. Exact Johnson lattice in O(N^3)

Let

```text
alpha_i=sum_(i->v)c_iv=F_T({i})-H.                        (17)
```

For distinct `i,j`, put `M_ij=V-{i,j}`, `k=m-1`, and

```text
b_u=c_ju-c_iu,                                           (18)
delta_ij(R)=alpha_j-alpha_i-sum_(u in R)b_u,
               R subseteq M_ij, |R|=k.                   (19)
```

Equation `(19)` is exactly the THM-4123 exchange increment.  If `k=0` or
`k=N-2`, its family has one member.  If `1<=k<=N-3`, fix one `k`-subset `R0`.
The Johnson graph on the `k`-subsets of `M_ij` is connected, and changing
`v` to `u` changes the subset sum by `b_u-b_v`.  Therefore

```text
g_ij,m=gcd(delta_ij(R):|R|=k)
      =gcd(delta_ij(R0), b_u-b_v : u,v in M_ij).          (20)
```

Taking

```text
d_m=gcd_(i<j) g_ij,m                                     (21)
```

is the exact Johnson-layer lattice, because the Johnson graph `J(N,m)` is
connected and its edge differences are precisely `(19)`.  A fixed standard
`m`-subset gives the anchor `a_m`; the exact support floor is

```text
L_m=a_m+d_m ceil((J_m-a_m)/d_m),                          (22)
```

with the usual constant-layer convention when `d_m=0`.

For each pair `i<j`, one prefix of the `b`-list and one gcd of its differences
serves every `m`.  Thus `(20)--(22)` compute all layers in `O(N^3)`, with no
`2^N` response vector.  This is the main exact-coset speedup.

### Hostile to a tempting extra shortcut

It is false that pair expansion forces every interior lattice to equal `2`.
All `202,012` rooted labelled presentations through quotient order six pass
that pattern, but the first strong order-seven generator witness is

```text
quotient label: 111100110011111111101
root:           2
child label:    1111001110011111111110100100
H(child):       387
layer gcds:     (6,6,6,6,6,6,6).                         (23)
```

So a fast order-eleven engine must retain `(20)--(22)` rather than replace
exact cosets by ordinary odd rounding.

## 6. Rational packet decomposition and analytic target

The pair structure gives a smaller exact carrier for the rational gate.
Work in the integer capacity scale.  For `u in U`, put

```text
p_uv=c_uv,
e_u=sum_(v in U-u)p_uv,
g_u=sum_(u->v in U)p_uv-sum_(v->u in U)p_uv,

x_u=c_au+c_bu,              y_u=c_au-c_bu,
s_u=+1 if r->u,             s_u=-1 if u->r,
kappa=c_ab=2H(Q).                                         (24)
```

Let

```text
P=sum_(u<v in U)p_uv,       Q=sum_(u<v in U)p_uv^2,
X=sum_u x_u,                Y=sum_u y_u,
Sx=sum_u s_u x_u,           Sy=sum_u s_u y_u,
X2=sum_u x_u^2,             Y2=sum_u y_u^2.              (25)
```

If `C4=4C_hd=sum_i h_i^(c)d_i^(c)` and
`Dcross=4D_4=sum_(disjoint e<f)c_ec_f`, direct expansion gives

```text
C4=sum_u (g_u-s_u x_u)(e_u+x_u)
   +(Sx X+Sy Y)/2+kappa(Sx+Y),                            (26)

Dcross=1/2[(P+X)^2+2kappa P+Q-sum_u(e_u+x_u)^2
           -(X^2+Y^2)/2+(X2+Y2)/2].                      (27)
```

The verifier checks `(26)--(27)` against the literal packet on every small
row.  They expose the exact analytic bottleneck.  The antisymmetric clone
data enter only through

```text
C4=C_symmetric+Y(kappa+Sy/2),
Dcross=D_symmetric+(Y2-Y^2)/4.                            (28)
```

For order eleven, every rational optimizer is central exactly when

```text
2|C4|<Dcross.                                             (29)
```

Thus a prospective analytic proof need not control the full mixed field: it
needs a strong-root inequality bounding the two antisymmetric terms in `(28)`
against the symmetric disjoint-edge reserve.  Exact-coset centrality still
requires the full capacity state and `(20)--(22)`.

## 7. Exact audit and scope

The primary verifier implements two disjoint routes: a literal endpoint DP
on the expanded child and the parent-once double-clone oracle `(5)--(9)`.
It compares the entire capacity matrix, cut field, layer moments, Johnson
lattices, and clone-decomposed packet on every rooted labelled quotient of
orders two through five:

| quotient order | rooted labelled rows |
|---:|---:|
| `2` | `4` |
| `3` | `24` |
| `4` | `256` |
| `5` | `5,120` |
| **total** | **`5,404`** |

It also freezes the directed-triangle hostile `H(P_0(C3))=5`, internal
capacity `6`, the order-seven gcd-six hostile `(23)`, the `10,240/4,608`
endpoint-state counts, and the `55+1` cut-tensor information budget.

Replay it with

```bash
python3 -B 04-computation/tournament_rooted_pair_mixed_twoear_tensor_thm4162.py
python3 -B -O 04-computation/tournament_rooted_pair_mixed_twoear_tensor_thm4162.py
PYTHONHASHSEED=271828 python3 -B 04-computation/tournament_rooted_pair_mixed_twoear_tensor_thm4162.py
```

The independent order-eleven cross-check uses a warning-clean C++ tensor
engine and a separate literal child engine. On nauty's pinned `0/1024`
strong order-ten shard they agree on all `20,580` rooted rows, including the
two semantic accumulators

```text
semantic_sum64=a3a1f225e6a370ee,
semantic_xor64=5ea64e1b98af2f79,                        (30)
```

the sharp-ratio and minimum-margin witnesses, zero failures, and fourteen
non-parity lattice rows. Their source hashes are bound in the front matter;
THM-4163 records the finite census and frozen shard outputs.

This theorem is an all-order algebraic/data-structural instrument. It does
not by itself prove order-eleven centrality, deduplicate pair presentations,
control prime tournaments, or identify actual response maximizers. **QED.**
