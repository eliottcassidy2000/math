---
id: THM-4099
title: "Squarefree gap transfer and the exact mixed-insertion boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For every finite tournament
  split into a nonempty base and an arbitrary finite inserted set, a squarefree
  gap polynomial has coefficient H(B union S) on every inserted face S. For
  two inserted vertices this is a four-state transfer matrix; degree r sees
  only base words with at most r bad old adjacencies. The sequential
  specialization is an exact derivative plus one-defect-repair formula.
  Proper-face and first-step count profiles do not determine the mixed
  coefficient, minimally at the transitive-triple/C3 boundary.
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness
related:
  - THM-012b-insertion-decomposition
  - THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
  - THM-4111-uniform-ear-average-and-recursive-selected-bank-growth
  - THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature
  - THM-4118-ear-response-lattice-and-stateful-unit-component-intervals
script: 04-computation/tournament_squarefree_gap_transfer_thm4099.py
output: 05-knowledge/results/tournament_squarefree_gap_transfer_thm4099.out
independent_audit_script: 04-computation/tournament_squarefree_gap_transfer_thm4099_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_squarefree_gap_transfer_thm4099_independent_audit.out
script_sha256: 533b821ba8401958565d5b76272315d554c96d36a08cff825ab62acbd5ee175a
output_sha256: 3aa4c9ce323fae69ac7fe415de48ffb126e8a4297fbaf88f994082de950480b1
independent_audit_script_sha256: e4e39b4c384cc7d19d2dc467650738f67f6c9ea60c823a98ae7c8da87d73cad2
independent_audit_output_sha256: 6c1acecef3f850efd10b2bc72e8904fc4be5596e90c2b36f9a93859e82dff19c
semantic_sha256: 405a8be93b12d59d1750e8caf4b90bef8555de797dc3badfa37c709de50a603b
independent_semantic_sha256: 2d9bb86849946e3402e6fd2896bbf953f1c786f6999f716c8384f5360d1990b8
hash_basis: raw LF bytes for files; canonical compact JSON for the semantic ledger
audit: >
  The primary referee compares closed local factors with direct block
  enumeration, matrix contraction with squarefree multiplication, and all
  four face coefficients with an independent endpoint-subset dynamic program
  on every labeled final tournament of orders 3 through 6. It also checks the
  sequential witness partition, deletion fibres, live defect layers, and the
  sharp hostile. The independent path checks every split face through order
  five, all 32,768 labelled two-insertion tournaments at order six, 3,145,728
  cutoff coefficients, the sequential identity, and both minimal boundaries.
  Normal and python -O transcripts are byte-identical. The audit repaired an
  overstatement about complete labelled fibres: only their count-profile
  quotient is proved insufficient.
---

# THM-4099 -- squarefree gap transfer and mixed insertion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem supplies the compositional sidecar requested by the boundary in
`THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness`. The
correct state is not a selected matching or a list of response counts. It is
the squarefree polynomial of actual inserted blocks in the gaps of actual
base words. For two insertions it has only four states.

The main bijection is proved for an arbitrary finite inserted set. The empty
base is excluded to avoid a separate boundaryless-gap convention; all
successive-insertion and hostile statements here have nonempty base.

## 1. Squarefree gap polynomial

Let `W` be a finite tournament with a fixed disjoint decomposition

```text
V(W)=B disjoint_union I,    B nonempty.                         (1)
```

For `S subseteq I`, put `X_S=product_(s in S) X_s` and work in the
squarefree commutative ring

```text
R_I=Z[X_i:i in I]/(X_i^2:i in I).                              (2)
```

Fix a permutation `P=(x_1,...,x_n)` of `B`. Its gaps are

```text
g_0=(-,x_1),  g_j=(x_j,x_(j+1)) for 1<=j<n,  g_n=(x_n,-).       (3)
```

For a gap `g` and `S subseteq I`, define `kappa_g(S)` as follows.

- At an internal gap `(a,b)`, it is the number of permutations
  `(s_1,...,s_r)` of `S` satisfying

  ```text
  a -> s_1 -> ... -> s_r -> b.                                 (4)
  ```

  For `S=emptyset`, this coefficient is `1[a->b]`.

- At the left boundary `(-,b)`, it is the number satisfying

  ```text
  s_1 -> ... -> s_r -> b.                                      (5)
  ```

- At the right boundary `(a,-)`, it is the number satisfying

  ```text
  a -> s_1 -> ... -> s_r.                                      (6)
  ```

  Both boundary coefficients are `1` for `S=emptyset`.

Define the local gap factor and global relative gap polynomial by

```text
G_g(X)=sum_(S subseteq I) kappa_g(S) X_S,

Z_(W/B)(X)=sum_(P in Sym(B)) product_(g in gaps(P)) G_g(X) in R_I.  (7)
```

> **Theorem 1.1 (squarefree gap-transfer bijection).** For every
> `S subseteq I`,
>
> ```text
> [X_S] Z_(W/B)=H(W[B union S]).                                (8)
> ```

### Proof

Expand one product in `(7)`. Since `X_i^2=0`, a term contributing to `X_S`
chooses pairwise disjoint subsets `S_g` from the labeled gaps, with

```text
S=disjoint_union_g S_g,                                        (9)
```

and chooses one locally directed ordering counted by `kappa_g(S_g)` in each
gap. Concatenating those ordered blocks with `P` produces a directed
Hamiltonian path of `W[B union S]`.

Conversely, delete `S` from any directed Hamiltonian path of `W[B union S]`.
The surviving vertices give one uniquely ordered base permutation `P`. The
maximal deleted blocks give a unique allocation `(S_g)_g` to its gaps and a
unique locally directed ordering in each gap. These constructions are
inverse. Thus the coefficient counts the Hamiltonian paths exactly. QED.

## 2. The four-state two-insertion transfer

Take `I={u,v}` and write `U=X_u`, `V=X_v`. Every local factor has the form

```text
G_g=a+bU+cV+dUV.                                               (10)
```

For an internal gap `(x,y)`, its coefficients are

```text
a=1[x->y],
b=1[x->u] 1[u->y],
c=1[x->v] 1[v->y],
d=1[x->u] 1[u->v] 1[v->y]
  +1[x->v] 1[v->u] 1[u->y].                                  (11)
```

At the left boundary `(-,y)`,

```text
a=1,
b=1[u->y],
c=1[v->y],
d=1[u->v] 1[v->y]+1[v->u] 1[u->y].                            (12)
```

At the right boundary `(x,-)`,

```text
a=1,
b=1[x->u],
c=1[x->v],
d=1[x->u] 1[u->v]+1[x->v] 1[v->u].                            (13)
```

Multiplication by `(10)` on coefficient columns in the ordered basis
`(1,U,V,UV)` is

```text
        [ a  0  0  0 ]
M_g  =  [ b  a  0  0 ].                                      (14)
        [ c  0  a  0 ]
        [ d  c  b  a ]
```

Consequently

```text
Z_(W/B)=H(W[B])
       +H(W[B union {u}]) U
       +H(W[B union {v}]) V
       +H(W) UV,                                               (15)
```

and the requested contraction is

```text
H(W)=sum_(P in Sym(B))
       e_UV^T (product_(g in gaps(P)) M_g) e_0.                 (16)
```

The matrices in `(16)` represent multiplication in a commutative ring, so
the coefficient is independent of the order chosen for the gap product.

## 3. Defect depth is polynomial degree

For a base permutation `P=(x_1,...,x_n)`, let

```text
delta(P)=#{j:x_(j+1)->x_j}.                                   (17)
```

> **Corollary 3.1.** The contribution of `P` to a monomial of total degree
> `r` is zero unless `delta(P)<=r`.

Indeed, each bad old adjacency is an internal gap whose local factor has zero
constant coefficient. A nonzero product term must assign at least one
inserted vertex to every such gap. The assigned subsets are disjoint, so their
total cardinality is at least `delta(P)`.

For two insertions, only the zero-, one-, and two-defect base layers can
contribute. A two-defect word can contribute only to `UV`, and every word with
at least three defects contributes zero. This is the exact finite filtration
behind the four-state transfer, not a heuristic truncation.

## 4. Sequential derivative and repaired failed cuts

Let `A` be a nonempty tournament and add a vertex `v`. For
`Q=(q_1,...,q_m) in Ham(A)`, set

```text
s_i=1[v->q_i],
tau_v(Q)=#{i:s_i=1 and s_(i+1)=0}.                             (18)
```

The legal insertion positions are the initial position when `s_1=1`, the
internal `0->1` transitions, and the final position when `s_m=0`. Binary
transition balance therefore gives

```text
number of legal v-insertions into Q = 1+tau_v(Q).               (19)
```

Let `Def_1(A)` be the permutations `R=(r_1,...,r_m)` of `V(A)` having exactly
one bad adjacency. If its unique bad adjacency is at `i`, put

```text
rho_v(R)=1[r_i->v and v->r_(i+1)],
F_(A,v)=sum_(R in Def_1(A)) rho_v(R),
L_(A,v)(z)=sum_(Q in Ham(A)) z^(tau_v(Q)).                      (20)
```

> **Theorem 4.1 (sequential composition).** One has
>
> ```text
> H(A+v)=L_(A,v)(1)+L'_(A,v)(1)+F_(A,v).                       (21)
> ```

### Proof

Deleting `v` from a Hamiltonian path of `A+v` gives either a Hamiltonian path
`Q` of `A`, or a permutation with exactly one bad adjacency. In the first
case, `(19)` counts the fibre. In the second, `v` was internal and its unique
cut repairs that bad adjacency, so the path is counted once by `rho_v`.
Conversely, every legal insertion and every repaired one-defect word produces
a Hamiltonian path. Their deletion words distinguish them, and the two cases
are disjoint. Summing gives `(21)`. QED.

Equivalently, if

```text
Phi_(A,v)(z,w)=L_(A,v)(z)+w F_(A,v),                            (22)
```

then

```text
H(A+v)=[(1+partial_z+partial_w) Phi_(A,v)] at (z,w)=(1,0).      (23)
```

Now take `A=T+u`. The full deletion partition of `THM-4094` gives

```text
L_(T+u,v)(z)
 =sum_(P in Ham(T)) sum_(Q in D_u(P)) z^(tau_v(Q))
  +sum_(Q in Orphan_u) z^(tau_v(Q)).                            (24)
```

Thus the actual first-step fibres and orphan path words, decorated with their
`v`-signatures, determine `L` exactly. Equation `(24)` accounts only for that
term. Direct evaluation of `(21)` may additionally retain `F`, for example
through the one-defect word bank `Def_1(T+u)` and its `v`-repair bits; no
insufficiency claim is made here about the complete labelled fibre system.
Formula `(7)` packages both layers without first scalarizing them.

## 5. Minimal mixed-coefficient hostile

Let the common base be the singleton `T={0}`. Compare two first insertions:

```text
W_tr:   0->u,
W_cyc:  u->0.                                                  (25)
```

In both, the first deletion count profile is

```text
(H(T), left-degree tuple, orphan count, H(T+u))
  =(1,(1),0,1).                                                (26)
```

Now give `v` the same labeled future incidences in both cases:

```text
0->v->u.                                                       (27)
```

Then `W_tr` is the transitive order `0,v,u`, whereas `W_cyc` is the directed
cycle `u,0,v,u`. Direct gap multiplication gives

```text
Z_tr =1+U+V+UV,
Z_cyc=1+U+V+3UV.                                               (28)
```

All three proper-face coefficients and the first-step count profile agree,
but the mixed coefficients are `1` and `3`. Sequentially, the transitive case
has `L(z)=1,F=0`; the cyclic case has `L(z)=z,F=1`. The lost coordinate is the
actual intermediate word and repaired cut, not another scalar response count.

This does not identify the two full labeled fibre systems: their unique
intermediate words are `(0,u)` and `(u,0)`. It proves exactly that the quotient
to count profiles and proper-face marginals is not compositional. The hostile
is minimal because every tournament on two vertices has one Hamiltonian path,
so an empty-base/two-insertion comparison cannot separate mixed counts.

## 6. Exact verification and scope

The primary referee exhausts every labeled tournament with the last two
labels designated `u,v`:

```text
base order       1       2       3        4
final order      3       4       5        6
tournaments      8      64    1024    32768.                    (29)
```

For every tournament it checks all four coefficients of `(15)` against an
independent endpoint-subset dynamic program. It independently enumerates the
local inserted blocks, checks the matrix orientation, checks the defect gate,
and verifies `(21)` as an equality of disjoint constructed witness sets with
the DP total. It also verifies the full first-deletion partition and the
minimal hostile. All gates use `require`; normal and optimized transcripts are
byte-identical.

The independent companion uses a separate permutation/assignment engine. It
checks every inserted face through order five, all order-six two-insertion
splits, `3,145,728` defect-cutoff coefficients, the sequential formula, and
the singleton hostile together with the empty-base minimality boundary.

Reproduce with

```bash
PYTHONHASHSEED=0 python3 \
  04-computation/tournament_squarefree_gap_transfer_thm4099.py
PYTHONHASHSEED=0 python3 -O \
  04-computation/tournament_squarefree_gap_transfer_thm4099.py
PYTHONHASHSEED=0 python3 \
  04-computation/tournament_squarefree_gap_transfer_thm4099_independent_audit.py
PYTHONHASHSEED=0 python3 -O \
  04-computation/tournament_squarefree_gap_transfer_thm4099_independent_audit.py
```

The proved polynomial is a Hamiltonian-path counting law. It does not turn a
nonzero mixed coefficient into chronological ordering, monodromy, or an LRC
gluing certificate; those targets require their own common-carrier and
endpoint/phase sidecars.
