# The SOP2--SOP3 typed-incidence compiler and the Rule 30 trace boundary

> **Status and scope.** This is a self-contained combinatorial repackaging of
> the proof architecture in Artem Chernikov, *SOP2 = SOP3*,
> [arXiv:2608.13291v1](https://arxiv.org/abs/2608.13291), not a claim of
> literature novelty.  The external equality `SOP2=SOP3` remains a **CITED
> VERY RECENT PREPRINT**.  The abstract theorem below is conditional only on
> its displayed set-system hypotheses and has its own elementary proof.  Its
> comparison with THM-3456 concerns the deliberately enriched free-input
> trace-incidence structure.  It proves no statement about Rule 30's named
> single seed and no LRC(14) statement.

## 1. Inheritance and the exact question

The closest proved repository mechanism is
[THM-3395](../01-canon/theorems/THM-3395-small-sheet-typed-cover-star-cochain.md):
an owner must retain its blocked coset and affine star gap until the common
source phase is glued.  Its canonical hostile is the `q=6` family
`(2,8,14)` with labels `(0,1,2)`: every pair gap fibre is nonempty, but the
compatible star fibre is empty.  The repaired near miss is therefore
quotienting a common-witness problem to bare vertices or pair feasibility.

Chernikov's proof uses the same discipline in a different category.  Its
vertices are witness--parameter pairs, and its load-bearing predicate is the
emptiness of an *oriented* shared-witness fibre.  The paper then splits on a
mixed partial type: either every such demand has a completion, or one failure
has a finite obstruction which tree homogeneity replicates.  Both branches
produce the same typed cross-ladder.

The question addressed here is whether this mechanism has a precise meaning
outside first-order model theory.  It does: it is a theorem about compact
homogeneous incidence set systems.  Compactness and homogeneity are genuine
hypotheses, however, and the theorem cannot be applied to a finite scalar
quotient or to one distinguished branch merely because the surrounding free
system has many branches.

## 2. The universal typed-fibre lemma

Let `W` be a witness set and `P` a packet set.  For every `p in P`, let
`A_p,B_p` be subsets of `W`.  Call a sequence

```text
z_i=(w_i,p_i),                   i<omega,              (1)
```

a **typed cross-ladder** when, for every `i<j`,

```text
w_i in A_(p_j),
w_j in B_(p_i),
A_(p_i) intersect B_(p_j) = empty.                    (2)
```

Define a directed relation on `W x P` by

```text
(w,p) R (w',q)
iff
  A_p intersect B_q = empty,
  w in A_q,
  w' in B_p.                                          (3)
```

> **Typed-fibre triangle lemma.** Every typed cross-ladder satisfies
> `z_i R z_j` for `i<j`, and `R` has no directed three-cycle.

The forward-edge assertion is exactly `(2)`.  If

```text
z_0 R z_1,       z_1 R z_2,       z_2 R z_0,          (4)
```

then the first edge puts `w_0 in A_(p_1)`, the second says
`A_(p_1) intersect B_(p_2)=empty`, and the third puts
`w_0 in B_(p_2)`, a contradiction.  This is the set-system content of
Chernikov's Lemma 2.3.  No compactness, tree, language, or model is used in
this final compiler.

The typing in `(1)` is essential.  As a smallest hostile, take

```text
W=P={0,1,2},                  A_i=B_i={i}.              (5)
```

The bare packet relation

```text
A_p intersect B_q = empty                             (6)
```

holds for every `p!=q`, so its directed version contains a three-cycle (and
both orientations of every edge).  The typed relation `(3)` remains
three-cycle-free.  The witness coordinate is not cosmetic decoration on the
packet graph.

## 3. Compact homogeneous incidence systems

Let `X,Y` be sets and let

```text
E subset X x Y                                        (7)
```

be an incidence relation.  Put

```text
N_x={y in Y:E(x,y)},
u # v  iff there is no x in X with E(x,u) and E(x,v),
D_u={v in Y:u # v}.                                   (8)
```

Thus `#` is the exact incompatibility relation induced by common
`X`-witnesses.

Suppose there are labels

```text
x_sigma in X,       sigma in omega^omega,
y_eta   in Y,       eta   in omega^(<omega),           (9)
```

satisfying the following four conditions.

1. **Branch incidence.** If `eta` is a proper initial segment of `sigma`,
   then `E(x_sigma,y_eta)`.

2. **Incomparable incompatibility.** If the finite nodes `eta,nu` are
   incomparable, then `y_eta # y_nu`.

3. **Mixed-fibre homogeneity.** For finite lists of leaves `sigma_t` and
   nodes `eta_l`, emptiness of

   ```text
   F((sigma_t);(eta_l))
     = intersection_t N_(x_(sigma_t))
       intersect intersection_l D_(y_(eta_l))          (10)
   ```

   depends only on the quantifier-free ordered meet-tree type of the combined
   leaf/node index tuple in
   `(omega^(<=omega); initial-segment, meet, lex-order, leaf-predicate)`.

4. **Countable compactness of the constraint family.** Every countable
   subfamily of the sets

   ```text
   {N_(x_sigma)} union {D_(y_eta)}                      (11)
   ```

   with the finite-intersection property has nonempty total intersection.

Condition 3 is deliberately only the fragment of treetop indiscernibility
used by the proof: it transports emptiness of a finite mixed fibre.  Condition
4 is the corresponding fragment of saturation/compactness: an empty
countable mixed fibre already has a finite empty subfibre.

## 4. The completion--finite-obstruction compiler

> **Compact homogeneous-tree compiler.** Every incidence system satisfying
> Conditions 1--4 produces a typed cross-ladder `(1)` and hence a relation
> `(3)` with an infinite complete forward chain and no directed three-cycle.

Set

```text
lambda_t = 0^t 1 0^omega,
c_t      = x_(lambda_t),
s_j      = y_(0^j),                                    (12)
```

and for `i<omega` consider the mixed fibre

```text
Gamma_i
  = intersection_(t<=i) N_(c_t)
    intersect intersection_(j>i) D_(s_j).              (13)
```

The proof splits exhaustively.

### 4.1 Every mixed fibre has a completion

Suppose `Gamma_i` is nonempty for every `i`, and choose `v_i in Gamma_i`.
Use packets

```text
p_i=(s_i,c_i),
A_(s,c)=D_s,
B_(s,c)=N_c.                                           (14)
```

If `i<j`, equation `(13)` gives

```text
v_i in D_(s_j)=A_(p_j),
v_j in N_(c_i)=B_(p_i).                                (15)
```

Also `0^i` is an ancestor of `lambda_j`, so branch incidence gives
`E(c_j,s_i)`.  No member of `D_(s_i)` can share the witness `c_j` with
`s_i`; hence

```text
A_(p_i) intersect B_(p_j)
  =D_(s_i) intersect N_(c_j)
  =empty.                                              (16)
```

Thus `(v_i,p_i)` is a typed cross-ladder.

### 4.2 One mixed fibre fails

Suppose `Gamma_m` is empty.  Countable compactness supplies a finite negative
subfamily witnessing the failure.  Adding all positive constraints and any
omitted intervening negative constraints preserves emptiness, so there is a
`k>=1` such that

```text
intersection_(t<=m) N_(c_t)
intersect
intersection_(l<k) D_(s_(m+1+l))
=empty.                                                (17)
```

The positive part alone is nonempty: `s_0` belongs to every `N_(c_t)`, so at
least one negative constraint is genuinely required.

For `a<omega`, define the purely combinatorial block indices

```text
mu_a       = 0^((m+1)a),
lambda^a_t = 0^((m+1)a+t) 1 0^omega,       t<=m,
rho^a_l    = 0^((m+1)a+m) 1 0^l,           l<k.        (18)
```

Let

```text
v_a=y_(mu_a),
p_a=((x_(lambda^a_t))_(t<=m),(y_(rho^a_l))_(l<k)),
A_(p_a)=intersection_(t<=m) N_(x_(lambda^a_t)),
B_(p_a)=intersection_(l<k) D_(y_(rho^a_l)).            (19)
```

For `a<b`, the node `mu_a` is an ancestor of every `lambda^b_t`, while every
`rho^a_l` is incomparable with `mu_b`.  Conditions 1 and 2 therefore give

```text
v_a in A_(p_b),                 v_b in B_(p_a).         (20)
```

The combined tuple of leaves `(lambda^a_t)` and nodes `(rho^b_l)` has the same
ordered meet-tree type as the base tuple of `(lambda_t)` and
`(0^(m+1+l))`.  Condition 3 transports `(17)` and gives

```text
A_(p_a) intersect B_(p_b)=empty.                       (21)
```

Again `(v_a,p_a)` is a typed cross-ladder, completing the proof.

This theorem is exactly the mixed-partial-type architecture of Chernikov's
Theorem 3.2 with the logical language removed.  It is not a second proof of
the paper's full model-theoretic theorem: obtaining Conditions 3 and 4 from an
arbitrary `SOP2` formula uses the cited treetop-indiscernible modeling and
monster-model machinery.

## 5. Audit of the cited proof and its sharp assumptions

The complete five-page v1 PDF was read, including the preliminaries, the
correction in Remark 2.7, Lemma 3.1, and both cases of Theorem 3.2.  The PDF
retrieved on 2026-08-15 has SHA-256

```text
c5588d89df72ff7bacab32434a85624ec7b46d9fe8a14a60f861c7e4079d7b56.  (22)
```

The proof's type flow is sound:

- Lemma 2.3 uses two polarized packet predicates `alpha` and `beta`; the
  relation retains both the witness and its packet.
- In Case 1, `Gamma_i` supplies the two cross-incidences, while the old branch
  witness `c_j` proves the required shared-fibre emptiness.
- In Case 2, compactness extracts a genuinely finite mixed obstruction;
  Lemma 3.1's explicit blocks preserve every meet, ancestor relation,
  lexicographic cut, and leaf designation needed to transport it.
- The final three-cycle contradiction uses the witness coordinate from the
  first and third edges against the empty fibre certified by the second edge.

Two boundaries are load-bearing.

First, countable compactness cannot be silently omitted.  The tails

```text
{n,n+1,n+2,...} subset N                              (23)
```

have nonempty finite intersections and empty total intersection.  Without
compactness there can be neither a global completion nor a finite obstruction
packet.

Second, homogeneity of the nonleaf skeleton is insufficient.  Chernikov's
Remark 2.7 records the gap in the earlier modeling proof: a copy of a
meet-closed nonleaf skeleton need not admit leaves in the prescribed
lexicographic cuts.  The paper repairs this by embedding the skeleton with
large gaps between successor directions before inserting leaves.  Condition
3 above is intentionally stated for the full mixed leaf/node tuple, so the
abstraction does not repeat that error.

Repository search found no prior standalone statement of `(1)`--`(21)` and
no direct `MISTAKE-*` entry for Chernikov's dichotomy.  This supports only the
description **new standalone repository abstraction**.  It is not a priority
or literature-novelty claim.

## 6. THM-3456: a fixed-formula SOP2 witness

THM-3456 studies a radius-one cellular automaton `F` on a finite alphabet
which is permutive in its left input.  Its triangular compiler gives the
homeomorphic coordinate system

```text
Theta_F(x)=(tau_F(x),x_(>0)),                          (24)
```

where `tau_F(x)` is the center trace.  In particular every infinite trace and
every finite trace cylinder occurs under free input.

The model-theoretic claim belongs to an explicitly enriched two-sorted
structure

```text
M_F=(configuration sort A^Z,
     finite-word sort A^(<omega),
     E(x;w) iff w is a prefix of tau_F(x)).             (25)
```

The word sort and the single incidence relation are essential.  If one used
one separate predicate for every concrete word, there would be a schema of
formulas rather than the one fixed formula required by `SOP2`.

In `(25)`, however, the fixed formula `E(x;y)` genuinely witnesses `SOP2`:

1. choose two alphabet letters and parameters `b_eta=eta` for
   `eta in 2^(<omega)`;
2. for every branch `sigma in 2^omega`, trace surjectivity from `(24)` gives a
   configuration `x_sigma` satisfying all
   `E(x_sigma;sigma restricted to n)`;
3. incomparable words cannot both be prefixes of one trace, so the
   corresponding pair of instances is inconsistent.

This argument occurs already in the standard structure and uses one fixed
formula.  It is stronger evidence than an analogy with a tree-shaped picture.

THM-3456 also has a direct `SOP3` certificate.  On the word sort define

```text
Q(u;v) iff
  every configuration incident to v is incident to u,
  and some configuration is incident to u but not v.  (26)
```

This first-order formula is strict inclusion of trace cylinders.  The words
`0^i` form an infinite forward chain under the appropriate orientation, and
strict inclusion has no directed cycle.  Thus the candidate theorem does not
need the external equality `SOP2=SOP3` for its own `SOP3` conclusion.
Chernikov supplies context and the more general compiler for situations where
no nested-cylinder order is definable.

### The saturation boundary

One should not claim that the concrete finite-word sort in `(25)` itself
satisfies Condition 4 of Section 3.  It does not.  Fix the trace `0^omega` and
let `c` be a configuration with that trace.  For

```text
eta_n=0^n 1,                                           (27)
```

every finite family

```text
N_c intersect intersection_(n<=N) D_(eta_n)           (28)
```

contains the longer prefix `0^(N+2)`, but the intersection over all `n` is
empty: every finite prefix `0^k` is comparable with `eta_n` once `n>=k`.

To run the compact homogeneous-tree compiler literally, pass to a sufficiently
saturated elementary extension of `M_F` and use the cited treetop-indiscernible
construction.  The direct fixed-formula `SOP2` and strict-cylinder `SOP3`
arguments in THM-3456 need no such passage.  Keeping these two proofs separate
prevents a false compactness claim about the standard word set.

## 7. Free-input theory versus the named Rule 30 seed

Equation `(24)` says that the complete trace together with the positive
initial boundary is a faithful coordinate system.  Projection to the trace
alone forgets which positive boundary was chosen and therefore forgets the
unique inverse-compiled left boundary

```text
beta_F(y,r).                                           (29)
```

The distinguished single seed is not selected by the free trace language.
For Rule 30 it is selected by the conjunction

```text
y_0=1,              r=0^omega,
beta_F(y,r)=0^omega.                                  (30)
```

All left-permutive rules on the same alphabet have isomorphic structures
`(25)`: the boundary recoding `Theta_G^(-1) Theta_F` preserves the complete
trace, and the word sort can be fixed.  This isomorphism does **not** preserve
the cellular-automaton map, the named seed, or its inverse-boundary predicate.
Rule 60 is the sharp hostile: it has the same free trace language, uniform
free-input cylinder law, and `SOP2`/`SOP3` incidence theory, while its
single-one seed has constant center trace.

There are three related but distinct statements:

```text
M_F has SOP2;                                           PROVED from free traces
(M_F,delta_0) has SOP2;                                true via its reduct
the branch incident to delta_0 witnesses SOP2;         FALSE/TYPED NO-GO.  (31)
```

Naming `delta_0` does not erase the tree elsewhere in the structure, but it
also does not make the one named trace branch into a tree.  Its incident
finite words form a single chain.  Consequently neither `SOP2`, `SOP3`, nor
uniform Bernoulli balance of free inputs implies absence of eventual
periodicity, limiting balance, or a computational lower bound for that branch.

A finite scalar universe has the same obstruction even earlier.  An infinite
forward sequence in a finite three-cycle-free target cannot repeat, because a
repetition would force a loop.  Thus the finite LRC row quotient and any fixed
finite trace table cannot themselves realize the compiler; an infinite lifted
witness/parameter fibre would be mandatory.

## 8. Connection contract and non-consequences

| field | exact content |
|---|---|
| source | compact homogeneous incidence tree, or after saturation the `SOP2` trace-cylinder tree |
| target | directed relation on typed `(witness, constraint packet)` vertices |
| map | completion packets `(s_i,c_i)` or replicated finite positive/negative obstruction packets, followed by `(3)` |
| preserved | two oriented cross-incidences and exact emptiness of the shared-witness fibre |
| destroyed | ancestry, meets, lexicographic position, internal cause of emptiness, and any distinguished-branch label |
| required sidecar | witness coordinate, packet identity, `A/B` polarity, full leaf/node homogeneity, compactness; for the Rule 30 seed, `r` and `beta_F` |
| cheapest abstract hostile | `(5)`: bare packet projection has directed triangles |
| cheapest trace hostile | Rule 30 versus Rule 60: identical free trace-incidence theory, radically different named-seed center behavior |
| closest repository hostile | THM-3395's `q=6` pair-feasible but nonclosing affine star |

This comparison creates no map from Rule 30 bits to an LRC phase, owner,
physical time, or current.  THM-3395 and THM-3456 share only the proved
methodological rule: existential fibres must remain typed until their common
witness or inverse boundary has been glued.  Neither the cited preprint nor
the abstract compiler gives a refined-ledger decrement or any LRC(14)
consequence.

The honest frontier after this abstraction is correspondingly sharp:

- the external `SOP2=SOP3` result remains cited pending the ordinary fate of a
  very recent v1 preprint;
- THM-3456's enriched trace-cylinder statement must be audited as precisely
  that structure, with the fixed formula and language declared;
- every Rule 30 prize question concerns the special solution of the boundary
  equation `(30)`, not the free-input trace theory.
