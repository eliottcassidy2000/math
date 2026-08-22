# The frozen U_full bridge needs exactly one joint coordinate in the minimum model

**Status: FINITE-EXACT pre-merge implementation sidecar and minimum-API
obstruction.**  The actual endpoint controls now retain factor and branch
lineage before summation.  The smallest full-marginal ambiguity is proved to
be one checkerboard coordinate, and one fixed joint address bit is necessary
and sufficient to recover it.  A calibrated four-record positive model
recovers the frozen `q_H-q_q5` bridge exactly, while an identical-marginal
hostile has bridge zero.  No map from the actual U_full cells to that joint
coordinate is constructed.  This is not a common-ancestry realization,
physical current, row exclusion, or LRC(14) result.

## 1. Inheritance and the narrowed question

The inherited frozen bridge is

```text
p = 572252886246508880869,

q_H  = (1,0,1),
q_q5 = (1,0,0),

B = A(q_H)-A(q_q5)
  = 389266878372286537904 mod p.                    (1)
```

The preceding exact probe rules out the three cheapest geometric relations:
equality of stored cut segments, equality of maximal cyclic `E` components,
and equality of the two circle points.  Their bridge values are respectively

```text
540653486701996040250,
540653486701996040250,
167726070588785644466,                              (2)
```

not (1).  The cross-component complement of the component relation is
nonzero and restores (1), but it has no ancestry meaning by itself.

The inheritance pass is therefore

```text
closest proved mechanism: THM-2471's finite Boolean ancestry fibre;
positive implementation:  audited THM-2594's linked-node, one-base product;
canonical hostile:        THM-2538's mixed-Haar transportation kernel;
corrected near miss:      MISTAKE-293, one fibre is not one circle point;
least-used sidecar:       E/Q factor, branch, wrap and boundary lineage.
                                                               (3)
```

The concept board has four objects: the actual U_full interval cells, the
THM-2471 sheet/horizon record, the `2x2` transportation kernel, and the
`F_13^3` address Fourier transform.  The anchor is the frozen bridge; the
niche is the minimum kernel theorem; the wildcard is exact boundary-lineage
retention.  Equation (2) means no further unlabelled geometric diagonal is a
good next move.

## 2. Actual factor-labelled cells before endpoint summation

The new extractor mirrors `fast_build_set`, but every surviving interval
retains

```text
(family, factor index, in/out/guard mode, branch, boundary side),
the parent positive branch,
the ordered cut lineage.                              (4)
```

Intersecting scaled `E` with periodic `Q` then retains, for every cell,

```text
(E component, Q component, periodic wrap turn,
 E lineage, Q lineage, left/right boundary provenance).       (5)
```

No adjacent-interval merge occurs before (4)--(5).  On the two cheap exact
controls the labelled projections reproduce the original interval lists
entry for entry and their endpoint contributions reproduce both public
factors:

| control | E cells | Q cells | labelled E--Q cells | overlap length | `AX` | `BY` |
|---|---:|---:|---:|---:|---:|---:|
| `ell=0` | 33,810 | 28,730 | 123,752 | 64,165,660,754,527,020 | 263,333,472,381,374,948,713 | 367,849,901,592,567,656,326 |
| `ell=v2` | 34,560 | 28,730 | 126,429 | 68,024,839,018,709,963 | 523,949,679,655,523,736,011 | 87,230,428,072,980,866,510 |

The corresponding left-cell digests are

```text
36b51bf77c363f9f283aa67925fb65e64b8b642b98c5fc023fd91695a9ee78ea,
76a586d218e659391e667a04765e426b192933330c64ddd5bb8cf3efb819a15f.
                                                               (6)
```

This is a real improvement to the computational object: factor labels are no
longer destroyed at the first interval projection.  It does **not** solve the
typing problem.  The available fields are still only

```text
circle interval, E/Q factor lineage, E/Q component,
periodic wrap turn, boundary provenance.                       (7)
```

There is no inherited function from (7) to

```text
base, root, owner sheet, word sheet, source sheet,
left horizon, right horizon, address.                          (8)
```

In particular, the circle interval in (7) is not an ancestry base, and two
linked THM-2471 nodes need not be one circle point.  Assigning (8) from (7)
without a new theorem would merely rename the geometric diagonal already
refuted in (2).

## 3. The exact `2x2` minimum-coordinate theorem

Let `J=(J_ij)` be a `2x2` joint table over a field of characteristic not two,
and let

```text
r_i = sum_j J_ij,       c_j = sum_i J_ij,       M=sum_i r_i. (9)
```

The row/column marginal map has rank three and exact kernel

```text
ker(marg) = span { C },

C = [[ 1,-1],
     [-1, 1]].                                             (10)
```

The complete `F_13` audit finds exactly the thirteen scalar multiples of
`C`.  Every observable depending only on the left index or only on the right
index vanishes on (10).  The mixed Haar coordinate

```text
h(J)=J_00-J_01-J_10+J_11                                (11)
```

pairs with `C` by `4`, so adjoining (11) raises the rank from three to four.
It reconstructs the table explicitly:

```text
J_00 = (h-M+2r_0+2c_0)/4,
J_01 = r_0-J_00,
J_10 = c_0-J_00,
J_11 = r_1-c_0+J_00.                                  (12)
```

Thus one joint scalar is necessary and sufficient in the smallest nontrivial
full-marginal problem.  This is the sharp local form of THM-2538's
`ker(e_U) tensor ker(e_V)` obstruction.

## 4. One ell-independent address map makes the bridge equal to the Haar coordinate

Use four complete formal records indexed by `(i,j) in {0,1}^2`.  Each record
has explicit fields

```text
(base,root,owner_sheet,word_sheet,source_sheet,
 left_horizon,right_horizon,address,left_factor,right_factor,measure).
                                                               (13)
```

In this minimum schema, `owner_sheet=i`, `word_sheet=j`, and
`source_sheet=i xor j`.  Define one map, before any character is selected,

```text
address(i,j) = q_H   if i xor j=0,
               q_q5 if i xor j=1.                         (14)
```

The implementation audits that the address function has only the record as
argument and no `ell` dependence.  The two labelled Boolean endpoint factors
are multiplied on each record before summing.  For
`ell=(alpha,beta,tau)`, form

```text
gamma_ell = sum_(i,j) J_ij zeta_13^(ell dot address(i,j)).  (15)
```

Normalized inverse `F_13^3` Fourier transform gives

```text
gamma^vee(q_H)-gamma^vee(q_q5)=h(J).                         (16)
```

The independent referee checks the four orthogonality sums directly:
diagonal address pairs give `13^3=2197`, off-diagonal pairs give zero.
Equation (16), unlike a relation selected separately for each `ell`, is the
Fourier transform of one fixed address map.

## 5. Positive and hostile bridge controls

The integer `B` in (1) is divisible by four.  Put

```text
lambda=B/4=97316719593071634476,

J_pos  = [[2lambda,       0],
          [      0,2lambda]],

J_flat = [[lambda,lambda],
          [lambda,lambda]].                                  (17)
```

Both are nonnegative integer measures.  They have identical row and column
margins

```text
(2lambda,2lambda),       total B.                             (18)
```

Every separate-left or separate-right scalar therefore agrees.  But

```text
h(J_pos)=B,             h(J_flat)=0.                          (19)
```

The full character-bank calculation gives

```text
                    q_H                         q_q5       bridge
J_pos   389266878372286537904                    0          B
J_flat  194633439186143268952   194633439186143268952        0.
                                                               (20)
```

The positive and hostile bank digests are respectively

```text
b5ec7b38ed4c44abffa1e760b053f4275281017d2f7667fb44532f02c4969f7f,
d9c2b1e9d5cba6a635c4c75dccd0d0684b17ec3ee05718fabc808a902a95a86c.
                                                               (21)
```

The positive in (17) is calibrated using the already known bridge `B`.
Consequently (20) proves sufficiency and minimality of the API coordinate;
it is not an independent derivation of `B` from the U_full endpoint geometry.
The hostile proves that even nonnegative integer couplings with all factor
margins fixed can change the bridge.

## 6. The minimum lawful record still missing from U_full

A genuine U_full implementation would have to refine (5), not replace it by
the formal records in (13).  Its minimum object is

```text
omega=(base,root,owner_sheet,word_sheet,source_sheet,
       left_horizon,right_horizon,address,
       E_lineage,Q_lineage,wrap_turn),                         (22)

pi_L(omega), pi_R(omega) = linked endpoint nodes,
r(omega) in F_13^3,                                           (23)
```

with one `ell`-independent `r`.  The endpoint response must have the form

```text
gamma_ell = sum_omega mu(omega)
  L(pi_L omega) R(pi_R omega) chi_ell(r(omega)),               (24)
```

where the product occurs before every marginalization.  THM-2594 proves that
this order of operations is meaningful on its own canonical row, with factors
at linked nodes `w_u`, `X_(u,a)`, and `Y_(q,e')` over one outer base.  It does
not supply (22)--(23) for U_full.

At bridge-only resolution, (10)--(16) sharpen the engineering target: it is
not necessary first to reconstruct an arbitrary full coupling.  It is enough
to construct one semantically typed mixed-Haar/address coordinate whose value
is `B`.  Conversely, adding any number of fields that depend only on one
endpoint factor cannot work, because all such observables annihilate (10).

## 7. Connection contract and strict boundary

| field | exact content |
|---|---|
| source | actual pre-merge U_full `E`/`Q` factor cells |
| visible refinement | factor/branch, component, wrap, and boundary lineage |
| target | frozen `q_H-q_q5` joint bridge only |
| quotient | row/column endpoint marginalization |
| exact kernel | one checkerboard line in the minimum `2x2` model |
| restoring sidecar | one joint parity/mixed-Haar address coordinate |
| positive | calibrated nonnegative four-record table, bridge `B` |
| hostile | same margins and total, bridge zero |
| destroyed by current API | pairing, THM-2471 sheets/horizons, common address |
| first unresolved map | actual labelled circle cells to (22)--(23) |
| next decisive test | construct that map and recover only (1); still no `K4` |

This result neither proves nor obstructs the existence of a lawful U_full
ancestry lift.  It excludes the claim that factor-labelled circle geometry or
separate endpoint fields already determine one.  It proves that the smallest
additional datum capable of resolving the local API kernel must be genuinely
joint.

No `K4` was evaluated.  No grouped exact-address coefficient `C(a;X,m)`,
all-unit projector `B(q)`, chronological arrival, physical current, scalar-row
exclusion, U_clock statement, or LRC(14) conclusion follows.

## 8. Reproduction and independent audit

Run from the repository root:

```text
python -B 04-computation/lrc_endpoint_ufull_minimal_joint_address_gate_20260816.py
python -B -O 04-computation/lrc_endpoint_ufull_minimal_joint_address_gate_20260816.py

python -B 04-computation/lrc_endpoint_ufull_minimal_joint_address_referee_20260816.py
python -B -O 04-computation/lrc_endpoint_ufull_minimal_joint_address_referee_20260816.py
```

Normal, optimized, and stored primary outputs agree byte-for-byte.  Normal and
optimized referee outputs agree byte-for-byte.  The referee does not import
the primary.  It exhausts the entire `F_13^4` margin kernel, all left-only and
right-only `F_13` functionals, the four address-orthogonality cases, and the
two full `2197`-character banks; it also audits the record schema,
`ell`-independent address signature, and absence of a `K4` call.

LF SHA-256 hashes are

```text
primary script:  fd5fcfc5f92385806f5ea5e77e854c823d2b3d6e34752c7aba777744201c8055
primary output:  3c68f4b3618abbcc2653eb1b2aee4729e59a57820749396ca757341ee7f91838
referee script:  e33714832d1638a7e0587d9d6f7eb40e1c8ec95b92958882ffa93220cf1c4212
referee output:  e83775507494f737c27fa564fe8a121730f7732cd6909330320597b414054a80
```

The primary and referee semantic digests are respectively

```text
adc0a12f3dbf8c3c626827e4df4c1d5e1deca97b7f5d4e5dfc17238cc2587294,
bce53a4d7f5845f7067c6d0c9ce6cb63f31d477c8818722b199378aee486db52.
```
