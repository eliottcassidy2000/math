---
id: THM-2753
title: "Six-edge parity erasure and three-matching resolvent restoration"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  The faithful S4 action on the six two-subsets of four sheets lands
  in A6, so ambient six-object permutation parity is identically even.  More
  sharply, a transposition and a nontrivial V4 double transposition have the
  same six-edge cycle type 1^2 2^2, although their actions on the three
  perfect matchings are respectively a reflection and the identity.  The
  matching action has kernel exactly V4, image S3, and its sign is the
  original quartic sign.  With a marked three-cycle the two binary choices
  generate A4 -> C3 and S4 -> S3.  The full labelled edge action is faithful,
  and one labelled edge or one mixed binary--ternary word separates the two;
  only unlabeled single-generator cycle/sign data are erased.  This is a
  finite resolvent-information theorem, not a Keller or LRC exclusion.
source: a4-resolvent-next-gate-scout-2026-07-28
depends_on: []
related:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2612-d4-deck-pole-tax-and-depressed-resolvent-gcd-gate
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2743-c2-c3-off-diagonal-projector-rank-and-s3-s4-compatibility-defect
  - THM-2746-c3-quotient-affine-v4-lifts-and-a4-projector-defect
script: 04-computation/s4_six_edge_matching_parity_thm2753.py
output: 05-knowledge/results/s4_six_edge_matching_parity_thm2753.out
script_sha256: 702111f146c2791f283b4c363613d37da8bfa04514174917cd5dfd25127a9d39
output_sha256: dfac9a8ba9248b1459fedb8479c719199567548a111ee5c7fc61d0ebb1c6f501
hash_basis: LF-normalized bytes
---

# THM-2753 -- six edges forget the parity that three matchings restore

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

The quartic `V4` torsor has two natural finite shadows.  Its six unordered
pairs are the edges of `K4`; its three nonzero translation directions are the
three perfect matchings

```text
12|34,                     13|24,                     14|23.       (1)
```

THM-2598 identifies `(1)` with the cubic-resolvent quotient
`S4/V4=S3`.  THM-2743 and THM-2746 identify the two live marked binary
generators as a transposition and a nonzero `V4` translation.  The theorem
below makes the intervening information loss exact:

```text
six-edge ambient sign and one-generator cycle type: parity lost;
three-matching quotient:                              parity restored. (2)
```

The qualification in `(2)` is load-bearing.  The full labelled six-edge
action is faithful.  What fails is the coarser unlabeled single-generator
cycle/sign statistic, not all six-edge data.

## 1. The two canonical actions

Let `X={1,2,3,4}` and put

```text
E=binom(X,2),                         |E|=6,
M={12|34,13|24,14|23},                |M|=3.             (3)
```

Relabelling the four sheets gives homomorphisms

```text
epsilon:S4 -> Sym(E)=S6,
mu:S4      -> Sym(M)=S3.                                  (4)
```

The edge action `epsilon` is faithful.  Indeed, if a permutation fixes every
two-subset, it fixes the three-edge star at every vertex and hence fixes every
vertex.  The matching action is the standard quotient

```text
ker(mu)=V4={1,(12)(34),(13)(24),(14)(23)},
im(mu)=S3.                                                 (5)
```

Thus `(4)` contains two different forgetful operations.  Passing from the
labelled edge action to its `S6` conjugacy class forgets how opposite edges
are paired.  Passing to `M` deliberately retains exactly that opposition
partition and forgets the `V4` translation.

## 2. Ambient six-edge parity is identically even

For a permutation `g` of an `n`-set, its action on the `k`-subsets has sign

```text
sgn(g on k-subsets)=sgn(g)^binom(n-2,k-1).                (6)
```

It is enough to check `(6)` on a transposition: precisely
`binom(n-2,k-1)` pairs of `k`-subsets are exchanged, and transpositions
generate the symmetric group.  At `(n,k)=(4,2)`, the exponent is `2`.
Therefore

```text
epsilon(S4) subset A6,
sgn_S6(epsilon(g))=+1                  for every g in S4. (7)
```

So the ambient sign of the six-edge permutation cannot be the quartic sign.
This does not contradict faithfulness: the original sign is still a
character of the labelled subgroup `epsilon(S4)`, but it is not the
restriction of the ambient `S6` sign.

The complete conjugacy-class table is

| sheet cycle | size | quartic sign | six-edge cycle | matching cycle |
|---|---:|---:|---|---|
| `1^4` | 1 | `+` | `1^6` | `1^3` |
| `2 1^2` | 6 | `-` | `2^2 1^2` | `2 1` |
| `2^2` | 3 | `+` | `2^2 1^2` | `1^3` |
| `3 1` | 8 | `+` | `3^2` | `3` |
| `4` | 6 | `-` | `4 2` | `2 1` |                 (8)

In particular a transposition and a double transposition become conjugate in
`S6`: both have cycle type `1^2 2^2` on `E`.  Hence every statistic of one
generator that factors through its unlabeled six-edge cycle partition has
the same value on these two quartic classes.  This includes fixed-point
counts, traces of all powers, the characteristic polynomial of the
permutation matrix, and the ambient permutation sign.

## 3. The matching quotient restores exactly the lost sign

The last column of `(8)` separates the collision:

```text
transposition          -> reflection of M, cycle type 1+2;
double transposition   -> identity on M.                     (9)
```

Moreover both maps

```text
sgn_S4:S4 -> C2,
sgn_S3 o mu:S4 -> C2                                      (10)
```

send a sheet transposition to `-1`.  Since transpositions generate `S4`,

```text
sgn_S4 = sgn_S3 o mu.                                    (11)
```

Equations `(7)` and `(11)` are the exact two-level signature:

```text
S4 --epsilon--> A6 --ambient sign--> 1,
 |
 mu
 v
S3 ------------sign----------------> C2.                 (12)
```

There is no function from the six-edge **cycle type** to the matching cycle
type, because the common edge type `1^2 2^2` in `(8)` has both matching
types `1^3` and `1+2`.  Retaining the opposition partition of the octahedron
`L(K4)` is precisely the ternary sidecar that makes `mu` available.

## 4. The marked binary--ternary branches are A4/C3 and S4/S3

Use zero-based labels and set

```text
delta=(01)(23),               tau=(03),               r=(012). (13)
```

Then exact closure gives

```text
<delta,r>=A4,                 mu(<delta,r>)=C3,
<tau,r>=S4,                   mu(<tau,r>)=S3.           (14)
```

Thus the same marked ternary rotation has two binary completions.  On the
six edges the binary generators alone both have type `1^2 2^2`; on the three
matching directions they are identity versus reflection.  This is the
finite representation behind the `A4/C3` row of THM-2746 and the `S4/S3`
row of THM-2743.

It also explains the discriminant character without overclaiming a transfer.
For a separable quartic, the discriminant square class is the sign character
of its sheet monodromy.  For the cubic resolvent, it is the sign character of
the matching action.  Equation `(11)` makes those characters equal.  The
stronger exact polynomial identity

```text
disc(quartic)=disc(integral resolvent cubic)             (15)
```

is THM-2598; it is not reproved from cycle types here.

## 5. Sharp boundary: labelled and joint edge data can distinguish

The edge representation itself loses nothing.  For the choices `(13)`, the
labelled edge `02` satisfies

```text
tau(02)=23,                        delta(02)=13.          (16)
```

Even without selecting one edge, a mixed word with the common ternary move
separates the branches:

```text
epsilon(tau r)   has cycle type 4+2,
epsilon(delta r) has cycle type 3+3.                     (17)
```

So the no-go is exactly one-generator and unlabeled.  A labelled edge image,
the octahedral opposition relation, or joint binary--ternary word data can
restore the distinction.  In particular:

- a six-vertex tournament orientation is additional data, not a canonical
  consequence of the edge set;
- no theorem says every tournament statistic factors through `(8)`;
- the matching quotient is a canonical small sidecar, not the only possible
  separating refinement.

This boundary is useful for the repo's octahedral `L(K4)` carriers.  Scalar
cycle data of one move cannot distinguish a `V4` translation from a quartic
reflection.  A proof-facing use must retain the opposite-edge/matching
colour, a labelled current, or a mixed-word observable.

## 6. Exact verification

Run

```bash
python 04-computation/s4_six_edge_matching_parity_thm2753.py
python -O 04-computation/s4_six_edge_matching_parity_thm2753.py
```

Both executions byte-match the stored `17`-line transcript
`05-knowledge/results/s4_six_edge_matching_parity_thm2753.out`.  The companion
uses explicit exceptions and no truth-bearing Python assertions.  It
enumerates all `24` sheet permutations, proves the edge action faithful,
checks its image lies in `A6`, computes the complete table `(8)`, identifies
the matching kernel and sign pullback, closes both marked groups in `(14)`,
and verifies the labelled and mixed-word hostiles `(16)--(17)`.

## 7. Boundary ledger

```text
PROVED HERE (candidate):  complete S4 edge/matching action table;
                          faithful six-edge action with trivial ambient sign;
                          transposition/double-transposition cycle collision;
                          matching kernel V4 and exact sign restoration;
                          marked A4/C3 versus S4/S3 images;
                          labelled and mixed-word sharp hostiles.

NOT PROVED:               a graph-quartic polynomial realization;
                          a new discriminant-divisor or Jelonek constraint;
                          exclusion of A4 or S4 Keller monodromy;
                          equivalence with a six-vertex tournament;
                          a physical LRC octahedral current or endpoint map;
                          JC(2), DC(2), or LRC(14).                      (18)
```

QED (candidate).
