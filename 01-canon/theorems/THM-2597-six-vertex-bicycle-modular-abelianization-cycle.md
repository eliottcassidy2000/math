---
id: THM-2597
title: "The first six-vertex bicycle is the modular abelianization cycle"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  In the tile graph K_6 minus its frozen base path, the
  unique nonzero bicycle is delta({1,4,5})=delta({2,3,6}), supported on the
  spanning cycle (1,3,5,2,4,6).  It is K_3,3 minus the frozen perfect
  matching and splits into the two derangement perfect matchings.  Its
  isometric embedding as Q_3 minus an antipodal pair has three Theta
  classes of two edges; a half-turn swaps each pair and a third-turn cycles
  the classes.  These commuting rotations realize C_2 x C_3, hence only
  the abelianization C_6 of PSL_2(Z).  Completing the support to K_3,3
  instead exposes all six S_3 perfect matchings and exactly the order-six
  Pfaffian obstruction of THM-2290.  The tile masks 873 and 150 form one of
  the two n=6 blue self-loops; the support rotations are not tournament
  automorphisms.  C_6 is not graceful, although deleting one edge gives a
  graceful P_6.  No LRC(14), graceful-tree, or general tournament
  consequence is claimed.
source: codex-2026-07-27-six-bicycle
depends_on:
  - THM-2467-bicycle-spaces-of-the-star-flip-split
  - THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete
related:
  - THM-648-blue-selfloops-only-at-even-n
  - THM-2587-deep-root-coshift-incidence-wall-and-theta-selector-no-go
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
external:
  - "Mario Krenn, Xuemei Gu, and Daniel Soltesz, Questions on the Structure of Perfect Matchings Inspired by Quantum Physics, arXiv:1902.06023v2."
script: 04-computation/six_vertex_bicycle_modular_abelianization_thm2597.py
output: 05-knowledge/results/six_vertex_bicycle_modular_abelianization_thm2597.out
script_sha256: 6a728dbb8e16cd84057b0cae4304628f3a66ff3436fb4b6b0ac727d81eebb66e
output_sha256: 9c64b359e40a8b90cb9d715a90b8f3c32cae0b9907cb3cac422fa7bfdc0babf1
hash_basis: normalized repository blobs (LF)
---

# THM-2597 -- the first bicycle is an order-six fork

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**  The bicycle, matching, partial-cube, modular, and graceful
claims have short proofs below.  The tile-mask, isomorphism, Hamiltonian,
automorphism, and blue-self-loop claims are finite-exact and reproduced by
the dependency-free companion.

This closes the open `n=6` identification requested in THM-2467.  It also
separates two mathematically different order-six objects:

```text
PSL_2(Z)_ab = C_2 x C_3 isomorphic to C_6,       commuting rotations;
PSL_2(F_2)  = S_3,                              matching permutations. (1)
```

The bicycle realizes the first.  Completing its support realizes the
second.  Collapsing both to “six states” destroys the decisive distinction.

## 1. The unique nonzero bicycle, explicitly

Let

```text
H=K_6 minus {(1,2),(2,3),(3,4),(4,5),(5,6)}.                (2)
```

This is the tile graph in THM-2467.  For a vertex subset `S`, write
`delta_H(S)` for its cut.  A cut is a bicycle exactly when every vertex has
even degree in it.

The recurrence proof of THM-2467 becomes completely explicit at `n=6`.
Put `x_i=1_[i in S]` and `sigma=sum_i x_i mod 2`.  Its boundary equation and
two-step recurrence give, from the seed `a=x_1`,

```text
(x_1,...,x_6)=(a,a+sigma,a+sigma,a,a,a+sigma).              (3)
```

The terminal and parity equations are automatic.  If `sigma=0`, (3) is
the empty or full cut and gives the zero vector.  If `sigma=1`, the two
complementary witnesses are

```text
S={1,4,5},                 V\S={2,3,6}.                     (4)
```

They define the same cut, so the unique nonzero bicycle is

```text
b=delta_H({1,4,5})=delta_H({2,3,6})
 ={(1,3),(3,5),(5,2),(2,4),(4,6),(6,1)}.                   (5)
```

Every vertex has degree two in (5); the displayed order is a spanning
`C_6`.  In the star-flip notation this is simultaneously

```text
b=star(1)+star(4)+star(5)
 =star(2)+star(3)+star(6)                                  (6)
```

over `F_2`.  Thus it is literally the one-dimensional failure of the
star/cycle transversality at the first nonzero residue of THM-2467.

## 2. A 3-by-2 partial cube

Embed the cycle vertices, in the cyclic order of (5), by

```text
1 -> 000,   3 -> 100,   5 -> 110,
2 -> 111,   4 -> 011,   6 -> 001.                          (7)
```

The image is `Q_3` with the antipodal vertices `010,101` deleted.  Along
the cycle the flipped-coordinate word is

```text
0,1,2,0,1,2.                                               (8)
```

For any two displayed vertices, the shorter cycle arc flips each differing
coordinate once and no equal coordinate unnecessarily; hence graph distance
equals Hamming distance.  This proves that (7) is an isometric embedding,
so `b` is a partial cube.  Its three Djoković--Winkler `Theta` classes are

```text
Theta_0={(1,3),(2,4)},
Theta_1={(3,5),(4,6)},
Theta_2={(5,2),(6,1)}.                                    (9)
```

Let `rho_2` rotate (5) by three positions and `rho_3` by two positions:

```text
rho_2=(1 2)(3 4)(5 6),
rho_3=(1 5 4)(3 2 6).                                    (10)
```

They commute, have orders two and three, and generate all six rotations of
the support.  More finely, `rho_2` swaps the two edges within every class
in (9), while `rho_3` cycles the three classes.  This is the exact `3 x 2`
co-occurrence frame suggested by the partial-cube picture.

The presentation

```text
PSL_2(Z)=<s,c | s^2=c^3=1>                                (11)
```

abelianizes to `C_2 x C_3`, which is cyclic of order six.  Sending `s` to
`rho_2` and `c` to `rho_3` therefore realizes the modular abelianization on
the bicycle.  It does **not** realize the free product: all commutators,
reduced-word history, and Bass--Serre/Farey position have been killed.

## 3. Two derangements, then the S3 completion

Use the bipartition

```text
L={1,4,5},                   R={2,3,6},
P_0={(1,2),(4,3),(5,6)}.                                  (12)
```

Then

```text
b=K_(3,3) minus P_0.                                      (13)
```

Relative to the identity matching `P_0`, its six edges are the two
derangement permutations in `S_3`.  Equivalently, the cycle has exactly
the two alternating perfect matchings

```text
M_+={(1,3),(5,2),(4,6)},
M_-={(3,5),(2,4),(6,1)},
b=M_+ disjoint_union M_-.                                 (14)
```

Each matching in (14) is a transversal of the partial-cube structure: it
uses one edge from every `Theta` class in (9).

Color every edge of `M_+` at both endpoints with one color and every edge
of `M_-` with a second color, all with unit weight.  Since an even cycle has
only its two alternating perfect matchings, its inherited-color response
has exactly the two monochromatic outputs.  This is a literal `d=2,n=6`
instance of the monochromatic perfect-matching construction studied by
Krenn--Gu--Soltész; THM-2290 identifies the selected pair kernel as its
faithful carrier.

Restoring `P_0` changes the answer qualitatively.  The support becomes all
of `K_(3,3)`, whose six perfect matchings are the six bijections `L->R`,
hence are literally `S_3`.  The two 3-cycles in (14), the identity, and the
three transpositions are all now present.  This is exactly the six-term
support used in THM-2290's first universal Pfaffian-gauge obstruction:
the product of the six permutation signs is `-1`, whereas every edge sign
appears twice and contributes product `+1`.  No tournament signing can make
all six hafnian monomials carry one Pfaffian sign.

Thus the exact order-six fork is

```text
C6 support:     two derangement matchings, commuting C2/C3 rotations;
K3,3 completion: all six S3 matchings, noncommuting parity obstruction. (15)
```

The abstract isomorphisms `PSL_2(F_2)=S_3` and `S_4/V_4=S_3` from candidate
THM-2598 point to the same group in (15), but not to a canonical action on
this support.  A common quotient group is not a transfer of selectors.

## 4. Exact tournament meaning

Use the repository tile order

```text
(6,1),(5,1),(4,1),(3,1),(6,2),
(5,2),(4,2),(6,3),(5,3),(6,4).                            (16)
```

The support (5) has bit mask

```text
873;                         tile complement 150.          (17)
```

Both masks are grid-symmetric and their fixed-base-path tournaments are
isomorphic.  Exhausting all `2^10` tilings and all `6!` relabelings gives
the smaller endpoints of the two `n=6` blue self-loop lines as

```text
19, 150.                                                   (18)
```

So `{150,873}` is exactly one of the two first blue self-loops anticipated
in THM-2467.  For the tournament at mask `873`,

```text
score vector=(2,3,2,3,2,3),
Hamiltonian paths=45,
Aut(T)=< (1 3 5)(2 4 6) > isomorphic to C3.                (19)
```

The five-cycle `(2 3 4 6 5)` relabels it to the complement-mask
tournament.  Crucially, neither support rotation `rho_2` nor `rho_3` from
(10) is a tournament automorphism.  The modular action belongs to the
**bicycle support**, not to the directed tournament.  This is the exact
sidecar that a tournament-only statement would lose.

## 5. Graceful boundary

The cycle support is not graceful.  Suppose a graceful labeling existed.
Its six edge differences would be `1,...,6`, with odd total

```text
1+2+3+4+5+6=21.                                          (20)
```

Modulo two, however, the sum of absolute edge differences equals

```text
sum_({u,v} in E) (f(u)+f(v))
 =sum_v deg(v) f(v)=0,                                   (21)
```

because every degree of `C_6` is two.  This contradicts (20).  The
companion also exhausts all `7P6=5,040` injective assignments.

Deleting any one cycle edge gives `P_6`.  Along its path order, the labels

```text
0,5,1,4,2,3                                               (22)
```

have differences `5,4,3,2,1`, so the resulting tree is graceful.  This is
a sharp scope boundary: the modular partial cube is nongraceful, while its
spanning path satisfies the graceful-tree prediction.  Nothing here
advances the graceful tree conjecture beyond that example.

## 6. What does not transfer to LRC

THM-2587 has six live wall pieces indexed by

```text
theta in C2  times  {low,middle,high},                     (23)
```

with each `theta` rail carrying exact masses `48,154,48`.  Its three-state
factor is an ordered interval/reflection wall, not a cyclic `C_3` action.
Therefore the shared cardinality six does not identify the LRC selector
wall with the bicycle `C_6` or with the `S_3` matching completion.

The precise information ledger is

```text
tile cut/cycle vector
  -> C6 support and its C2 x C3 rotations
  -> two selected perfect matchings
  -> [add the frozen matching and four new matching fibres]
  -> K3,3 / S3 parity obstruction.                         (24)
```

The first two arrows forget tournament orientation and free modular words.
The completion arrow changes the matching family; it is not a harmless
relabelling.  No LRC owner, clock, phase, or terminal current occurs in
this chain.

## 7. Exact reproduction

```bash
python 04-computation/six_vertex_bicycle_modular_abelianization_thm2597.py
python -O 04-computation/six_vertex_bicycle_modular_abelianization_thm2597.py
```

Both modes reproduce the stored transcript byte-for-byte after LF
normalization.  The companion independently checks the cut/cycle census,
star words, all partial-cube distances and `Theta` actions, perfect
matchings, support rotations, all `1,024` tile masks with exact `6!`
canonicalization of every grid-symmetric blue-loop candidate, (19), the
explicit complement isomorphism, and the graceful hostile.
Every truth-bearing check remains active under `python -O`.
