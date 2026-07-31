---
id: THM-2968
title: "Quartic edge and oriented-cycle S4 complements"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over a
  characteristic-zero field, the six-edge and oriented-four-cycle lifts of a
  full-S4 quartic, on their common separability open and in the exact common
  matching labelling frozen below, are the two signed-pair complements H0 and
  H1.  Pointwise rho_or(g)=z^sgn(g)rho_edge(g) for all 24 sheet permutations.
  Their permutation/fifth-exterior characters are respectively
  (1+2+3,1+2+3) and (1+2+3sgn,sgn+2+3).  THM-2769 supplies a generic-S4
  physical witness, but its V4 normalization ramifies: it is an affine
  boundary hostile, not a quasi-etale Keller layer or a Keller map.
source: codex-thm2968-quartic-edge-orientation-s4-complements-2026-07-30
audit: >
  Independent hostile canonicalization fixed the characteristic-zero scope,
  the exact common orientation-pair gauge, the 24 pointwise central-twist
  tests, the explicit transposition control, the ramified/non-Keller boundary,
  and the requirement that the diagnostic loop already be known to have
  transposition monodromy.  Normal and optimized replays byte-match the
  LF-stored transcript.
depends_on:
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-2862-modular-level-three-four-congruence-ladder-and-inequivalent-six-lifts
  - THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product
  - THM-2965-modular-signed-pair-complements-and-v4-compound-shadow
  - THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile
related:
  - THM-2887-quaternionic-arf-lift-of-the-semantic-v4-and-global-carry-no-go
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary
script: 04-computation/quartic_edge_orientation_s4_complements_thm2968.py
output: 05-knowledge/results/quartic_edge_orientation_s4_complements_thm2968.out
script_sha256: 6b3bb5dbe92824c99e267fd7c7da3d3582701b742345ac6160e16f93ffb79d66
output_sha256: 0d11a1bd5d7709d477a0d87a1112d0dd728ada2b3e5a3cf4771faa5851315692
hash_basis: LF-normalized bytes
---

# THM-2968 -- quartic edge and oriented-cycle `S4` complements

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and statement

The nearest proved mechanism is
[THM-2753, six-edge parity erasure and matching restoration](THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration.md):
the faithful action of `S4` on the six edges of `K4` has identically even
ambient six-point sign, while the action on the three perfect matchings is the
quotient `S4/V4=S3` and recovers the quartic sign.  [THM-2862, the two
six-point lifts](THM-2862-modular-level-three-four-congruence-ladder-and-inequivalent-six-lifts.md)
constructs the oriented-four-cycle action, whose ambient sign is the quartic
sign.  [THM-2864, the polynomial sextic
resolvents](THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product.md)
makes both actions polynomial and computes their discriminants.

The finite target is the signed-pair group from
[THM-2965, signed-pair complements and compound
shadow](THM-2965-modular-signed-pair-complements-and-v4-compound-shadow.md):

```text
W=F2^3 semidirect S3,
(delta,sigma): e_(i,b) -> e_(sigma(i),b+delta_i).
```

Put

```text
p(delta)=delta_1+delta_2+delta_3 mod 2,
s(sigma)=parity(sigma),
z=((1,1,1),(0,1,2)),
H0=ker p,                    H1=ker(p+s).                (1)
```

Work over a field of characteristic zero.  Let a separable quartic have full
`S4` Galois/monodromy action on its four labelled roots, and restrict to the
open on which the edge and oriented-cycle sextic generators are both
separable.  Label the three opposite-edge pairs by

```text
(((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))             (2)
```

and, in exactly the same matching order and gauge, label the three inverse
orientation pairs by

```text
(((0,2,1,3),(0,3,1,2)),
 ((0,3,2,1),(0,1,2,3)),
 ((0,1,3,2),(0,2,3,1))).                                (3)
```

Here a four-cycle is written modulo cyclic rotation, normalized to begin at
`0`.  If `epsilon(g)` is the parity bit of `g`, then the induced signed-pair
embeddings satisfy the pointwise identity

```text
rho_or(g)=z^epsilon(g) rho_edge(g)       for every g in S4, (4)

rho_edge(S4)=H0,                         rho_or(S4)=H1.    (5)
```

Thus every edge move makes an even number of within-pair flips.  An
oriented-cycle move differs in all three within-pair gauges precisely for an
odd sheet permutation; equivalently its within-pair flip parity is the
matching-permutation parity.

## 2. Proof of the two complement identities

For any permutation of three two-element blocks, exchanging two whole blocks
is a product of two point transpositions.  Hence only within-block flips
contribute to its sign on six points:

```text
sign_6(delta,sigma)=(-1)^p(delta).                        (6)
```

Both physical actions have the same three-block permutation `sigma=mu(g)`.
This common quotient is the quartic normal-`V4` frame

```text
S4 -> S4/V4=S3,                                           (7)
```

equivalently the three roots of the resolvent cubic.  THM-2864 proves the
coefficient identity

```text
Disc(resolvent cubic)=Disc(quartic)=Delta,                (8)
```

while THM-2753 gives `s(mu(g))=epsilon(g)`.

For edges, THM-2753 says the six-point sign is always `+1`; `(6)` gives
`p=0`, so `rho_edge(S4)` lies in `H0`.  The action is faithful, and both
groups have order `24`, hence its image is all of `H0`.

For oriented four-cycles, THM-2862 says the six-point sign is the quartic
sign.  Equations `(6)--(8)` give `p=s`, so the image lies in `H1`.  Its point
stabilizer is `C4`, with trivial core in `S4`; the action is faithful, and the
common order `24` gives equality with `H1`.

The exact gauge `(3)` strengthens these image equalities to `(4)`.  The
sheet transposition and ternary generator

```text
a=(0,1,3,2),                         b=(2,0,1,3)          (9)
```

generate `S4`, and direct relabelling gives

```text
rho_edge(a)=(000,(0,2,1)),       rho_or(a)=(111,(0,2,1)),
rho_edge(b)=rho_or(b)=(011,(1,2,0)).                     (10)
```

Since `z` is central, `(10)` proves `(4)` on every word in `a,b`.  The exact
companion additionally checks `(4)` separately for all `24` permutations,
so no generator-convention or tuple-order assumption is hidden.

This identifies the common three-block quotient more sharply: the edge and
orientation sextics are two lifts of one quartic-`V4`/resolvent-cubic `S3`
frame, and the lifts differ by the central `C2` sign graph.  This is also the
physical interpretation of THM-2965's `C2*C3` calculation.  The binary
and ternary signed-pair lifts are
`s0=rho_edge(a)` and `t=rho_edge(b)`.  The first can be replaced by `z s0`
without changing order because `(z s0)^2=1`, whereas twisting the second gives
`(z t)^3=z`, not `1`.  This is an interpretation of the already-proved
signed-pair statement, not a new Keller conclusion.

## 3. The exact permutation and fifth-exterior characters

Let `U_edge` and `U_or` be the two six-dimensional permutation modules over
characteristic zero.  Their cycle tables are

| sheet type | edge type | orientation type |
|---|---|---|
| `1^4` | `1^6` | `1^6` |
| `2 1^2` | `2^2 1^2` | `2^3` |
| `2^2` | `2^2 1^2` | `2^2 1^2` |
| `3 1` | `3^2` | `3^2` |
| `4` | `4 2` | `4 1^2` |

Writing the characteristic-zero irreducibles of `S4` as
`1,sgn,2,3,3sgn`, exact character inner products give

```text
U_edge =1+2+3,                       Lambda^5 U_edge=1+2+3;
U_or   =1+2+3sgn,                    Lambda^5 U_or=sgn+2+3. (11)
```

Indeed `Lambda^5 U=det(U) tensor U*`.  The edge determinant is trivial and
the orientation determinant is `sgn`; permutation modules are self-dual.
These are exactly the `H0` and `H1` source decompositions of THM-2965.  The
determinant twist in the second fifth exterior power is the retained quartic
sign.

## 4. Polynomial discriminants and one physical full-`S4` witness

For a depressed quartic

```text
f(X)=X^4+pX^2+qX+r,
```

THM-2864 writes the edge sextic as `E(Y)=S(Y^2)`, where `S` is the matching
cubic, and constructs the orientation sextic `O`.  On the common separability
open their discriminants are

```text
Disc(E)=64 q^2 Delta^2=(8qDelta)^2,                       (12)
Disc(O)=2^66 q^12 J_or^4 Delta^3
       =Delta (2^33 q^6 J_or^2 Delta)^2.                 (13)
```

Thus the coefficient square classes are respectively `1` and `[Delta]`, the
polynomial shadows of `H0` and `H1`.

[THM-2769, the full-`S4` affine-boundary
hostile](THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile.md) proves
that over `C(t)`

```text
F_t(Y)=Y^4-2Y^2-8tY+1-4t                               (14)
```

has generic Galois group `S4`.  Here `p=-2`, `q=-8t`, `r=1-4t`, and exact
substitution into THM-2864 gives

```text
Delta=-4096 t^2(27t^2-14t+3),                            (15)
J_or=4096(t-1)(27t^3-25t^2+8t-1).                       (16)
```

Since `q Delta J_or` is not the zero polynomial, both sextic covers are
separable on a common dense open of this actual generic-`S4` family.  Hence
`(14)` is a physical witness for both complement images in `(5)`.

Its affine meaning is strictly negative.  THM-2769 proves that at `t=0` the
three Kummer root valuations have parity word `110`, so the normalization in
the `V4` layer ramifies.  It therefore does **not** supply the quasi-etale
`V4` normalization required of a Keller layer.  The family is an affine
boundary hostile, not a Keller map; no Jelonek divisor or Keller source is
constructed here.

## 5. Cheapest diagnostic and exact boundary

Once full `S4` quotient monodromy and faithfulness are known, a loop already
known independently to have **transposition monodromy** decides the complement
bit.  For the transposition `a` in `(9)`, `(10)` gives the exact controls

```text
edge:          (000,(0,2,1)),          cycle type 2^2 1^2;
orientation:   (111,(0,2,1)),          cycle type 2^3.    (17)
```

The provenance of the transposition is load-bearing: this theorem does not
infer a local branching model or assert that an arbitrary odd loop is a
transposition.  In characteristic zero the same complement bit can instead
be read from `(12)--(13)`, provided the full-`S4`, faithfulness, and
separability hypotheses have already been established.  A square sextic
discriminant without those hypotheses need not define an `S4` complement,
and on an `A4` lane the two complements restrict identically.

The action-level bridge is also not a multiplication bridge.  THM-2965's
incidence map retains only constant and first-Walsh layers, while
[THM-2951, fifth-compound reconstruction and its
boundary](THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary.md)
and THM-2965 give the diagonal hostile

```text
fifth-compound weight 450 != balanced-third weight 30.   (18)
```

Thus the physical `H1` orientation section does not turn that incidence map
into a multiplication-compound intertwiner.  It also does not identify the
orientation sextic with the reduced real six-point scheme of
[THM-2950, the quartic resolvent
frame](THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame.md).
No SFC(4), degree-four Keller, `JC(2)`, `DC(2)`, GMC, or LRC conclusion
follows.

## 6. Exact evidence

Run

```text
python 04-computation/quartic_edge_orientation_s4_complements_thm2968.py
python -O 04-computation/quartic_edge_orientation_s4_complements_thm2968.py
```

Both modes byte-match the LF-stored transcript

```text
05-knowledge/results/quartic_edge_orientation_s4_complements_thm2968.out.
```

The companion enumerates all `24` sheet permutations, all `48` signed-pair
elements, both `24`-element complements, both six-point actions, the exact
common labelling `(2)--(3)`, the pointwise identity `(4)`, both homomorphism
laws, the binary/ternary controls, both conjugacy-class tables, and all
character inner products in `(11)`.  It then substitutes `(14)` into the
THM-2864 discriminant formulas and checks `(15)--(16)` and both square-class
identities.  Every truth-bearing gate is an explicit `require`/exception;
there is no Python `assert` or floating-point decision.

LF-normalized SHA256:

```text
script  6b3bb5dbe92824c99e267fd7c7da3d3582701b742345ac6160e16f93ffb79d66
output  0d11a1bd5d7709d477a0d87a1112d0dd728ada2b3e5a3cf4771faa5851315692
```

**QED.**
