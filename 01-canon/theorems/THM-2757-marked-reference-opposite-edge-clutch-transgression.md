---
id: THM-2757
title: "Marked-reference opposite-edge clutch transgression"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  One marked vertex
  and three arms canonically
  orient the three opposite-edge pairs of K4 into a 2-by-3 star/opposite
  table.  The opposite-edge difference is a scaled Hadamard transform with
  D D^T=4I and D^T D=4I-J.  Every column total is automatically constant, so
  a formal two-row gain (alpha,beta) has charged C3 output exactly through
  (alpha-beta) times the marked-reference standard component.  Swapping an
  unequal reference into one equal arm yields both charged characters, but
  a transposition and V4 double transposition remain value-indistinguishable
  on the equal corolla.  This is a finite representation/design theorem; it
  constructs no physical LRC four-point carrier or wing-diagonal operator.
source: root/marked-reference-opposite-edge-transgression-2026-07-28
audit: >
  root-zero-clutch-audit/2026-07-28 (independent marked table, Hadamard,
  general transgression, F13 charge/census, A4/S4 hostile, representation
  scope, and normal/-O/hash replay: ACCEPT after field/source repairs)
depends_on: []
related:
  - THM-2721-semantic-inner-triangle-equal-following-amplitude-and-current-reanchoring-no-go
  - THM-2750-arm-blind-clutch-no-go-and-minimal-marked-leakage
  - THM-2751-root-zero-clutch-mayer-vietoris-wing-shear
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-2756-opposite-edge-projectors-parity-cancellation-and-integral-clutch
script: 04-computation/lrc14_marked_reference_opposite_edge_clutch_thm2757.py
output: 05-knowledge/results/lrc14_marked_reference_opposite_edge_clutch_thm2757.out
script_sha256: 0a86047b240e66d643a5a99693c800fabcce9331944c961dcbf2b97473d7a666
output_sha256: b0c6c207762dde15a5e0d9f5ab13035b014dfd4f3e58c871c01ab663b9d4d8e0
hash_basis: LF-normalized bytes
---

# THM-2757 -- marked-reference opposite-edge clutch transgression

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2750 isolates the smallest positive way to charge an equal three-arm
corolla: add one fixed reference and move it into an arm.  THM-2753 identifies
the three perfect matchings of four points as the canonical ternary resolvent,
while THM-2756 splits the six edges into matching and standard channels.  The
present candidate combines those observations at the marked-coordinate level.

The mark does something precise: it removes the sign gauge in each opposite
edge pair.  The six edges become a literal `2 x 3` table, and its column sums
are automatically equal.  Thus the conditional invariant-total hypothesis
which appeared in the now-refuted first THM-2751 candidate has an exact finite
model.  MISTAKE-313 is load-bearing, however: the `(12,2)` specialization below
is only a formal algebraic control, not a physical root-zero wing operator.

## 1. A marked cone gives a canonical two-by-three table

Let `k` be a commutative ring, let

```text
X={0,1,2,3},                         r=(1 2 3),          (1)
```

and mark vertex `0`.  For a vertex field `v=(v0,v1,v2,v3)`, order the three
perfect matchings by the arm incident to the mark:

```text
m1=01|23,                 m2=02|13,                 m3=03|12. (2)
```

Define the star and opposite rows

```text
M=(v0+v1, v0+v2, v0+v3),
L=(v2+v3, v1+v3, v1+v2).                              (3)
```

The arm cycle permutes the three columns and preserves the two row names.
Every column is a decomposition of the same four-vertex total:

```text
M_i+L_i=sigma(v):=v0+v1+v2+v3.                       (4)
```

So a marked four-point carrier supplies the constant-column-total condition
without an additional hypothesis.

Put `D=M-L`.  In the ordered vertex basis,

```text
D = [1  1 -1 -1]
    [1 -1  1 -1]
    [1 -1 -1  1].                                     (5)
```

Direct integral multiplication gives

```text
D D^T=4 I_3,                       D^T D=4 I_4-J_4.    (6)
```

Consequently over a field of characteristic different from two, `D` induces
a scaled isomorphism from the four-point standard quotient
`k^4/k(1,1,1,1)` to the three marked opposite-edge differences.  Formula
`(6)` is the marked-coordinate form of THM-2756's negative edge block.  More
precisely, the vertex-sum edge tables `(3)` form only the four-dimensional
`1 direct-sum [31]` submodule of the full six-edge module: their constant
matching totals kill the `[22]` summand.  This theorem does not identify its
marked chart with all of THM-2756's
`1 direct-sum [22] direct-sum [31]` carrier.

## 2. Exact two-row clutch transgression

For two formal row gains `alpha,beta`, fold the two rows after applying them:

```text
T_i=alpha M_i+beta L_i.                                (7)
```

The division-free identity is

```text
2T_i=(alpha+beta)sigma(v)+(alpha-beta)D_i.             (8)
```

If two **and three** are invertible and `Q_3` denotes the mean-zero projector
away from the equal-column line, then

```text
Q_3 T=((alpha-beta)/2) Q_3 D
     =(alpha-beta) Q_3 M.                              (9)
```

In particular, if `alpha-beta` is a unit,

```text
Q_3T=0  iff  v1=v2=v3.                                (10)
```

Indeed `D_i-D_j=2(v_i-v_j)`.  Thus the clutch does not create charge from a
reference merely because it is present; it detects exactly the failure of the
three arm values to remain equal.  The mark supplies a coordinate for the
standard component, not the component itself.

The corresponding `2 x 3` mixed rectangle is

```text
Delta_ij=M_i+L_j-M_j-L_i=D_i-D_j,                     (11)
```

and `(8)` gives

```text
T_i-T_j=((alpha-beta)/2) Delta_ij.                     (12)
```

This is the exact ANOVA/transgression dictionary: the external charged mode,
the marked opposite-edge standard component, and the nonzero mixed rectangle
are the same obstruction after multiplication by units.

## 3. The marked swap and the formal F13 control

Now work over `F_13`, choose the primitive cube root `omega=3`, and retain the
formal control

```text
(alpha,beta)=(12,2).                                   (13)
```

Start with a fixed reference plus equal arms,

```text
v=(x,A,A,A).                                           (14)
```

Its table and every output `(7)` are column-flat.  Apply the positive
permutation `s=(0 1)` while keeping vertex `0` as the marked coordinate.  The
new value field is

```text
s v=(A,x,A,A),                   delta=x-A,             (15)
```

and exact substitution gives

```text
D=delta(1,-1,-1),
Q_3T=delta(11,1,1).                                    (16)
```

With normalized three-point Fourier transform, both charged coefficients are

```text
That(1)=That(2)=12 delta.                              (17)
```

The rectangles and target differences are

```text
(Delta_12,Delta_13,Delta_23)=delta(2,2,0),
T_1-T_2=T_1-T_3=5 Delta_12.                            (18)
```

Hence `delta!=0` fires both nontrivial arm characters.  Exhaustion of all
`13^4=28561` vertex fields gives exactly `13^2=169` flat tables and `28392`
charged tables, in agreement with `(10)`.

The word **formal** in `(13)` is essential.  MISTAKE-313 retracts the physical
root-zero interpretation of the old `(12,2)` common/wing diagonal.  The
correct frozen-clock carrier has a coefficient-null left wing and nonzero
right wing, so no linear row-diagonal operator maps one to the other.  This
section proves what a lawful unequal two-row gain would do if a future
physical carrier supplied it; it does not assert that THM-2749 or repaired
THM-2751 does.

## 4. Sharp A4/S4 value collision

The marked table retains more than an unlabeled edge cycle type, but an equal
corolla still has a sharp blind spot.  Let

```text
s=(0 1),                         d=(0 1)(2 3).          (19)
```

On `(14)`,

```text
s v=d v=(A,x,A,A),                                      (20)
```

so their complete star/opposite value tables, rectangles, and charged outputs
are identical.  Yet their actions on the three matchings are

```text
s: (m1,m2,m3) -> (m1,m3,m2),
d: (m1,m2,m3) -> (m1,m2,m3).                           (21)
```

Thus the value transgression detects that the marked reference moved, but not
whether the mover was a transposition generating the `S4` branch or a `V4`
double transposition generating the `A4` branch.  The matching-label action or
a mixed group word remains necessary, exactly as in THM-2753.  This is the
minimal hostile to upgrading `(16)` into an `A4/S4` discriminator.

## 5. Holotopy and physical boundary

The construction may be viewed as a marked cone on the three arms.  Star
edges and opposite triangle edges form the two rows; `D` is their signed
opposite-edge boundary, and `(9)` is its transgression into the charged arm
sector.  This is an exact finite linear model, not a claim that the LRC
current/following nerve contains that cone.

The closest physical candidates each miss one datum:

| source | retained | missing |
|---|---|---|
| THM-2721 | positive equal three-arm corolla and four changed-source controls | no four-point `S4` action, edge-incidence table, or target reanchoring |
| THM-2749 | marked root-zero common clutch and primitive `C13` target characters | its target colours are not external `C3` arms or four vertices |
| MISTAKE-313 frozen-clock `FINITE-EXACT` repair | a genuine one-sided target cofiber in one frozen clock chart | no linear left-to-right wing gain or marked four-point carrier; THM-2751 remains unavailable until separately re-promoted |
| THM-2756 | canonical rational/integral `2/3` edge splitting | no physical LRC six-edge packet or endpoint current |

Therefore the theorem solves the **finite representation design** but not its
realization.  A physical application needs one literal four-state carrier on
common support, a lawful marked-reference move, the edge-sum/incidence map,
and endpoint/current typing preserved through `(7)`.  Merely observing four
addresses or three target characters is insufficient.

## 6. Exact companion and ledger

Run

```bash
python 04-computation/lrc14_marked_reference_opposite_edge_clutch_thm2757.py
python -O 04-computation/lrc14_marked_reference_opposite_edge_clutch_thm2757.py
```

Both modes byte-match the stored transcript.  The dependency-free companion
uses integral arithmetic and exhaustive `F_13` enumeration, with explicit
exceptions and no truth-bearing assertions.  It checks `(4)--(12)`, all
`28561` vertex fields, `(15)--(18)` for all `169` pairs `(x,A)`, and the
distinct matching actions in `(21)`.

```text
PROVED HERE:             marked-reference 2 x 3 edge table;
                         scaled Hadamard/opposite-edge identities;
                         exact general clutch transgression and rectangles;
                         complete F13 charge criterion and census;
                         minimal marked-swap charged profile;
                         sharp equal-corolla A4/S4 value collision.

NOT CONSTRUCTED:         a physical unequal root-zero wing operator;
                         an LRC four-point or six-edge carrier;
                         semantic endpoint/current or owner transport;
                         row exclusion or LRC(14).                       (22)
```

An independent hostile audit rederived the integral Hadamard identities,
general clutch and rectangle laws, exact field hypotheses, complete `F_13`
census, marked-swap characters, and the equal-corolla `A4/S4` collision.  It
also checked that the vertex-sum carrier omits THM-2756's `[22]` block and
that MISTAKE-313 supplies only finite-exact frozen-clock context.  Normal,
optimized and stored outputs byte-match after LF normalization, and both
declared hashes are exact.  No physical LRC realization was inferred.

QED.
