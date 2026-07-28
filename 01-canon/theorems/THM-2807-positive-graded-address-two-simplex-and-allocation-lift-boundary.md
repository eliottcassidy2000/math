---
id: THM-2807
title: "Positive graded-address two-simplex and allocation-lift boundary"
status: >
  RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
  INDEPENDENT AUDIT.  Three literal positive semantic cylinders form a
  commuting address-translation two-simplex with edges Z1, Z2^4079, and
  Z1 Z2^4079.  The THM-2791 tau-three unit is exactly its diagonal after
  deleting the pure-Z1 vertex from the positive tau-twelve filler.  Modulo
  169 the vertical edge collapses but retains transgression 10.  A unique
  affine odometer element fixes the initial address and sends the pure-Z1
  edge to the graded diagonal.  This is collapsed weighted-carrier/address
  holotopy, not factor covariance, endpoint allocation, a root/Cech map,
  row exclusion, or LRC(14).  Until independent promotion, no proved result
  may depend on this candidate.
source: root/positive-graded-address-simplex-2026-07-28
depends_on:
  - THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2791-full-arm-orbit-transfer-and-lower-central-chord
  - THM-2792-cyclic-unit-intertwiner-and-positive-naturality-boundary
related:
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2763-carrier-equivariant-endpoint-address-extension-and-gauge-obstruction
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
script: 04-computation/lrc14_positive_graded_address_two_simplex_thm2807.py
output: 05-knowledge/results/lrc14_positive_graded_address_two_simplex_thm2807.out
script_sha256: 15e32522b5ee913cb4223c037836d27244ee1743f1e14dbc628d4e6f659708cb
output_sha256: c80c1243fd1570c37f6ba444dc91bce8e0ce4bee7ebec284b206cbcd84058894
hash_basis: LF-normalized bytes
---

# THM-2807 -- the graded semantic chord bounds a positive address triangle

**RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
INDEPENDENT AUDIT.**

THM-2782 supplies a positive physical central arm and THM-2791 extends it to
a two-point full-depth chord whose first nontrivial carry digit is `10`.
That chord could still have been an isolated positive edge.  It is not.

The exceptional target column `tau=12` fills a third positive cylinder
between the two THM-2791 endpoints.  All three arrows are exact translations
of whole weighted cylinders.  In the nerve of the physical address action
they form a commuting two-simplex:

```text
                         Z2^4079
             n_+ ----------------------> n_a
              ^                            ^
              |                            |
           Z1 |                            | Z1 Z2^4079
              |                            |
              +------------ n_0 -----------+
```

Here the slanted arrow is read from `n_0` to `n_a`.  The `tau=3` packet is
exactly the diagonal face: it agrees with `tau=12` at `n_0,n_a` and is empty
at `n_+`.

This resolves the bare address-level homotopy behind the graded chord.  It
does not lift the triangle to THM-2772's endpoint-origin/allocation pullback.
That distinction is load-bearing: equality of weighted cylinders and an
address-group conjugacy do not give factorwise covariance or a common
allocated endpoint atom.

## 1. The three exact physical addresses

Use

```text
p=13,                         M=p^6=4826809,
n_0=3454614,
n_+=n_0+13=3454627,
n_a=n_0+689364=4143978.                              (1)
```

The addresses lie in THM-2782's clock-one target coordinate, with source
rail eight, root-one target chart, predecessor carry six, and semantic word

```text
E3 -> D^6 -> Q_(3,{1,2}).                              (2)
```

Let `C_tau(n)` be the exact narrow target cylinder at address `n`, restricted
inside THM-2782's full lawful `(sigma,tau)=(0,tau)` present packet, both
relative-present safeties, rail-eight weight, root-one half, and delayed
carry-six coefficient.

Put

```text
mu=60781651775958960/371293,
c =790161473087466480,
w =27581135604.                                       (3)
```

Direct exact restriction gives the complete selected table

```text
                n_0             n_+             n_a
tau=3          (mu,c)           (0,0)           (mu,c)
tau=12         (mu,c)           (mu,c)          (mu,c). (4)
```

Every nonzero cell in `(4)` is one whole cylinder with constant raw rail
weight `w`.  More strongly, as literal weighted intervals,

```text
C_3(n_0)=C_12(n_0),             C_3(n_a)=C_12(n_a).     (5)
```

Thus the THM-2791 two-point packet is obtained from this three-vertex
`tau=12` packet by deleting precisely `n_+`.  This is not an equality of
separately integrated marginals.

All three centres have semantic record `(3,{1,2})`, predecessor carry six,
source target-label `sigma=0`, and target label `tau=12`.  The two diagonal
centres additionally retain target label `tau=3`.  The exact positive
controls therefore live on the same typed semantic sheet, not merely in the
same numerical address orbit.

## 2. The lower-central triangle

The three address gaps are

```text
n_+-n_0 =13,
n_a-n_+=689351=169*4079,
n_a-n_0 =689364=13*53028,
53028=1+13*4079.                                       (6)
```

In THM-2788's depth-six odometer notation

```text
Z1=O^13,                         Z2=O^169,              (7)
```

where `O:n |-> n+1`.  Hence the three edges in `(6)` are exactly

```text
Z1,                    Z2^4079,                    Z1 Z2^4079. (8)
```

Their first vertical residue is

```text
4079=10 mod13,                                          (9)
```

the first transgression computed independently for the diagonal in
THM-2791.

The endpoint-origin quotient sees only addresses modulo `169`.  There

```text
(n_0,n_+,n_a)=(85,98,98) mod169.                        (10)
```

Thus the vertical edge collapses, and the pure edge and graded diagonal are
two physical lifts of the same quotient edge `85 -> 98`.  Their difference
is not lost at full depth: it is exactly `Z2^4079`.

This is the precise holotopy statement.  The quotient triangle is degenerate,
while its full-address lift has a nontrivial vertical kernel edge.  It is not
a new cohomology class: the three arrows already form a boundary in the
translation-action nerve.  What survives quotienting is the need to choose
between two lifts, with transition `(9)`.

## 3. Whole-cylinder composition

The three physical circle shifts are

```text
epsilon_1 =7(n_+-n_0)/13^6 =7/371293,
epsilon_2 =7(n_a-n_+)/13^6 =28553/28561,
epsilon_d =7(n_a-n_0)/13^6 =371196/371293.              (11)
```

They satisfy

```text
epsilon_1+epsilon_2=epsilon_d mod1.                     (12)
```

Exact endpoint arithmetic, including wraparound, gives

```text
T_(epsilon_1) C_12(n_0)=C_12(n_+),
T_(epsilon_2) C_12(n_+)=C_12(n_a),
T_(epsilon_d) C_12(n_0)=C_12(n_a).                     (13)
```

Each equality in `(13)` preserves the entire interval and its constant
weight `w`.  Thus the triangle commutes before integration.

The phrase **weighted-carrier triangle** is intentional.  The companion
reconstructs the Boolean product and proves that its restriction is exactly
one interval of weight `w`; it does not exhibit an identity on every hidden
factor label along all three arrows.  THM-2791 separately proves literal
rail-sheet ancestry on the diagonal.  No corresponding factorwise statement
for both `tau=12` edges is asserted here.

## 4. A unique bare affine conjugacy

THM-2788 also has the address dilation

```text
X:n |-> 14n.                                             (14)
```

The element `14=1+13` has order `13^4` modulo `13^5`.  Exhaustive exact
discrete logarithm, equivalently repeated use of the odd-prime lifting law,
gives the unique exponent

```text
k=23098 mod13^4,
14^k=53028 mod13^5.                                    (15)
```

At depth six,

```text
m=14^k mod13^6=2652079.                                (16)
```

Set

```text
v=(1-m)n_0 mod13^6=352469,
g=O^v X^k.                                              (17)
```

Then

```text
g(n_0)=n_0,                    g(n_+)=n_a,              (18)
```

because

```text
m*13=689364 mod13^6.                                   (19)
```

So the pure `Z1` edge and the graded diagonal are conjugate in the bare
affine address group while their initial vertex is fixed.  The base-`13`
digits of `k` are

```text
k=(10,8,6,10)_13.                                      (20)
```

Equations `(15)--(20)` are an address statement only.  They do not say that
the physical factor packet is covariant under `X^k`, and they do not upgrade
`g` to a map of endpoint currents.  Indeed the existence of the translation
triangle `(13)` makes the missing datum sharper: the issue is naturality of
the retained factors and endpoint allocation, not reachability of the
addresses.

## 5. What the simplex repairs, and what it does not

The source-to-target contract is now:

| item | exact content |
|---|---|
| source | the three `tau=12` whole cylinders in `(4)` |
| target face | THM-2791's `tau=3` diagonal |
| map | delete `n_+`; identities `(5)` on the retained vertices |
| address homotopy | the commuting translation triangle `(13)` |
| quotient loss | `n_+=n_a mod169` |
| retained transition | `Z2^4079`, first residue `10` |
| affine comparison | unique bare `O^vX^k` in `(17)` |
| first unconstructed lift | one common endpoint origin and Boolean allocation atom |

Consequently the next exact object is not another cyclic Fourier invariant.
THM-2792 already supplies the unique abstract cyclic-module isomorphism and
THM-2802 reconstructs its projective coefficient orbit.  This theorem now
supplies the missing address homotopy as well.

The cheapest physical test is instead to pull the explicit three-cylinder
simplex through THM-2772's

```text
A_tilde x_(G_full) P                                   (21)
```

while retaining one endpoint origin, `(L,R,q,Delta)`, all four allocation
states, semantic word/owner, triangle, and common Abel boundary.  One should
test, in order:

1. a nonempty common allocation atom over `n_0`;
2. lawful transport along the pure and graded edges;
3. nonzero mixed face `(P_0-P_1)(Q_0-Q_1)` before Fourier;
4. equality of the two endpoint allocations after `(10)`, with transition
   `Z2^4079` retained;
5. the determinant/root-Cech invoice; and only then
6. the thirteen Fourier values of the resulting same-ancestry cycle.

A failure at the first step is an owner/allocation obstruction.  A failure
at steps two or four is genuine carry holonomy.  Success would construct a
new physical endpoint cycle; it would not identify that cycle with the
canonical THM-2625/2790 bank without a second common-current cospan.

## 6. Sharp boundaries

1. **No factor covariance.**  Whole weighted cylinders translate as in
   `(13)`; hidden factor labels need not.
2. **No endpoint origin.**  The mod-`169` coincidence `(10)` is not a choice
   of the missing `13`-element endpoint-origin fibre.
3. **No allocation.**  The triangle contains no bare/source/target/both
   `K4` and no common mixed face.
4. **No automatic physical intertwiner.**  The affine conjugacy `(17)` lives
   in the address group.  It is not THM-2792's pointwise torsor map, which
   that theorem explicitly withholds.
5. **No nontrivial quotient class.**  The full triangle is already a
   translation-nerve boundary; only the vertical lift is forgotten modulo
   `169`.
6. **No ledger consequence.**  No endpoint current, root/Cech correction,
   row exclusion, or LRC(14) conclusion follows.

## 7. Exact companion

Run

```bash
python 04-computation/lrc14_positive_graded_address_two_simplex_thm2807.py
python -O 04-computation/lrc14_positive_graded_address_two_simplex_thm2807.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_positive_graded_address_two_simplex_thm2807.out.
```

The companion rebuilds the inherited THM-2782 interval carrier once and uses
exact rational arithmetic throughout.  It verifies all six selected cells
in `(4)`, the literal interval equalities `(5)`, semantic/carry labels, all
three translations, lower-central and quotient arithmetic, the complete
`13^4` discrete-log search in `(15)`, and the affine equations `(16)--(19)`.
It contains explicit exception gates, no truth-bearing Python assertions,
no floating point, and no scratch dependency.

**Awaiting independent audit; not QED.**
