---
id: THM-2807
title: "Positive graded-address two-simplex and allocation-lift boundary"
status: >
  RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
  INDEPENDENT AUDIT.  Three literal positive semantic cylinders form a
  commuting address-translation two-simplex with edges Z1, Z2^4079, and
  Z1 Z2^4079.  The THM-2791 tau-three unit is exactly its diagonal after
  deleting the pure-Z1 vertex from the selected tau-twelve three-vertex
  restriction.  The quotient modulo 169 collapses the vertical edge, while
  its full-depth lift retains transgression 10.  One unique depth-five
  exponent class has thirteen full-depth affine lifts; all agree on the
  selected low-digit-seven sheet.  This is collapsed weighted-carrier/address
  holotopy, not factor covariance, endpoint allocation, a root/Cech map, row
  exclusion, or LRC(14).  Until independent promotion, no proved result may
  depend on this candidate.
source: root/positive-graded-address-simplex-2026-07-28
depends_on:
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2791-full-arm-orbit-transfer-and-lower-central-chord
related:
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2763-carrier-equivariant-endpoint-address-extension-and-gauge-obstruction
  - THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2792-cyclic-unit-intertwiner-and-positive-naturality-boundary
script: 04-computation/lrc14_positive_graded_address_two_simplex_thm2807.py
output: 05-knowledge/results/lrc14_positive_graded_address_two_simplex_thm2807.out
script_sha256: 11cdbe3c6cc7f9d5b6b24863ced71eb91cc84adc67fe38a3f8a3e637362453fb
output_sha256: a6fbc42c5a9fa7b84fa42a6fb625228851385777155e22653e360e190ea765d9
hash_basis: LF-normalized bytes
---

# THM-2807 -- the graded semantic chord bounds a positive address triangle

**RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
INDEPENDENT AUDIT.**

THM-2782 supplies a positive physical central arm and THM-2791 extends it to
a two-point full-depth chord whose first nontrivial carry digit is `10`.
That chord could still have been an isolated positive edge.  It is not.

The exceptional target column `tau=12` contains a selected third positive
cylinder between the two THM-2791 endpoints.  On this three-vertex
restriction, all three arrows are exact translations of whole weighted
cylinders.  In the nerve of the physical address action they form a
commuting two-simplex:

```text
                         Z2^4079
             n_+ ----------------------> n_a
              ^                            ^
              |                            |
           Z1 |                            | Z1 Z2^4079
              |                            |
              +------------ n_0 -----------+
```

Here the slanted arrow is read from `n_0` to `n_a`.  On these three vertices,
the `tau=3` packet is exactly the diagonal face: it agrees with `tau=12` at
`n_0,n_a` and is empty at `n_+`.  The other `118` positive `tau=12`
cylinders are not deleted, classified, or claimed to lie in this simplex.

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
E3 -> Q_(3,{1,2}).                                    (2)
```

The delayed coefficient is separately the `D^6` carry-six coefficient; that
delay is not part of the semantic record in `(2)`.

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

Thus the THM-2791 two-point restriction on these addresses is obtained from
the selected three-vertex `tau=12` restriction by deleting precisely `n_+`.
The global `tau=12` support has `121` positive cylinders, so no global packet
deletion statement is intended.  Equation `(5)` is not an equality of
separately integrated marginals.

All three centres have semantic record `(3,{1,2})`, predecessor carry six,
source target-label `sigma=0`, and target label `tau=12`.  The two diagonal
centres additionally retain target label `tau=3`.  At all three centres the
complete semantic stability radius equals the open cylinder radius, so these
fixed labels and the semantic record are constant throughout the cylinders.
The full `sigma`-label banks differ, and no equality of those banks is
claimed.  The exact positive controls therefore live on one fixed typed
semantic sheet, not merely in the same numerical address orbit.

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

## 4. One depth-five class and thirteen full-depth affine lifts

THM-2788 also has the address dilation

```text
X:n |-> 14n.                                             (14)
```

The element `14=1+13` has order `13^4` modulo `13^5` and order `13^5`
modulo `13^6`.  Exhaustive exact discrete logarithm, equivalently repeated
use of the odd-prime lifting law, gives the unique depth-five exponent class

```text
k_0=23098 mod13^4,
14^k_0=53028 mod13^5.                                  (15)
```

Choose the least representative `k_0=23098`.  Its depth-six multiplier and
fixed-vertex translation are

```text
m_0=14^k_0 mod13^6=2652079,
v_0=(1-m_0)n_0 mod13^6=352469.                         (16)
```

There are not unique full-depth lifts.  For `t in F_13`, put

```text
k_t=k_0+t*13^4,
m_t=14^k_t mod13^6=m_0+t*13^5,
v_t=(1-m_t)n_0 mod13^6=v_0+6t*13^5,
g_t=O^v_t X^k_t.                                       (17)
```

The thirteen pairs `(m_t,v_t)` are distinct and all satisfy

```text
g_t(n_0)=n_0,                    g_t(n_+)=n_a,
g_t O^13 g_t^(-1)=O^689364.                            (18)
```

Indeed `m_t*13=689364 mod13^6` for every `t`.  Moreover, if `n=7 mod13`,

```text
g_t(n)-g_0(n)=t*13^5(n+6)=0 mod13^6.                  (19)
```

Thus all thirteen affine lifts agree on the entire low-digit-seven sheet
containing `n_0,n_+,n_a`, but they differ off that sheet; at `n=0`, for
example, consecutive lifts differ by `6*13^5`.  The induced conjugacy on the
selected sheet is unique, while its extension to the full depth-six address
space is a thirteen-element torsor.

For clarity, the least exponent has expansion

```text
k_0=10+8*13+6*13^2+10*13^3=(10,6,8,10)_13             (20)
```

in conventional high-to-low digit order.  Equations `(15)--(20)` are address
statements only.  They do not say that the physical factor packet is
covariant under any `X^k_t`, and they do not upgrade any `g_t` to a map of
endpoint currents.  The thirteen-lift ambiguity is structurally adjacent to
the missing endpoint-origin fibre, but this theorem does not identify those
torsors.  The issue is naturality of the retained factors and endpoint
allocation, not reachability of the addresses.

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
| affine comparison | unique depth-five class; thirteen full-depth lifts `(17)` agreeing on the selected sheet |
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
4. **No automatic physical intertwiner.**  The thirteen affine lifts `(17)`
   live in the address group.  Their common restriction to one residue sheet
   is not THM-2792's pointwise torsor map, which that theorem explicitly
   withholds.
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
in `(4)`, the literal interval equalities `(5)`, semantic/carry labels and
stability radii, all three translations, lower-central and quotient
arithmetic, the complete `13^4` discrete-log search in `(15)`, all thirteen
full-depth affine lifts, their common residue-seven restriction, and a
hostile point off that sheet.  It contains explicit exception gates, no
truth-bearing Python assertions, no floating point, and no scratch
dependency.

**Awaiting independent audit; not QED.**
