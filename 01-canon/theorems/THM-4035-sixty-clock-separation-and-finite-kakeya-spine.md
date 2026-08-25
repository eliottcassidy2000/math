---
id: THM-4035
title: "Sixty-clock separation and the finite four-dimensional Kakeya spine"
status: >
  PROVED arithmetic clock separation + FINITE-EXACT AP-selector fibres and
  finite-field transversality + VERIFIED-EXACT.  The THM-4029 phase-template
  labels, Fibonacci states modulo 10, and triangular states modulo 30 are
  pointed 60-cycles with different mechanisms.  Their scalar shadows do not determine
  the AP law.  Over F_61, three clocks chart the nonzero torus of P^3 and one
  clock has both a four-wise transverse twisted-cubic embedding and a wholly
  narrow planar embedding.  This proves no Euclidean or finite-field Kakeya
  conjecture and no new LRC case.
source: root + ap60_mechanism + period60_audit + kakeya4_bridge, 2026-08-24
depends_on:
  - THM-4029-lrc14-ap-cover-twelve-owner-rational-tail
related:
  - THM-536-lrc-seven-sector-sturmian-partial-sum-reframe
  - THM-1425
  - THM-4055-sixty-dyadic-response-fibre-law
  - THM-4026-sun-two-four-six-eight-binomial-counterexample
  - THM-4027-sun-two-four-six-eight-universal-modular-solubility
  - THM-4028-sun-two-four-six-eight-average-order-criticality
  - HYP-2235
  - MISTAKE-489
script: 04-computation/sturmian_fibonacci_triangular_kakeya_sixty_clock_thm4035.py
output: 05-knowledge/results/sturmian_fibonacci_triangular_kakeya_sixty_clock_thm4035.out
source_reference: 05-knowledge/reference/CORE-PAPERS-KAKEYA-4D-2026-08-24.md
script_sha256: 80e578c500f8d9fd8b007e5755d6ed457b1aa0d0b96eaaff76f20fa23805ea99
output_sha256: e9d29aa6d2683952605c123411e774ee9d40854ae19db5070f0f97d610300735
---

# THM-4035 -- sixty-clock separation and a finite Kakeya spine

**PROVED + FINITE-EXACT + VERIFIED-EXACT.**  The common number `60` gives an
exact common clock, but not a common cause.  A full state can relabel the clock;
a scalar shadow generally cannot.  The four-dimensional Kakeya connection is
an exact finite-field transversality control, not a transfer theorem.

## 1. Inherited object and the source of the AP clock

Let `a(m)` and `D(m)=1-a(m)` be the seven-sector AP-cover sequence of
THM-4029, and write `n=m-1`.  For `m>=13`, that theorem proves

```text
D(n+1)=(1/7) sum_j C_(r,j)/(n-c_(r,j)),
r=n mod 60,                 0<=c_(r,j)<=5.             (1)
```

The twelve persistent owners are the reduced rationals `p/q` with `q<=6`.
Their last occupied tracks are

```text
E_s(n)=n-((n-s) mod q).
```

Consequently the phase is the simultaneous track state for denominators
`1,...,6`, and

```text
60=lcm(1,2,3,4,5,6).                                  (2)
```

This `60` is created by rational-owner localization.  It is not created by
the aperiodicity of the underlying Sturmian/mechanical word.
THM-4029's onset boundary is retained: `m=12` is the sharp hostile predecessor,
and the twelve-owner formula begins at `m=13`.  Separately, MISTAKE-489 warns
that a finite-table endpoint cannot certify a span crossing; its corrected
row is `(7,8,10,13,26,infinity)`.

Aggregate equal shifts in `(1)` and call the resulting rational function
`L_r(n)`.  **FINITE-EXACT:** the sixty functions `L_r` are pairwise distinct;
there is no proper phase period.  Expanding at infinity gives

```text
L_r(n)=127/(35n)+A_1(r)/n^2+A_2(r)/n^3+O(n^-4),        (3)
A_k(r)=(1/7) sum_j C_(r,j)c_(r,j)^k.
```

The scalar `A_1` has 25 values and `A_2` has 42 values, but the ordered pair

```text
r |--> (A_1(r),A_2(r))                                 (4)
```

is injective on `Z/60`.  Thus two exact tail coefficients form a compact
phase fingerprint.  They still forget the individual owner, side, missing
sector, and winning track.

## 2. Fibonacci and triangular clocks

Put `F_0=0`, `F_1=1`, and `T_r=r(r+1)/2`.

### Fibonacci modulo ten

For

```text
A=[[0,1],[1,1]],
```

the state `(F_r,F_(r+1))` advances by `A`.  Modulo two, `A` has order three.
Modulo five,

```text
A^5=3I,
```

so `A^20=I`.  The proper divisors `1,2,4,5,10` of 20 do not return the
identity: the first three follow by direct multiplication, while
`A^5=3I` and `A^10=4I`.  Hence `A` has order 20 modulo five.
The Chinese remainder theorem gives

```text
pi(10)=lcm(pi(2),pi(5))=lcm(3,20)=60.                 (5)
```

The starting column `v=(0,1)^T` is cyclic: `v` and `Av=(1,1)^T` form a basis
modulo ten because their determinant is `-1`, a unit.  If `A^P v=v`, then
`A^P(Av)=A(A^P v)=Av`, so `A^P` fixes a basis and equals the identity.
The state orbit therefore has the full matrix order 60, and

```text
r |--> (F_r,F_(r+1)) mod 10                            (6)
```

is a bijection from `Z/60` to its pointed orbit and conjugates `r|->r+1` to
the Fibonacci state transition.

### Triangular numbers modulo thirty

For every modulus `M`, a period `P` must satisfy

```text
T_(r+P)-T_r=P r+P(P+1)/2 = 0 mod M                    (7)
```

for every `r`.  Subtracting consecutive instances forces `M|P`.  If `M` is
odd, `P=M` works.  If `M` is even, an odd multiple of `M` leaves the constant
term congruent to `M/2`, whereas `P=2M` works.  Thus the minimal period is

```text
tau(M)=M  (M odd),             tau(M)=2M  (M even).   (8)
```

In particular `tau(30)=60`.  The pointed affine state

```text
(u,t)=(r mod 30,T_r mod 30),
(u,t) |--> (u+1,t+u+1) mod 30                         (9)
```

has an injective 60-cycle.  Equivalently,
`r|->(T_r,T_(r+1)) mod 30` is injective: the difference determines `r mod 30`,
and the two possible lifts differ in the first coordinate by 15.

The modulus is load-bearing.  Triangular last digits have period 20, not 60;
`T_r mod 60` has period 120.

## 3. Exact quotient-loss audit

The AP phase-template labels `r|->L_r`, the full Fibonacci states in `(6)`,
and the full triangular modular states in `(9)` are conjugate to the same
pointed successor cycle `Z/60`.  This is an exact address conjugacy, not an
identity of their evaluators or a derivation of one law from another.  The
evaluated AP value `D(n+1)=L_(n mod 60)(n)` also retains the unbounded height
`n` and is not itself periodic.

The scalar quotients expose the distinction:

| observable on `r mod 60` | image size | largest fibre | determines `L_r`? |
|---|---:|---:|---|
| `F_r mod 10` | 10 | 8 | no |
| `T_r mod 30` | 12 | 8 | no |
| `(F_r mod 10,T_r mod 30)` | 48 | 2 | no |
| `(F_r mod 10,T_r mod 30,r mod 2)` | 60 | 1 | yes, as an address |
| `(A_1(r),A_2(r))` | 60 | 1 | yes, from the AP law itself |

There are exactly twelve doubleton fibres in the combined arithmetic scalar.
Every doubleton joins opposite parities, so one bit is sufficient and, for
that quotient, necessary.  The cheapest hostile is

```text
r=0,15: (F_r mod 10,T_r mod 30)=(0,0),
A_1(0)=37/14 != 31/14=A_1(15).                         (10)
```

The full Fibonacci state is therefore a reversible *address* for the AP
phase.  It does not explain the AP coefficients: any pointed 60-cycle would
provide such an address.  This is the precise survivor of the apparent
periodicity bridge.

## 4. A four-dimensional finite-field transversality spine

In `F_61`, take `sqrt(5)=26` and

```text
phi=(1+sqrt(5))/2=44.
```

Then `phi^2=phi+1` and `phi` has multiplicative order 60.  Thus the nonzero
affine torus in projective direction space has the exact chart

```text
(a,b,c) |--> [1:phi^a:phi^b:phi^c],
(Z/60)^3  -->  P^3(F_61),                              (11)
```

whose image has `60^3=216,000` directions.  The full projective space has

```text
|P^3(F_61)|=61^3+61^2+61+1=230,764.                   (12)
```

This gives a dimensionally typed statement: a generic nonzero chart in
four-dimensional direction space needs **three independent** 60-clocks.  One
clock supplies only a curve.

The direction parameter is itself a Fibonacci-state observable:

```text
t_r=phi^r=F_(r+1)+43F_r mod 61.                       (13)
```

There are two sharply different equivariant embeddings of one clock.  With
`t_r=phi^r`, define

```text
B(r)=[1:t_r:t_r^2:t_r^3],
N(r)=[1:t_r:1:t_r].                                   (14)
```

For distinct `r_1,...,r_4`, the four row representatives of `B` have

```text
det(B(r_i))=product_(i<j)(t_(r_j)-t_(r_i)) != 0.       (15)
```

Hence every one of the `binom(60,4)=487,635` phase quartets is four-wise
transverse.  By contrast, every `N(r)` lies in the fixed vector plane
`x_0=x_2`, `x_1=x_3`, so all `487,635` four-determinants vanish.  The same
clock therefore supports both maximally broad and wholly narrow geometry.

There is a further projective warning.  The Fibonacci matrix itself has order
60 modulo 61, but

```text
A^15=11I mod 61.                                      (16)
```

The two eigenvalues are `phi=44` and `psi=18`; their ratio is 16, whose order
is 15 (`16^3=9`, `16^5=47`, `16^15=1`).  Hence 15 is the least positive
scalar power.  As above, the starting vector and its image form a basis, so a
projective return of that vector forces a scalar matrix power.  Its 60-state
vector orbit therefore has exactly 15 projective directions.  Scalar state
advance disappears under projectivization: a recurrence clock is not a
direction atlas unless this gauge loss is audited.

The triangular shadow supplies a second hostile quotient:

```text
u_r=phi^(2T_r)=phi^(r(r+1)).                           (17)
```

This sequence has minimal period 60 but only 12 distinct values.  Exact
periodicity therefore does not imply 60 distinct directions, much less a
Kakeya direction cover.

The same multiplicative clock parameterizes the nonzero input residues for
the four binomial roles in THM-4026--4028.  If

```text
S_k={C(t,k) mod 61:t in F_61},       k=2,4,6,8,
```

then exact enumeration gives image sizes `(31,24,26,24)`, and already
`S_2+S_4=F_61`.  The Sun counterexample is `21 mod 61`, so this clock gives
no congruence obstruction even before the degree-six and degree-eight roles
enter.  It is an input atlas, not an explanation of the height-sensitive
empty integer fibre.

## 5. Typed Kakeya relevance and boundaries

The Euclidean four-dimensional Kakeya set and maximal-function conjectures
remain **OPEN**.  The cited multilinear theory detects broad configurations
through wedge volumes, while the four-dimensional planebrush handles a
plany/narrow regime.  Equation `(15)` is the exact finite-field analogue of
the broad determinant predicate; equation `(14)` supplies its hostile narrow
control.  The current external status and imports are pinned in
[`CORE-PAPERS-KAKEYA-4D-2026-08-24.md`](../../05-knowledge/reference/CORE-PAPERS-KAKEYA-4D-2026-08-24.md).

| source | target/map | preserved predicate | destroyed information | required sidecar |
|---|---|---|---|---|
| AP phase `r` | full Fibonacci or triangular state | successor and phase address | AP owners and coefficients | owner/track signature or `(A_1,A_2)` |
| `(Z/60)^3` | torus chart `(11)` | nonzero projective direction | zero-coordinate/boundary charts | full projective atlas |
| `Z/60` | twisted cubic `B` | four-wise determinant nonvanishing | almost every direction and every basepoint | direction completion and affine lines |
| `Z/60` | planar control `N` | phase injectivity | four-fold transversality | determinant/plane label |
| directions | Kakeya tube family | nothing without another map | placement, overlap, shading and scale | basepoints, multiplicity, two-ends and multiscale nonconcentration |

In particular, the 60-point spine is not a finite-field Kakeya set: it covers
only 60 of 230,764 directions and specifies no line in any direction.  A
finite-field calculation does not import Euclidean angles, tube volume, or
induction on scale.

## 6. Frontier consequences

For LRC(14), `(4)` supplies a compact exact coordinate for the AP tail, but it
does not repair the AP-container loss.  A sparse set still needs occupied
tracks, gaps, anchoring and owner placement.  The first decisive extension is
to replace the consecutive horizon `E_s(n)` by a sparse occupied-track ledger
and test exactly when the no-skip propagation survives.

The finite Kakeya control opens three bounded experiments without claiming a
transfer:

1. Fourier-decompose the six owner-denominator blocks on `C_60` and ask
   whether sparse-gap operations preserve any nontrivial character packet.
2. Generalize the seven-sector theorem: for `K` sectors, test whether the
   eventual owners are exactly `q<K`, whether the candidate phase is
   `lcm(1,...,K-1)`, and whether that period is minimal.
3. Extend HYP-2235 at `p=61` by phase-coloring a genuine one-line-per-direction
   family using the three-clock torus atlas, while retaining concurrency,
   pinned multiplicity, boundary charts, and determinant rank.
4. For the Sun 2-4-6-8 problem, retain the four `C_60` input phases at `p=61`
   and measure bounded-lift cost and exact fibre multiplicity after the output
   residue has already saturated.  THM-4027 forbids treating any fixed
   residue clock as a pointwise obstruction.

None of `(1)--(17)` proves AP extremality for sparse sets, LRC(14), a Kakeya
volume or dimension bound, the Euclidean Kakeya conjecture, or the repo's
distinct Arithmetic-Kakeya certificate target.

## 7. Replay

From the repository root:

```text
python3 -B 04-computation/sturmian_fibonacci_triangular_kakeya_sixty_clock_thm4035.py
python3 -B -O 04-computation/sturmian_fibonacci_triangular_kakeya_sixty_clock_thm4035.py
```

The two output streams are byte-identical.  **QED.**
