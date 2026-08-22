---
id: THM-3679
title: "LRC p-adic scalar seam and total-speed checksum"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  At every
  depth 13^n, four
  packet charts for each of two distinct blocker sources force a jointly
  target-zero address to be globally scalar.  On the all-coordinate-unit
  relation fibre this scalar hostile exists exactly when 13^n divides the
  total speed sum.  Hence for every positive integer row it dies by depth
  nu_13(sum w_i)+1.  For one source alone, its all-unit blind tower persists
  at every depth exactly when the source speed and the sum of the other eight
  speeds have equal 13-adic valuation; otherwise it dies at the first unequal
  layer.  In the equal-valuation case one explicit primitive exact relation
  realizes the hostile tower, whereas the combined two-source chart map is
  injective on the exact relation lattice.  THM-2334 realizes each higher
  difference pushforward by exact
  coordinate twists at a sufficiently delayed word clock.  At any depth where
  the all-unit blind fibre is empty, nonzero total mass of the all-unit
  restricted current forces a nonzero joint-chart target without a phase-cone
  premise.  Reduction compatibility at N_n=7*13^n identifies the all-91-unit
  total at every depth with the base scalar sum_q B(q), while the higher
  character twists remain septimally charged.
  No one-packet marginal, covering-row exclusion, or LRC(14) conclusion is
  proved.
source: kps-s195 / THM-3676 Bockstein continuation, 2026-08-21
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-3676-lrc-cross-source-blind-intersection-and-scalar-seam
related:
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-3658-lrc-mod169-carry-fourier-block-intertwiner
  - THM-3676-lrc-cross-source-blind-intersection-and-scalar-seam
script: 04-computation/lrc_padic_scalar_seam_total_speed_checksum_thm3679.py
output: 05-knowledge/results/lrc_padic_scalar_seam_total_speed_checksum_thm3679.out
script_sha256: cb9e80f688254bf6c83927c20306e33a3e716ddd16e1e70febf577cc96baecde
output_sha256: c357dc4643e251b42807d1c086b27928d2c3b1f2ffa9b1c3c9875f402d4c2c8d
semantic_sha256: 7afb7edb68614deebafcf4d544ffa6c4e28be94ed1b656ddb269aead70b5fb77
hash_basis: raw LF bytes
audit: >
  lovelace-2026-08-21 independently reconstructed the finite-ring kernels,
  checksum criterion, one-source valuation law, exact hostile relation, and
  higher-depth reduction compatibility.  The hostile audit also supplied a
  pure four-way ANOVA measure showing that every one-chart marginal may
  cancel although the joint pushforward is nonzero, and it isolated the
  remaining role-to-packet intertwiner.  No theorem conclusion relies on a
  deepest blocker carrying a nonzero semantic current.
---

# THM-3679 -- every positive row eventually fails the scalar checksum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
lifts the exact
equality graph of THM-3676 from `F_13` to every finite ring
`Z/13^n Z`.  It identifies a finite depth at which the last mod-thirteen
scalar hostile must disappear.  It does not prove that the semantic LRC
current is transported to that depth.

The independent audit reconstructed the ring argument without field division,
the exact integer kernel, and the reduction tower.  It then attacked the
strongest tempting marginal conclusion: an explicit signed all-unit measure
has nonzero total and nonzero joint target but zero pushforward through each
of the four individual charts.  Section 6 records this sharp boundary.

## 1. Higher packet differences

Retain the nine labelled coordinates

```text
U={0,1,2,3,4,5},                 B={1',2',3'},       (1)
```

and let `w in Z_(>0)^9` be the fixed speed row.  The six labels in `U` are
units modulo thirteen and the three blocker speeds are divisible by thirteen,
as in THM-3676.  For `n>=1` put

```text
A_n=Z/13^n Z,
K_n={r in A_n^9:sum_i r_i w_i=0 in A_n}.            (2)
```

Throughout this theorem, an address coordinate called a **unit** means a unit
of `A_n`, equivalently nonzero modulo thirteen.  The all-`91`-unit bank implies
this condition, but the converse need not hold modulo seven; no septimal
support claim is made here.

Fix a selected source `j in B`, write `{a,b}=B minus {j}`, and for distinct
`k,l in U` define the full depth-`n` packet difference

```text
pi^j_(k,l,n)(r)=(r_a-r_k,r_b-r_l) in A_n^2.         (3)
```

At `n=1`, this is the target-difference map underlying THM-3671/3676.  At
larger `n`, equation (3) is first a relation-address refinement; Section 6
records its exact coordinate-twist realization and the remaining semantic
qualification.

Use the four charts

```text
C={(0,1),(2,3),(4,5),(1,0)}.                       (4)
```

Their equality graph is connected on `U union {a,b}` over every coefficient
ring.  Consequently

```text
pi^j_(k,l,n)(r)=0 for every (k,l) in C

iff

r_i=c_j for every i!=j,                            (5)
```

where the source coordinate `r_j` remains free.

## 2. Two sources force the global scalar diagonal

Let `j,j' in B` be distinct and impose (5) for both source-dependent packet
families.  A unit coordinate belongs to both outside-source sets, so their
two common values agree.  The coordinate `j` belongs to the outside-source
set for `j'`, and conversely for `j'`.  Hence every coordinate agrees:

```text
intersection_(s in {j,j'}) intersection_((k,l) in C)
  ker(pi^s_(k,l,n))

 ={c 1:c in A_n}.                                  (6)
```

This identity is exact over `A_n`; no field rank or division is used.  Adding
the third source family changes neither side.

## 3. The relation equation becomes one checksum

Intersect (6) with `K_n`.  A scalar address `c 1` is a relation exactly when

```text
c W=0 mod 13^n,                 W=sum_i w_i.         (7)
```

Now impose the all-coordinate-unit condition inherited from the all-`91`-unit
address bank:

```text
r_i not congruent 0 mod 13 for every i.             (8)
```

On the scalar diagonal, (8) says that `c` is a unit of `A_n`.  Multiplication
by `c` is then invertible, so (7) is equivalent to

```text
13^n divides W.                                    (9)
```

Therefore the exact jointly blind all-unit fibre is

```text
{c 1:c in A_n^times},       if 13^n divides W,

empty,                      otherwise.             (10)
```

At `n=1`, (10) recovers precisely THM-3676's zero-sum scalar seam.  At
`n=2`, it says that the first divided-difference layer kills that seam unless

```text
169 divides sum_i w_i.                              (11)
```

## 4. Bockstein tower and finite death depth

Suppose a packet difference is zero modulo `13^h`.  Its next divided
difference

```text
(pi^j_(k,l,h+1)(r)/13^h) mod 13                    (12)
```

is well-defined, and vanishes exactly when the packet difference is zero
modulo `13^(h+1)`.  Thus simultaneous vanishing of the ordinary target and
the first `n-1` divided-difference layers is exactly the left side of (6) at
depth `n`.  Equations (9)--(10) are therefore a complete scalar Bockstein
tower, not merely a mod-169 coincidence.

Because every LRC speed is positive,

```text
0<W<infinity.
```

Put `v=nu_13(W)`.  Then the scalar all-unit hostile exists through depth `v`
and is empty at depth

```text
n_*=v+1.                                            (13)
```

No positive speed row can remain jointly scalar-blind at every 13-adic
depth.

## 5. One source often kills its own blind tower

The two-source checksum is universal, but many rows do not need an owner
handoff.  Fix one source `j`, put

```text
S_j=W-w_j=sum_(i!=j)w_i,
a=nu_13(w_j),                    b=nu_13(S_j).       (14)
```

Both valuations are finite because all speeds are positive.  By (5), a
one-source blind all-unit address has one unit value `c` off the source and
another unit value `d` at the source.  Its relation equation is

```text
c S_j+d w_j=0 mod 13^n.                            (15)
```

After division by the unit `c`, put `u=d/c`.  If `n<=min(a,b)`, both terms in
`S_j+u w_j` vanish modulo `13^n`, so a blind lift exists.  If
`n>min(a,b)` and `a!=b`, exactly one term has the smaller valuation and
cancellation is impossible.  Finally, if `a=b=lambda<n`, divide by
`13^lambda`; the two quotients are units, and

```text
u congruent -(S_j/13^lambda)(w_j/13^lambda)^(-1)
  mod 13^(n-lambda)
```

has a unit lift.  Therefore the exact criterion is

```text
one-source blind all-unit fibre at depth n is nonempty

iff

n<=min(a,b) or a=b.                                 (16)
```

When `a!=b`, the first empty depth is `min(a,b)+1`, and every deeper fibre is
empty.  A particularly useful corollary follows from the non-Archimedean
triangle law.  If

```text
nu_13(W)<nu_13(w_j),                                (17)
```

then `nu_13(S_j)=nu_13(W)`, so this one fixed source loses its all-unit blind
fibre already at depth `nu_13(W)+1`.  Taking `j` to be the deepest blocker
therefore gives the cheapest *algebraic* one-source split on every row whose
total-speed valuation is smaller than that blocker depth.  Its physical use
additionally requires that this source own a nonzero marked current.  The
canonical semantic construction universally guarantees a shallow owner, not
a deepest owner, so no deepest-source handoff is inferred here.

The criterion is sharp because it is the reduction of an exact integer
relation.  Put

```text
g_j=gcd(w_j,S_j).
```

All integer solutions of the one-source blind equation are multiples of the
primitive vector whose off-source coordinate is `w_j/g_j` and whose source
coordinate is `-S_j/g_j`.  Both coordinate values are nonzero modulo thirteen
exactly when

```text
nu_13(w_j)=nu_13(S_j).
```

Thus the persistent branch of (16) is realized at every depth by one exact
all-unit relation, rather than by unrelated finite-ring lifts.

Across two sources, the integer equality graph forces `r=c 1`.  The exact
relation equation is then `cW=0` in `Z`.  Since `W>0`, it forces `c=0`.
Consequently the combined two-source chart map is injective on the exact
relation lattice.  Equations (9)--(13) quantify how many 13-adic digits are
needed before that characteristic-zero injectivity becomes visible.

## 6. Exact twist realization and semantic boundary

THM-2334 permits reduction of the exact relation-address current modulo every
integer `N`.  At `N=13^n`, push that current through (3).  Fourier inversion
identifies its two-dimensional transform with the coordinate twists

```text
ell=s(e_a-e_k)+t(e_b-e_l),          s,t in A_n,      (18)
```

which translate the nine base coordinates by `ell_i/13^n`.  If the delayed
word clock is `R=13^K` with `K>=n`, its translated coordinate is shifted by

```text
R ell_i/13^n in Z,
```

so every transported word factor is unchanged.  THM-2365 permits the fixed
positive word clock to be chosen arbitrarily large.  Thus each higher packet
difference has an exact physical coordinate-twist profile on any one fixed
owner current; no new carry-convolution identification is needed.

The analytic nonvanishing premise does not grow with the depth.  Let
`C_(rho,13^n)` be THM-2334's Abel residue current and let `U_n` be its
all-coordinate-unit locus.  Absolute regrouping at `rho<1` gives

```text
M_(rho,n)=sum_(r in U_n)C_(rho,13^n)(r)

 =sum_(lambda in Lambda;
       lambda_i not congruent 0 mod 13 for every i)
    C_rho(lambda),                                  (19)
```

because being a unit modulo `13^n` depends only on the first digit.  The
right side is independent of `n`.  Each finite residue sum has a boundary
limit by THM-2334, so

```text
M_n=M_1                    for every n>=1.           (20)
```

Thus, for one fixed delayed current, its nonzero mod-thirteen all-unit total
mass can be regrouped directly at the finite checksum-killing depth; no
separate nonvanishing theorem is needed at each Bockstein layer.  Equation
(20) does not compare currents built with two different delayed clocks.

There is a sharper version which retains the septimal all-unit projector.
Put

```text
N_n=7*13^n.
```

Reduce the same exact relation current modulo `N_n`.  Reduction
`K_(N_(n+1))->K_(N_n)` pushes `C_(rho,N_(n+1))` to `C_(rho,N_n)`.  Moreover,
an exact relation coordinate is a unit modulo `N_n` if and only if it is
coprime to `91`; this condition is independent of `n`.  Hence

```text
T_n^91
 :=sum_(y in K_(N_n); every y_i in (Z/N_n Z)^times)
      C_(N_n)(y;1-)

 =sum_(lambda in Lambda; gcd(lambda_i,91)=1 for every i)
      C(lambda)

 =T_1^91.                                           (20a)
```

At depth one, THM-2334 partitions this same all-`91`-unit current by the
mod-thirteen packet target, so

```text
T_n^91=T_1^91=sum_(q in (Z/13 Z)^2) B(q).           (20b)
```

Thus the p-adic checksum argument needs only the single base scalar premise

```text
sum_q B(q)!=0.                                      (20c)
```

This is strictly a compatibility statement for finite residue regroupings of
one fixed current.  It does **not** make the higher physical twists neutral
modulo seven: THM-2334 shows that their word charge remains septimally active.
Consequently (20b) neither supplies a role-to-packet intertwiner nor identifies
currents belonging to two different owner words.

There is also an exact joint-chart consequence which avoids the phase-cone
premise of THM-3676.  Let `Pi_(j,n)` collect the eight coordinates of the four
maps (3), and let `Pi_(j,j',n)` collect all sixteen coordinates for two
sources.  Let `mu` be any complex residue measure on `K_n`, supported on the
all-coordinate-unit locus, and assume

```text
M=sum_(r in K_n)mu(r)!=0.                            (21)
```

At the first empty depth in (16), the zero fibre of `Pi_(j,n)` does not meet
the support.  At depth `n_*` in (13), the same is true for
`Pi_(j,j',n_*)`.  In either case the joint pushforward has value zero at its
zero label, while the sum of all its fibre aggregates is `M`.  Hence some
nonzero joint label has nonzero aggregate.  No sign or common half-plane is
used.  The complete character group of this joint pushforward consists of
arbitrary linear combinations of the dipoles in (18), so THM-2334 realizes
the conclusion by finitely many physical coordinate twists.

What remains missing is the cross-source semantics.  THM-2305 does not yet
transport literally one nonzero grouped all-unit residue current, in the same
labelled coordinates and phase gauge, through two distinct legitimate
owner-source packet families.  Merely applying the second set of coordinate
twists to the first owner's current does not prove that it is the second
owner-aligned packet current.  The new exact interface is

```text
valuation-mismatch route:
  one nonzero all-unit current at the first empty depth in (16)
  + one-source noncancellation
  -> some higher target is nonzero;

valuation-aligned route:
  common all-unit current at depth n_*
  + sourcewise noncancellation
  + two source packet actions
  -> some higher target is nonzero.                              (22)
```

Equation (13) proves that there is no algebraic scalar hostile left at that
depth, equation (16) gives the cheaper one-source split, and (18) realizes
each action separately.  The theorem does not supply a nonzero all-unit
total mass (21), nor the common-current and one-packet marginal premises needed
in the aligned route (22).  A nonzero joint aggregate can cancel after
marginalizing away the other chart coordinates.  No scalar row is excluded
and LRC(14) remains open.

The marginal warning is sharp even before analytic semantics are imposed.
The independent hostile audit uses

```text
w=(1,2,3,4,5,6,0,0,0) mod 13
```

and constructs four relation vectors `v_i`, one visible in each chart and
zero in the other three.  Every nonempty subset sum `r_S=sum_(i in S)v_i` can
be normalized to an all-coordinate-unit relation.  The signed measure

```text
mu(r_S)=(-1)^(|S|+1),              empty!=S subset {1,2,3,4}
```

has total mass one and a nonzero joint pushforward, but every individual
chart pushforward is `delta_0`.  This is the pure four-way ANOVA interaction:
joint noncancellation cannot be marginalized chart by chart without a new
structural property of the canonical interval-factor coefficients.

## 7. Exact finite-ring controls

The companion computes Smith normal forms over `Z`, not ranks over a field.
Each one-source chart matrix has seven unit invariant factors and nullity two;
every two-source matrix has eight unit invariant factors and nullity one, and
the third source adds no rank.  Thus (5)--(6) hold after base change to every
coefficient ring.

Independently, it exhausts the unit equation (15), the scalar checksum (7),
and the divided-difference zero upgrade over `Z/3^n Z` through `n=4`.  The
finite ledger contains `12,651` hostile cases and `13,007` active gates.  The
odd prime three is deliberately distinct from the LRC prime thirteen; these
are finite-ring controls for the proof mechanism, not an exhaustive proof of
the all-`n` statement.

```bash
python -B 04-computation/lrc_padic_scalar_seam_total_speed_checksum_thm3679.py
python -B -O 04-computation/lrc_padic_scalar_seam_total_speed_checksum_thm3679.py
```

Normal and optimized transcripts are byte-identical to the pinned output.
**QED.**
