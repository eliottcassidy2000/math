---
id: THM-3048
title: "Prime-holonomy parallel-correspondence arrival and norm descent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For a nonnegative correspondence between two C_p root torsors, the product
  of the pointed determinants over all p equivariant pointings is a canonical
  unpointed norm whose first log ghost is total root mass.  If the whole
  correspondence is parallel around a loop with nonzero relative holonomy
  and p is prime, holonomy acts transitively on the pointing torsor: all
  pointed determinants coincide, the norm is their pth power, and every
  pointing has arrival exactly total root mass/p.  At p=13, if THM-2591's
  root-deck holonomy is lawfully inherited as the endpoint relative holonomy,
  its 7a class would force arrival from a parallel physical table.
  Current canon has vertexwise selectors but no such whole-table common
  carrier, so this is a sharp conditional bridge and excludes no LRC row.
source: kind-pasteur-2026-08-01-holotopy-arrival
audit: >
  An independent immutable-file hostile audit ACCEPTED the correspondence and
  pointing operator types, relative-holonomy sign, complete determinant
  conjugacy, norm and first-ghost normalizations, prime/generator boundary,
  conditional THM-2591 application, Reynolds physicality tax, and the zero,
  composite, nonparallel, affine-slope, cemetery, and semantic-horn hostiles.
  It replayed normal, optimized, and stored output byte-for-byte, matched both
  LF hashes, and passed the documentation checker after startup compaction.
depends_on:
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-3044-pointed-correspondence-determinant-hall-dual-and-cycle-boundary
related:
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2565-target-active-self-return-and-future-root-overlap-audit
  - THM-2604-unshifted-future-root-accessibility-and-selector-cross-mixing-boundary
  - THM-2613-canonical-root-diagonal-opposite-shift-section
  - THM-2835-q11-semantic-word-horn-and-bockstein-blind-support-no-go
script: 04-computation/lrc14_prime_holonomy_parallel_correspondence_arrival_thm3048.py
output: 05-knowledge/results/lrc14_prime_holonomy_parallel_correspondence_arrival_thm3048.out
script_sha256: a42bca68fcb43ccc088d4b75bdda22046602b2e4c4d2872ef7d3f91f1ba65eed
output_sha256: dfa203df7c842027f4d3d76a3038b569ffe6c5cd196d470e9e67e69cb5880f2d
hash_basis: LF-normalized bytes
---

# THM-3048 -- prime holonomy can force pointed arrival

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The root-holotopy frontier has treated nonzero chart holonomy as an obstacle
to choosing one global root pointing.  For a whole correspondence table the
same holonomy has an opposite use.  If the table is parallel around the loop,
holonomy permutes its possible pointings while preserving their determinants.
At prime order, any nonzero holonomy visits every pointing.  The ambiguity is
then harmless because every choice has the same arrival, strictly positive
whenever `M_root>0`.

This mechanism requires a transported **table**, not one selected root in
each chart.  That distinction is exactly where current LRC canon stops.

## 1. The canonical unpointed pointing norm

Let `K=C_p`, and let `H,B` be fixed-action `K`-torsors.  Let

```text
C:R^B -> R^H                                             (1)
```

be a nonnegative root-to-root correspondence, written as a matrix with rows
indexed by `H` and columns by `B`.  Cemetery columns, if present in a larger
THM-2545 table, are omitted from `C`; put

```text
M_root=sum_(h in H,b in B) C(h,b).                       (2)
```

By THM-2611, the fixed-action equivariant bijections
`iota:H->B` form a `K`-torsor of cardinality `p`.  Let
`J_iota:R^H->R^B` be the corresponding permutation operator.  Following
THM-3044, define

```text
D_iota(t)=det(I+t C J_iota),
N_C(t)=product_(iota in Iso_K(H,B)) D_iota(t).           (3)
```

The product `N_C` is canonical without choosing a pointing.  Independent
changes of the `H` and `B` gauges permute `Iso_K(H,B)` and conjugate the
endomorphisms `C J_iota`, so `(3)` is unchanged.  Its first formal log ghost
is

```text
[t]log N_C(t)
 =sum_iota tr(CJ_iota)
 =sum_(h,b) C(h,b)
 =M_root.                                                (4)
```

For each pair `(h,b)`, exactly one equivariant pointing sends `h` to `b`,
which proves the middle equality.  Equation `(4)` is a total-mass identity,
not yet semantic localization: the product has deliberately averaged over
all pointings.

## 2. Parallel transport turns holonomy into arrival

Choose temporary origins and write the endpoint translations accumulated
around one chart loop as

```text
P_H=P_u on H,              P_B=P_v on B.                (5)
```

The relative pointing holonomy is

```text
g=v-u in C_p.                                             (6)
```

Assume the entire physical correspondence is parallel around the loop:

```text
C=P_H C P_B^(-1).                                       (7)
```

Condition `(7)` is the closed-loop form of an edgewise family transported on
one common carrier and returning to the same semantic table.  It is much
stronger than choosing one positive vertex in every chart.

Index the pointings by `c in C_p`, with `iota_c(h)=h+c`.  Their permutation
operators obey

```text
J_(c+g)=P_B J_c P_H^(-1).                               (8)
```

Using `(7)` gives the exact conjugacy

```text
CJ_(c+g)
 =C P_B J_c P_H^(-1)
 =P_H(CJ_c)P_H^(-1).                                   (9)
```

Consequently

```text
D_(c+g)(t)=D_c(t)                                      (10)
```

as complete polynomials, and every formal log ghost is constant on the
subgroup orbit `c+<g>`.

If `p` is prime and `g!=0`, that orbit is all of `C_p`.  Hence there is one
polynomial `D(t)` such that

```text
D_c(t)=D(t) for every c,
N_C(t)=D(t)^p.                                          (11)
```

The formal `p`th root in `(11)` is unique because its constant term is `1`.
Combining `(4)` and `(11)` yields

```text
[t]log D_c(t)=M_root/p              for every c.        (12)
```

Thus if `M_root>0`, every equivariant pointing has positive semantic arrival.
One does not trivialize the nonzero holonomy; one uses its transitivity to
make the pointed answer independent of the missing section.

## 3. The THM-2591 alternative

THM-2591 proves that its seven-chart root-deck system has holonomy

```text
g=7a !=0 in F_13                                      (13)
```

for every marker `a in F_13^*`.  Suppose in addition that the head and later
endpoint torsors of one THM-3044 correspondence are lawfully placed on those
same physical charts so that their induced **relative** pointing holonomy is
`g=7a` (or `-7a` after reversing orientation).  Then any nonzero nonnegative
root correspondence which is parallel around that loop satisfies

```text
semantic arrival=M_root/13>0.                           (14)
```

This is an alternative to THM-2591's mixed-square flattening invoice.  Either
construct a correction which kills `(13)`, or retain `(13)` and construct one
whole-table parallel correspondence.  In the latter route the thirteen-fold
cover is unnecessary for the determinant: the unpointed norm descends as a
perfect thirteenth power.

Current canon does not supply the premise.  THM-2591 explicitly grants only
vertexwise positive selectors and not an identification of the two physical
chart carriers.  THM-2565, THM-2604, and THM-2613 provide genuine pointings on
scoped packets, but none proves that one full THM-2545 semantic/root table is
parallel around the THM-2591 loop while retaining the required word,
endpoint, address, and repair fields.  Equation `(7)`, not another marginal
or selected root, is the cheapest decisive test.

## 4. Reynolds construction and its physicality tax

Algebraically, any seed correspondence `C_0` can be made parallel by the
orbit sum

```text
Cbar=sum_(j=0)^(p-1) P_H^j C_0 P_B^(-j).               (15)
```

Indeed `P_H Cbar P_B^(-1)=Cbar`.  When `g!=0` and `p` is prime,
`M_root(Cbar)=p M_root(C_0)`, so `(12)` becomes

```text
[t]log D_c(Cbar;t)=M_root(C_0)              for all c. (16)
```

This is a lawful LRC construction only if all `p` transported copies in
`(15)` are separately physical on the same semantic stratum and may be added
before the determinant is formed.  If `P_H,P_B` are merely chart gauges, the
sum is artificial.  It is then the exact analogue of a Reynolds duty term:
formally covariant, but not a common-carrier current.  The theorem does not
infer physicality from the algebraic average.  Sharply, with `P_H=I`,
`P_B=P_1`, and `C_0=P_1`, the seed has zero identity-pointed arrival while
`Cbar` is the all-ones matrix and has arrival `13` at every pointing.  The
orbit sum has enlarged one displacement to all thirteen, and the determinant
must be taken only after that support-changing sum.

## 5. Sharp failure boundaries

### Zero holonomy

At `g=0`, `(10)` compares each pointing only with itself.  On `C_13`, the
positive parallel table `C=P_1` has total root mass `13` but zero arrival at
the identity pointing.  Thus nonzero holonomy is load-bearing.

### Failure of whole-table parallelism

Keep `g=1`, but take `P_H=I`, `P_B=P_1`, and `C=P_1`.  The identity-pointed
arrival is again zero, while `(7)` fails.  Vertexwise compatible labels do not
replace table covariance.

### Composite-order holonomy

On `C_6`, take `g=2` and the circulant nonnegative table whose displacement
ledger is the indicator of `{1,3,5}`.  It is parallel, but its six pointed
first ghosts are

```text
(0,6,0,6,0,6).                                         (17)
```

Holonomy is nonzero yet has two orbits.  Prime order, or more generally
`<g>=K`, is the exact transitivity hypothesis.

### Fixed-action pointings

The fixed-generator clause from THM-2611 is also load-bearing.  On `F_13`,
take `P_H=P_1`, `P_B=P_2`, and

```text
C(h,b)=1_(b=2h).                                        (17a)
```

This table is parallel and every translation pointing has arrival `1`, as
the theorem predicts.  If affine pointings of slope `2` are admitted, their
arrival ledger becomes `(13,0,...,0)`.  The theorem compares the thirteen
fixed-action identifications, not the larger affine group.

### Cemetery mass

If the larger correspondence has positive mass only in a cemetery column,
then `M_root=0` and `(12)` has no positive conclusion.  Root-root mass is
load-bearing.

### The unpointed norm is displacement-blind

For `C=uP_d` on `C_p`, with odd prime `p`, the norm is independent of `d`:

```text
N_C(t)=(1+tu)^p(1+(tu)^p)^(p-1).                       (18)
```

One pointing gives the identity permutation and the other `p-1` pointings
give `p`-cycles.  Formula `(18)` proves that `(4)` alone does not locate the
semantic factor.  Parallel nonzero holonomy is what turns the blind norm into
the pointed root `(11)`.

### A root pointing is not a semantic return

THM-2835's `q=11` word horn has a square-zero semantic operator: the two
source words point into `QA`, with no outgoing `QA` edge.  Tensoring that
operator with any perfectly pointed root endomorphism remains square-zero, so
`det(I+tN_11 tensor A)=1`.  The correspondence in `(7)` must itself be the
semantic head-to-later-root table.  A separate root pointing cannot create a
missing semantic return arrow.

## 6. Exact companion

The dependency-free exact companion checks:

- `14` nonnegative parallel prime-order tables at `p=3,5,13`, including full
  polynomial conjugacy, equality of all pointed determinants, the norm power,
  and the exact `M_root/p` first ghost;
- `9` arbitrary-table norm first-ghost identities;
- `21` pure-shift norms and the exact `p=13` formula `(18)`;
- the zero-holonomy, nonparallel, composite `C_6`, affine-slope, and Reynolds
  support hostiles, and records the cemetery/root-mass boundary;
- transitivity of every THM-2591 holonomy `7a`, `a in F_13^*`.

Reproduce with

```text
python 04-computation/lrc14_prime_holonomy_parallel_correspondence_arrival_thm3048.py
python -O 04-computation/lrc14_prime_holonomy_parallel_correspondence_arrival_thm3048.py
```

Both runs equal the stored six-line transcript byte-for-byte.  This theorem
does not construct `(7)`, identify the chart carriers, preserve a full
terminal packet, or exclude an LRC row.

**QED.**
