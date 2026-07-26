---
id: THM-2382
title: "Saturated septimal seven-bin root-fibre closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In
  THM-2367's only-c3-dominant k=2 alternative, the saturated profile
  (t,b)=(5,2), W=7 is impossible. On every generic 7^(M+1)-orbit on
  which c3 is safe, the five top q masks and the two top blocker masks
  occupy the seven top bins exactly once each. After writing cj=13Cj,
  choose a generic base phase at which all three quotient blockers are
  safe. All thirteen inverse roots then have all three original
  blockers safe, so the septimal partition would require one of the
  five q masks at every root. This demands thirteen q incidences, while
  each 13-unit q mask has at most two, for total capacity at most ten.
  The proof is an almost-everywhere scalar-cover exclusion of this one
  septimal alternative only. It excludes no individual thirteen-adic
  row, lands no target, settles none of the other k=2 alternatives, and
  does not prove LRC(14).
source: codex-2026-07-25-saturated-septimal-root-fibre
depends_on:
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
related:
  - THM-2377-septimal-valuation-collision-and-bockstein-carry-gate
  - THM-2378-hard-septimal-root-stalk-closure
script: 04-computation/lrc14_saturated_septimal_root_fibre_closure_thm2382.py
output: 05-knowledge/results/lrc14_saturated_septimal_root_fibre_closure_thm2382.out
script_sha256: e3a29874372f227aadafda065904b7930949795dae1d448c3d4c36643d1e8f1d
output_sha256: 781ef09462a9ac2b29838196734a62cabbba2ebc52fa38b12637ef5cb2016667
hash_basis: working-tree bytes (LF)
---

# THM-2382 -- saturated septimal seven-bin root-fibre closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2367 leaves four alternatives when both of the first two blockers
are at or below the septimal top layer. This theorem excludes only the
saturated one:

```text
(t,b)=(5,2),                 W=7.                    (1)
```

The mechanism is a composition of two finite fibres with deliberately
different jobs:

```text
7^(M+1)-orbit:
  scalar coverage + saturated top weight
  -> the seven top masks partition the seven bins;

13-root fibre:
  all quotient blockers safe
  -> five 13-unit q masks would have to cover thirteen roots
  -> 13 <= 5*2, contradiction.                       (2)
```

No base point is confused with one of its thirteen inverse roots. That
typing is load-bearing: quotient-blocker safety at `y` controls the
original blockers at the roots `(y+h)/13`, not generally at `y`
itself.

## 1. Scalar-cover and valuation hypotheses

Retain a canonical first-depth-one scalar cover

```text
C_H subset
  union_(i=1)^5 D_(q_i)
  union D_(c_1) union D_(c_2) union D_(c_3)          (3)
```

almost everywhere, where

```text
D_v={x in T:||vx||<1/14},

C_H={x in T:||Hx||>1/7},

E_H={x in T:||Hx||<1/7}.                             (4)
```

Away from endpoints, (3) is equivalently the cover of `T` by

```text
E_H,D_(q_1),...,D_(q_5),D_(c_1),D_(c_2),D_(c_3).    (5)
```

The inherited thirteen-adic typing is

```text
13 does not divide Hq_1...q_5,

13 divides c_1,c_2,c_3.                              (6)
```

Put

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5)).
```

Assume the only-`c_3`-dominant saturated alternative of THM-2367:

```text
M>0,

nu_7(H)<M,

nu_7(q_i)=M                    for i=1,...,5,

nu_7(c_1)=nu_7(c_2)=M,

nu_7(c_3)>M.                                         (7)
```

Thus the seven ordinary danger masks

```text
Q_top={D_(q_1),...,D_(q_5),D_(c_1),D_(c_2)}         (8)
```

are precisely the weight-seven top layer. The guard `E_H` is strictly
below it and `c_3` is strictly above it.

## 2. Reconstruction of the saturated seven-bin partition

Let

```text
N=7^(M+1).                                          (9)
```

Because `N` divides `c_3`, the truth value of `D_(c_3)` is constant on
each additive `N`-orbit

```text
O_z={z+j/N:j in Z/NZ}.                              (10)
```

Consider an orbit on which `c_3` is safe and which avoids all mask
endpoints and the exceptional null set in the almost-everywhere cover.
Such generic orbits fill the `c_3`-safe set up to a null set.

Partition its indices into the seven top bins

```text
B_r={j in Z/NZ:j=r mod 7},             r in Z/7Z.   (11)
```

Each speed `v` in (8) has `nu_7(v)=M`. Its ordinary danger mask is
therefore constant on each `B_r` and occupies exactly one of the seven
bins. Let

```text
m_r
 =#{v in {q_1,...,q_5,c_1,c_2}:D_v occupies B_r}.   (12)
```

There are seven top masks, each occupying one bin, so

```text
sum_(r in Z/7Z)m_r=7.                               (13)
```

The lower guard has an exact, binwise count. If

```text
nu_7(H)=h<M,
```

then on a fixed bin the induced `H`-grid has order

```text
7^(M-h),
```

with multiplicity `7^h`. Its danger arc has length `2/7`; away from
the aligned endpoints it therefore contributes

```text
7^h * 2*7^(M-h)/7
 =2*7^(M-1)
 =2N/49                                             (14)
```

points to every `B_r`.

On the present orbit `c_3` contributes nothing. Counting incidences in
one bin, the scalar cover (5) consequently gives

```text
m_r N/7+2N/49 >= N/7,
```

or

```text
7m_r+2>=7.                                          (15)
```

Since `m_r` is an integer, (15) gives `m_r>=1` for every `r`.
Together with (13), this forces

```text
m_r=1                         for every r in Z/7Z.  (16)
```

Hence, for almost every physical point `x` at which `c_3` is safe,
exactly one of

```text
D_(q_1),...,D_(q_5),D_(c_1),D_(c_2)                (17)
```

contains `x`.

This is only an almost-everywhere statement on the `c_3`-safe set. At
an aligned strict endpoint a top mask can occupy no bin, and nothing
here asserts a partition where `c_3` is dangerous. Those exclusions
are harmless for the positive-measure choice in the next section.

## 3. A generic quotient-blocker-safe base exists

Write

```text
c_j=13C_j,                         j=1,2,3.          (18)
```

Every ordinary danger comb has Haar measure `1/7`, including its
closure. Therefore

```text
mu(
  T minus
  union_(j=1)^3 closure(D_(C_j))
 )
 >=1-3/7
 =4/7.                                               (19)
```

Let `G` be the full-measure set of physical points at which the
partition conclusion of Section 2 is valid whenever `c_3` is safe.
Choose representatives `y in [0,1)` and, for `h in {0,...,12}`, use
the `h`-th inverse root of multiplication by thirteen,

```text
R_h(y)=(y+h)/13.                                    (20)
```

The set of `y` for which some `R_h(y)` lies outside `G` is null: it is
a finite union of inverse images of a null set under these affine
branches on `[0,1)`.
Consequently (19) contains a base phase `y` such that

```text
y notin closure(D_(C_1))
          union closure(D_(C_2))
          union closure(D_(C_3)),                   (21)

R_h(y) in G                       for every h.       (22)
```

This is the only generic choice in the proof.

## 4. The thirteen-root capacity contradiction

Put

```text
x_h=R_h(y)=(y+h)/13.                                (23)
```

For every blocker and every root,

```text
c_j x_h
 =13C_j (y+h)/13
 =C_j y+C_j h
 ==C_j y                         mod 1.             (24)
```

Thus (21) implies that all three original blocker masks are strictly
safe at every `x_h`. In particular `c_3` is safe, so (22) and (17)
apply. Since `c_1` and `c_2` are safe as well, the unique top danger at
each root must be one of the five `q_i`. Therefore

```text
sum_(i=1)^5 #{h:x_h in D_(q_i)}=13.                 (25)
```

Fix one `q=q_i`. By (6), multiplication by `q` permutes the thirteen
root labels. The points

```text
q x_h,                         h=0,...,12,           (26)
```

form a translate of the equally spaced thirteen-grid. An open arc of
length `1/7` contains at most two points of this grid: three distinct
grid points would span at least `2/13`, while

```text
2/13>1/7.                                           (27)
```

Hence

```text
#{h:x_h in D_(q_i)}<=2                              (28)
```

for every `i`. Summing (28) contradicts (25):

```text
13<=5*2=10.                                         (29)
```

This proves that no scalar cover satisfying (3), (6), and (7) exists.
Equivalently, THM-2367's saturated

```text
k=2,                  (t,b)=(5,2),                  (30)
```

alternative is empty.

## 5. Why each hypothesis is visible

The proof is sharp at each of its two fibre interfaces.

1. **Genericity is real.** On an aligned seven-grid, the strict
   radius-`1/14` arc can have zero rather than one top point. Likewise
   an aligned guard grid can lose an endpoint. The argument removes
   these null phases; it does not silently change strict masks to
   closed ones.

2. **Saturated weight is real.** Equation (16) uses both `m_r>=1` and
   `sum m_r=7`. In the unsaturated `W=2,k=2` alternatives, THM-2367's
   count instead has `sum m_r=W=2`, and its inequality is only
   `7m_r>=0`. Thus a collision vector such as
   `(2,0,0,0,0,0,0)` survives the bin count. The partition conclusion
   cannot be exported to those branches.

3. **All quotient blockers must be safe.** If one `D_(C_j)` is
   dangerous at the base, (24) makes its original blocker dangerous
   on all thirteen roots, and the q-capacity contradiction disappears.
   The three-comb measure bound (19) is what supplies a simultaneous
   safe base.

4. **Thirteen-unit q labels are load-bearing.** If `13|q`, its mask can
   be constant and dangerous on all thirteen roots. The bound (28)
   would then be false.

5. **Five q labels are well below the threshold.** The same argument
   contradicts even six unit q masks, since `6*2<13`; with seven masks,
   the raw capacity becomes `14` and no contradiction follows.

6. **No thirteen-adic row is singled out.** The proof uses only the
   common first-depth-one facts (6). It empties a septimal alternative
   uniformly across every row in which that alternative could occur;
   it does not decrement the `165`-row ledger.

## 6. Relation to the collision gate

THM-2377 says that nonzero high-target modes need a repeated lower
septimal valuation and a Bockstein carry. The saturated branch has the
largest possible repeated layer: all seven ordinary top masks share
valuation `M`. Section 2 supplies the missing rooted sidecar for that
coarse collision datum—their seven bin positions form a transversal,
not an arbitrary same-layer packet. There are `7!=5,040` possible
labelled transversals; their `21` collision pairs split intrinsically as
ten `q-q`, ten `q-c`, and one `c-c` pair.

The leading `F_7` cancellation data do not choose a Bockstein carry. If
`U,V` are units modulo `49` and

```text
U+b_0 V==0 mod 7,
```

then the seven lifts `b=b_0+7t` preserve that leading cancellation,
while

```text
(U+bV)/7 mod 7
```

runs through every residue because `V` is a seven-unit. Thus the
collision gate alone cannot exclude the saturated atlas. The
thirteen-root fibre kills its rooted transversals before any
target-mode landing is required.

Thus THM-2377 is explanatory rather than a logical dependency here.
Its Bockstein carrier remains relevant to the three unsaturated `k=2`
alternatives, where the bin transversal fails.

## 7. Exact companion

The dependency-free exact companion:

- checks generic lower-guard counts on septimal grids and exhibits the
  aligned-endpoint failure;
- checks that a generic top-layer ordinary tooth occupies one
  seven-bin and exhibits the strict-endpoint zero-bin hostile;
- exhausts all `1,716` weak compositions of seven into seven bins and
  finds `(1,1,1,1,1,1,1)` as the unique vector satisfying (13)--(15);
- enumerates all `7^7` labelled top-bin assignments, retaining exactly
  the `7!=5,040` transversals, and records their `10+10+1` pair types;
- verifies that the `W=2,k=2` collision vector survives its weaker
  inequality;
- checks the blocker pullback identity (24) and a nonmultiple hostile;
- computes exact unions for every ordered quotient-blocker triple
  `1<=C_j<=12`, verifying the uniform `4/7` safe-measure floor;
- exhausts the exact phase cells of the thirteen-root word, obtaining
  ordinary counts `{1,2}` and guard counts `{3,4}`;
- exhibits a `13|q` mask with all thirteen roots dangerous; and
- checks the decisive capacities `10<13`, `12<13`, and the sharp raw
  failure `14>=13`;
- verifies on all `1,764` ordered unit pairs modulo `49` that fixed
  leading cancellation leaves every Bockstein carry possible; and
- records a physical root-only hostile with

  ```text
  y=1/97,

  H=1,

  (q_1,...,q_5)=(112,126,175,238,301),

  (c_1,c_2,c_3)=(91,2366,107653).
  ```

  It has exactly the valuation/divisibility profile (7), all blockers
  are safe on the thirteen roots, the guard word is `{0,1,12}`, and
  the five disjoint two-point q words cover the other ten roots. Thus
  guard plus q support covers the root fibre, while the rooted
  septimal top-owner partition fails at precisely three roots. This
  shows why the Section 2 partition, rather than root marginals alone,
  is the missing coordinate.

Run

```bash
python3 04-computation/lrc14_saturated_septimal_root_fibre_closure_thm2382.py
python3 -O 04-computation/lrc14_saturated_septimal_root_fibre_closure_thm2382.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_saturated_septimal_root_fibre_closure_thm2382.out
```

byte for byte.

## 8. Independent audit and conclusion

An independent agent reconstructed the proof from the hypotheses,
checked the quotient/base-versus-root typing, reran the exact companion
under both normal and optimized Python, and verified the stored
transcript and hashes. The audit also found and repaired the hostile's
initially insufficient thirteen-adic blocker profile; the final packet
has

```text
nu_13(c_1,c_2,c_3)=(1,2,3),

nu_7(c_1,c_2,c_3)=(1,1,2).
```

All forty companion checks use an explicit optimization-safe
`require` function. No Python `assert` is used.

The saturated alternative's emptiness is therefore proved. **QED.**
