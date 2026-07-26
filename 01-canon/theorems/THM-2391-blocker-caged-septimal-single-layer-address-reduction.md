---
id: THM-2391
title: "Blocker-caged septimal single-layer address reduction"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. In the
  last k=2,(t,b)=(1,0) septimal lane, THM-2388's collision cage forces,
  on every generic C_3-safe top-q fibre, the two divided low-blocker
  words to partition the guard word exactly and disjointly, while all
  four lower q words lie inside that guard word. A sharp cyclic-period
  lemma then forces H, the four lower q labels, C_1, and C_2 to have one
  common septimal valuation a<M. Primitivity forces a=0, so THM-2390's
  weight-seven alternative disappears and only its weight-eight
  one-double word remains. At M>=2 every contained ordinary word has
  guard-gauge step +/-1 or +/-2; the blocker partition is either two
  contiguous halves or an even/odd interleaving. At M=1 only a binary
  address survives. Pullback by thirteen turns the quotient-blocker
  partition into D_c1 plus D_c2 equals E_(13H) on the 13q top comb
  outside c_3. This is a structural reduction, not a branch exclusion,
  row decrement, target landing, or proof of LRC(14).
source: codex-2026-07-26-blocker-caged-septimal-address
depends_on:
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2388-thirteen-root-multiplicity-reflection-and-blocker-caged-toothpick-law
  - THM-2390-septimal-layer-kraft-peeling-and-heavy-word-reduction
related:
  - THM-2377-septimal-valuation-collision-and-bockstein-carry-gate
  - THM-2385-two-top-septimal-blocker-collision-reduction
  - THM-2390-septimal-layer-kraft-peeling-and-heavy-word-reduction
script: 04-computation/lrc14_blocker_caged_single_layer_address_thm2391.py
output: 05-knowledge/results/lrc14_blocker_caged_single_layer_address_thm2391.out
script_sha256: 8b39faf0114f995c1aba45baa80d901f146d8c90783d5e4a710f27fc64925512
output_sha256: fb27e0f6d543c37584357d2900dbfb353c6eab88053fcc64befdbd987570a08c
hash_basis: working-tree bytes (LF)
---

# THM-2391 -- blocker-caged septimal single-layer address reduction

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2390 proves that the only remaining septimal lane has a lower layer
of weight seven or eight.  THM-2388 supplies a different object: every
collision of two unit masks must land in one of three **divided**
blockers.  Combining them does more than count weight.  It puts all eight
lower units into one layer and leaves a labelled finite address cage:

```text
two divided blockers partition the two-address guard word;

each of four lower q words chooses addresses inside that guard word;

the original blockers are the thirteen-drifted versions of the
two partition words.                                      (1)
```

The word is binary when the top and lower layers are adjacent.  Across a
larger gap it has only four possible slopes.  This is the precise
Fano/`chi_7` remnant; no tournament on the runner labels alone retains
the needed slopes, translates, and thirteen-drift.

## 1. The last branch and the collision cage

Retain THM-2367's primitive scalar cover

```text
C_H subset
  union_(i=1)^5 D_(q_i)
  union D_(c_1) union D_(c_2) union D_(c_3),        (2)

D_v={x:||vx||<1/14},

E_H={x:||Hx||<1/7}.                                 (3)
```

The six labels `H,q_1,...,q_5` are units modulo thirteen,
`c_j=13C_j`, and the five `q_i` are pairwise distinct.  In the sole
remaining septimal alternative,

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5))>0,

nu_7(q_*)=M,

nu_7(H),nu_7(q_i),nu_7(c_1),nu_7(c_2)<M
                                      for q_i!=q_*,

nu_7(c_3)>M.                                        (4)
```

Division by thirteen preserves every septimal valuation.  Define the
six-unit multiplicity

```text
K=1_(E_H)+sum_(i=1)^5 1_(D_(q_i))                  (5)
```

and the quotient-blocker cage

```text
B=D_(C_1) union D_(C_2) union D_(C_3).              (6)
```

THM-2388 proves

```text
{K>=2} subset B                                    (7)
```

almost everywhere.  This direction is essential: it is a statement at
the physical collision point about the **divided** blockers.

## 2. Exact partition on a top bin

Put

```text
N=7^(M+1)
```

and take a generic additive orbit

```text
O_z={z+j/N:j in Z/NZ}.                              (8)
```

Because `nu_7(C_3)>M`, the `C_3` word is constant on this orbit.
The set of `C_3`-safe orbit bases has measure `6/7`; deleting the scalar
cover's null exceptional set, its `N` translates, and the finitely many
strict endpoints still leaves full measure inside that safe set.

On such an orbit the top word `D_(q_*)` occupies one complete residue bin

```text
Q={z+(r+7s)/N:s in Z/(N/7)Z},                      (9)
```

of `N/7` points.  Every ordinary subtop word has exactly `N/49` points
in `Q`, and the subtop guard has exactly `2N/49`.

If `x in Q intersection E_H`, then the top `q_*` and the guard both
contribute to (5), so `K(x)>=2`.  Equations (6)--(7), together with
`C_3`-safety, give

```text
Q intersection E_H
 subset
 Q intersection (D_(C_1) union D_(C_2)).           (10)
```

The left side has `2N/49` points.  Each set on the right has `N/49`
points.  Thus there is no capacity slack:

```text
Q intersection E_H
 =
 (Q intersection D_(C_1))
 disjoint_union
 (Q intersection D_(C_2)).                         (11)
```

Likewise, for every lower `q_i`, a point of
`Q intersection D_(q_i)` has `K>=2`.  Hence

```text
Q intersection D_(q_i)
 subset Q intersection E_H,             q_i!=q_*.  (12)
```

Equations (11)--(12) hold on every generic `C_3`-safe orbit, not merely
after averaging over those orbits.

## 3. A cyclic period lemma

Let `L=7^M`.  Reindex the bin (9) by `s in Z/LZ`.  If

```text
v=7^a u,                     0<=a<M,       7 does not divide u,
```

then the slice of `D_v` in `Q` is the lift from

```text
p_v=7^(M-a)                                      (13)
```

of a unit-gauge cyclic interval of `p_v/7` residues.  The guard slice is
the corresponding lift of an interval of `2p_H/7` residues.  Arbitrary
base phase changes only the cyclic translates.

> **Period-containment lemma.**  Let `A` be a lifted width-one word of
> least period `p`, and let `G` be a lifted width-two word of least
> period `q`, where `p,q` are positive powers of seven in a common
> cyclic group.  If
>
> ```text
> A subset G,                                           (14)
> ```
>
> then `p=q`.

**Proof.**  A nonempty proper cyclic interval has no nonzero translation
stabilizer, so the displayed periods are indeed least periods.

If `p<q`, put the guard in its unit gauge.  In every `q`-period its
complement contains a consecutive run of `5q/7` residues.  Since
`p<=q/7`, that run contains `p` consecutive residues and hence a full
residue system modulo `p`.  The nonempty `p`-periodic word `A` meets it,
contradicting (14).

If `p>q`, put `A` in its unit gauge.  Its interval of `p/7` consecutive
residues contains `q` consecutive residues because `p>=7q`.  The
nonfull `q`-periodic guard misses at least one of them, again
contradicting (14).  Therefore `p=q`. QED.

Apply the lemma to each inclusion in (11)--(12).  It gives one common
lower valuation

```text
nu_7(H)
 =nu_7(q_i)                         for q_i!=q_*
 =nu_7(C_1)
 =nu_7(C_2)
 =:a<M.                                             (15)
```

Since `c_j=13C_j`,

```text
nu_7(c_1)=nu_7(c_2)=a.                              (16)
```

If `a>0`, then every coefficient in the nine-factor scalar row is
divisible by seven: the seven labels in (15)--(16) have depth `a`, while
`q_*` and `c_3` have still larger depth.  Primitivity therefore forces

```text
a=0.                                                (17)
```

This strengthens THM-2390: its entire lower weight eight lies in the
primitive layer zero.  The weight-seven heavy word is impossible.
Consequently the terminal seven-root word supplied there has
multiplicity one at six addresses and multiplicity two at exactly one.

## 4. The four slopes across a nonadjacent gap

The remaining finite word sharpens further when

```text
M>=2.
```

Put

```text
m=7^M,                         k=m/7.
```

Multiply the bin address by the unit `H` and translate it so the guard
word is

```text
I={0,1,...,2k-1} subset Z/mZ.                       (18)
```

An ordinary lower word of speed `v` becomes a `k`-term arithmetic
progression with unit step

```text
delta_v=H v^(-1) mod m.                             (19)
```

If that progression lies in `I`, take the balanced representative of
its step.  Consecutive progression points both lie in an interval of
diameter `2k-1`; because

```text
2(2k-1)<m=7k,
```

their integer difference is the unique lift of `delta_v` in that
range.  The progression therefore has no hidden wrap and

```text
(k-1)|delta_v|<=2k-1.
```

Since `k>=7`,

```text
delta_v in {+1,-1,+2,-2}.                           (20)
```

If two such progressions partition `I`, there are only two set-types.
A step of absolute value two fills one parity class, forcing the other
word to fill the other parity class.  A step of absolute value one is a
contiguous `k`-block; its complement is another allowed progression
only when the two blocks are the two end halves.  Thus (11) is either

```text
contiguous type:
  {0,...,k-1} disjoint_union {k,...,2k-1},

parity type:
  {0,2,...,2k-2} disjoint_union {1,3,...,2k-1}.     (21)
```

In ratio notation every lower ordinary label satisfies

```text
v/H in {+1,-1,+1/2,-1/2} mod 7^M.                  (22)
```

For the four pairwise-distinct lower `q` labels, either two repeat one
class in (22), in which case

```text
7^M divides q_i-q_j,                                (23)
```

or all four classes occur.  Modulo seven the latter case is the balanced
cross

```text
{1,6,4,3},                                          (24)
```

with two quadratic residues and two nonresidues.  This is the exact
`chi_7`/Fano probe left by the branch.

When `M=1`, one has `k=1`.  Directions disappear: the two blocker
singletons are simply the two guard addresses, and each of the four
lower `q` singletons chooses one of them.  This adjacent-layer binary
cage is a genuine separate boundary, not a failed version of (20).

## 5. The thirteen-drift identity

Equation (11) has a useful physical pullback.  In circle notation it is
the almost-everywhere identity

```text
1_(E_H)
 =1_(D_(C_1))+1_(D_(C_2))

on D_(q_*) intersection D_(C_3)^c.                 (25)
```

Substitute `y=13x`.  Since `c_j=13C_j`,

```text
1_(D_(c_1))+1_(D_(c_2))
 =1_(E_(13H))

on D_(13q_*) intersection D_(c_3)^c.               (26)
```

Thus the two original low blockers are not arbitrary occupants of the
weight-eight word: on the thirteen-drifted top comb they are exactly a
width-two guard at speed `13H`.  The congruence

```text
13=-1+2*7                                           (27)
```

shows why the adjacent binary address reverses while the next septimal
digit carries a nontrivial translation.  Any recursive or tournament
encoding must retain the blocker label, the chosen guard address, and
this carry; quotienting to a two-colouring loses the preserved
predicate.

## 6. Sharpness and residual

The cyclic conclusions are sharp as finite words.

- At `M=1`, any two distinct singleton addresses partition a two-point
  guard word, and four labelled lower singletons can choose those two
  addresses with repetition.
- At every `M>=2`, both the contiguous and parity partitions in (21)
  exist, and all four slopes in (20) occur.
- Individual containment has thin additive hostile families such as
  `v=H+q_*`, because

  ```text
  ||Hx||=||(v-q_*)x||
       <=||vx||+||q_*x||<1/7
  ```

  on `D_v intersection D_(q_*)`.  Therefore no proof may replace the
  two-blocker partition by an unsupported one-mask anti-shield.

What remains is global complementarity across base phase: two centred
low combs must keep partitioning the moving guard word on the
`C_3`-safe top comb, while their thirteen multiples realize (26) and the
actual lower word has one double.  Exact finite searches find no global
two-comb partition in the first canonical boxes, but that no-go is not
proved here.

No thirteen-adic row is removed, the ledger remains `165`, and LRC(14)
remains open.

## 7. Exact companion

The dependency-free companion:

- exhausts lifted septimal words through modulus `343`, checking their
  cardinalities and least periods;
- verifies the period-containment lemma on every distinct-depth word
  pair through modulus `49`;
- exhausts all contained progressions at moduli `49` and `343`, obtaining
  exactly the four slopes in (20);
- classifies every complementary pair as contiguous or parity type;
- checks the binary `M=1` boundary, the four-class pigeonhole, and the
  thirteen-address drift; and
- retains explicit positive controls for the two partition types and
  for the additive one-mask hostile.

Run

```bash
python3 04-computation/lrc14_blocker_caged_single_layer_address_thm2391.py
python3 -O 04-computation/lrc14_blocker_caged_single_layer_address_thm2391.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_blocker_caged_single_layer_address_thm2391.out
```

after LF normalization.  Every executable check raises explicitly under
optimized Python.
