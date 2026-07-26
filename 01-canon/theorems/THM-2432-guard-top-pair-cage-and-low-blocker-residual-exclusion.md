---
id: THM-2432
title: "Guard-top pair-cage and low-blocker residual exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The two
  positive-depth THM-2427 residuals (k,t,b,W)=(2,0,0,2) and
  (2,5,0,7) are empty. In the first type, the seven lower ordinary
  words partition every c_3-safe bin outside the top guard; a
  thirteen-root descent cages the two low blockers under c_3, and
  THM-2385 excludes both unequal and equal lower depths. In the second
  type, the guard and five top ordinary words form an exact one-fold
  seven-bin partition, so every ordinary pair is c_3-caged. Scaling
  gives two 91-unit sources under one 91-divisible target, contradicting
  a transverse two-source comb-escape lemma. Five valuation types
  remain, no scalar profile row is removed, and LRC(14) remains open.
source: mac-mini-2026-07-26-guard-top-pair-cage
depends_on:
  - THM-2385-two-top-septimal-blocker-collision-reduction
  - THM-2427-guard-top-thirteen-root-capacity-and-residual-types
related:
  - THM-2430-guard-top-common-ninety-one-root-tiling-spectrum
  - THM-2431-repeated-step-rounding-exclusion-of-guard-top-zero-blocker-types
script: 04-computation/lrc14_guard_top_pair_cage_thm2432.py
output: 05-knowledge/results/lrc14_guard_top_pair_cage_thm2432.out
script_sha256: 58bf8c52c8ab5ffb4ae0db252813a61ba099b82d18c9d1ecfe5e238eb8c6c006
output_sha256: 3cd9180a2ac17de96b9444f99ed9a66ae498bc1a60a63b239e93d58176ac095d
hash_basis: working-tree bytes (LF)
---

# THM-2432 -- exact partitions cage two pairs that cannot be caged

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2427 reduces the positive-depth deep-`c_3` branch to four valuation
types. Two apparently different types carry the same hidden object:

```text
an exact one-fold partition away from c_3
  -> a pair of source combs can overlap only under c_3
  -> scale out their common top layer
  -> two source combs cannot fit inside one transverse deeper comb. (1)
```

For `(2,0,0,2)` the exact partition is made by the lower words outside
the top guard. For `(2,5,0,7)` it is made by the top guard and all five
ordinary words. The first branch invokes THM-2385 directly; the second
requires its `7 x 13` transverse analogue.

## 1. Typed setup

Retain THM-2427's primitive almost-everywhere scalar cover

```text
C_H subset
  union_(i=1)^5 D_(q_i)
  union D_(c_1) union D_(c_2) union D_(c_3),             (2)

D_v={x:||vx||<1/14},             C_H=E_H^c,

E_H={x:||Hx||<1/7},              c_j=13C_j.
```

Put

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5))>0
```

and work on the remaining deep-blocker/top-guard side

```text
nu_7(H)=M<nu_7(c_3).                                  (3)
```

The labels `H,q_1,...,q_5` are units modulo thirteen. The tuple

```text
(k,t,b,W)
```

has the meaning fixed in THM-2427. Endpoints and inherited exceptional
sets are always removed before a finite root or orbit disintegration;
their finitely many affine pullbacks remain null.

## 2. The low-blocker type has an exact lower partition

Assume first

```text
(k,t,b,W)=(2,0,0,2).                                    (4)
```

Thus the guard is top, while every `q_i,c_1,c_2` is strictly below
depth `M`, and `c_3` is strictly above it. Put

```text
N=7^(M+1),

L=sum_(i=1)^5 1_(D_(q_i))+1_(D_(c_1))+1_(D_(c_2)).       (5)
```

On every generic `c_3`-safe `N`-orbit, each of these seven lower
ordinary words contributes exactly `N/49` incidences to each of the
seven top residue bins. The top guard occupies two complete bins. In
each of the five other bins the scalar cover forces lower coverage,
while the total lower incidence is exactly

```text
7N/49=N/7,
```

the size of the bin. Coverage is therefore one-fold:

```text
L=1                  almost everywhere on
                     (E_H union D_(c_3))^c.            (6)
```

No exactness is asserted inside the guard or the high blocker.

## 3. Thirteen roots cage the two lower blockers

Suppose there were a base point

```text
y in D_(C_1) intersection D_(C_2)
       minus closure(D_(C_3)).                         (7)
```

Choose it away from the relevant null set and form

```text
x_h=(y+h)/13,                  h in F_13.              (8)
```

All three physical blocker words are constant on this root fibre:

```text
1_(D_(c_j))(x_h)=1_(D_(C_j))(y).                      (9)
```

Hence both lower blocker summands in (5) equal one at every root and
`c_3` is safe. The guard speed is a thirteen-unit, so its arc of length
`2/7` meets at most four of the thirteen roots. At least nine roots lie
in the set from (6), where (6) says `L=1` although the two blocker
summands already give `L>=2`. This contradiction proves

```text
D_(C_1) intersection D_(C_2)
  subset closure(D_(C_3))                 almost everywhere. (10)
```

The difference of the left side and the closed target is open.
Consequently (10) is a literal containment: a nonempty difference
would have positive Haar measure.

## 4. THM-2385 closes both lower-depth possibilities

The depths in (4) give

```text
nu_7(C_1),nu_7(C_2)<M<nu_7(C_3).                     (11)
```

If the two source depths differ, THM-2385's unequal-depth intersection
lemma gives

```text
mu(D_(C_1) intersection D_(C_2)
      minus closure(D_(C_3)))=6/343,                  (12)
```

contradicting (10).

If the depths agree at `r<M`, write

```text
C_1=7^r a,             C_2=7^r b,             C_3=7^r c.
```

Then

```text
7 does not divide ab,                  49 divides c.   (13)
```

Surjectivity of multiplication by `7^r` transports (10) to

```text
D_a intersection D_b subset closure(D_c),
```

contradicting THM-2385's two-layer single-target lemma. Thus

```text
(2,0,0,2) is empty.                                    (14)
```

No strict thirteen-adic blocker profile was used.

## 5. The five-top type has an exact top partition

Now assume

```text
(k,t,b,W)=(2,5,0,7).                                    (15)
```

All five ordinary labels and the guard have exact depth `M`; both
low blockers lie below `M`; and `c_3` lies above it. On a generic
`c_3`-safe `N`-orbit, let `m_r` be the top-word multiplicity in bin
`r`. THM-2427 gives

```text
7m_r>=W-k=5,
```

so every `m_r>=1`. But the guard contributes two top bins and the five
ordinary words one each, whence

```text
sum_(r in F_7)m_r=W=7.
```

Therefore

```text
m_r=1                         for every r in F_7.       (16)
```

The top words form an exact one-fold partition off the high blocker:

```text
1_(E_H)+sum_(i=1)^5 1_(D_(q_i))=1
                         almost everywhere on D_(c_3)^c. (17)
```

In particular, for every pair `i!=j`,

```text
D_(q_i) intersection D_(q_j)
  subset closure(D_(c_3)).                              (18)
```

Again the possible difference is open, so the almost-everywhere
containment is literal.

Scale an arbitrary pair and the target by their common top layer:

```text
q_i=7^M a,               q_j=7^M b,               c_3=7^M c.
```

Exact top depth and the thirteen-unit labels give

```text
gcd(ab,91)=1,
```

while `13|c_3` and `nu_7(c_3)>M` give

```text
91|c.                                                    (19)
```

Surjectivity of multiplication by `7^M` turns (18) into

```text
D_a intersection D_b subset closure(D_c).              (20)
```

The next lemma contradicts (20).

## 6. Transverse two-source comb escape

> **Lemma.** Let `a,b,c` be positive integers satisfying
>
> ```text
> gcd(ab,91)=1,                         91|c.            (21)
> ```
>
> Then
>
> ```text
> D_a intersection D_b minus closure(D_c)
> ```
>
> contains a nonempty open interval.

**Proof.** Suppose the containment held and put `g=gcd(a,b)`. At every
common zero `x=k/g`, both source inequalities hold strictly. Hence the
cyclic subgroup generated by `c/g` would lie in the closed centred arc
of radius `1/14`. Every nontrivial finite cyclic subgroup has an element
of centred norm at least `1/3`. Therefore

```text
g|c.                                                     (22)
```

Because `g` is a `91`-unit, division by `g` preserves (21). The
surjective change of variable `y=gx` reduces to

```text
gcd(a,b)=1.
```

After interchanging the sources, assume `0<a<b`; the equal case would
give `a=b=1` and is covered by the same central argument below.

If `c>b`, choose

```text
1/(14c)<x<min(1/(14b),13/(14c)).                      (23)
```

Both sources are dangerous and the target is strictly safe. Thus
containment forces `c<=b`. Equality is impossible because `7|c` and
`7` does not divide `b`, so

```text
0<c<b,                         b>=92.                  (24)
```

At the centres `x=n/b` of the `b`-teeth, multiplication by `a`
permutes the residues modulo `b`. Let `lambda` be the balanced
representative of

```text
lambda==c a^(-1) mod b,                 |lambda|<=b/2. (25)
```

Containment would imply, for every residue `j`,

```text
||j/b||<1/14       =>       ||lambda j/b||<=1/14.      (26)
```

If `A=|lambda|>=2`, choose

```text
j=floor(b/(14A))+1.
```

Since `b>=92`,

```text
0<j<b/14,                 b/14<Aj<=4b/7.              (27)
```

If `Aj<=b/2`, its centred residue already exceeds `b/14`; otherwise it
is at least `b-4b/7=3b/7`. Both alternatives contradict (26). Hence

```text
lambda in {0,+1,-1}.                                   (28)
```

The zero case contradicts `0<c<b`. The positive case gives `c=a`,
contrary to their different septimal divisibility. Thus

```text
lambda=-1,                         b=a+c.              (29)
```

Write

```text
b=14m+r_0,                         1<=r_0<=13,
```

and choose `n` with `an==m mod b`. Since `c>=91>r_0`, put

```text
delta=(c-r_0)/(14b)
```

and choose

```text
0<eta<
 min(1/(7b),(1/14-delta)/a,delta/c).                   (30)
```

At

```text
x=n/b-1/(14b)+eta
```

the three phases modulo integers are

```text
b x=-1/14+b eta,

a x=delta+a eta,

c x=-(1/14+delta)+c eta.                               (31)
```

The first two have norm strictly below `1/14`; the last has norm
strictly between `1/14` and `1/7`. All inequalities persist on a
neighbourhood. This is the required open escape and the final
contradiction. QED.

Applying the lemma to (19)--(20) proves

```text
(2,5,0,7) is empty at positive depth.                  (32)
```

## 7. Exact residual list

Before this theorem, THM-2427 left

```text
M=0:
  (0,5,0,7), (1,5,1,8), (2,5,2,9);

M>0:
  (1,5,0,7), (2,0,0,2), (2,5,0,7), (2,5,1,8).        (33)
```

Equations (14) and (32) reduce the positive-depth list to

```text
M>0:
  (1,5,0,7), (2,5,1,8).                               (34)
```

Exactly five valuation types remain in total. This is a valuation-type
reduction, not a thirteen-adic row exclusion: each of the `165` scalar
profile rows may still realize another valuation type.

## 8. Sharp boundary and concurrent route

Pure `91`-root geometry cannot replace the pair-cage argument. At

```text
y=12345/65537,

H=1,

(q_1,...,q_5)=(731,1046,318,1775,1047),
```

the guard has `26` roots, each ordinary word has `13`, and the six
words partition the entire common root stalk. The exact minimum
distance from every tested boundary is

```text
19503/11927734>0,                                      (35)
```

so the tiling persists on an open chamber. THM-2430 classifies all such
finite tilings and likewise shows that they are locally physical.
Neither object is a global scalar-cover packet.

Reserved THM-2431 pursues a different repeated-step rounding mechanism
for the broader `t=5,b=0` family. This theorem neither assumes nor
promotes that route; its exclusion of the single positive-depth
`(2,5,0,7)` type is independently supplied by (17)--(31).

## 9. Exact companion

Run

```bash
python3 04-computation/lrc14_guard_top_pair_cage_thm2432.py
python3 -O 04-computation/lrc14_guard_top_pair_cage_thm2432.py
```

The dependency-free companion:

- checks `54` top-guard fibres, `108` lower words, and all `756`
  lower per-bin counts;
- reconstructs all generic thirteen-root phase chambers and the exact
  ordinary/guard size histograms `{1,2}` and `{3,4}`;
- verifies `1,170` blocker pullbacks and `432` unequal-depth
  intersection fibres, including the exact `6/343` wall measure;
- independently replays all `72,600` old `49`-layer escape cases;
- finds the unique admissible multiplicity vector among all `1,716`
  weak compositions of seven;
- checks the subgroup and multiplier inequalities and constructs exact
  open-escape witnesses for all `62,160` transverse `91`-layer cases,
  exercising the subgroup, central, multiplier, and boundary-perturbation
  mechanisms;
- verifies the strict physical hostile and its exact margin (35); and
- reproduces the seven-to-five residual type reduction.

Every truth-bearing check raises explicitly. Normal and optimized
transcripts must byte-match

```text
05-knowledge/results/lrc14_guard_top_pair_cage_thm2432.out
```

after LF normalization.

## 10. Independent hostile audit

Two independent audits reconstructed both partitions, the
almost-everywhere-to-literal openness upgrades, the scaling maps, and
the gcd normalization in the transverse lemma. One audit exhausted
`6,048` exact chamber triples with smaller sources and found no
containment; the other reconstructed all terminal `lambda=-1` cases
through source size `240`. Both confirmed that `91|c`, rather than the
older literal hypothesis `49|c`, is sufficient because gcd division
preserves the full transverse layer and leaves `c>=91`. No endpoint,
owner, valuation-zero, or primitivity exception was found.
