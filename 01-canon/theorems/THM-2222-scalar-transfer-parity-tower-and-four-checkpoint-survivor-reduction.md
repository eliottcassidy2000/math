---
id: THM-2222
title: "Scalar transfer-parity tower and four-checkpoint survivor reduction"
status: >
  PROVED + VERIFIED-EXACT. In the scalar five-unit/three-blocker branch,
  the signed unit residual is a -1 eigenfunction of the unnormalized
  thirteen-fold transfer operator. A cover therefore forces its positive
  set into every even divided-blocker union and its negative set into every
  odd one. In particular, lambda_1>=6 forces a union of three danger combs,
  with a primitive coefficient triple bounded by
  B_6=3,877,322,523,365,316, to survive four consecutive 169-checkpoints
  with measure at least 961/6930. Thus one explicit finite extremal
  inequality S_4(B_6)<961/6930 would force lambda_1<=5 and remove 455
  valuation profiles. That inequality remains OPEN. It is enough to bound
  an explicit adaptive upper LP using only singleton, pair, and triple
  overlaps of the four checkpoint unions. A geometric-chain control has
  three-checkpoint mass 916159/4826809 above the target, while both its
  four-checkpoint mass and its cubic-moment upper certificate lie below
  the target. Thus four is the first plausible checkpoint depth; no
  universal extremality statement or proof of LRC(14) is claimed.
source: codex-gauss-2026-07-24-scalar-transfer-parity-tower
depends_on:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
related:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2210-nested-binomial-minorant-and-adaptive-moment-lp-hierarchy
  - THM-2215-scalar-depth-234-affine-needle-capacity-exclusion
  - THM-2216-residual-capacity-hinge-gram-law
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
script: 04-computation/lrc14_scalar_parity_tower_four_checkpoint_thm2222.py
output: 05-knowledge/results/lrc14_scalar_parity_tower_four_checkpoint_thm2222.out
script_sha256: 1ee54ce28580d198c5f4919ec6c50e001f3300a59c2e804bb256ad46d3afa900
output_sha256: 88621cf4bd723d2d6338231fa808d120c6869d1d497eba49d86e2f0d2846acfd
hash_basis: working-tree bytes (LF)
---

# THM-2222 -- the scalar transfer-parity tower

This theorem does not empty the deeper scalar ledger. It replaces 455 of
its profiles by one exact finite four-checkpoint extremal problem and
identifies the first checkpoint count that is not defeated by the
canonical geometric chain.

## 1. Setup and the signed residual

On `T=R/Z`, put

```text
D_a={x:||ax||<1/14},
C_H={x:||Hx||>1/7}.                                  (1)
```

Work in the scalar `5+3` branch of THM-2198. Thus `H,q_1,...,q_5` are
positive thirteen-units, the three actual blocker coefficients `c_j` have

```text
1<=lambda_1<=lambda_2<lambda_3,
lambda_j=nu_13(c_j),                                 (2)
```

and a putative counterexample satisfies

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j)             (3)
```

almost everywhere.

Define

```text
R=1_(C_H)-sum_(i=1)^5 1_(D_(q_i)),                  (4)

A_+={R>0}
   =C_H\union_(i=1)^5 D_(q_i).                       (5)
```

THM-2198 proves the uniform mass floor

```text
measure(A_+)>=delta_5=961/6930.                      (6)
```

Let

```text
T(x)=13x mod 1,
(Lf)(y)=sum_(x:T(x)=y)f(x)                           (7)
```

be the multiplication map and its unnormalized transfer operator.

## 2. Exact transfer eigenline

For `13` not dividing `a`, direct root counting gives, away from finitely
many endpoints,

```text
L 1_(D_a)=2-1_(D_a).                                 (8)
```

Indeed, write the thirteen preimages as `(y+k)/13`. Multiplication by `a`
permutes the root digits, and the danger count is the number of integers
in an interval of radius `13/14` centered at `ay`: it is one when
`||ay||<1/14` and two otherwise.

The analogous guard count is

```text
L 1_(C_H)=10-1_(C_H).                                (9)
```

The complement count is the number of integers in an interval of radius
`13/7`: it is three when `||Hy||<=1/7` and four otherwise, so the safe
count is ten minus the guard bit. Endpoint choices are null.

Equations (4), (8), and (9) telescope:

```text
LR
 =(10-1_(C_H))-sum_i(2-1_(D_(q_i)))
 =-R.                                                (10)
```

Consequently

```text
L^kR=(-1)^kR                         for every k>=0. (11)
```

The useful object is therefore not the scalar measure of the unit
residual but its one-dimensional sign representation under transfer.

## 3. Divided blockers and parity inclusions

For every integer `0<=k<=lambda_1`, define

```text
B_k=sum_(j=1)^3 1_(D_(c_j/13^k)),
U_k={B_k>0}=union_(j=1)^3D_(c_j/13^k).              (12)
```

Outside the null set in (3),

```text
R<=B_k composed T^k.                                 (13)
```

If `R=1`, the point lies in the guard and avoids all five unit dangers, so
the cover supplies an actual blocker. Divisibility gives

```text
1_(D_(c_j))(x)
 =1_(D_(c_j/13^k))(T^k x).
```

If `R<=0`, inequality (13) is automatic because `B_k>=0`. This argument
also includes `k=0`; retaining that original-cover checkpoint is
load-bearing.

Apply the positive operator `L^k` to (13). Since every one of the
`13^k` preimages has the same `T^k` image,

```text
(-1)^kR<=13^kB_k.                                    (14)
```

Images of the discarded null set are null. Hence

```text
A_+ subset U_k                         for even k,   (15)
{R<0} subset U_k                       for odd k.    (16)
```

This is the transfer-parity tower.

## 4. One three-comb shift register

Put `s=floor(lambda_1/2)` and set

```text
U=U_(2s).
```

For `0<=r<=s`,

```text
U_(2r)=T^(-2(s-r))U.                                 (17)
```

Using every even inclusion in (15), including `k=0`, gives

```text
A_+
 subset intersection_(j=0)^s T^(-2j)U.              (18)
```

Thus the same three-comb set must accept the positive residual word at
`s+1` consecutive checkpoints of the map `T^2(x)=169x`.

In particular, if `lambda_1>=6`, take

```text
d_j=c_j/13^6.
```

Then (18) contains the four-checkpoint consequence

```text
measure K_4(d_1,d_2,d_3)>=delta_5,                  (19)

K_4(d)
 =intersection_(r=0)^3
    union_(j=1)^3 D_(169^r d_j).                    (20)
```

No capacity averaging or root-family quotient is used in this reduction.

## 5. The exact finite target

In THM-2203's fixed scalar section, write `c_j=13s_j`. The corresponding
original speed coordinate is

```text
208s_j=16c_j.
```

THM-2203's largest-speed tooth gives

```text
208s_j<=w<=(13/14)91^12.                             (21)
```

If `lambda_1>=6`, then `d_j=s_j/13^5` is a positive integer. Equations
(21) give

```text
d_j<=B_6
 :=floor(91^12/(224*13^5))
 =3,877,322,523,365,316.                             (22)
```

The three blockers are distinct. Sort the `d_j` and divide by their common
gcd. This does not change (20)'s measure: multiplication by the common
factor pulls the set back under a Haar-preserving circle endomorphism.

Define the finite extremal quantity

```text
S_4(B)
 =max measure K_4(d_1,d_2,d_3),                     (23)

where
 1<=d_1<d_2<d_3<=B,
 gcd(d_1,d_2,d_3)=1.
```

Every row is exactly rational. One may construct `D_a` from the endpoints

```text
(14m-1)/(14a), (14m+1)/(14a)
```

on the circle, merge the three unions in (20), intersect the four interval
lists, and sum their rational lengths.

Equations (6), (19), and (22) prove the finite implication

```text
S_4(B_6)<961/6930
          implies lambda_1<=5 in every scalar survivor.        (24)
```

The strict inequality in (24) is **OPEN**.

Among the `1,136` profiles left after the four depth-at-most-three
exclusions, exactly

```text
sum_(c=7)^19 binom(c-5,2)=binom(15,3)=455           (25)
```

have `lambda_1>=6`. Thus (24) would leave `681` of that ledger. THM-2213
and THM-2215 already exclude `(3,3,4)` and `(2,3,4)`, respectively, both
among those `681`. Relative to current proved canon, (24) would therefore
reduce `1,134` live profiles to `679`. Later independent profile
exclusions should be subtracted separately.

## 6. Why four checkpoints are the first plausible target

Let `A=169` and take the geometric triple

```text
d=(1,A,A^2).                                         (26)
```

Put

```text
X_j(x)=1_(D_1)(A^j x).
```

Then membership in the shifted three-comb union at checkpoint `r` says

```text
X_r or X_(r+1) or X_(r+2)=1.                         (27)
```

The three-checkpoint intersection contains `D_(A^2)` outright, so its
measure is at least `1/7>delta_5`. Exact transfer gives more. Squaring
(8) yields

```text
L^2 1_(D_1)=24+1_(D_1).                             (28)
```

Therefore `(X_j)` is the stationary reversible Markov chain with

```text
stationary law (6/7,1/7),

transition matrix
 [[145,24],
  [144,25]]/169,                                    (29)
```

where states are `0,1`. The event (27) at consecutive checkpoints is the
binary language with no word `000`. Exact finite-state multiplication
gives

```text
measure K_3(1,A,A^2)
 =916159/4826809
 =delta_5+1710418421/33449786370,                   (30)

measure K_4(1,A,A^2)
 =3385513/33787663
 =delta_5-1286905579/33449786370.                   (31)
```

Thus no universal three-checkpoint estimate can reach the target. The
canonical chain crosses below it at exactly four checkpoints. Nothing here
proves that this chain maximizes (23).

## 7. A cubic-moment upper certificate

There is no need to construct the fourfold intersection in (20) in order
to certify (24). For a fixed triple `d`, put

```text
Y_r=1_(union_j D_(169^r d_j)),  0<=r<=3,
K=Y_0+Y_1+Y_2+Y_3,
p_k=measure{K=k},
M_r=E binom(K,r).                                    (32)
```

Then `measure K_4(d)=p_4`. If `t=p_4`, binomial inversion through degree
three gives

```text
p_3=M_3-4t,
p_2=M_2-3M_3+6t,
p_1=M_1-2M_2+3M_3-4t,
p_0=1-M_1+M_2-M_3+t.                                (33)
```

As `t` increases, only `p_1` and `p_3` decrease. Hence the exact largest
value of `p_4` compatible with `M_0,...,M_3` is

```text
U_3(d)
 =1/4 min(M_3, M_1-2M_2+3M_3).                     (34)
```

Equivalently, the two sharp pointwise majorants of `1_(K=4)` are

```text
binom(K,3)/4,
[K-2binom(K,2)+3binom(K,3)]/4.                      (35)
```

The second has values `(0,1/4,0,0,1)` on `K=0,...,4`; choosing the
smaller majorant after seeing the moment packet is the adaptive step.
Each `M_r` with `r<=3` is a sum of measures of `r`-wise intersections
among the four checkpoint unions. Therefore the sufficient finite target

```text
max_d U_3(d)<961/6930,                              (36)
```

over the same primitive bounded triples as (23), needs only singleton,
pair, and triple overlaps and implies (24). It is stronger than the
desired `S_4` inequality, not equivalent to it.

For the geometric triple (26), the exact checkpoint-count histogram has

```text
p_1=127310580000/965009442943,
p_4=3385513/33787663.                               (37)
```

The first majorant alone fails:

```text
M_3/4=864315097/5710115047
     =0.151365... > delta_5.                        (38)
```

But the adaptive second majorant is binding and succeeds:

```text
U_3
 =p_4+p_1/4
 =128521281793/965009442943
 =delta_5-5245941691819/955359348513570.            (39)
```

Thus degree three contains enough information on the canonical hostile
chain, but only after retaining the nonmonotone binomial-basis
combination. This is the upper-LP dual of THM-2210's adaptive lower
hierarchy and a decisive positive/hostile control for any finite search.

## 8. Representations and the remaining obstruction

The binary shift-register view is literal, not metaphorical: (20) records
four observations of one three-comb state under a fixed endomorphism.
Its transition word is analogous to THM-2195's cyclic run ledger, while
the eigenline (10) is the useful stopping coordinate suggested by the
biased-coin problem. Neither analogy supplies the open maximum in (23).

THM-2218 gives the correct labelled Fourier/group-algebra carrier for
individual lift capacities. It also explains why a coarser recursion is
dangerous: family sums retain only the zero mode, Fourier magnitudes lose
relative phase, and reduced Hasse jets lose integer order. At the marked
root-sheet level, translation cost is a genuine cyclic Hamming derivative;
summing phases before taking absolute values permits cancellations just as
forgetting tournament block markers permits cheaper transport.

The present theorem deliberately avoids that loss by tracking set
membership through time. Its remaining obstacle is now precise:

```text
prove the finite four-checkpoint extremal inequality (24),
or find a triple witnessing its failure and retain the extra coordinate
that distinguishes that witness from the geometric chain.              (40)
```

A flat enumeration through `B_6` is finite but enormous. The likely next
tools are a ratio/valuation quotient, exact interval uncrossing, or a
Fourier cross-phase bound compatible with all four dilates.

## 9. Exact referee

Run

```bash
python3 04-computation/lrc14_scalar_parity_tower_four_checkpoint_thm2222.py
python3 -O 04-computation/lrc14_scalar_parity_tower_four_checkpoint_thm2222.py
```

The companion checks (6), (22), (25), the stationary law (29), and every
fraction and gap in (30)--(31) using exact rational arithmetic. Ordinary
and optimized outputs are byte-identical. It explicitly labels (24) open
and does not enumerate `S_4(B_6)`.

The transfer identities and cover implication prove the all-depth parity
tower; the computation is only an arithmetic and finite-state referee.
This theorem is a reduction, not a proof of LRC(14). QED.
