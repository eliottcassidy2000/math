# Four-comb nonuniformity is false; the four-residue chart gives a stronger cone

This note attacks the replacement proposed in THM-1141 and records the
research ledger promoted to THM-1148.

## 1. Exact refutation of the `4/3` actual-mean lemma

Write

```text
D_k={t in R/Z : ||kt||<1/14}.
```

The core

```text
P={1,2,3,4,5,6,7,8}
```

has the maximal safe component

```text
J=[1/14,13/112].
```

Indeed its midpoint is safe for all eight core speeds, the left endpoint is
the `p=1` wall, and the right endpoint is the `p=8` wall.  The four distinct
killers

```text
(k1,k2,k3,k4)=(108,109,110,111)
```

are legal because `108>13 max(P)=104`.  Exact endpoint subtraction gives five
survivor components in `J`, with lengths

```text
319/55944, 305/55944, 97/18648, 277/55944, 13/3024.
```

Consequently

```text
survivor mass       A = 955/37296,
component count     C = 5,
actual mean       A/C = 191/37296,
longest gap          L = 319/55944,
L/(A/C)                = 638/573 = 1.113438... < 4/3.       (1)
```

Thus the THM-1141 proposal is false if “mean gap” has its ordinary meaning,
even on a **maximal** core-safe component, with distinct integer moduli, at a
legal LRC scale.  Arbitrary truncation makes the obstruction stronger: inside
any region where the endpoint word has no wall event, shrink around a regular
tooth cell and the finitely many gap lengths become arbitrarily close, so no
universal multiplicative factor `c>1` can follow from distinctness alone.

This counterexample does not refute the sharp four-comb conclusion.  In fact

```text
7 k4 L = 319/72 = 4.4305... > 1.                         (2)
```

The point is structural: low gap variance can coexist with a very large mean.
Here the nearby teeth coalesce, reducing the survivor component count.  The
missing carrier in THM-1141 is therefore not “irregularity of gaps” alone; it
is the pair

```text
(tooth-overlap cluster count, metric gap word).          (3)
```

The phrase “`4/3` times the mean” is also ambiguous in THM-1141.  If it means
`4/3` times the crude union-bound baseline `3/(7 sum k_i)`, it is the absolute
inequality `L>=4/(7 sum k_i)`, not a nonuniformity theorem about the actual
gaps.  Equation (1) refutes the latter and leaves the former as essentially
the four-comb problem itself.

## 2. The corrected local dichotomy

There is a simple exact statement behind (3).  Consider a cyclic tooth chart
of circumference `1/k1` containing one full tooth from each of four combs.
If their union has `q` connected clusters, then it leaves `q` cyclic safe
gaps of total length at least

```text
1/k1 - (1/7) sum_i 1/k_i.
```

For `q<=3`, the longest gap is therefore at least

```text
[1/k1-(1/7)sum_i 1/k_i]/q > 1/(7k4).                   (4)
```

For the strict inequality at `q=3`, multiply by `7k4` and use

```text
7k4/k1-k4 sum_i(1/k_i)
 =6k4/k1-k4/k2-k4/k3-1
 >=4k4/k1-1 >3.
```

So nonuniformity is needed, if at all, only in a chart with **four disjoint
tooth clusters**.  Applying a blanket `4/3` demand to connected or partially
connected configurations discards exactly the overlap surplus which makes
(2) easy.  Boundary fragments and charts containing two teeth from one
modulus must still be ledgered before (4) becomes a global theorem; the
statement here isolates that boundary obligation rather than hiding it in an
unqualified interval claim.

## 3. A new four-residue multiplier lemma

The genuinely useful replacement is already latent in THM-1134's cardinality
ledger:

> For every subset `A` of `Z/13Z` with at most four elements, some nonzero
> multiplier `u` makes the largest cyclic gap of `uA` at least seven grid
> units.

For four-element sets, the `715` labelled subsets split into seven affine
orbits.  Representatives, orbit sizes, witnesses, and witness gap words are

```text
representative  size  u   gap word
0 1 2 3           78  1   1 1 1 10
0 1 2 4          156  1   1 1 2 9
0 1 2 5          156  1   1 1 3 8
0 1 2 6          156  2   2 2 8 1
0 1 3 4           78  1   1 2 1 9
0 1 3 9           52  2   2 3 1 7
0 1 3 11          39  1   1 2 8 2.
```

The sizes sum to `715`; the sixth row proves sharpness.  Sets of smaller
cardinality can either be extended to four points and then pruned, or checked
directly; their exact minimum best gaps for cardinalities `1,2,3,4` are

```text
13,12,9,7.                                                (5)
```

## 4. Multiplier Kakeya cone

Let `P subset {1,...,12}` and let

```text
k1<k2<k3<k4,
B in [k1,k4] integer,
a_i=k_i-B,
A=max_i |a_i|,
M=max(A,84).
```

Apply (5) to `{a_i mod 13}` and choose `t0=u/13` so that the offset centres
have a seven-unit cyclic gap.  Put

```text
epsilon=14/(365M),             I=[t0-epsilon,t0+epsilon].
```

Every core speed remains safe on `I`.  At the worst value `M=84`,

```text
epsilon=1/2190<1/2184,
1/13-12 epsilon=339/4745>1/14.                         (6)
```

The initial vertical safe width in the seven-unit centre gap is

```text
7/13-1/7=36/91.
```

Both walls together drift by at most

```text
2A epsilon<=28/365<1/13,                               (7)
```

so grid clusters do not cross.  Their common safe gap has width at least

```text
36/91-28/365=10592/33215>11/73.                        (8)
```

Fix an arc `X` of length `11/73` in that gap.  If

```text
B>=15M,                                                  (9)
```

then

```text
B|I|>=15*(28/365)=84/73=1+|X|.                         (10)
```

The slope-`B` needle therefore has a complete preimage of `X` inside `I`.
It is safe for the core and all four killers, and its length is

```text
11/(73B)>1/(7k4),                                      (11)
```

because `B<=k4` and `77>73`.

Thus (9) is an arbitrary-core, arbitrary-shape four-comb cone.  With the
upper midpoint

```text
Delta=k4-k1,
A=ceil(Delta/2),
B=k1+A,
```

a transparent corollary is

```text
k1>=max(1260,7(Delta+1)).                              (12)
```

This strictly strengthens THM-1128's
`max(1272,26(Delta+1))`: outside the new cone one must have approximately
`k4/k1>8/7`, rather than merely `k4/k1>27/26`.  It therefore removes most of
the clustered region for which THM-1141 sought a global nonuniformity law.

The constants `15` and `84` are Pareto-minimal for this fixed phase-blind
scalar ledger: a universal chart, complete-preimage winding, and no
grid-cluster crossing.  This is not a global optimality statement about other
Kakeya constructions.  Any safe needle longer than `1/(7k4)` needs a
vertical arc longer than `1/7`, so a complete preimage costs more than `8/7`
windings.  Suppose first that `B<=14M`.  When `M=A`, order preservation gives

```text
B|I|=(B/A)(2A epsilon)<14/13<8/7.
```

When the floor is active with `M<=84`, core safety gives
`2 epsilon<=1/1092`, and hence

```text
B|I|<=14M/1092<=14/13<8/7.
```

Thus ratio `14` cannot work with any floor at most `84` in this ledger.  At
ratio `15`, the complete-preimage winding budget requires
`15(2M epsilon)>8/7`, or `2M epsilon>8/105`.  Core safety
requires `epsilon<=1/2184`, hence

```text
M>1092*(8/105)=83.2.
```

The least integer floor is therefore `84`, exactly as used above.

The best integer chart can also be written without a midpoint loss.  If
`Delta=k4-k1`, maximizing `B/max(A,84)` over integer `B in [k1,k4]` gives

```text
Delta<=84:        (B,M)=(k4,84),
85<=Delta<=168:   (B,M)=(k1+84,84),
Delta>=169:       (B,M)=(k1+ceil(Delta/2),ceil(Delta/2)). (13)
```

Thus the exact cone gate is respectively `k4>=1260`, `k1>=1176`, or
`k1>=14 ceil(Delta/2)`.  The script exhaustively compares (13) with the
brute optimizer through all transition values and both midpoint parities.

## 5. The four-comb separated-ratio gate

The mass/component argument which THM-1134 wrote for five combs has a useful
four-comb form.  Let `J` be a nonwrapping core-safe interval of length `ell`,
let `k1<k2<k3<k4=K`, and put

```text
Q4=ell(21K-7 sum_i k_i)-6K sum_i(1/k_i)-39.            (14)
```

Then `Q4>0` forces a survivor component longer than `1/(7K)`.

Indeed, THM-1097's sharp one-comb discrepancy gives survivor mass

```text
A0>=3ell/7-(6/49)sum_i(1/k_i).                         (15)
```

A tooth of `D_ki` meeting `J` has its centre in an interval of scaled length
`ki ell+1/7`, hence there are at most `ki ell+8/7` such teeth.  Four combs
therefore contribute at most `ell sum_i ki+32/7` teeth.  The complement of
their union in connected `J` has at most one more component, so

```text
C0<=ell sum_i ki+39/7.                                 (16)
```

The constants are not mnemonic approximations: direct expansion gives

```text
7K*(right side of (15))-(right side of (16))=Q4/7.     (17)
```

Thus `Q4>0` means `7K A0>C0`, proving the claim.  Isolating the self terms
gives the equivalent and often clearer form

```text
Q4=ell(14K-7(k1+k2+k3))
      -6K(1/k1+1/k2+1/k3)-45.                         (18)
```

Along an integer ratio ray `(k1,k2,k3,K)=m(a,b,c,d)`, (18) becomes

```text
7m ell(2d-a-b-c)-6d(1/a+1/b+1/c)-45.                  (19)
```

Consequently every ray with `2d>a+b+c` enters the Q4 tail at an explicit
finite scale.  Rays with `2d<=a+b+c` can never be certified by this gate;
at equality its value is a fixed negative rational.

## 6. Exact residual after the four current gates

For a fixed eight-speed core `P`, let `ell(P)` be its longest safe-component
length.  The four proved gates certify a legal quadruple whenever at least one
of the following holds:

```text
(i)   Delta<=30                                      [THM-1133],
(ii)  the exact piecewise chart condition (13) holds [new cone],
(iii) Q4(ell(P);k1,k2,k3,k4)>0                       [equation (14)],
(iv)  the exact three-step Phi transfer below succeeds [THM-1137].
```

The live mainline correction THM-1137 is essential here: the former THM-1140
claim that a one-period window always leaves a `6/(7k)` component is false.
The sharp transfer is

```text
Phi(x)=min(6/7,(x-1/7)/2),             x>=1.           (20)
```

An independent exact atlas over all `495` cores gives

```text
min_P ell(P)(13 max(P)+1)=72/35>13/7,                 (21)
```

attained by `P={1,2,6,7,8,9,10,11}` with `ell=1/70`.  Hence every legal
`k1` puts the first transfer in its saturated regime.  Set

```text
c1=6/7,
x_i=(k_i/k_(i-1))c_(i-1),
c_i=Phi(x_i),                         i=2,3,4.         (22)
```

Gate (iv) succeeds precisely when every `x_i>=1`; then `c4>=3/7>1/7` and
the final survivor component has length greater than `1/(7k4)`.

As a clean corollary, three adjacent ratios at least `9/5` suffice.  At the
boundary `r=9/5`, the exact ledger is

```text
inputs:   54/35,       63/50,        3519/3500,
outputs:   7/10,      391/700,       3019/7000>1/7.    (23)
```

The optimal common-ratio threshold within (22) is the positive root
`r_*=1.797111878...` of

```text
6r^3-r^2-2r-28=0;                                     (24)
```

`9/5` is a transparent rational version.  This strengthens THM-1140's
corrected coarse `7/3` cone and replaces the false historical `7/6` claim.

These four gates do **not** close all shapes.  The first primitive
proof-method residual in the bounded ray census is

```text
(k1,k2,k3,k4)=m(3,4,5,6),             m>=53.          (25)
```

It is legal for every core, its span is `3m>30`, and the cone fails.  Its Q4
linear coefficient vanishes exactly and

```text
Q4=-6*6*(1/3+1/4+1/5)-45=-366/5.                      (26)
```

The exact transfer starts at `c1=6/7` and gives

```text
x2=(4/3)(6/7)=8/7,       c2=1/2,
x3=(5/4)(1/2)=5/8<1,                                  (27)
```

so (22) cannot continue.  Thus the whole ray (25) is an infinite proof-method
residual, not a counterexample and not a finite-bank omission.
It sits exactly on the Q4 slope wall and one unit of ratio arithmetic beyond
the transfer wall, making it a sharper next target than an unstructured
“clustered majority.”

For normalized ratios

```text
(x,y,z)=(k1,k2,k3)/k4,
```

every asymptotic residual obeys the necessary conditions

```text
x<=7/8,              x+y+z>=2,
and the exact three-step Phi path (22) fails.           (28)
```

The boundary `x=7/8` retains an integer-midpoint parity sidecar: even offset
span is certified at equality, odd span is not.  A census of all `3,646,069`
primitive rays with `1<=a<b<c<d<=100` finds `3,054,412` in the strict cone,
strict Q4 tail, or exact-transfer union and `591,657` residual (`16.2273%`).  This
is a diagnostic census of ratio space, while (i)--(iv) are the exact finite
predicate.

## 7. Tournament/carrier audit

Possible vertices were challenged as runners, combs, individual teeth,
survivor gaps, section boundaries, wall events, residues, and proof
obligations.  Ordering exact endpoints gives a transitive tournament with
score multiset `{0,1,...,N-1}`, singleton SCCs, no directed cycles, and one
Hamiltonian path after ties are coalesced.  That tournament forgets interval
length, endpoint owner, overlap cluster, and wall slope, so it does not
preserve either (1) or (11).

The useful finite tournament is instead the cyclic residue word after the
multiplier switch `a -> ua`: its Hamiltonian cycle carries integer gap labels.
The naked orientation again loses the LRC predicate; the faithful quotient is

```text
(cyclic residue word, integer gap labels, wall-slope budget, needle slope).
```

That carrier explains both outcomes: tooth-cluster count kills the false
nonuniformity claim, while a labelled seven-unit residue gap proves the new
cone.

## Verification

Run

```bash
python3 04-computation/lrc14_r5_nonuniformity_refutation_multiplier_cone_codex_20260718.py
```

The script uses `Fraction` throughout, reconstructs the exact counterexample,
enumerates and disjointly covers all seven affine orbits, checks sharpness on
all `715` four-sets, guards the piecewise chart formula and named cone
identities, proves the core-atlas minimum, and performs the primitive-ray
census.  The short algebra deriving Q4 from the displayed mass/component
bounds remains the paper calculation in THM-1148.
