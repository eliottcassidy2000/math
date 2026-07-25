---
id: THM-2190
title: "Basis-safe floor and height-500 rank-six harvest"
status: >
  PROVED + VERIFIED-EXACT. Let L be a saturated rank-r lattice of integer
  relations on a thirteen-coordinate row with no zero coordinate. For
  0<=r<=5, the Haar mass of the 1/14-safe cube in K_L is at least
  (6-r)6^(12-r)/7^(13-r). In rank six the same basis-cover argument, with
  a strict facet-intersection audit, still proves positive instancewise Haar
  mass, but this argument supplies no coefficient-independent floor.
  Combining the uniform floors with THM-2185's degree-500 Jackson ledger
  proves that every zero-safe distinct positive thirteen-speed row has
  dim_Q W_500 at least six.
  Optimizing the same exact Jackson ledger rank by rank sharpens this to
  dim_Q W_178>=3, dim_Q W_204>=4, dim_Q W_262>=5, and dim_Q W_450>=6.
  Every twelve-speed deletion then has five independent relations of height
  at most 235800. Finiteness of the possible bounded six-planes and the
  instancewise rank-six positivity also imply a universal finite, presently
  uncomputed height H_* with dim_Q W_(H_*)>=7. The theorem does not compute
  H_*, give a rank-seven safe-cube floor, produce rank eight, or prove
  LRC(14).
source: codex-2026-07-24-basis-safe-floor
depends_on:
  - THM-2185-rank-two-safe-cube-floor-and-height-500-rank-three-harvest
  - THM-2164-relative-packet-rank-harvesting
related:
  - THM-2054-relative-fejer-whole-product-decorrelation
  - THM-2164-relative-packet-rank-harvesting
  - THM-2178-mod14-transverse-code-rank-harvest
script: 04-computation/lrc14_basis_safe_floor_rank_six_harvest_thm2190.py
output: 05-knowledge/results/lrc14_basis_safe_floor_rank_six_harvest_thm2190.out
script_sha256: c5f0e317abbb152a2723315a5dcde06658e6c798df2da9f45796291b3c96e8dc
output_sha256: 9f7e38267266e9964385cc64670afa66c0f811a93da257f7a2166b17f7248ae9
hash_basis: working-tree bytes (LF)
---

# THM-2190 -- basis-safe floor and height-500 rank-six harvest

Put

```text
I={x in R/Z:||x||<1/14},
J=(R/Z)\I=[1/14,13/14],
p=measure(I)=1/7,
q=measure(J)=6/7.                                    (1)
```

For `v in Z^13`, write

```text
Lambda(v)={m in Z^13:m.v=0},
W_H(v)=span_Q(Lambda(v) intersection [-H,H]^13).      (2)
```

If `L subset Z^13` is saturated, define

```text
K_L={x in (R/Z)^13:l.x=0 mod 1 for every l in L}.     (3)
```

The main change of viewpoint is to choose independent **coordinate
characters**, rather than classify every short relation among danger events.
After a finite torus cover, those characters become separate one-coordinate
clocks. Only the `r` characters outside that basis need a union bound.

## 1. The finite-cover coordinate basis

Assume

```text
rank(L)=r,
L subset Lambda(v),
v_i!=0 for every i,
d=13-r.                                               (4)
```

Saturation makes

```text
X=Z^13/L                                              (5)
```

a free abelian group of rank `d`, canonically the character group of the
connected torus `K_L`. Let

```text
c_i=e_i+L in X                                       (6)
```

be the thirteen coordinate characters. Every `c_i` is nonzero: if `c_i=0`,
then `e_i in L subset Lambda(v)`, contrary to `v_i!=0`.

The `c_i` generate `X`, so choose and relabel `d` of them,

```text
c_1,...,c_d,                                          (7)
```

which form a basis of `X tensor Q`. Their integer span `X_0` has finite
index

```text
D=[X:X_0].                                            (8)
```

Multiplication by `D` kills `X/X_0`. Hence there is an injective
finite-cokernel map

```text
iota:X -> Z^d,
D x=sum_j iota(x)_j c_j.                              (9)
```

In particular,

```text
iota(c_j)=D e_j                 (1<=j<=d).             (10)
```

Pontryagin duality reverses (9) into a surjective finite-cover homomorphism

```text
rho:(R/Z)^d -> K_L.                                   (11)
```

The pushforward of Haar probability under a continuous surjective compact
group homomorphism is Haar probability. Therefore every Haar-mass statement
on `K_L` can be computed on the cover. Equations (9)--(10) say that

```text
c_j(rho(y))=D y_j                    (1<=j<=d),        (12)
```

while each remaining coordinate pulls back to

```text
c_i(rho(y))=a_i.y
```

for a nonzero vector `a_i=iota(c_i) in Z^d`. Injectivity of `iota` and
`c_i!=0` are exactly what exclude a constant extra character.

## 2. The basis-safe floor

Let

```text
E_D={y in R/Z:D y in J},
B=E_D^d.                                              (13)
```

Multiplication by the nonzero integer `D` preserves Haar pushforward, so

```text
measure(E_D)=q,
measure(B)=q^d.                                       (14)
```

For an extra character `a_i`, put

```text
A_i={y:a_i.y in I}.                                   (15)
```

Choose an index `j` with `(a_i)_j!=0`. For every fixed choice of the other
coordinates, multiplication by `(a_i)_j` pushes Haar on the `j`-th circle
to Haar. Thus

```text
measure(A_i intersection B)
 <=p q^(d-1).                                         (16)
```

Indeed, integrate first in `y_j`: the danger fibre has total mass `p`, and
intersecting it with `E_D` can only decrease that mass. The other `d-1`
basis-safe coordinates contribute `q^(d-1)`.

The preimage of `K_L intersection J^13` is

```text
B \ union_(i=d+1)^13 A_i.                             (17)
```

There are `r` extra characters. The union bound and (14)--(16) give

```text
measure_(K_L)(K_L intersection J^13)
 >=q^d-r p q^(d-1)
 =q^(d-1)(q-rp)
 =(6-r)6^(12-r)/7^(13-r).                            (18)
```

For the ranks used below, the exact floors are

```text
r=2: 241864704/1977326743,
r=3:   30233088/282475249,
r=4:    3359232/40353607,
r=5:      279936/5764801.                             (19)
```

They decrease across this table and remain positive through `r=5`.
This proof neither assumes positivity nor distinctness of the coordinates;
it uses only (4).

## 3. What really happens at rank six

At `r=6`, `d=7` and the right side of (18) is zero. The apparent equality
case is nevertheless impossible for any fixed lattice. The reason is that
the fibre estimate (16) is always strict.

> **Strict facet lemma.** Let `a in Z^d\{0}` and choose `j` with `a_j!=0`.
> Then
>
> ```text
> {y:
>    D y_j in I,
>    D y_k in interior(J) for k!=j,
>    a.y in I}
> ```
>
> contains a nonempty open set.                              (20)

To prove this, put

```text
g=gcd(D,a_1,...,a_d),       h=D/g,       b=a/g.        (21)
```

Write

```text
y_k=(m_k+z_k)/D,
m_k in Z/DZ.                                          (22)
```

Choose the centre `z_j=0` and `z_k=1/2` for `k!=j`.
Because

```text
gcd(h,b_1,...,b_d)=1,
```

the residues `b.m mod h` run through all of `Z/hZ`. At the chosen centre,
the numerator of `a.y` is congruent modulo `h` to

```text
b.m+(1/2) sum_(k!=j)b_k.                              (23)
```

If the displayed sum is even, choose `m` so that (23) is zero. If it is
odd, choose `m` so that its circular distance from zero is `1/2`. The
available numerator variation inside the base danger/safe intervals has
half-width

```text
|b_j|/14+(3/7)sum_(k!=j)|b_k|.                        (24)
```

In the odd case, the second sum contains at least one nonzero term and
`b_j!=0`, so (24) is at least

```text
1/14+3/7=1/2.                                        (25)
```

The target danger interval contributes the additional positive numerator
radius `h/14`. Thus in either parity case the centres can be perturbed
strictly inside all intervals so that `a.y in I`. This gives the open set
in (20).

Now integrate `A={a.y in I}` first in `y_j`. Up to null endpoints, the
deficit in (16) is exactly

```text
p q^(d-1)-measure(A intersection B)
 =measure(A intersection {D y_j in I}
             intersection_(k!=j){D y_k in J}).        (26)
```

The strict facet lemma makes (26) positive. Consequently, when `r=6`,

```text
measure_(K_L)(K_L intersection J^13)
 >=q^7-sum_(i=8)^13 measure(A_i intersection B)
 >q^7-6p q^6
 =0.                                                  (27)
```

This is instancewise strict positivity, not a uniform numerical floor.
The proof gives no coefficient-independent lower bound for the six positive
deficits in (26). The next open analytic interface is therefore

```text
inf_L measure_(K_L)(K_L intersection J^13)>0 ?        (28)
```

over saturated rank-six relation lattices satisfying (4), or an effective
lower bound in terms of a bounded lattice presentation. Such a bound would
give a direct explicit seventh-relation height. Section 6 instead obtains an
existential universal height from finite bounded-seed compactness.

## 4. One-height rank-six harvest

Assume now that `v` is a distinct positive thirteen-speed row with

```text
measure(S(v))=0,
S(v)={t:v_i t in J for every i}.                      (29)
```

THM-2185 proves

```text
dim_Q W_500(v)>=3.                                    (30)
```

Suppose for contradiction that

```text
r=dim_Q W_500(v) in {3,4,5},                          (31)
```

and define the full saturated bounded-relation lattice

```text
L=W_500(v) intersection Z^13.                         (32)
```

This choice is load-bearing. Taking merely some selected relations would
not identify every bounded Fourier resonance. The lattice in (32) is
saturated, has rank `r`, and is contained in `Lambda(v)`.

Use THM-2185's normalized squared-Fejer approximant

```text
q_251=J_251*1_J,
degree(q_251)<=500,
||q_251-1_J||_1<=eta_251,
eta_251<21/12500.                                     (33)
```

A Fourier tuple in `[-500,500]^13` survives the line average exactly when it
belongs to `Lambda(v)`. Every such tuple belongs to `W_500(v)` by definition,
and hence to `L`. Conversely, the annihilator of `K_L` is `L`. Therefore the
finite Fourier expansions give

```text
integral_(R/Z) product_i q_251(v_i t) dt
 =integral_(K_L) product_i q_251(x_i) dx.              (34)
```

Every coordinate character is nonzero on both averaging groups. Product
telescoping, (29), (33), and (34) imply

```text
measure_(K_L)(K_L intersection J^13)
 <=26 eta_251
 <546/12500.                                          (35)
```

But (18)--(19) show that the smallest floor for `r in {3,4,5}` is the
rank-five value, and

```text
279936/5764801-546/12500
 =175809327/36030006250
 >0.                                                  (36)
```

This contradicts (35). Together with (30), it proves:

> **Height-500 rank-six harvest.** Every zero-safe distinct positive
> thirteen-speed row satisfies
>
> ```text
> dim_Q W_500(v)>=6.                                  (37)
> ```

The strict rank-six statement (27) does not automatically continue this
argument to rank seven: (34)--(35) require a lower bound larger than the
fixed smoothing error, not merely an instancewise positive mass.

## 5. The degree-optimized nested harvest

The height `500` proof has the advantage of using one inherited approximant.
The same squared-Fejer/Jackson family gives a stronger nested result when its
degree is optimized separately at each relation rank.

For `N>=2`, let `bar_eta_N` be THM-2185's upper bound obtained by replacing
`pi` with `355/113` in the exact `L1` error formula for the normalized
squared-Fejer approximant, whose degree is

```text
H_N=2N-2.                                             (38)
```

The companion exhausts the finite rational sums and finds the first `N` for
which

```text
26 bar_eta_N < alpha_r,
alpha_r=(6-r)6^(12-r)/7^(13-r):                       (39)
```

```text
current rank r       2       3       4       5
first N              90      103     132     226
degree H_N           178     204     262     450.      (40)
```

At `N=89,102,131,225`, respectively, the same exact rational-`pi` ledger
does not clear (39). Thus (40) records the first passing values for this
specific kernel and cap, not an optimality theorem over all approximants.
Convenient strict rational wrappers at the passing values are

```text
bar_eta_90  <47/10000,
bar_eta_103 <41/10000,
bar_eta_132 <2/625,
bar_eta_226 <373/200000.                              (41)
```

THM-2164 starts the induction with

```text
dim_Q W_105(v)>=2.                                    (42)
```

Suppose first that `dim_Q W_178(v)=2` and put
`L=W_178(v) intersection Z^13`. Exactly as in Section 4, every degree-`178`
line resonance belongs to this full saturated lattice, so the line and
`K_L` averages of the `N=90` tensor agree. The rank-two floor and the first
column of (40) contradict zero safe measure. Hence

```text
dim_Q W_178(v)>=3.                                    (43)
```

Repeat this argument at the increasing heights in (40). At each step, if
the desired next rank has not yet appeared, the full lattice

```text
L=W_(H_N)(v) intersection Z^13                        (44)
```

has exactly the current rank `r`; its annihilator therefore has the floor
`alpha_r`, and (39) gives the same smoothing contradiction. This proves

```text
dim_Q W_178(v)>=3,
dim_Q W_204(v)>=4,
dim_Q W_262(v)>=5,
dim_Q W_450(v)>=6.                                    (45)
```

In particular, (37) is an immediate weaker corollary of (45). The sequential
use of the full lattice (44), rather than the saturation of an arbitrarily
selected subfamily, preserves exact agreement of all bounded Fourier
resonances.

## 6. A non-effective uniform seventh relation

The strict rank-six positivity in (27) has a further consequence once it is
joined to the bounded six-tuple in (45).

Let `F` be the finite, nonempty collection of rational six-planes spanned
by independent integer tuples bounded respectively by the height profile

```text
(105,105,178,204,262,450)                             (46)
```

and containing no coordinate vector `e_i`. For `S in F`, put

```text
L_S=S intersection Z^13,
delta(S)=measure_(K_(L_S))(K_(L_S) intersection J^13). (47)
```

Every `L_S` is saturated of rank six and has nonzero coordinate characters.
The hypothesis `L subset Lambda(v)` in Sections 1--3 was used only to force
this last nonvanishing; equivalently, because `e_i notin S`, a generic vector
of `S^perp` has every coordinate nonzero. Section 3 therefore proves
`delta(S)>0`. Since `F` is finite,

```text
delta_* = min_(S in F) delta(S)>0.                    (48)
```

Let `epsilon_N` denote the actual `L1` approximation error of the normalized
squared-Fejer approximant. The approximate-identity theorem gives
`epsilon_N -> 0`. Choose `N_*>=226` so that

```text
26 epsilon_(N_*)<delta_*,
H_*=2N_*-2.                                           (49)
```

This `H_*` is universal, although no useful numerical value is extracted
here.

Suppose a zero-safe row had `dim_Q W_(H_*)(v)=6`. By (45), it has a nested
six-tuple with profile (46); call its span `S`. Since
`W_450(v) subset W_(H_*)(v)` and both `S` and `W_(H_*)(v)` have dimension
six, they agree. Thus

```text
W_(H_*)(v) intersection Z^13=L_S.                    (50)
```

The degree-`H_*` line and `K_(L_S)` Fourier averages consequently have the
same complete bounded resonance set. Product telescoping would give

```text
delta_*<=delta(S)<=26 epsilon_(N_*),
```

contrary to (49). Therefore:

> **Finite-height rank-seven theorem.** There is a universal finite integer
> `H_*` such that every zero-safe distinct positive thirteen-speed row
> satisfies
>
> ```text
> dim_Q W_(H_*)(v)>=7.                                (51)
> ```

The proof is finite and constructive in principle: enumerate `F`, compute
the rational-polyhedral Haar masses in (47), and then choose `N_*` using
either the exact error formula with successively sharper certified rational
bounds for `pi`, or any explicit approximate-identity rate. It is
non-effective in the present theorem only in the practical sense that this
finite enumeration and a numerical `H_*` have not been carried out.

## 7. Five explicitly bounded relations on every deletion

The nested spaces in (42)--(45) admit independent relations

```text
r_1,...,r_6
```

with respective height bounds

```text
(h_1,...,h_6)=(105,105,178,204,262,450).              (52)
```

Fix a coordinate `i`, and let `s` be the first index for which
`(r_s)_i!=0`. If there is no such index, any five relations already live on
the deletion. If `s=6`, the first five do. Otherwise retain
`r_1,...,r_(s-1)` and put

```text
u_j=(r_j)_i r_s-(r_s)_i r_j,          s<j<=6.         (53)
```

These are five relations in total, all vanishing in coordinate `i`. A
triangular coefficient comparison in the basis `r_1,...,r_6` proves their
independence. Their heights obey

```text
||u_j||_infinity<=2 h_s h_j<=2 h_s h_6.              (54)
```

For `s=1,2,3,4,5`, the worst bounds are

```text
94500, 94500, 160200, 183600, 235800,                 (55)
```

while at `s=6` the retained first five have height at most `262`. Hence:

> **Deletion rank-five corollary.** Every twelve-speed deletion of a
> zero-safe distinct positive thirteen-speed row has relation rank at least
> five inside coefficient height `235800`.                    (56)

Using only (37) and six arbitrary height-`500` relations would give the
simpler bound `2*500^2=500000`; the nested flag pays less.

The non-effective seventh relation in (51) similarly gives six independent
relations on every deletion at a universal finite height. More precisely,
choose

```text
r_7 in Lambda(v)\span_Q{r_1,...,r_6},
||r_7||_infinity<=H_*.
```

Apply the same first-nonzero-pivot elimination to the ordered seven-tuple.
If the pivot is among `r_1,...,r_6`, every resulting relation has height at
most

```text
2*450*H_*=900H_*.
```

If only `r_7` is nonzero in the deleted coordinate, retain the first six.
Thus every deletion has relation rank at least six by height `900H_*`.
This is an absolute bound, but it remains nonnumeric because `H_*` does.

## 8. Exact referee and boundaries

The companion uses only exact integer and `Fraction` arithmetic. It:

1. reconstructs (18)--(19) and verifies monotonicity through rank five;
2. checks the exact one-height Jackson comparison gap (36);
3. independently reconstructs the squared-Fejer coefficients at `N=251`
   and verifies THM-2185's rational-`pi` error cap;
4. exhausts the exact rational-`pi` ledgers to verify the first-passing
   thresholds and adjacent failures in (40);
5. checks the rational constant identity in the parity-radius estimate used
   by the strict facet lemma; and
6. verifies every explicit height in (52)--(55).

Normal and optimized executions are required to agree with the frozen
transcript.

The theorem does not compute `H_*`, prove a rank-seven safe-cube floor,
produce an eighth independent relation, give a finite speed cap, or prove
LRC(14). QED.
