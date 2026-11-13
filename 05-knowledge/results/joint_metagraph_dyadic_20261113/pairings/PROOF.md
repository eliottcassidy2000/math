# The dyadic clock of bounded joint Hamilton fibres

**Status: PROVED IN NOTES, P1–P4; independently audited.** The arguments
are elementary derivations from the stated inherited mechanisms. The
executables give **FINITE-EXACT** controls. See the root [audit](../AUDIT.md)
and byte-level [promotion record](../PROMOTIONS.json). Chan–Pak completeness
remains **OPEN**.

This packet extracts a modular consequence of the existing fixed-H rational
series. It does not claim the SCC recurrence, those rational series, the
minimum pairing-defect lemma, or the abstract free binary-action rank as new.
No external-priority claim is made for polynomial sequences modulo powers of
two. The new local deliverable is their exact application to **joint ordinary
and self-converse fibres**, with a computable least parity clock and an
all-order low-H optimal-pairing defect.

## Inheritance and the live board

The inherited source is revision
`f90cd879d156a37f0027d1719a25b8731918db6a`.

* Closest proved mechanism: [Session24, §§1–2](../../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md)
  gives finite strong cores, rational fixed-H generating functions, exact
  polynomial degree and converse-fixed quasipolynomials. The joint star and
  square-central recurrences are recovered in the
  [November2 spectrum packet](../../tournament_joint_parity_20261102/spectrum.md).
* Pairing input: [November3 P3](../../paired_recursions_20261103/PAIRINGS.md)
  proves that an H-preserving involution commuting with converse has minimum
  fixed-point count equal to the number of odd self-converse H fibres.
* Canonical hostiles: the unique H=1 class; H=19, which is realized but never
  self-converse; and H=9, with **two** strong self-converse classes at order5.
  [THM-4922, first creators](../../../../01-canon/theorems/THM-4922-first-creators-of-hamiltonian-path-counts-converse-fibers-endpoint-genealogy.md)
  and MISTAKE-735 explicitly prohibit the one-class-or-one-converse-pair guess.
* Corrected near miss: ordinary multiplication adds H-exponents, whereas
  ordered SCC concatenation multiplies H. The proposed literal toothpick
  lift already fails at physical-port incidence and at root-forgetting:
  `3z+z^3` is rooted insertion mass, while `F_3=z+z^3`.
* Least-used sidecar: the **valuation at x=1 of the core numerator after
  reduction modulo2**. Real leading degree alone loses its cancellations.

The board is: (anchor) joint H multiplicities; (niche) converse-fixed fibres;
(wildcard) dyadic residue clocks; plus strong-core numerator cancellation and
physical insertion incidence. The resulting connection is between the first
four. Toothpick's physical dyadic decomposition is not transported.

## Definitions and inherited finite-core identities

Let `a(n,h)=[z^h]F_n(z)` count ordinary order-n tournament classes, and let
`c(n,h)` count the self-converse ones. Empty order has H=1. Put

```text
A_h(x)=sum_n a(n,h)x^n,       C_h(x)=sum_n c(n,h)x^n.
```

Delete singleton SCCs from an ordinary tournament. The surviving ordered
word is its nontrivial strong core. Let `b(h,r,m)` count these words with
H-product h, length r and total vertex order m. These inventories are finite.
The inherited exact expression is

```text
a(n,h)=sum_(r,m) b(h,r,m) binom(n-m+r,r),              (1)
```

with the summand zero for n<m. Equivalently

```text
A_h(x)=sum_(r,m) b(h,r,m) x^m/(1-x)^(r+1).           (2)
```

Each nontrivial strong component has H>=3. Thus r<=floor(log_3 h).
Camion's classical Hamilton-cycle theorem gives |S|<=H(S) for a nontrivial
strong S, by cutting the cycle at each vertex. Consequently every such core
has m<=h: the sum of factors >=3 is at most their product. For h=1 the core
is empty and (1) gives a(n,1)=1 for every n.

A self-converse SCC word has the unique form `U -> K -> U^op`, where K is
empty, a singleton, or one nontrivial strong self-converse component. If the
nontrivial core of U has product u, length t and order m, and K has product v
and order d, then h=u^2 v. Its singleton-padding contribution is

```text
[n>=2m+d and n=d (mod2)] binom((n-d)/2-m+t,t).        (3)
```

The multiplicity is the number of ordinary half-core words times the number
of actual allowed centers. In particular `2t<=floor(log_3 h)`.
The bare skeleton order `2m+d` is at most h: if u=1 then m=0; otherwise
m<=u, d<=v and `2u+v<=u^2v` for u>=3,v>=1. This includes the empty and
singleton center conventions. Therefore n>=h is a safe common onset for
all the literal binomial expressions in (1) and (3).

The sharper logarithmic onset is already in Session24. The conservative
linear onset here makes the uniform modular statement transparent and
requires no floating-point logarithm or unproved strong-spectrum input.

## P1. A finite joint cutoff has an explicit dyadic modular clock

Fix integers B>=3 and k>=1. Let R be the unique integer with
`3^R<=B<3^(R+1)`, and set

```text
P(B,k)=2^(k-1+ceil(log_2(R+1))).                     (4)
```

For **every n>=B and every h<=B**, simultaneously,

```text
a(n+P(B,k),h)=a(n,h) (mod 2^k),
c(n+P(B,k),h)=c(n,h) (mod 2^k).                     (5)
```

For B=1 or2, only h=1 is realized, so period1 works from order0 for every
precision. The empty cutoff has the same trivial convention. Even h are
identically zero by Rédei, and nonrealized odd h are permitted.

**Proof.** Write `L=k-1+ceil(log_2(R+1))`, so
`L=k+floor(log_2 R)`. For 1<=j<=r<=R, the elementary identity

```text
v_2 binom(2^L,j)=L-v_2(j)>=k                         (6)
```

holds: in `(2^L/j) binom(2^L-1,j-1)` the latter factor is odd, since every
binary digit of `2^L-1` is1. Vandermonde's formula therefore gives
`binom(N+2^L,r)=binom(N,r) mod2^k`. Apply this to each term of (1).

For (3), the even order increment P preserves the center parity and shifts
the binomial's upper argument by P/2. If t=0 there is no change. Otherwise
`t<=floor(R/2)` and

```text
v_2 binom(2^(L-1),j)>=L-1-floor(log_2 floor(R/2))=k
```

for 1<=j<=t. Vandermonde applies again. The onset n>=B puts every skeleton
in range, so no zero-truncation boundary is crossed. Summing proves (5). ∎

This period is an upper bound. The integer multiplicities need not be
periodic: already a(n,3)=n-2 for n>=3. No bound on the full, untruncated
polynomial F_n is obtained by allowing B to vary with n.

## P2. The least eventual parity period is exactly a reduced pole order

For a realized h let `r_h=max{r:b(h,r,m)>0}`, and construct the integral
polynomial

```text
P_h(x)=sum_(r,m) b(h,r,m)x^m(1-x)^(r_h-r),
A_h(x)=P_h(x)/(1-x)^(r_h+1).                         (7)
```

Reduce P_h modulo2. Let e_h be the reduced denominator exponent:

```text
e_h=max(0,r_h+1-ord_(x=1)(P_h mod2)),                (8)
```

with e_h=0 if the numerator is zero. For an unrealized h also set e_h=0.
Then the sequence `a(n,h) mod2` has **least eventual period**

```text
p_h=2^ceil(log_2(max(1,e_h))).                       (9)
```

The same sequence, hence the same least period, is `c(n,h) mod2`.

**Proof.** Converse is an involution preserving H, so its nonfixed classes
occur in pairs and `a(n,h)=c(n,h) mod2` at every order. After cancelling the
power of (1-x), the generating series over F2 has denominator `(1-x)^e_h`.
If e_h=0 it is a polynomial; its tail is zero and has period1. This case
includes nonzero polynomial quotients, not only the zero series.

If e_h>0, multiplying by `1-x^(2^j)=(1-x)^(2^j)` makes the series polynomial
exactly when `2^j>=e_h`; therefore that power of two is an eventual period.
Conversely, any eventual period p makes `(1-x^p)A_h` polynomial. Write
p=2^j q with q odd. The multiplicity of x=1 in `1-x^p` over F2 is exactly
2^j, forcing e_h<=2^j. This proves minimality in (9), among **all** positive
integer periods. Polynomial parts only affect the finite initial segment. ∎

For a finite collection of h, the least eventual period of the entire
parity **vector** is max_h p_h, because all periods are powers of two.
Its sum as an integer-valued statistic can have a smaller period; no
minimality transfer from vector to sum is assumed.

The hostile H9 shows why reduction must come first:

```text
A_9(x)=2x^5/(1-x)^2+x^6/(1-x)^3.
```

The two strong H9 atoms cancel modulo2, leaving a pole of order3 and
least parity period4. An inventory of mere H-support cannot reconstruct
this answer. In contrast the H15 inventory gives an even numerator, so
the whole parity sequence vanishes although the H15 stratum is nonempty.

## P3. Exact pairing defect through H=21 at every order

Let X_(n,B) contain the ordinary order-n classes with H<=B. Among involutions
J on this set that preserve H and commute with converse R, the inherited
minimum fixed-point formula is

```text
Delta_B(n)=sum_(h<=B) (c(n,h) mod2).                 (10)
```

For completeness, use R itself on its nonfixed pairs. Within each
self-converse H fibre pair consecutive elements in any fixed ordering,
leaving one exactly when that fibre is odd. This commutes with R and
attains (10); parity on the R-fixed locus gives the matching lower bound.
Thus P1 makes Delta_B eventually periodic with a dyadic period bound.
This is a sharp count of fixed classes, not a native metagraph matching.

The independent complete labelled-orbit census through order6 below gives
the following nontrivial strong inventory. Each row lists
`(order, ordinary multiplicity, self-converse multiplicity)`:

| h | complete strong inventory |
|---:|---|
|3|(3,1,1)|
|5|(4,1,1)|
|7|empty|
|9|(5,2,2)|
|11|(5,1,1)|
|13|(5,1,1)|
|15|(5,2,2), (6,2,0)|
|17|(6,2,2)|
|19|(6,2,0)|
|21|empty|

Its **all-order completeness** uses the **CITED** Busch bound, not an
extrapolation: [Theorem1 and Corollary1, EJC13 (2006), N3](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
give minimum25 at order7 and larger minima thereafter. The primary PDF was
reopened 2026-11-13. The same low-H inventory and the necessary cutoff are
**inherited**, independently replayed here; they are not a new census claim.

Equations (1) and (3) now give, for n>=6, the odd self-converse fibres

```text
h=1: always;          h=3: n odd;        h=5: n even;
h=9: floor(n/2) odd;  h=11,13: n odd;
h=7,15,17,19,21: never odd.
```

The same total formula holds also at n=4,5 by their direct census. Hence

```text
Delta_21(n)=2+2(n mod2)+(floor(n/2) mod2),  n>=4.     (11)
```

The residue classes n=0,1,2,3 modulo4 give respectively `2,4,3,5` fixed
classes. Its least eventual period is exactly4. The onset4 is sharp:
Delta_21(3)=2, while (11) would give5. At n=6 the complete unrestricted
minimum defect is6, whereas this H<=21 defect is3; truncation is essential.

## P4. An exact unbounded pairing-rank clock inside one H fibre

This is an explicit extraction from the low-H inventory and the inherited
[finite-action classification P1](../../paired_recursions_20261103/PAIRINGS.md),
not a new abstract action theorem. At every n>=6,

```text
a(n,19)=2(n-5),        c(n,19)=0.                    (12)
```

Indeed19 is prime, so its core has just one of the two strong order-six
classes, and singleton padding has n-5 positions. Those two classes are a
converse pair. No self-converse center can have H19, and19 has no nontrivial
square divisor, so no order admits a self-converse H19 class.

On the nonempty order-n H19 fibre, the largest rank of a free elementary
abelian 2-group action commuting with converse is therefore exactly

```text
d_19(n)=1+v_2(n-5).                                 (13)
```

The action preserves H because its whole carrier is a single H fibre. Its
construction remains the inherited choice of numbered class blocks; no
native metagraph move is asserted.

The integer sequence d_19 is unbounded and **not eventually periodic**.
For any proposed positive period P choose j>v_2(P) and n=5+2^j beyond its
onset. Then d_19(n)=j+1, but d_19(n+P)=v_2(P)+1. In contrast, for every k>=1
the capped observable min(k,d_19(n)) has exact period `2^(k-1)` (period1
for k=1). For k>=2, divisibility by `2^(k-1)` determines the maximal cap;
the two residue classes separated by `2^(k-2)` distinguish the last bit.
This is a concrete warning against removing the fixed-precision qualifier
from P1, even on one fixed-H fibre with a complete exact inventory.

## Connection contract, failure boundaries, and checks

**Source:** ordered nontrivial strong SCC words plus singleton-gap lengths.
**Target:** the vector of joint ordinary/self-converse fibre counts modulo
2^k, and at k=1 the optimal H-preserving converse-compatible pairing defect.
**Map:** exact SCC padding binomials, then reduction modulo2^k; at k=1
cancel the numerator at x=1 to obtain the least clock. **Preserved:**
order, actual H-product, ordinary class multiplicity, converse-fixed
parity, and exact optimum fixed-point count on the specified cutoff.
**Destroyed:** individual tournament identities, native adjacencies,
insertion/root incidence, integer multiplicities after modular reduction,
and every H>B fibre. **Needed sidecars:** actual strong inventories,
ordered core lengths/orders, center types, and numerator multiplicities.
**Cheapest decisive tests:** H9 parity cancellation, H15's nonempty even
fibres, H19's empty self-converse locus, and n=3 failure of the eventual
Delta21 formula.

The toothpick comparison stops at this point: the present binary clock is
a theorem about residues of finite-core padding counts. It supplies neither
a physical endpoint decomposition nor a complete class-growth recurrence,
and no map to the toothpick segment population is asserted. The unrestricted
Hamilton spectrum and a compact full A000568 recurrence are untouched.

Reproduce from the repository root:

```sh
python3 -I -B 05-knowledge/results/joint_metagraph_dyadic_20261113/pairings/check.py
python3 -I -B -O 05-knowledge/results/joint_metagraph_dyadic_20261113/pairings/check.py
```

The checker imports no repository computation or class catalogue. It
partitions **every labelled orientation** at orders0..6 into full S_n
orbits, checks the orbit mass, computes H independently by memoized path
extension and literal permutations, and tests converse on actual orbits.
It then compares two independent counting paths: ordered finite-core
padding and first-SCC dynamic convolution through n=160. It checks the
joint dyadic law at precisions1..4, least parity periods, every low-H
coefficient against the actual census, and an explicit optimal commuting
involution on every class through6. Normal and optimized receipts are
byte-identical; [receipt.json](receipt.json) records the code hash and
consequence objects. [Source pins](source-pins.json) freeze dependencies.
