# Ordered even/odd phase and a physical switch between count ranges

**Status:** exact two-unit carry, exact continuation jump, and an exact
eventual spectrum of range differences for one shared-suffix grammar.
No new prime-counting formula or universal Chan–Pak coverage is claimed.

## Möbius inversion and an order-sensitive alternative

The ordinary Möbius function is zero whenever **any** prime square divides
n, not just at n itself equal to a prime square. On squarefree n it is
(−1)^k, with k the number of distinct prime factors. Thus
mu(2^v m)=mu(m), −mu(m), 0 for odd m and v=0,1,≥2, respectively.
The exact inversion uses
```
J(x)=sum_(p^k <= x) 1/k,
pi(x)=sum_(n <= log_2 x) mu(n) J(x^(1/n))/n.
```
Substitution and the identity sum_(d|k) mu(d)=[k=1] prove it,
with the convention that endpoints are included. Replacing J by
li gives the smooth Riemann-type approximation, not the exact
prime count: the actual explicit formula has zero terms and other
corrections. See [Riemann's manuscript](https://www.claymath.org/collections/riemanns-1859-manuscript/).

The proposed even-then-odd vs odd-then-even signs cannot be a
scalar function of the product: 2·3 and 3·2 are the same n, and
their phases would have to be opposite. Also, the product of two even
integers is divisible by 4 and would receive zero from ordinary mu,
not +1. In a weighted ordered sum
with symmetric weights, they cancel in pairs. To retain the sign,
keep the ordered inputs, or use two noncommuting channels. A
minimal explicit implementation labels the even and odd channel
by the [Pauli matrices](https://quantum.cloud.ibm.com/learning/en/courses/foundations-of-quantum-error-correction/stabilizer-formalism/pauli-operations-and-observables):
```
E = X, O = Y,       E²=O²=I,       EO=iZ, OE=−iZ.
```
The odd-odd phase +1 is a choice here, not a deduction from the
prompt. Reading the Z=+1 coordinate gives the scalar two-letter
rule; to concatenate more letters, preserve the operator rather
than replacing iZ by i. This is an associative bookkeeping model,
not a replacement for mu in the prime formula. In tournament
terms, projecting a witness to its count is precisely the kind
of projection that can erase the two ordered channels.

## Exact same-height alternatives with distinct tape products

This sharpens the [one-unit carry](DISTINCT_TAPE_PRODUCTS.md).
Use any physical two-tape core with deletion count g, two-cover
count B≥g, and reads (h,U,V,Ubar,Vbar) as in
[the folded-center transfer](FIGURATE_OVERLAP.md). Pick any
positive r with rg≡1 mod4, and set
```
a=8rB−(rg+1)/2,     b=8rB+rg.
E=(4a,4b),          F=(4a+2,4b−2).
```
These positive pairs are ordered and balanced (larger ≤ twice
smaller), since B≥g. The sum is fixed; the product rises by 12rg.
Each address is even, so both signed boundary terms vanish.
Write b=4z+1: a is odd, and the binary checksum recurrence gives
```
C(4a+2)−C(4a)=1−a,
C(4b−2)−C(4b)=−2z−2.
```
Their sum is −12rB by the definitions of a,b. Thus the product
term B·12rg and the checksum term g·(−12rB) cancel identically.
The **entire** center change is 2(Ubar−U), with no search.
For symmetric unsigned reads U=Ubar this is exact equality.
It works at every depth and also with the singleton core
(g=B=1, all four endpoint reads zero).

In the folded construction, the repair palette is the same
for E and F. At depth 65, all 14,725 positions of its certified
central interval coincide exactly. This is not an increase in
the number of distinct heights: the gain is alternative marked
physical continuations at the same height.

For the singleton core, r=1 gives tapes (28,36) and (30,34).
The product/checksum pairs are (1008,54) and (1020,42), so both
counts are 1063. Literal subset DP gives rooted profiles
(1063,1,1847,1) and (1063,1,1967,1) on order 12. Trial division
confirms 1063 prime. The distinct physical tapes (64,114) give
7441=7·1063, independently recounted, with no claim of a paired
interface.

## Exactly which continuation bit changes

Both tape widths stay fixed during the move: 4a→4a+2 cannot
cross a power of two, and 4b→4b−2 cannot do so since b≡1 mod4
and b>1. Write w for the right width. The change in right tape
incontacts is at its low positions: 4b ends in 100, while 4b−2
ends in 010, with all higher bits the same. Hence the width-w
bit reversal increases by 2^(w−2)−2^(w−3)=2^(w−3).
The right tape's outward weighted start is the complement
of that reversal. Since the left address grows by two, the
complete weighted root start t changes by exactly
```
Delta t = 2^w·2B − g·2^(w−3) = 2^(w−3)(16B−g).
```
Because g is odd, the two continuations agree in their lower
w−3 bits and disagree at bit w−3. This is a physical even/odd
channel assignment: strip the shared lower bits and read parity.
In the small example both t values are odd and total H is the
same, so neither H nor the parity of t identifies the channel.
The first differing bit is bit 3, as the formula predicts.

A common ordinal sink resets the ordinary root end count to
zero without changing the height, doubles the weighted starts,
and leaves the deletion count unchanged. Now follow both with
the same [retained-root word](ANS_INVERSE.md). If W is the word
read least significant bit first, the final counts differ by
```
D·W,     D=2^(w−2)(16B−g).
```
This follows also from the linear transfer: the common checksum
depends on the deletion and word, not on the initial t. For odd
W the first differing bit of the returned counts is bit w−2;
more generally it is w−2+v2(W) for W>0. The small experiment
literally recounts 4757/4997 at W=1 and 12143/12863 at W=3.

This is a rule between **absolute** ranges, not just a parity
label. Given two requested outputs of one word and a fixed common
source height, subtract them and divide by D to get the sole
possible W. Then compute one exact returned count to check
the absolute baseline. The phase information alone is insufficient.

## Exact interval-difference test, and why it does not suffice

Suppose both initial representations have the full repaired
interval H0+P·[l,u]. Under one common suffix W each interval
is translated by a fixed amount. Its width remains P(u−l),
and the two lower endpoints are separated by DW. Allow different
repair digits at the two ends. The exact differences over this
language are therefore
```
DW + Pq,      W >= 0,      −R <= q <= R,         R=u−l.
```
For the folded core, P=2M with M odd, g divisible by every prime
factor of M, and B coprime to M. Thus gcd(D,P)=2. The relative
displacement of the two windows can realize **every even
residue modulo P**, by choosing W mod M. Conversely, for a
prescribed gap Delta and a fixed q the exact W is forced.
Equivalent target algorithm: first solve W mod M, then check
whether some nonnegative W in that progression leaves
|Delta−DW|≤PR.

In fact **all sufficiently large even gaps** occur in this grammar
if and only if 2R+1≥D/2. Proof: as q runs through consecutive
integers it visits min(2R+1,D/2) residues of Pq mod D, and each
one yields the entire sufficiently large tail of that class as W
increases. Otherwise the fraction of even gaps covered in the
tail is exactly min(1,(2R+1)/(D/2)). In the present family R grows
linearly in depth but D grows much faster, so not even the necessary
difference coverage is universal. There is an even sharper limit: after any **nonempty numerical**
suffix W>0, the two returned windows are disjoint. To share a
count they need DW divisible by P and DW≤PR. The first forces
W≥M, so the second would require D/2≤R, which fails here.
Absolute selection would still require putting one of the ends
at the target, even if the difference threshold were met. A constructive next target is
to enlarge R faster than the first differing-bit scale, or to
choose several compatible marked suffixes independently rather
than requiring the same W.

## Bonus: overlaps at every congruent product gap

For the folded core let g=GM, P=2M and G=1105.
Fix any L≥2 with 2L≡1 mod M. In the same two-unit move take
b=4z+2eta+1, a=b−L, eps=(1−L) mod2, eta=0 or 1. The
product change is 4(2L−1) and the checksum change is
−6z+L−1−eps−2eta, directly by the two recurrence evaluations.
The normalized center change is K−3Gz, with
```
K = 2B(2L−1)/M + G(L−1−eps−2eta)/2.
```
Choose z nearest to K/(3G). Since 3G=3315, the separation is
at most 1657 repair steps. The bound B/g≥2 suffices to give
z≥L/2 (use the displayed K and subtract 1/2 for rounding),
so both tape pairs are positive and balanced. Their common
interval therefore contains at least 13,068 positions at depth
65, for *every* such product gap. This is a proof of many
overlaps, not a proof that centers occur at every desired height.

Replay: `python3 exact_product_carry.py` from this packet.
It includes independent subset counts, both continuation channels,
and exact bigint controls for arbitrary-scale and rounded carries.
