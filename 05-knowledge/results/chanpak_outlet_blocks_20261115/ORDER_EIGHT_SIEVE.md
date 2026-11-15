# An order-eight palette obstruction on the fixed progression 1 mod 390

The unconditional presentation bound of the transitive-spine outlet-block
family crosses one half when arbitrary strong order-eight cores are added.
This does not mean the new cores have supplied the missing targets.
A target-first divisibility sieve and the inherited modular checksum
mixing recover a strict upper bound in an explicit fixed progression.

For the inherited complement-star address e and an arbitrary outlet block
B, the exact count is N = h + 2P e − hJ(e), where (h,P) is its ordinary
and two-optional-slot profile. Let M be odd. A profile with
gcd(h,P,M)>1 cannot hit any prescribed progression N=1 mod 2M, since
every count in its chart is divisible by that gcd.

## A fixed-admissibility row has exactly its predicted progression mass

If gcd(h,P,M)=1, then in any odd progression mod 2M the presentation
mass is exactly 1/M of the unrestricted row mass:
```
lim_(X→∞) [M/X] #{e≥1 : N(e)≤X, N(e)≡b mod 2M}
 = (1/h) log[2P/(2P−h)]       (b odd).
```
The claim concerns addresses, not distinct values. Here is the argument
connecting the physical row shape to the modular programme.

Set q=2^k near √X and e=q a+x, 0≤x<q as in OUTLET_BLOCKS.
The split identity has J(e)=J(x)+J(a)+a·rev_k(x)+epsilon, where
epsilon = [top bit of x]·[low bit of a]. The real cutoff in a
is one initial interval up to O(1) count errors per x, and its
endpoint is X/[q(2P−h rev_k(x)/q)] + O(1). On each low-word x,
the modular condition is the one in
```
Pq a + P x − h J(e)/2 = (b−h)/2 mod M .
```
The discrepancy of the high digits (a mod 2M,J(a)/2 mod M) from
the uniform law is sublinear on any long interval. This is the
inherited prefix-programming and interval-mixing theorem: the pair
(a mod 2M,p(a)/2 mod M) has a uniform stationary digit law because
each state has two equally weighted incoming digit steps, and the
fixed-length suffix theorem gives a strictly positive common
transition power, hence contraction on aligned dyadic intervals.
Every interval is a sum of at most two aligned dyadic intervals per
scale. Replacing p(a)/2 by J(a)/2 = floor(a/2)−p(a)/2 preserves uniformity.
The discrepancy estimate applies irrespective of the low-word x
parameters, since for fixed M only finitely many different tests on
the state are possible.

Write d=gcd(h,M). Reducing the congruence mod d, exactly a fraction
1/d of the a coordinate classes satisfy it, since gcd(Pq,d)=1.
For each such class, the function J(e)/2 differs from J(a)/2 by
a constant depending only on (x,a mod 2M); division by 2 is legitimate
as J(x) and a·rev_k(x)+epsilon are both even. Thus exactly a fraction
d/M of the checksum classes satisfy the full congruence. Summing the
sublinear errors over q low words gives o(X), and the same real
Riemann sum proves the stated mass.

## A fixed small modulus reaches the strong order-eight boundary

The inherited complete class census for order eight gives 1,567
distinct strong profiles, and the strong palette through order eight
has 1,721 distinct profiles after including the independent labelled
checks through seven. In each composition the ordinary and optional-slot
counts multiply on the SCC chain. A component having gcd(h_i,P_i,M)>1
therefore eliminates the product chart from progression 1 mod 2M.
Some remaining products can still be eliminated because a factor from
h_i matches a factor from P_j of a different core. Retaining these
invalid mixed profiles and charging the same positive mass to them
only makes the upper bound weaker.

Choose **M=195=3·5·13**. Exactly 1,034 of the 1,721 profiles survive
the single-core sieve gcd(h_i,P_i,M)=1. Use this list in (4) of
OUTLET_BLOCKS. For fixed M the same majorant applies to normalized
presentation counts, with only the finite multiplicative factor M;
hence the limit may still pass through the unbounded SCC-multiplicity
sum. The exact rational certificate is less than .4975. Since
there are X/(2M)+O(1) targets at most X in this progression, after
renormalizing this bounds the attainable fraction by strictly less
than .995: **more than half a percent of all 1 mod 390 targets are
missed in lower density** even allowing arbitrarily many SCCs of
order up to eight in B.

One can exclude all 739 nonprimitive individual profiles by taking
M as the product of all their odd gcd prime factors (23 primes).
Then 982 profiles survive; the certified bound is .485574, giving
roughly 2.88 percent missed in that more dilute progression. The
smaller explicit progression is the stronger practical placement
diagnostic. No integer is excluded from tournaments as a whole.
Nor does local non-solvability explain this restriction: the
singleton outlet has all odd-compatible residues, including 1 mod
390. The issue is that real presentation mass in that progression
still leaves a positive gap.

## Reproduction and provenance

The inherited order-eight canonical class file is
`05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/n8_classes.csv`
at commit e42efb712986613af94481e9f57253f1f6bef4b7. Its SHA-256 is
5226fb6d6cb71e2c99dc50536744f044c31af55e173e3be9f5f722a39edf3bc2.
The order-eight completeness is inherited, not claimed from the small
orientation replay. To recheck exact profiles and rational bounds:
```sh
python check_order8.py /path/to/n8_classes.csv
```
Omit the argument to use the extracted profiles8.txt directly.
The earlier n≤7 profile census and the small literal verifiers are
independent of the order-eight class source.
