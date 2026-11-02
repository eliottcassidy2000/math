# Session141: source-preserving recovery of twenty and the 1–2–4 chart

Status: **INHERITED / PROVED elementary deductions / FINITE-EXACT controls**.
Source pin: `c190e11b7795ed4e2d2f3aeb667e86d08238cdbf` in
`eliottcassidy2000/math`. No novelty or priority claim. The all-size extension
in section 5 is proved here; the preceding classifications are recovered,
not rediscovered. A transcript phrase is evidence of the question, not a
mathematical premise. The actual user lists are preserved in `check.py`.

## 1. The normalization and the two meanings of “three terms”

For positive n put sigma(n)=sum_{d|n}d, s(n)=sigma(n)-n, and
Delta(n)=s(n)-n. Proper divisors include1. Thus

    s(20)=1+2+4+5+10=22,    s(20)-1=21=20+1.

The displayed `2+4+5+10` is the sum of the nonunit proper divisors, not
s(20). The older `aliquot_divisor_overlap_20261028/THEORY.md` explicitly
preserved this distinction. Its predicate is

    C(n): there exist positive a,b with ab=n and s(n)=1+a+b.

**Inherited all-size theorem:** C(n) iff tau(n)=4 OR Delta(n)=2; the branches
are disjoint. In the first branch n=p^3 or pq for distinct primes. In the
second the pair is (1,n), so its compressed sum repeats1 and includes n;
it is not the literal list of proper divisors.

Proof for reference: if 1<a<b, the three distinct proper divisors1,a,b
already sum to s(n), hence no others exist. If a=1 the equation is exactly
s(n)=n+2. If a=b is prime, s(a^2)=1+a<1+2a. If a=b is composite, a has
least prime factor p with a>=2p, and the distinct proper divisors1,a,a^2/p
have sum at least1+3a>1+2a. Conversely tau4 has exactly the two stated
shapes. Prime cubes are deficient; pq has excess2-(p-1)(q-1)<=0, so the
branches do not overlap.

At tau6, shapes are p^5 or p^2q. Prime powers are deficient; the second
shape has excess2 iff q=1+2p/(p^2-p-1). Since q>=2, p is2 or3; these give
q5 or11/5. Therefore20 is the unique source with FIVE proper divisors
satisfying C. This was already proved in the pinned overlap packet.

## 2. What the user's lists actually support

The first supplied list is `16,10,15,14,21,35,22,33`; every entry except16
is a product of two DISTINCT primes, hence has exactly three proper divisors.
The second list `12,18,20,28,24` consists of two-prime sources with at least
five proper divisors, not one fixed arity:24 has seven. Its first four
have five. The negative list `4,9,2,8,49,121,27,169,16` consists entirely
of prime powers. In particular16 occurs in both the first and negative
lists. No Boolean predicate of n alone can simultaneously accept and
reject that occurrence. Moreover8 and27 satisfy C although in the negative
list, and20 satisfies C although in the second list. C is an overlap
coordinate, not a literal reconstruction of all list memberships.

A consistent recovered organization gives prime powers priority, then
separates distinct semiprimes from other mixed-prime exponent shapes;
retain C as a second flag. This does not assert the user's intended rule.
A plausible reason16 bridges the lists is its inverse fiber:

    s^{-1}(15)={16,33};  s^{-1}(16)={12,26}; s^{-1}(22)={20,38}.

These exact fibers connect different exponent shapes without identifying
them. They cannot establish why a numeral was put in a particular list.

## 3. The meaningful roles of15,20,24,26

The map must be named before the numbers can be associated:

| source n | factorization | proper terms | s(n) | sigma(n) | Delta(n) |
|---:|---|---:|---:|---:|---:|
|15|3*5|3|9|24|-6|
|20|2^2*5|5|22|42|2|
|24|2^3*3|7|36|60|12|
|26|2*13|3|16|42|-10|

Complete small fibers, certified by explicit finite bounds in `check.py`:

    s^{-1}(7)={8}; s^{-1}(15)={16,33}; s^{-1}(16)={12,26};
    s^{-1}(20)={34,361}; s^{-1}(21)={18,51,91};
    s^{-1}(22)={20,38}; s^{-1}(24)={529};
    s^{-1}(26)={46}; s^{-1}(28)={28}; s^{-1}(36)={24};
    sigma^{-1}(24)={14,15,23}; sigma^{-1}(42)={20,26,41};
    sigma^{-1}(60)={24,38,59}.

The sigma42 fiber also has a short exhaustive factorization proof. Each
prime-power block p^e exactly dividing n must have sigma(p^e) dividing42.
For e1, the possible blocks are2,5,13,41, with sums3,6,14,42. For e2,
the sum is odd and divides21; only p2 gives a block, namely4 with sum7.
For e3 only p2 or3 could have sum<=42, and their sums15 and40 do not
divide42. For e4 only p2 is small enough and gives31; e>=5 has sum>=63.
The allowable blocks can multiply to total42, respecting distinct prime
supports, only as4*5,2*13,or41. Hence the displayed fiber is complete,
and its aliquot images are22,16,1 respectively. This proves a genuine
shared-sigma relation, NOT any isomorphism of the sources' divisor sets.

In particular20 and26 share TOTAL divisor sum42, not aliquot sum. At a
fixed sigma value, s(n)=sigma-n is an anti-diagonal in the (n,s(n)) plane;
keeping only the sigma value loses the source and excess. Inverse s fibers
and inverse sigma fibers are different partitions. Twenty reaches the
1–2–4 decomposition through its literal aliquot tail:

    20 ->22 ->14 ->10 ->8 ->7 ->1 ->0.

The other relevant exact tails are

    26 ->16 ->15 ->9 ->4 ->3 ->1 ->0,
    24 ->36 ->55 ->17 ->1 ->0.

We stop at0; s(0) is not defined here. This is a finite orbit calculation,
not a theorem about all aliquot dynamics.

## 4. The geometric-series map, and the precise prime guards

For prime p and e>=1,

    s(p^e)=1+p+...+p^(e-1),  sigma(p^e)=1+p+...+p^e.

Thus7=1+2+4=s(8),15=1+2+4+8=s(16)=sigma(8), and
26=2(1+3+9)=2s(27). The latter is a scalar multiple, not s(26).
Also24=sigma(15)=(1+3)(1+5), whereas s(24)=36. These identities
have different carriers and are not interchangeable.

For a>=1 and an ODD PRIME q,

    s(2^a q)=(2^a-1)q+2^(a+1)-1,
    Delta(2^a q)=2^(a+1)-1-q.

The a2 chart gives12->16,20->22,28->28 at q3,5,7; the a3,q3 chart gives
24->36. The semiprime a1 chart gives26->16. In general its output agrees
with the a2 chart exactly when the second odd prime r=3q+4 is prime:
` s(4q)=s(2r) iff r=3q+4 `. Thus q3 produces12/26 and q5 produces20/38;
q7 produces r25 and fails the prime guard. Likewise
` sigma(4q)=sigma(2r) iff 3r=7q+4 `; q5 gives20/26. These are guarded
fiber correspondences, not conjugacies of aliquot iterations.

For q composite, even the displayed affine formulas can fail: treating
q9 as prime gives s(18)=12 instead of its actual21. The same issue applies
to geometric-series expressions at nonprime bases.

The pinned 1–2–4 packet additionally proves a MARKED bridge: on the ordinary
Collatz cycle{1,2,4}, Collatz is multiplication by4 modulo7, inverse to the
Paley7 triangle traversal1->2->4->1. In F8=F2[a]/(a^3+a+1), the map
r->a^r turns residue doubling into Frobenius; {a,a^2,a^4} has rank2 and is
a Fano line, NOT a normal basis. No global Collatz descent modulo7 follows.
The source snapshot contains its proof and counterexamples. Cardinality
coincidences involving15,24 or26 do not supply any missing transfer map.

## 5. Extension of the recovered low-arity excess-two classification

**PROVED.** The complete set of positive integers n with tau(n)<=17 and
s(n)=n+2 is

    {20,104,464,650,1952}.

Equivalently, these are all excess-two sources with at most SIXTEEN proper
divisors. The pinned overlap packet already proves the claim for tau<=12;
we retain that attribution and give the remaining all-size case split.

The inherited two-prime theorem says that every excess-two integer with
at most two distinct prime factors is

    n=2^a(2^(a+1)-3), a>=2, with2^(a+1)-3 prime.

Its short proof uses Delta(2^a p^b)+1=(2^(a+1)-p)R_b(p), where
R_b=1+p+...+p^(b-1). At excess2 the product is3; b>=2 would give
R_b>=p+1>=4, impossible. An odd two-prime source is deficient since its
abundancy is strictly below(3/2)(5/4)<2, and a prime power is deficient.
For tau<=17 the family permits a<=7. The guarded parameters are
5,13,29,61,125,253 at a2..7; the last two are composite. This accounts
for20,104,464,1952 and no others with at most two prime factors.

The additional tau values13..17 have the following exponent patterns
(up to permutation, generated independently in `check.py`):

    13:(12); 14:(13),(6,1); 15:(14),(4,2);
    16:(15),(7,1),(3,3),(3,1,1),(1,1,1,1); 17:(16).

Only tau16 introduces shapes with at least three primes. For a squarefree
four-prime source, an EVEN n has4|sigma(n), whereas2n+2 is2 modulo4. For
an ODD n the maximum abundancy occurs at3*5*7*11 and is
(4/3)(6/5)(8/7)(12/11)=2304/1155<2. Thus that shape is excluded.

It remains n=p^3qr with distinct primes. If n is even, either p2 and
(q+1)(r+1) is divisible by4, or q2/r2 and sigma(p^3) is divisible by4.
Again4|sigma(n) contradicts2n+2. Suppose n is odd. For p>=5,

    sigma(n)/n <=(156/125)(4/3)(6/5)=1248/625<2,

because sigma(p^3)/p^3 decreases with p and the two remaining distinct
odd primes have reciprocal-product bound(4/3)(6/5). This upper bound
may admit an illegal repeated prime, which only makes it a weaker valid
upper bound. For p3, write n=27qr with distinct odd q,r>=5. Excess2 is

    40(q+1)(r+1)=54qr+2,
    (7q-20)(7r-20)=533=13*41.

Both factors are at least15. The positive factor pairs of533 are(1,533)
and(13,41), neither admissible. Therefore no such source exists. This
exhausts the exponent patterns and proves the all-size extension.

The boundary is genuine:130304=2^8*509 has tau18 and excess2, so extending
the same five-element list to seventeen proper divisors would be false.
No classification of the entire tau18 layer, all excess-two sources,
odd excess-two sources, or infinitude of the guarded prime family is made.

## 6. Inheritance and verification contract

Closest mechanism: the source-preserving C classification and exact inverse
bound; hostile:16 at conflicting list memberships,8/27 for C versus
prime-power priority,650 for a two-prime extrapolation,130304 for our arity
boundary. Least-used sidecar: source exponent shape at a common inverse
fiber. The live board is divisor arity / excess / sigma fiber / s fiber /
geometric-series address / marked finite-field action.

The source map sends a labeled divisor list to its sum; it preserves that
sum and loses the prime labels and multiplicities. A fiber connection owes
those labels as sidecars. The marked finite-field map preserves Frobenius
orbits, but needs marking and linear rank; it is not an aliquot-map quotient.

Run `python3 -B check.py > normal.json`, `python3 -B -O check.py > optimized.json`,
then `cmp normal.json optimized.json`; use `python3 -B verify.py --replay`
for manifest/source-pin and replay verification. All code is standard-library
Python, all failures are explicit (active under -O), and no network or original
checkout is required for replay. A sieve and an independent multiplicative
oracle agree on every n<=200000. Finite checking does not replace section5's
all-size proof. Source snapshots are deliberately limited to consulted files;
the giant historical mistakes ledger was not read or copied.
