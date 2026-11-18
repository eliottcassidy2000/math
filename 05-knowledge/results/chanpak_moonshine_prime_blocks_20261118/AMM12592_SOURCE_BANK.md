# A source can add capacity; a change of prefix cannot

**Status:** quantitative all-height theorem for the specified retained-root
alphabet; exact out-of-training-block benchmark. No universal Chan–Pak claim.

## Inheritance from AMM 12592

[Elliot Glazer's AMM 12592](https://sites.math.rutgers.edu/~zeilberg/EM26/AMM326.pdf)
asks for exact fair extraction with an unknown coin bias under a
deadline tied to the initial constant run. The relevant prior research
is [orientation splitting and the dyadic donor](../../../01-canon/theorems/THM-3344-orientation-splitting-saves-exactly-one-dyadic-donor-bit.md),
and the [audit of unrestricted tail routing](../../../01-canon/theorems/THM-3024-degree-resolved-hall-closes-the-general-class-amm-floor.md).
The former uses the initial orientation to divide supply that a
coarser ledger hides. The latter warns that apparent failure at a
truncation boundary need not be true infeasibility, and that aggregate
supply is useless unless the needed class can receive it.

Here the deadline is the exact upper bound on how long a word
can first land below the target. The supply is physical construction
words, and the receiver is a particular absolute integer. The
analogue of splitting the initial orientation is to retain the first
construction bit. There is no random coin here, and the argument
below does not require a transport conjecture for AMM.

## Freeze the sources, then move forward

The baseline is the three rooted profiles (29,11,31,13), (23,2,50,9),
(25,2,50,9) from the [exact inverse decoder](ANS_INVERSE.md).
The four entries are (Hamiltonian count, paths ending at the root,
weighted starts at the root, deletion count). Candidate sources
were the distinct profiles with h>1 and t≥2d from the archived
[456 order-seven representatives](../actual_hamilton_multiplicities_20261021/session24-multiplicities/n7_representatives.txt).
There were 2,806 candidate profiles before excluding any baseline
profiles. Choose eight, greedily maximizing the number of *additional
distinct odd values* in [1024,4096), with lexicographic profile ties.
The bank and all thresholds stayed fixed above 4096. This is a
single retrospective experiment with an explicitly separated
training slice, not a claim of statistical independence from wider
repository history.

| Untouched half-open target block | Old odds | With new bank (odds) | Old primes | With new bank (primes) | Total primes | Old both p,7p | With new bank both |
|---|---:|---:|---:|---:|---:|---:|---:|
| [4096,16384) | 962 | 5451 | 238 | 1184 | 1336 | 133 | 1062 |
| [16384,65536) | 3923 | 21931 | 714 | 4098 | 4642 | 437 | 3679 |
| [65536,262144) | 15463 | 87585 | 2546 | 14458 | 16458 | 1529 | 13041 |

The last two columns allow distinct sources and words for p and 7p.
The upper baseline includes the known literal sevenfold copy of
(29,11,31,13), so any old native pair is credited. The upper bank
is likewise frozen before scoring. "Additional" means added to
this baseline, not globally novel counts. For example the newly
selected (11,4,12,5) source has a five-vertex realization with
masks (6,28,24,17,1) and root 4. Words 001111101 and 100100100011
give the prime 4111 and its sevenfold partner 28777, checked
by literal independent Hamiltonian-path counts.

The first-hit deadline permits exact labels for this restricted
language. Strip all trailing zeros: they do not affect the Hamiltonian
count. After appending a 1, H'=h+t−e≥t. Both t updates strictly
increase t because t>d. A partial word whose t is already ≥ the
upper endpoint cannot contribute any new total in the block.
Thus the negative entries are not artifacts of word-length truncation.
They also do not rule out other tournament constructions.

## The all-height capacity theorem

Fix one physical seed S=(h,e,t,d) with t>d and d>0. Count **words**
rather than distinct counts, using only nonempty words ending in 1
so that padding by 0 is never counted. For any fixed 0<a<b, the
number whose outputs lie in [a·2^m,b·2^m) is
```
(b−a) 2^(m−1) · kappa(S) + O(2^(m/2)),
kappa(S) = (2/d) log(t/(t−d)).
```
The logarithm is natural. The implied constant may depend on the seed
and the interval. In particular this is an upper bound on the
capacity to fill *distinct* odd totals in every such interval.
It is not a prime-density bound.

**Proof.** For a word of length k split after p=floor(k/2), with
q=k−p. The normalized t value at the cut, A=t_p/2^p, is uniform
on the equally spaced grid t−d·V/2^p, 0≤V<2^p; here V reads the
first bits most significant first. Independently, B=W/2^q is
uniform on the grid in [1/2,1), because W reads the last q bits
least significant first and the last bit is 1. The exact suffix
formula gives
```
H/2^k = A B + O(2^(−p)+2^(−q)).
```
Indeed the initial H at the cut is at most t_p; the suffix debit
is at most d(2^q−2)+d once p≥1. Couple each grid to a corresponding
uniform variable on its interval, by jittering in the adjacent
grid cell. The change in AB is O(2^(−p)+2^(−q)). Its continuous
density f is bounded because A≥t−d>0. Thus the empirical mass in
any interval differs from the product-uniform mass by O(2^(−k/2)),
uniformly over endpoints. This is a deterministic grid count; the
random variables describe its empirical distribution.

At every eligible length k=m+j, the t bounds just after and
before the final 1 imply (t−d)2^(k−1)≤H≤t·2^k. Only boundedly
many j can enter the height block. Length k carries 2^(k−1)
words; summing the density on the rescaled interval gives
```
2^(m−1) integral_[a,b] sum_j f(z 2^(−j)) dz + O(2^(m/2)).
```
Conditional on A=u, the product density is 2/u on [u/2,u].
For any z>0 and almost every u, **exactly one** j has
z 2^(−j) in [u/2,u], so the sum of rescaled densities integrates
to (1/d) integral_[t−d,t] 2 du/u = kappa(S), independent of z.
This is where the exact dyadic stacking occurs. The seed itself
contributes at most one extra distinct count. QED.

## Two consequences for a selector

1. **Conservation under refining the source.** The two children
   have t values 2t and 2t−d, with the same d. Their log capacities
   add exactly to the parent's:
   the split point 2t−d cancels in the sum of the logarithms.
   A prefix-free set of descendant seeds can spend but cannot
   increase this mass. This also follows from the disjoint
   normalized t subintervals corresponding to the prefixes.
   A new seed with the same t,d may still shift h and avoid
   existing collisions; if it is literally in the parent language
   it cannot add new outputs. Refining a decoder is not the same
   thing as adding a new physical source.
2. **Initial-bit cuts.** Put g=gcd(t,d). The suffix expression
   over the full word gives H≡h−e b_0 (mod g). If only one of those
   classes is needed, the other first-bit branch is not a donor.
   The available word capacities per unit of the full odd stream
   in the two branches are respectively
   (2/d) log(t/(t−d/2)) and (2/d) log((t−d/2)/(t−d)).
   For example (5,2,6,3) has substantial capacity, but can never
   land on 1 mod 3 at any height. The other orientation actually
   supplies 0 mod 3, not the missing 1 mod 3. This is a physical
   cut, not an estimated density. The same reasoning refines any
   finite fixed prefix, but does not claim equidistribution modulo
   an arbitrary odd modulus.

A concrete obstruction is already stronger than the coarse binary
length bound: the baseline plus the top new source (3,1,3,1)
has word capacity
```
(2/13)log(31/18) + (4/9)log(50/41) + 2log(3/2) < 0.983.
```
Thus that four-source bank must miss a positive fraction of odd
totals in every sufficiently large fixed-ratio block, even before
collisions. The replay certifies the strict numerical comparison
with rational bounds on the atanh series. The full frozen bank
has more than unit word capacity; it still misses targets because
supply is not injective or aligned with the demanded counts. In
the largest holdout it has 267,394 terminal words for 98,304 odd
slots, but their projection hits just 87,585 distinct slots.

**Next proof test.** For a target rejected by the inverse debit
check, find an independent source that changes the exact debit
or shifts the count while preserving the continuation reads.
Spend existing certificates only when they actually reach the
defective class at the same height. A global surplus, a change of
prefix, and a congruence-only pass each fail to prove this. An
AMM-style replacement would require an exact physical map from an
overrepresented target to a missing one that retains the necessary
continuation profile. No such universal map is asserted here.

Replay: `python3 benchmark.py` checks the frozen bank and literal
witnesses. Pass `n7_representatives.txt` as an argument to also rerun
training; its SHA-256 is pinned in the script. `fixtures.py` is the
prior packet's independent subset counter.
