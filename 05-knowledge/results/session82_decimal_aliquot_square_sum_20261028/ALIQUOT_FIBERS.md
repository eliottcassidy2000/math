# Exact inverse aliquot fibers, prime pairs, and the square-deficit coordinate

Status: **PROVED** elementary reductions; **FINITE-EXACT** complete computed fibers; **CONDITIONAL/OPEN** global untouchability assertions. No literature-priority claim for the algorithms or identities.

## Inheritance and boundaries

We use `s(n)=sigma(n)-n` on positive integers, so `s(1)=0`; `s(0)` is not defined. Session 78 `ALIQUOT.md` already established finite fibers for targets `t>=2`, the bound `n<=(t-1)^2`, its prime-square exception, the distinct-prime formula, and prime-dilation reconstruction. Session 81 `ALIQUOT_OBSERVERS.md` already established prime-power reconstruction, valuation-profile inversion, and the multiplier 40 collision. None is presented as new here.

The new deliverable is a complete, reproducible computation of the requested million-scale fibers without a trillion-source scan, with a second implementation and a precise bridge from fibers to centered prime pairs. Live board: scalar fibers / factor-prefix inversion / prime-pair gap coordinate / parity and squares / directed forward chains. Closest inherited mechanism: factored divisor sums. Hostile: repeated primes and the enormous prime-square source. Least-used sidecar: largest prime and its exponent.

## 1. Complete factor-prefix inverse algorithm

Let `t>=2`. Every source is composite. Write uniquely

`n=m q^a`, where q is the largest prime divisor, `a>=1`, and every prime factor of m is smaller than q.

Coprimality gives the exact formula

**`t=s(m)q^a + sigma(m)(1+q+...+q^(a-1))`.**

If `m>1`, then `s(m)>=1`, so `sigma(m)<t` and hence `m<t`. For exponent 1 the candidate is uniquely

**`q=(t-sigma(m))/s(m)`.**

It must be integral, prime, and strictly larger than every prime divisor of m. For `a>=2`, the displayed polynomial is strictly increasing in the integer q>=2. Thus a sieve of sigma and largest-prime data only through t−1, followed by guarded polynomial inversion, suffices. Exponents can stop once the polynomial at `q=largestprime(m)+1` exceeds t. A binary search upper bound for all a>=2 is `floor(sqrt(t/s(m)))`. The m=1 branch is the geometric-sum equation `t=1+q+...+q^(a-1)`; a=1 is impossible, a=2 gives q=t−1, and a>=3 gives q<=sqrt(t). This proves termination, soundness and completeness for `independent_prefix_audit.py`. It does not assert optimal complexity or novelty.

The primary `check_aliquot_fibers.py` instead builds increasing-prime prefixes recursively. At state `(m,S=sigma(m))`, it handles a final exponent 1 using the linear formula. It loops over the next prime p and its exponent. Two lower bounds stop the p loop once both exceed t:

* final exponent>=2: `s(m)p²+S(p+1)`;
* non-final exponent>=1: `s(m)p²+2Sp+S` (replace the later prime by p to obtain a lower bound).

For a chosen full power p^a, state n=m p^a has exact divisor sum `Sn=S(1+p+...+p^a)`. Recurse only when `s(n)q+Sn<=t`, where q is the next prime after p. Any valid continuation contains such an nq as a divisor, and s is nondecreasing under divisibility; hence this pruning cannot remove a solution. Direct terminal equality is tested for exponent>=2. Each source has its unique increasing factorization path. The finite prime list ends at t: any prime divisor q of a composite source is proper, so `1+q<=t`. These facts prove the DFS complete, independently of empirical validation.

## 2. Exact requested fibers

Run `python3 check_aliquot_fibers.py` to generate the full source lists in `results.json` beside the script; generated JSON is not part of the published source package. This compact description is also exact. Let

`G(t)={p q : p<q are primes and p+q=t−1}`.

| Target | Complete fiber | Cardinality |
|---|---|---:|
| 1001 | `G(1001) union {11191}` | 28+1=29 |
| 1000001 | `G(1000001) union {2660799,13011919,16390159,30001459,77519461,190124401}` | 5402+6=5408 |
| 1001000 | `{718720,1001998998001}` | 2 |

Non-semiprime factorizations:

* `11191=19²*31`;
* `2660799=3*79*103*109`;
* `13011919=17*59*12973`;
* `16390159=17*829*1163`;
* `30001459=31²*31219`;
* `77519461=79²*12421`;
* `190124401=199²*4801`;
* `718720=2^7*5*1123`;
* `1001998998001=1000999²`, with 1000999 prime.

The last fiber exposes the inherited prime-square exception at full scale: its large source is exactly `(t−1)²`, while the other source is only 718720. Checking sources only up to t, or only through the ordinary `(t−1)²/4` branch, would miss it.

Direct target data (not inverse fibers):

* `1001=7*11*13`, `s(1001)=343`;
* `1000001=101*9901`, `s(1000001)=10003`;
* `1001000=2³*5³*7*11*13`, `s(1001000)=2143960`.

## 3. Prime-pair and square-deficit bridge

For distinct primes p<q, **`s(pq)=1+p+q`**. Unique factorization makes the map from unordered distinct prime representations of t−1 to the squarefree-semiprime part of F_t a bijection. This is not a bijection onto the whole fiber: t=13 already has the non-semiprime source27, and the new exact target lists exhibit further failures.

For odd t>=9 the relevant primes are odd. Put `h=(t−1)/2` and `r=(q−p)/2`. Then

**`n=pq=h²−r²`, with `r>=1` and both h−r,h+r prime.**

Conversely every such r gives exactly one semiprime source. Equivalently `(t−1)²−4n=(q−p)²`. This is a genuine square-deficit representation with a primality sidecar, not a square-sum theorem. Larger gaps give smaller sources at fixed t. The gap is between the two chosen primes, which need not be consecutive; no uniqueness or conjecture about consecutive prime gaps is involved. Dropping r/primality loses the source or the validity guard. The r=0 diagonal is forbidden: `s(p²)=p+1`, not2p+1. The hostile p=q=2 would incorrectly make target 5 touchable.

The separate square parity connection is also exact. In the factorization of sigma, an odd-prime geometric sum is odd iff its exponent is even; the power-of-two factor has odd sigma regardless of exponent. Therefore sigma(n) is odd iff n is a square or twice a square. Consequently odd n has even s(n) iff n is a square; even n has odd s(n) iff n is a square or twice a square. This supplies a valid filter, not a classification of fibers or a transfer of a sum-of-two-squares theorem.

## 4. Untouchable primes: what follows, and what remains open

Targets 2 and 5 have empty fibers (direct exhaustive small-target checks also certify this). One has `s(4)=3`, `s(8)=7`, and in general **`s(2^a)=2^a−1`** for every a>=1. Thus every Mersenne number is touchable, whether prime or composite. Every target p+1 for prime p is touchable through p². For odd prime p, p+3 is touchable through2p; p=2 is excluded because2p then repeats a prime.

**Conditional theorem.** If every even integer E>=8 is a sum of two distinct primes, then 5 is the only odd untouchable number, hence 2 and 5 are the only prime untouchable numbers. Proof: t>=9 odd has E=t−1 and the distinct-semiprime source; the remaining odd targets 1,3,7 have sources 2,4,8. The prime-only conclusion is weaker than the all-odd conclusion; neither converse is claimed. Ordinary Goldbach permitting equal primes cannot be inserted directly into this proof; its repeated-prime case requires additional work, and no such implication is claimed here.

Primary literature: [Pomerance–Yang, *On untouchable numbers and related problems*, 14 June 2012, introduction](https://math.dartmouth.edu/~carlp/uupaper3.pdf) explicitly states this distinct-prime conditional and the conjectural absence of further odd untouchables. Retrieved 2026-10-28; the search did not establish a later resolution. This session does not prove the conjecture, exhaust all primes, or treat a finite computation as a global proof. The paper also discusses prior fast untouchability algorithms, so no literature novelty is asserted for our implementation.

## Verification and stop conditions

Requires Python 3 and SymPy (tested with the version pinned in `requirements.txt`). Keep the two fiber Python scripts in one directory. All generated JSON paths are relative to the scripts, not the shell working directory. From this directory:

```
python3 -m pip install -r requirements.txt
python3 check_aliquot_fibers.py
python3 independent_prefix_audit.py
python3 -O check_aliquot_fibers.py
```

The DFS compares all complete fibers t=2..300 against an independent divisor sieve through 299², using the inherited global bound. Million-target sources are individually re-encoded by independent factorized sigma. The second implementation uses a full prefix divisor sieve and monotone polynomial binary search, not DFS, and compares all three full large lists. Checks raise explicit errors rather than removable assertions. The prefix-audit command generates `independent_audit.json` as a local verification record; it is not a required input or a published artifact.

Boundaries: complete enumeration for precisely these three targets; no full classification of odd fibers or prime untouchables; no proof of Goldbach or prime-gap distribution; no claim that all sources are semiprimes; no forward-chain termination theorem. Generated JSON contains no timestamps or timing fields; normal and `-O` runs produce identical bytes.

## Additional independent audit

A separate implementation using only integer arithmetic, a divisor/prime sieve, and monotone polynomial inversion reproduced every source in all three requested fibers. It separately verified all factorizations and every target 2 through 300 against a direct proper-divisor sieve over the theorem-complete source universe. This audit did not import either public inversion implementation. Its scope is the stated exact fibers, not all million-scale odd targets.
