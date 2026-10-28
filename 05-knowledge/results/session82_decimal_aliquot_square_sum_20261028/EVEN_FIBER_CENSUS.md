# Complete even-target fiber census from two source families

Status: **PROVED** elementary completeness; **FINITE-EXACT** census; **CITED / REPRODUCED**, not a new literature algorithm. Companion to [ALIQUOT_FIBERS.md](ALIQUOT_FIBERS.md).

## Exact source universe

For any integer V>=2, all sources of even targets t in 2..V belong to the disjoint union

1. even integers n in 2..2(V−1), and
2. odd squares m² with 3<=m<=V−1.

Proof. If n is even and s(n)>=2, then n>=4, and1,n/2 are distinct proper divisors. Thus `t>=1+n/2`. If n is odd and t is even, sigma(n)=n+t is odd. Factoring sigma shows every prime exponent of the odd n is even, so n=m². Since n=1 has target0, m>=3. The proper divisors 1,m give `t>=1+m`. This proves the bounds and parity split. Conversely every enumerated source is retained exactly when its actual s-value is an even target in range, so the resulting grouped lists are full inverse fibers, not truncated source samples.

`even_range_census.py` computes ordinary sigma by a divisor sieve through 2(V−1), computes smallest-prime factors only throughV−1, and evaluates sigma(m²) multiplicatively from m. It never sieves toV². These loops use O(V log V) elementary sieve additions and O(V) stored indices, in the usual unit-cost integer-operation model; no bit-complexity claim is needed. The program returns the actual grouped source lists; the generated JSON retains aggregate statistics rather than a bulky half-million-list dump. Generated JSON is not included in the published package.

## Exact result at V=1,001,000

* 500,500 even target values;
* 150,390 empty fibers;
* 642,193 total sources across all those fibers;
* maximum fiber cardinality 18, attained exactly at targets 560554,640634,895994;
* target 1001000 has exactly `{718720,1001998998001}`.

Running the commands below generates the complete cardinality histogram and counts at each 100000 threshold in `even_range_results.json` beside the scripts. The even-empty count through 1000000 is 150231. Its sum with the one known odd untouchable 5 equals the published total 150232; this comparison is **not** an independent verification that no additional odd untouchables occur in that interval, since this program enumerates only even targets.

## Recovery from primary literature

[Pomerance–Yang 2012, §5](https://math.dartmouth.edu/~carlp/uupaper3.pdf) already uses the even-source/odd-square split for near-linear untouchability enumeration. It also gives a sharper bound for composite odd roots: m<V^(2/3). Indeed choose a proper divisor g of composite m with sqrt(m)<=g<m; then gm is a proper divisor of m², so `s(m²)>gm>=m^(3/2)`. Prime roots need separate handling via `s(p²)=p+1`. Our simpler implementation keeps all odd roots below V and retains full fibers, not only touchability flags. The published table's ten total counts at 100000,200000,...,1000000 equal our even-empty counts plus 1 in every case. Source accessed 2026-10-28. No new density assertion follows.

## Verification

Requires the adjacent `check_aliquot_fibers.py` and SymPy as pinned in `requirements.txt`; run the target checker first to generate its input `results.json`.

```
python3 -m pip install -r requirements.txt
python3 check_aliquot_fibers.py
python3 even_range_census.py
python3 -O even_range_census.py
```

Both use explicit exceptions. An independent full source sieve checks every even target through 120 under the inherited universal source bound. The million-scale target 1001000 agrees with both the factor-prefix DFS and separate polynomial inversion. Additional exact DFS comparisons check 52,124,560554,640634,895994,1000000, including all maximal-cardinality targets. The first 30 empty even targets match the published initial sequence (excluding odd 5); ten aggregate counts match the primary table as described above.

Stopping boundary: this finite census closes all even fibers only throughV=1001000. It does not enumerate all odd fibers, prove any global untouchability conjecture, or produce an asymptotic density. The giant square source demonstrates why an unqualified linear source cutoff would be false; it is the two-family proof, not the size of the computation, that certifies completeness.

## Sequence identification and additional audit

The inverse-fiber count for targets t≥2 is [OEIS A048138](https://oeis.org/A048138); its offset excludes target1, whose fiber is all primes and is infinite. Zero entries are untouchables; singleton fibers are listed by [A057709](https://oeis.org/A057709). We use the explicit singleton-fiber definition rather than relying on historical variations in the word hermit.

An additional independent multiplicative-recurrence implementation reproduced every one of the500,500 fiber lists, not just the aggregate histogram. It checked the complete small-source oracle through120, all boundary cutoffs V=2..150, and the three maximal fibers by separate polynomial inversion. No infinite-density conclusion follows from this finite result.
