# Three certified aliquot orbits: square sums, primitive triples, and a frontier

**Status: FINITE-EXACT** for the traces; **PROVED** for the elementary bridge
identities. Write `s(n)=sigma(n)-n` for positive integers. Zero is a terminal
marker, not a positive vertex or an evaluated `s(0)`.

## 1. Complete short orbits and the long certified prefix

```
1000001 → 10003 → 1437 → 483 → 285 → 195 → 141 → 51 → 21 → 11 → 1 → 0
1001 → 343 → 57 → 23 → 1 → 0
```

These contain 11 and 5 evaluated edges, including `s(1)=0`. Their first common
positive vertex is 1. The first explicitly joins the `21→11→1→0` tail;
the second joins `57→23→1→0`.

The third orbit begins

```
1001000 → 2143960 → 4307240 → 6769240 → 8548760 → 12435640.
```

There are **1144 certified edges**, with these checked endpoints (seed index 0):

```
n_1107 = 34824938264087318322868834739622267021533597146401245964
n_1137 = 120743803148050878826315520222770879154553267874534033029280572
n_1144 = 3102001427119207440325481174042778229464844199474079305061910724
```

All evaluated sources in this finite prefix are even. None repeats, and no
positive consecutive pair has square sum. This does **not** establish eventual
termination, divergence, permanent evenness, or absence of future cycles.
The last term is not evaluated.

## 2. A genuine square-sum bridge and a normalization firewall

For prime `p`, `s(p^k)=1+p+...+p^(k−1)`. In the second orbit,

```
s(7³)=1+7+49=57,       343+57=400=20².
```

Thus `343→57` is both an aliquot edge and a square-sum adjacency between
distinct positive vertices. It is the only such consecutive edge in the two
short orbits. The identity `1+0=1` does not enter the positive square-sum graph.
For prime cubes generally, the test is whether `(p+1)(p²+1)` is square;
the complete prime-power classification is given separately in [SQUARE_SIGMA_EDGES.md](SQUARE_SIGMA_EDGES.md). The prime hypothesis
matters: `s(4³)=s(64)=63`, not `1+4+16=21`.

For `n>1`, divide `n+s(n)=sigma(n)` by `g=gcd(n,s(n))`. This yields a positive,
pairwise-coprime triple `(A,B,C)=(n/g,s(n)/g,sigma(n)/g)` with `A+B=C`:

| Source | Primitive triple | Removed gcd |
|---:|---|---:|
|1000001|(1000001,10003,1010004)|1|
|1001|(143,49,192)|7|
|1001000|(275,589,864)|3640|

The original edge has square sum iff **`gC` is square**, equivalently iff `g`
and `C` have the same squarefree kernel. A square `C` alone is insufficient:

- `483→285` gives `161+95=256`, but the original sum is `768=3·256`.
- `141→51` gives `47+17=64`, but the original sum is `192=3·64`.

For `343→57`, `g=1`, so both tests agree. Keeping `g` reconstructs the original
source and target. This aliquot normalization is not the separate binary-power
construction `(u−2^a)+2^a=u`; it does not identify their dynamics or prove ABC.

## 3. A source-preserving valuation rule—with an observed escape

If `v2(n)=2` and `7|n`, then `7|sigma(n)` because `sigma(4)=7`, hence
`7|s(n)`. Moreover,

```
v2(s(n))=2  iff  8|sigma(n).
```

Indeed, writing `n=4m` with `m` odd reduces the criterion to subtraction
modulo 8. Odd `b=v7(n)` is sufficient for the divisibility, since
`sigma(7^b)=(1+7)(1+7²+...+7^(b−1))`. This is a **one-step** condition,
not a permanence theorem.

The third orbit retains `v2=2` and divisibility by 7 at every source index
709–824, but then escapes:

```
n_824 = 154723849679327693492 = 4·7²·13·60723645871007729
n_825 = 184478436156121486288 = 16·7²·13·18100317519242689.
```

At index 824, `v2(sigma(n))=2`; subtraction raises the successor's valuation
to 4. A later run covers indices 852–1143 (292 certified sources), with no
claim about its future. The verifier checks both runs and the escape.

## 4. Exact factor-split evidence and the final resource frontier

Two initially unresolved cofactors were split by bounded, seeded runs of the
already-installed SymPy 1.14 ECM implementation
([official documentation](https://docs.sympy.org/latest/modules/ntheory.html#sympy.ntheory.ecm.ecm)):

```
83940594934600503096994848435730837699779204259589
 = 5511013429899318414019 · 15231426307036603838278075031

107583730854949977748239831120755574285637006670563966403
 = 7612737786472819497221147 · 14132068366536545816669370486649.
```

All four factors have complete recursive primality certificates. Their exact
products are independently verified. The first successful search used
`B1=200000, B2=2000000, max_curve=40, seed=82003` (30.06 seconds), after
smaller unsuccessful stages. The second used
`B1=1000000, B2=10000000, max_curve=80, seed=82005`, succeeding on curve 15
(101.47 seconds), after all 80 curves at bounds `200000/2000000`, seed 82004,
finished without a factor (99.95 seconds). Times are observed, not guarantees.

The session stops at

```
n_1144 = 2²·7·11·R,
R = 10071433204932491689368445370268760485275468180110647094356853.
```

This 62-digit `R` is composite: the verifier checks that `2^(R−1) mod R`
equals `1059352591669555434863290548381996608006601416003032122563405`,
not 1. Its unresolved factorization prevents a certified next divisor sum.
A single factorization exhausted its 60-second budget; no further ECM
escalation was performed. A resource limit is not mathematical impossibility.

## 5. Regenerate and independently verify—without bundled data

The two companion scripts generate their JSON locally; no preexisting certificate bundle is required. The generator requires Python, installed SymPy 1.14.0, and POSIX
`SIGALRM` support; native Windows without that support is rejected. The
verifier requires only Python's standard library and exact integer arithmetic.

```
python3 generate_aliquot_orbits.py --output /tmp/aliquot-session82.json
python3 verify_aliquot_orbits.py /tmp/aliquot-session82.json > /tmp/aliquot-check.json
python3 -O verify_aliquot_orbits.py /tmp/aliquot-session82.json > /tmp/aliquot-check-O.json
cmp /tmp/aliquot-check.json /tmp/aliquot-check-O.json
```

Generation takes several minutes. Each factorization has a 60-second wall
budget; `--seconds` must be positive. A timeout raises an exception rather
than emitting a completed packet. `--seconds 120` is an explicit optional
budget change, never automatic escalation. The four exact factors above are
hints used only when division succeeds; they are separately prime-certified.
No repetition of the two explicit ECM discovery searches or external factoring service is required. Ordinary SymPy factorint may still invoke ECM internally.

These deliberate corruption controls must each exit nonzero; repeat with
`python3 -O` to check optimized mode:

```
python3 verify_aliquot_orbits.py /tmp/aliquot-session82.json --negative-control edge
python3 verify_aliquot_orbits.py /tmp/aliquot-session82.json --negative-control prime
```

The first changes a successor; the second changes a primality witness to 1.
All validation gates use explicit exceptions, not removable assertions.

### Why the prime certificates are conclusive

For each `p>2`, recursively factor `p−1` into certified primes and supply `a`
with `a^(p−1)=1 mod p` and, for every prime `q|p−1`,
`gcd(a^((p−1)/q)−1,p)=1`. For any prime divisor `r` of `p`, these conditions
force the multiplicative order of `a mod r` to contain every full prime-power
part of `p−1`. Thus `p−1` divides `r−1`, forcing `r≥p`, hence `p=r` is prime.
The recursion strictly decreases to 2. SymPy's probable-prime screening is
never accepted in place of this independently verified Lucas certificate.

### Fresh replay verification

A fresh from-seed regeneration reproduced every row of all three frozen traces,
including all 1144 C edges. The combined packet verifies **5806 distinct
primes**. Normal and optimized verifier outputs agree byte-for-byte; successor
and prime-witness corruption controls fail in both modes. Zero and negative
`--seconds` arguments are rejected before computation.

## Independent audit

A separate topological prime-certificate traversal and Horner divisor-sum implementation verified all1160 combined edges and all5806 primes, both ECM splits, the exact unevaluated frontier, unchanged earlier prefixes, source/endpoint novelty, and the valuation escape. It also checked the fresh from-seed generated bundle against all frozen rows. Its verification did not rely on probable-prime screening.
