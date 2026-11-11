# Summand growth: exact closure, shrinking frontier, and binary inverse

Date: 2026-11-11. **PROVED** statements are in the two proof files;
**FINITE-EXACT** controls are separately scoped in the receipt. No claim
of literature priority or a completed tournament-spectrum proof is made.

[CLOSURE.md](CLOSURE.md) supplies four new operations on the strict
paired-parent graph:

1. For any finite positive seed set S, its strict closure is exactly
   `S union union_(a<b in S)(a+b+H(S))`. The same mixed-seed ideal theorem
   holds in positively weighted commutative monoids, including strict
   multiplication of seeds at least 2.
2. A retained seed m gives a fixed finite state of 2m residue minima,
   plus seed exceptions, with exact updates when any new seed is adjoined.
   The state decides all future integer targets, not only a truncation.
3. Every finite seed set of at least two values has normalized synchronous
   depth `log2 n+O_S(1)` on large attainable targets.
4. Every such set eventually has top three values spaced by d=gcd(S),
   after which `M_(t+1)=2M_t−d`. A consecutive three-block catches the
   maximum through an exact shrinking deficit-window compiler. The proof
   gives a finite computable transition bound, not its minimum onset.

[TOP_SUFFIX.md](TOP_SUFFIX.md) resolves the remaining lower-hole boundary:
eventually `S_t=C(S) intersection [1,M_t]` at every later generation. It
also gives exact eventual cardinality and target-arrival formulas. The
sharp minimum onset remains open; its controls and hashes are separate.

[INVERSE.md](INVERSE.md) retains the full labelled shadow-incidence
observer rather than only its birth degree. Its inverse follows binary
ancestors, has every coefficient in `{−1,0,1}`, and has an exact logarithmic
row size. The multiplicand inverse is a synchronized tensor of the same
binary kernels across prime exponents: synchronizing the correction scale
before summing is essential. The witness `6→36` gives coefficient 2 even
though each separate exponent-coordinate inverse gives 0.

## Recovered connections and the retained sidecars

| Source | Operational transfer | Preserved | Information that must remain attached |
|---|---|---|---|
| May 3 strict summand cospans; S259 translated two-seed semigroup | Mixed-support ideal for arbitrary seeds | Exact attainable values | Original seed types; ordinary monoid alone permits illegal doubling |
| May 31 natural sum/product graphs | Positively weighted commutative-monoid closure | Strict two-parent generation | Actual monoid element and expression, not only its real logarithm |
| S259 long parity transient | Consecutive-block chase in moving deficits | Exact upper window and eventual maximum | Initial block and its window; lower holes are discarded |
| July involution-parity reflection; September doubling complement | Binary correction inverse | Every input weight | Full labelled observation sequence, not fibre counts |
| Multiplicand square-deletion inverse | Common-scale tensor on prime exponents | Full inverse coefficient | One synchronized scale j; multiplying separately summed inverses fails |

The inheritance, source paths, quantifiers, failure boundaries, proofs,
and exact connection contracts are in the proof files. In particular,
ordinary reachability, paired-parent generation, and weighted incidence
inversion remain distinct predicates.

## Reproduction

From this directory:

```sh
python3 -B check.py
python3 -O -B check.py
```

The receipt checks 65,280 closure memberships, 6,606 finite seed updates,
4,097 moving-window updates, 1,023 persistent top-three controls, four
complete frontier chases, four multiplicative seed sets through target
1,000, and direct integer-matrix inversion through order 128. Explicit
exceptions remain active under `-O`. `--write` regenerates the receipt.
