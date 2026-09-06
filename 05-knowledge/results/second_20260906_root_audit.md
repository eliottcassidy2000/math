# Independent audit of the complete mixed Smith law and balanced decoder obstruction

**Verdict: PASS in the stated scopes.** The root reviewer accepts the analytic
proofs and exact witnesses in [the mixed Smith classification](second_20260906_smith.md)
and [the balanced decoder obstruction](second_20260906_decoder.md).
This audit imports neither producer. It does not claim a general higher-jet
classification or a solution of LRC(14).

## Complete mixed (2,2,1) observer

Clearing the value and derivative at the first doubled affine node gives
two unit factors and the stated residual three-by-three matrix. Its first
determinantal ideal has valuation `min(2A,A+v_p(2),2B)`. Its determinant
has valuation `4A+2B+2C`. The reciprocal Hermite formulas give the largest
factor, including both potentially canceling numerators `a+2b` and `3a-2b`.
Because exactly three factors remain, these three data determine the entire
partition. This inference would fail for a larger residual observer.

The case split uses the ultrametric constraint that the minimum of A,B,C
occurs at least twice. When the doubled pair is closest, cancellation can
alter the partition only in the dyadic adjacent-depth case. Its active range
is `e>=2`; at `e=0,1` the other entries fix the largest factor before that
unit bit can matter. The two doubled-node swaps preserve the bit because
`tau -> tau/(2tau-1)` fixes every odd residue modulo four. Changing the
unit-separated reference changes tau by a multiple of `2^e`, using the
Pluecker identity; this is enough exactly at the stated active depth.
Rescaling primitive directions and applying a unimodular projective map
cancel in the bracket ratio.

The projective transfer uses complete local jets and an invertible local
parameter change from the [audited projective result](creative_20260906_smith_bridge.md).
A common integral affine chart exists unless all three directions occupy
the three classes of `P^1(F_2)`. In that exception every bracket is a unit;
the determinant is a unit, so every Smith factor is a unit. Thus neither
infinity nor the dyadic chart exception is silently discarded.

The [independent program](../../04-computation/second_20260906_root_audit.py)
builds full five-by-five ordinary monomial value/derivative matrices, and
uses SymPy's integer Smith normal form. This differs from the producer's
rational DVR pivot algorithm and residual-minor calculations. The universe
contains all distinct nonzero signed affine nodes a,b in `[-12,12]`, plus
labelled closest-pair arrangements at primes 2,3,5,7, outer depths 0,...,5,
inner increments 1,...,4, and both unit signs: **1,059 matrices**, with
**4,236 full prime partitions**. Every determinant is also checked directly.
Both inherited dyadic twins occur in the bank. The full capped-exponent
kernel ratio is independently checked at depths 2,...,15 and all integer
precisions 1,...,4e+3. It is exactly two on `3e+2<=N<=4e`, and one elsewhere.
The source proof, rather than the finite bank, supplies arbitrary primes,
depths, units and projective charts.

## Actual balanced decoder obstruction

The literal five-denominator construction is primitive: each denominator
is missing from one leaf and the denominators are pairwise coprime. Its
five star edges have primitive coefficient sum 356. The seven-component
links have allowed inert sums 50,110,332. The physical scales retain the
box and coprimality, so these are actual primitive thirteen-speed entries.

For each of the two scales, the independent program reconstructs every
pair edge from prime factors of its primitive coefficient sum, and checks
the resulting relation matrix has rank eleven with no cross edge. It
excludes all two-small/one-large mixed relations by their cleared
distinguished coefficient being greater than Q; every one-small/two-large
relation is excluded by the stronger dominance inequality with the two
largest U labels. Thus the full bounded relation span, not just selected
pair rows, equals the decoder span.

Every one of the 1,716 seven-subsets of each row has gcd one. This independently
checks **3,432 subsets**; every larger subset contains one of them and hence
also has gcd one. The producer additionally checks all retained joint-shadow
profiles and complement words. This is why the failed earlier 355-chain
cannot substitute for the repaired 331-chain.

Both rows fail the balanced minimum comparison and both simple phase-grid
comparisons. For every one of the six maximum-endpoint pairs, the new native
relation gate is paid (`delta=1,c<=Q`) but its phase inequality fails.
The proof obstruction therefore persists after trying all those pairs.
The program independently evaluates the full-row safe phases: `1/2` for
the odd row and `11/23` with clearance `3/23` for the mixed row. It also
checks every rational numerator at denominators 2,...,22 and all odd
sixteenth phases in the mixed row. These are failures of specified
sufficient certificates, not counterexamples to safety.

The incoming cross-divisor criterion in `8e560f2142` asks for
`lcm(D,v)<=3Q/28` at this split. Since every smaller label is already above
that threshold, `lcm(D,v)>=v` rules out its score for every choice, including
shared-factor cancellation. The independently verified actual box, full
entry and subset profiles therefore refute the proposed implication that
those data always force a qualifying score. This additional consequence
uses no unverified search or assumption about which pair is selected.

## Reproduction

```bash
python3 -B 04-computation/second_20260906_root_audit.py
python3 -B -O 04-computation/second_20260906_root_audit.py
```

Both runs produce the same [frozen output](second_20260906_root_audit.out).
The only nonstandard dependency is SymPy (version 1.9 in this session).
Every check uses an explicit exception, so optimization does not disable
the audit. Raw source and output hashes are recorded in the
[session manifest](second_20260906_manifest.json).
