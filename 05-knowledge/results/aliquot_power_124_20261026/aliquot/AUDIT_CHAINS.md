# Independent audit of the doubling, power and abc interfaces

2026-10-26. **PASS.** No substantive correction required. This is a proof
and source audit, with computational replay and an additional direct census.
The abc conclusion remains conditional.

Reviewed `../chains/RESULT.md`, SHA256
`2485fa87ed6ad33be1af70008048cfc2bbfd3fcd51dfd1ee1ccb5e0a29b1aa71`,
and `../chains/check.py`, SHA256
`63956a289445c8ade9386ed1453a9d390cf733534a5733432a4f40867809b540`.

The final pre-promotion proof SHA256 is
`ddbd208d5188c3eaa673dd855bcf6a37e453e01645c9e5f33abdf268e8b8b354`.
The bounded delta adds the explicit empty power/combined fibres at N=1
and the 268,651-gate receipt count; both were already checked above.
PASS applies to that final version as well. A status-only audit promotion
does not alter the mathematical scope of this report.

## Exact carriers and addresses

P1 counts all positive row bases in the doubling table, so its first row
is the odd core and its multiplicity is v2(N)+1. Restricting to odd roots
would change the table and remove those duplicates. P2 counts power
addresses by common divisors of the prime exponents. Exponent1 belongs to
the declared table; the separate nontrivial-power statement correctly
subtracts it. The first distinct-base duplicate with exponents at least2
is16; numbers below16 have no second such base.

In P3 with odd core U>1, k must divide every odd exponent and the
two-adic exponent of the base ranges from0 to floor(A/k). These conditions
are both necessary and sufficient. The pure-two case needs its own
formula because the odd exponent list is empty and base1 is forbidden:
the pairs (b,k) are exactly the positive lattice points bk≤A, giving
sum_(k=1)^A floor(A/k). For example, N=8 has five combined addresses and
N=16 has eight. N=1 has none in the power or combined table, as the checker
explicitly verifies. No infinite multiplicity is silently retained.

The semidirect operation law has the correct order: applying a kth power
after multiplication by2^j multiplies that exponent by k. It gives36 in
one order and18 in the other when n=3,j=1,k=2. The first-row distinction
at36 therefore reflects two different tables, not a contradictory count.

## Delayed rows and the square boundary

The delay criterion follows exactly from firstrow(n^k)=odd(n)^k and
n=2^a u. Its strict inequality is essential. The first even base is6,
but the unrestricted first base is3. The integer-log expression for the
first delayed exponent retains that strict inequality and is correct
also for odd n>1. Powers of two never enter the delayed set.

Formula(9) partitions bases by their exact two-adic valuation. Subtracting
the odd-cofactor count up to the exact integer root excludes precisely
u^(k−1)≤2^a. The nondelayed bound(10) splits at2^a of size
X^((k−1)/k); both geometric tails are O_k(X^(1/k)). This proves the
claimed density without extrapolating the finite scan.

At X=4^m, the low-valuation contribution is2^m−1, the intermediate
contribution is2^(m−1)−1 and the final pure power contributes1. Their sum
is3·2^(m−1)−1. This includes m=1, whose intermediate sum is empty.
An independent direct enumeration of every even n≤4^m for1≤m≤7 gives
2,5,11,23,47,95,191 nondelayed bases. Here the direct predicate was
odd(n)^2≤n, not the author's integer-root census formula.

## Perfect-number and abc maps

The square/odd-core square commutes, while the subsequent arrows use the
prime guard M=2^(a+1)−1. Since M>2^a, the row M² exceeds N=2^a M.
Its aliquot chain is M²→2^(a+1)→M→1→0. It does not enter the fixed point
N: s(M²)=2^(a+1) could equal N only if M=2, excluded here. The example
M=15 correctly breaks the prime-square divisor formula. The explicit
failure of odd-core projection to conjugate s keeps the lost coordinate
visible. I read [King's Theorem3 and proof](https://www.jonathanking.org/PDF/mersenne-perfect.pdf);
the Euclid–Euler scope used here is correct and agrees with Conrad's
independently checked Theorem1.1 in the sibling aliquot note.

For the primitive aliquot triple, division by gcd(n,s(n)) gives positive
pairwise-coprime entries for n≥2. Abundance order and perfectness survive,
but every perfect integer becomes(1,1,2); the scale recovers the original
edge. The prime-power triple is a different construction: its second term
is p^k−1, hence is coprime to p. The radical identity follows, and the
displayed conditional lower bound is the correct rearrangement of abc.
Its constant is uniform in p and k for each fixed epsilon, but is not a
supplied effective numerical cutoff. I independently read the stated
quantifiers in [Granville–Tucker, printed page1227](https://dms.umontreal.ca/~andrew/PDF/abc.pdf).

## Replay and scope

Normal and optimized executions of the author's checker pass268,651 gates, are byte-identical
and match its frozen receipt, SHA256
`4847fb13f19b3ce377478e4e77f6fd3c46f853efa07fc8438d5c70565bc6aa75`.
The reviewer copies are `chains.replay.normal.json` and
`chains.replay.optimized.json`. Run:

```sh
python3 -B ../chains/check.py > chains.replay.normal.json
python3 -B -O ../chains/check.py > chains.replay.optimized.json
cmp chains.replay.normal.json chains.replay.optimized.json
cmp chains.replay.normal.json ../chains/normal.json
```

The literal table enumeration is independent of the factor-exponent
formulas, and the divisor sieve is independent of multiplicative sigma.
The declared finite universes and invalid-input controls match the proof.
No count is used to prove global aliquot termination, odd-perfect-number
nonexistence, or abc. All such boundaries remain intact.
