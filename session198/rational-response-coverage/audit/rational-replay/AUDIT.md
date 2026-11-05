# Independent audit: exact rational signed-response coverage

Verdict: PASS, including the arithmetic g(n) formulation and C2 offset corollary.

- PROOF.md SHA256: b33f84878bb69367993f8043c75aff0c4fc8d6dc323399ef87127ba069cc2cb3
- verify.py SHA256: f9e2ee15265704469505b8b15953a12fd5cb1b0b34c7f91fd19e9507d349d2e8
- Both frozen receipts SHA256: 5a574f62381b977581728a64bea05eff937c3c1e9da884c32c9bc8e08615efa2

The periodic-word construction works for every stated coprime even-a/odd-b
slope. A is even and lies between B/2 and B-1, so the chronological word starts
with O and ends with Z. This removes the initial E boundary term and ensures a
strong final source. The exact raw transfer and K=aR-bC are correct.

The positivity lemma is an all-order proof, not a finite check. The formal
cover orbit has repeating digits in the required chronological order, stays
strictly between0 and1, and returns to its initial value. Since the first bit
is zero, the cyclic ending-count convention adds no missing boundary term.
For each maximal run of ones, backward iteration from its following zero gives
the displayed sum of deficits, and hence a strictly positive run contribution.
Therefore K/b>0. K is integral by its formula and even by the actual unit-
deletion source parity check. No formal periodic state is claimed to be physical.

The modular-to-exact step is sound. Coprimality of2a and b provides odd seed
residues for every odd target modulo K, including K=2. The inherited interface
compiler supplies those residues with unbounded actual D and exact s=E=d=1,
and retains transitive deletion. A single O makes the response positive and
arbitrarily large because2a-b>0. Its difference from the target is a positive
multiple of K; the stated number of genuine word repetitions reaches the target
exactly and finishes strong. The cone bounds prove necessity outside the open
slope interval. No efficient-size claim is made.

The arithmetic formulation is exact: starting at cover n=2, the maps n->2n
and n->2n-1 give profile (g(n),n mod2,n,1), and inverse reduction reaches every
n>=2. Odd n>=3 are strong. The modular seed also obeys this formula by the g
recurrence. The3/4 hostile control is valid: for even m its response is1 mod4,
and for odd m>=3 the strong cone makes it at most-3, excluding target3.

The C2 identity gives every odd offset with a variable physical parent. Fixed
offset sources can have unbounded D, so the parent counts are unbounded.
Eliminating U proves the explicit parent restriction n=19-2delta modulo34.
Thus this corollary does not supply either child for an arbitrary given parent.
The two literal controls157->313 and187->375 are correct C2/C2 pairs.

Final isolated normal/optimized replays reproduce the frozen receipts. Every
proof/helper dependency is unchanged from previously reviewed packets.
Additional controls in ../rational-independent.py and
../RATIONAL_INDEPENDENT.json check:
- 528 admissible slopes through denominator101, using C=A-g(A+1)+1 rather
  than the packet's bit-pair summation, with18 arithmetic inputs per slope;
- all550 compressed seed recipes against the original right-to-left g
  recurrence, rather than the constructor's left-to-right bit evaluator;
- 60 fully materialized exact arithmetic n witnesses via the closed cover
  recurrence and direct g(n), rather than repeated four-coordinate transfers;
- 17 literal strong unit-deletion source graphs recovered from their cover
  integers, plus independent graph assembly/counting of both C2 offset pairs;
- 9,999 odd inputs for the3/4 hostile control.

All controls pass. The proof establishes exact signed-response surjectivity
under its precise hypotheses. It does not establish positive Hamilton-count,
full C2/C3, cofinite, or universal tournament-spectrum coverage.
