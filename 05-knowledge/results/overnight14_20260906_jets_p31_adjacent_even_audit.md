# Independent audit of the adjacent even ideals and actual p31 factor pair

**Status: PASS, independent analytic audit; FINITE-EXACT verification.**
Root read the complete [producer proof](overnight14_20260906_jets_p31_adjacent_even.md)
after deriving the adjacent-ideal mechanism independently of its code. The
universal even-minor formula, both complete determinantal ideals, the two
ordered factors and their precision window are accepted without repair.
The note may be promoted to PROVED + INDEPENDENTLY AUDITED. A complete
all-parameter p31 partition is not claimed.

## 1. Exact scope and proof review

The object is the full 48-column Hasse observer at three nodes
0,31^e,31^e a, with derivative orders0 through15 at each and a,a-1 units.
Clearing the node-zero identity block is integral and invertible. It leaves
16 unit factors and a32-by32 residual matrix. For positive e, residual
ideal E_r is precisely the full ideal of order16+r, and adjacent differences
are actual ordered Smith exponents. At e0 every factor is a unit, a separate
boundary case. Selected-minor values alone would not suffice.

For rank2s at general multiplicity m, the highest possible row-order sum
uses both copies of m-s through m-1, uniquely; the lowest possible column
sum uses m through m+2s-1, uniquely. Their difference is3s². Lowering all
selected degrees and derivative orders by h=m-s gives the complete residual
multiplicity-s bank. The falling-factorial identity factors by columns and
rows and retains the exponent c-j. Its determinant is therefore exactly
the stated rational scalar times [a(a-1)]^(s²), with the positive sign for
the declared node/derivative order. Gauss's lemma makes the scalar integral.
This is an identity of determinants; it does not make these row scalings
unimodular or preserve an entire Smith form.

At rank28, the scalar is43!42!/(15!³14!³), of31-valuation2. The monomial
is a unit, so the unique minimum attains588e+2 at every lift. Every other
minor has a nonnegative integral column excess plus row deficit. At total
excess one there are four row replacements and one column replacement.
Every retained derivative order is at least one and column31 is present.
That column is identically zero over F31[a], because binom(31,j) is divisible
by31 for1<=j<=15. Thus all five normalized polynomials have31-content at
least one. They have valuations at least589e+1; higher bands have at least
590e. Both meet the required lower bound for e>=1, including e1. The
minimum supplies attainment, proving the entire E28=588e+2 ideal.

At rank30 the unique minimum has weight675 and scalar45!/15!³ of valuation
one. The next possible weight is676, already sufficient for all e>=1.
No rank28 zero-column claim is imported into row patterns that now permit
j0. The unique minimum attains675e+1, proving the whole E30 ideal.

Subtracting the already independently audited E29=631e+1+kappa proves
lambda29=43e-1+kappa and lambda30=44e-kappa. Their gap e+1-2kappa is
nonnegative and vanishes exactly at e1,kappa1. They occupy positions45,46
of the full48-factor list. Both residue classes occur, for example a4/a3.
Neither a generic residue fit nor a finite-depth extrapolation supplies
these conclusions.

For the pair's kernel exponent, cap each of its two valuations at the
integer precision b. Moving one valuation unit from the larger factor to
the smaller changes the sum by one exactly for43e<=b<=44e-1. The interval
has e integer precisions. Outside it, the sum is unchanged. The total pair
valuation87e-1 is fixed. This statement concerns the identified pair;
other factors might change between arbitrary inputs in the two classes.

## 2. Independent exact path and universe

The [independent program](../../04-computation/overnight14_20260906_jets_p31_adjacent_even_audit.py)
imports no mathematical producer. It checks45 signed even minors by rational
Gaussian determinants at s1..5, m=s,s+2,s+5, and a2,-2,4. A falling-factorial
product separately reconstructs the shift scalar and both factorial forms.
Every row and column complement at ranks28 and30 is enumerated separately,
verifying the unique minima, all five rank28 next-band cases, retained
column31, and its polynomial divisibility. No enumeration of the Cartesian
product of all minor choices is claimed.

For full matrix controls the referee materializes the residual bank directly.
It repeatedly eliminates arbitrary unit pivots and divides common31-adic
layers. This differs from the producer's full48-row, globally least-valued
pivot algorithm. All operations preserve the Smith exponents; visibility is
checked at finite precision48e+5. The resulting32 exponents have sum768e,
and their first28,29,30 sums and positions29,30 agree with all formulas.

The complete finite universe has every29 admissible residue at depths0,1,2,3,
with shifted higher-digit representatives, plus e5,a879; e4,a-20;
e7,a=3+31^4; and e6,a=4-31^3:120 matrices. These include the individual-
packet cancellation hostile, negative lifts, and depths beyond the producer
bank. Every integer precision through44e+2 for e1..50 checks the exact pair
kernel window, the unchanged sum, and the equality boundary.

Both modes pass132,821 always-active gates with byte-identical LF output.
The finite controls corroborate the universal proof; they do not extend it
to other primes or to the unclassified remaining factors.

```
python -B 04-computation/overnight14_20260906_jets_p31_adjacent_even_audit.py
python -B -O 04-computation/overnight14_20260906_jets_p31_adjacent_even_audit.py
```

The source and transcript are frozen with this report in the checkpoint
manifest. The reviewed producer source is SHA256
`e0eb5d1b403187706ccb9902747a58060e157885d5e8cfbefea12d757eb67bfb`;
its LF output is `9c7cd02d21bfea5ffdb80dac65cf572fd9904909b18e35a2d7a1d6ce575e1288`.

**Filing:** root read these proofs and audits and integrated the fourteenth
checkpoint. Reproduction commands are relative to the repository root.
