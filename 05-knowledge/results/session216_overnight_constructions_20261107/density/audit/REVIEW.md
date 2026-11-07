# Independent audit: Session216 logarithmic refinement

Verdict: PASS WITH ONE LOCAL CORRECTION. Both headline inequalities are valid as stated. The optional optimization-barrier proof needs two related algebraic corrections, detailed below. No Git publication performed; parent owns release.

## Frozen dependency and scope

Audited `session216-density/RESULT.md`, `check.py`, `check.json` against repository HEAD `0f3dad8c75d02db9e4060746280fb58f4d69633e`.

The dependency is `01-canon/theorems/THM-4887-one-apex-energy-at-scale-two-to-the-n-over-three.md`, statement (D1), specifically its Delta=0 constants. Its quantifiers are m>=1, C=2^m, P=2^(k-1)>=8C, n=k+m, and 0<=Delta<=C. It expressly applies at every admissible split, not only the m>=floor(n/3) range used in its other corollaries. The new choice is therefore licensed. The range and tournament realization are inherited from that file, and the symmetric finite-check identity is inherited from `01-canon/theorems/THM-4890-one-apex-collision-balance.md` (A).

Targeted correction search found MISTAKE-730 about random-model scope; the present proof does not use that model. This is an optimization of a proved bound, not a new collision estimate or power exponent. Positive density, linear energy, density one, and Chan–Pak remain open. No external literature-priority assertion was checked or is required.

## Algebra, quantifiers, rounding

The split m=floor((n-log2(n))/3)+1 gives 1<c<=2 exactly, including equality c=2 when the inner quantity is integral. The integer implementation agrees. For n>=12, m>=1 and B=N/(2C^2)>=2^(n/3+(2/3)log2(n)-3)>8. Thus P, C and B are appropriate integer powers of two and all D1 hypotheses hold.

The four normalized terms printed in the main proof expand correctly. Convexity bounds the first line by the larger endpoint 2+2 ln2. The second-line replacement drops a nonpositive term and uses c/2<=1 against a positive bracket; its derivative is negative for n>=12. The bound lnB<=n ln2/2 follows from log2(n)<=n/4+3/2, whose derivative check is valid starting at 12. Both summands in the resulting last-line bound decrease there.

Cauchy–Schwarz uses exactly |X_n|=N/2, explaining the factor 29.2=4(7.3). Every represented value is <=N since the largest word is N-1. For real X>=4096, n=floor(log2 X)>=12, N<=X<2N, and n<=ln X/ln2. Both replacements have the correct lower-bound direction, including X a power of two. No floor or boundary loss is missing.

## Elementary exact constant certificate

One may replace all informal decimal evaluations by

    0.6931 < ln2 < 0.6932,
    ln12 < 2.485,
    12^(2/3) > 5.24.

Then the normalized energy is less than

    (2+2*0.6932)
      +(33.74-12*0.6931+8*2.485)/12
      +(2+9*0.6932)/(16*5.24)
    =114125651/15720000 <7.3.

The spectrum constant is certified by the entirely rational inequality

    29.2^3*4 <60^3*(0.6931)^2.

The logarithm bounds are rigorously reproducible using the positive series
ln x=2 sum_(j>=0) z^(2j+1)/(2j+1), z=(x-1)/(x+1), at x=2,3 and ln12=ln3+2ln2. After 30 terms, the omitted tail is at most 2 z^61/[61(1-z^2)]. `audit.py` proves the rational comparisons using `Fraction`, not floating-point assertions. The root bound follows from 5.24^3<144.

## Required local correction: optimization barrier

The paragraph retaining the logarithmic term currently loses one power of C. The correct formulas are

    12PC^2 ln(P/C)=6NC ln(N/(2C^2)),
    d/dC [C ln(N/(2C^2))]=ln B-2.

The derivative remains positive since B>=8 and ln8>2 (certified already by ln2>0.6931). The evaluation at C0=(N/n)^(1/3) in the draft already uses the correct squared denominator. Also C0 lies below the admissible upper endpoint sqrt(N)/4 for n>=12. Hence the remaining proof works unchanged: C<=C0 gives 8P^2>=2Q, while C>=C0 gives at least (3 ln2/2)Q. This is only a floor on the displayed numerical D1 right-hand side, never on actual energy. With the algebra repaired, its narrow claim about optimizing this same bound is justified.

## Independent finite controls

Run from the common workspace root:

    python3 session216-density-audit/audit.py > session216-density-audit/audit.json

Universe: every word in [2^(n-1),2^n), n=12..18, with the optimized split at each n. An independent high-to-low bit-append evaluation replaces the author's memoized recursion table; sorting into runs replaces the histogram-based support/energy computation. All exact energies and supports match the supplied JSON. At n=12, all 4,194,304 ordered pairs are literally compared, independently confirming E_12=6000. The symmetric identity and range are checked word by word throughout. A deliberate omission of the parity term fails on 512,1024,...,32768 words at n=12..18, so the identity test detects its intended hostile mutation. These finite checks are implementation controls, not the all-n proof.

## Administrative boundary

Source checkout contained the parent's unrelated Session216 edits; none were modified. Startup packet was read. Targeted processor status could not run because `.machine-id` is absent; no registration was fabricated. No new theorem namespace was reserved, and no source checkout changes, commits, or pushes were made by this auditor.
