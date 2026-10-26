# Independent audit of the growing twin sieve

**PASS, 2026-10-26.** Reviewed the all-parameter proofs in the frozen producer
`../sieve/README.md`, SHA-256
`5c825bdc81fedc3fd920012396ef5dda0acf2558d61ec0dde0b8e8d5a7f7cf06`.
No mathematical correction is required. This audits W1–W3, their equality
boundaries, and the stated finite experiment; it does not establish twin-prime
infinitude or an improved asymptotic sieve estimate.

## Proof review

- **W1, domains and lifts.** Prime 2 has one forbidden residue; an odd tested
  prime has two. A newly adjoined prime is coprime to the prior modulus, so
  the two displayed lift indices are distinct. Each allowed class has exactly
  ell−2 children. The integer D-depth changes by tM/2 because M is even.
  This is residue refinement, not a Berggren child operation.
- **Exact floors.** The lower-end subtraction at 2 counts precisely the
  integers from 3 through H, including negative floor arguments. The
  squarefree Möbius sum is an exact finite indicator of coprimality with
  the tested primorial. The root multiplicity is multiplicative and the
  prime-2 factor is one. Both stated error bounds follow by summing absolute
  residue errors; retaining the minimum of the two bounds is valid. The
  exponential term count is acknowledged rather than suppressed.
- **W2, prime self-hits.** An accepted non-strict point must have one entry
  equal to a tested odd prime ell. A composite other entry, at most ell+2,
  would have a proper prime factor at most sqrt(ell+2)<ell; this is tested
  and cannot use an equality exception. Thus both entries are prime and the
  lower one is at most B. Conversely those twin pairs pass every proper-factor
  test and fail strict avoidance. This proves the disjoint partition for
  every B>=2, not merely below a height threshold. An accepted even x>=3 is
  impossible because the prime-2 test fails without equality.
- **W3, strict boundary.** If H+2<q², any composite candidate entry has a
  proper prime factor below the next untested prime q, so it is removed.
  Equality is deliberately excluded. At x=q²−2 the upper entry is the
  untested square q²; the lower entry survives precisely when prime, since
  a composite lower entry has a factor below q. This proves the stated iff
  and explains why q²−2 need not attain the first false positive.
- **Descending recurrence.** With z=floor(sqrt(H+2)), the next prime is
  strictly larger than sqrt(H+2), and z<H for H>=3. Thus the patch is
  pi_2(z), and the formula descends to the declared zero base below 3.
  Counting by the lower twin member permits its upper member to exceed H;
  no endpoint is silently discarded.
- **Limit and height.** The coherent class −1 has signed representative −1,
  centre zero, and triple (−1,0,1). Its nonnegative representatives M_B−1
  grow unboundedly. Every fixed positive x>=3 eventually fails a strict
  wheel because some prime divisor of x is eventually tested. The restored
  intersection is exactly the actual twin-prime set, which does not imply
  that this intersection is infinite. The quantifiers are kept distinct.

The numerical hostile at B=61 is correctly typed: 4487=7*641 is removed;
4757=67*71 survives, and 4759 is prime. All smaller candidate x were inspected
independently at this stage. Among prime cutoffs through 61, this is the first
failure of the prime-lower / upper-prime-square pattern. No global pattern is
inferred.

## Independent computation

The producer's 380,043 explicit gates were rerun normally and with `-O`; both
outputs match its frozen normal/optimized receipts byte for byte. Producer
receipt SHA-256:
`feaf5d1c3496b1eb5975c1396e9c1d3e06aeaf0e9dab04437c0deb6eeeaaf6ce`.

Our `check_audit.py` imports no producer module. It uses literal trial division
for primality, tests divisibility directly, and enumerates forbidden roots
literally instead of inheriting the producer's CRT root construction. It checks
integer cutoffs B=2,...,101 (including nonprime cutoffs), all self-hit points
3<=x<=1000, the descending recurrence at every 3<=H<=10000, and full threshold
boxes within that height. It independently checks wheel/floor/Möbius formulas
through cutoff 11 and searches every candidate through 5000 for all prime
cutoffs through 61.

**490,523 independent gates pass**, with normal and optimized receipts identical:
`2e759ce60ef82017ec7365de70ff3b0af42a55315c9abe2e863184d384cbb90c`.
The triangle transport formulas are inherited and producer-tested; the
independent sieve audit does not reclassify the entire Berggren tree.

Reproduce from this directory:

```
python3 -B check_audit.py
python3 -B -O check_audit.py
```

## Source and promotion boundary

The auditor independently opened Adam J. Harper's author-hosted
[Elementary Methods notes, Section 1, Definition 1.1 and Proposition 1.2](https://warwick.ac.uk/fac/sci/maths/people/staff/harper/elementarymethodsnotes1.pdf)
and read the indicator/floor proof and remainder discussion on pages 1–3.
This supports the classical inclusion-exclusion attribution; the specific
two-form self-hit and height assertions were checked directly above.

After PASS, the owner promoted only the status/link text. The auditor
reversed those edits and recovered the exact reviewed hashes, confirming
that all mathematics and computation stayed unchanged. Final README SHA-256:
`8f6e97daa359cd1e24eb70f0be8565436213fa6facfa4120e9762b768785766c`;
final SOURCES SHA-256:
`f473f7bd73b23516b1ab49721b208fca5116f138bb75df5e6ef960f4693620a7`.
The reviewed README hash at the start of this report remains the proof-audit
lineage. The finite computations remain finite; the all-parameter conclusions
rest on the proofs.
