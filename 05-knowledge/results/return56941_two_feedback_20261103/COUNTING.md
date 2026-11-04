# A000568: an exact batched tail, and the benchmark that rejects it as default

**PROVED reformulation of an inherited kernel; FINITE-EXACT full-prefix
controls and local timings. No new enumeration record or fastest algorithm.**

The ordinary A000568 counter forgets all tournament representatives. It
cannot classify Hamilton values or actual metagraph lines. Preserve that
scope when combining it with the named-source work in this packet.

## C1. Turn the inherited unit tail into integer polynomial products

Inherit the odd-cycle Burnside sum and
[THM-514 unit-tail kernel](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md).
For a non-singleton odd core mu let d be its mass, r its number of cycles,
E its orientation-orbit exponent and z its permutation centralizer factor.
The contribution after adding k singleton cycles is

    (2^E/z) * 2^(binom(k,2)+kr)/k! .                     (C1)

For a requested whole prefix through N put D=odd(N!). Let

    U_r(x)=sum_(cores mu, mass<=N, length=r) D*(2^E/z)*x^d,
    V_r(x)=sum_(0<=k<=N) D*2^(binom(k,2)+kr)/k!*x^k.

Then

    D^2 sum_(n=0)^N A000568(n)x^n
      = sum_r U_r(x)V_r(x) mod x^(N+1).                  (C2)

This is direct coefficient extraction from(C1), not a new Burnside formula
or state quotient. Both polynomials have integral coefficients: their odd
denominators divide D, and the binary part of the factorial denominators
is cancelled by the respective orbit powers. `fmpz_poly.mul_low` evaluates
only the requested coefficients. Every final division by D^2 is checked.

`count_prefix.py` computes the core using the inherited packed future-gcd
response DP with common denominator D. `bucket_prefix.py` instead preserves
Session87's degree buckets and local normalization O(d)=odd(d!), then
converts its terminal weights by D/O(d) before the same convolution.
Both also expose their inherited direct terminal observer as a control.

## C2. Recover the best relevant baseline before measuring

The first baseline was
[October21 packed DP](../tournament_enumeration_20261021a/check.py).
Three alternating full-prefix trials gave local medians:

|N|older packed DP|common-denominator FLINT convolution|
|---:|---:|---:|
|200|7.109 s|7.168 s|
|256|31.646 s|27.055 s|

That does **not** establish a frontier improvement. Incoming work routed
back to the already faster
[Session87 engine](../session87_four_vertices_enumeration_20261028/FAST_COUNTS.md),
which retains degree-local odd factorials, loops over multiplicity outside
response states, and folds the terminal polynomial before multiplication.
The subsequent same-workload three-trial comparison gave:

|N|Session87 Python buckets|FLINT buckets, folded tail|FLINT buckets, convolution|
|---:|---:|---:|---:|
|200|1.946 s|2.419 s|3.237 s|
|256|9.981 s|12.599 s|13.805 s|

Therefore **retain Session87 as the recommended full-prefix implementation**.
The convolution is an exact alternative, not a demonstrated speedup over
the best recovered baseline. Larger integer objects, denominator scaling
and Python/FLINT transitions outweigh its batching here. These measurements
come from a concurrent workstation, not controlled hardware performance
claims; all sample timings and exact outputs are retained.

`count_comparison.json` records the older baseline; `frontier_comparison.json`
records the corrected performance comparison. No false numerical result
was found. The initial apparent benefit at256 was relative to an older
implementation, and must not be quoted as a repository-frontier speedup.

## C3. Exact outputs, incoming independent certificate, and future scope

Every coefficient through256 agrees between the common-denominator engine,
the independent older integer implementation and Session87. Small orders
0..30 also agree in seven arithmetic/layout variants, including the
rational control retained in `count_comparison.json`. Full a(256) has
30,957 bits,9,319 decimal digits and2-adic valuation131.

Incoming commit83d4954e15 independently publishes
[the full256 certificate](../orbit_spectrum_session_20261103/counting/full256.json).
Its whole-prefix hexadecimal-line SHA256 and full-integer byte SHA256 agree
with our output, checked in `count_controls.json`. That packet records prior
repository scalar computations through300. Thus neither256 nor this prefix
is presented as a new global or repository enumeration record.

Incoming
[Session175's projective quotient](../session175_tournament_spectrum_20261104/a000568/README.md)
retains S_j=j*r-R_j for one exact target mass and compensates by2^((n-d)r).
It is a stronger equivalence for that **single coefficient**. It cannot be
silently reused for all degrees in(C2), where the remaining mass varies.
The distinction between proportional responses at one target and a whole
polynomial of responses is the same practical issue that forces extra
ordered-run information in the source-return problem.

The incoming
[certified leading-bit observer](../orbit_spectrum_session_20261103/counting/NOTE.md)
has a different output contract: a finite positive remainder proves a
requested high-bit prefix without computing the entire integer. It is
the appropriate improvement when only leading digits are needed. Neither
it nor(C2) supplies a Hamilton witness or an automorphism transporter.

Next counting work should reduce core-state arithmetic or support genuinely
shared queries. Retuning a soluble terminal tail alone is now a measured
low-value route at200/256. The constructor work instead needs actual
nonuniform incidence and joint path laws, which this scalar enumerator loses.

The closing fetch also recovered
[Session105's integral projective lineage through IC1--IC2](../second_feedback_18979_20261103/PROOFS.md).
Its latest bucket comparisons independently retain Session87 as the faster
full-prefix choice. Integral projective normalization and degree buckets
must therefore remain credited to that earlier work, even when rederived
in recent packets. This does not change the measured outcome above.
