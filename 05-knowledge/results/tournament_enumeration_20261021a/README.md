# Faster exact tournament counting with packed responses and odd denominators

2026-10-21. **PROVED:** the uniform odd-denominator lemma, its exact integer
transition, and the lossless packed future-response representation for every
finite order bound. **FINITE-EXACT:** the declared runs and benchmarks below.
**INHERITED:** Burnside enumeration, the future-response quotient, and the
pair-free-core cycle-index substitution. **OPEN:** a polynomial-complexity
all-order counting algorithm or a smaller general mathematical state quotient.
No historical-priority claim is made for factorial clearing or bit packing.

This packet accelerates the existing algorithm without changing the counted
objects or its state equivalence. On this host, the order-200 A000568 run fell
from 50.30 seconds to 7.94 seconds; both execute 2,994,814 transition
contributions and reach a maximum of 96,836 states. The improvement is chiefly
exact integer arithmetic, with a smaller contribution from packed profiles.
There is no new order-11 tournament census.

## 1. Inheritance, target and connection contract

The closest proved mechanisms are the
[future-response recurrence](../tournament_future_boundary_20261020/README.md),
the [live-divisor and unit-tail kernel](../defect_boundary_20261020/boundary.md),
and the [canonical weighted-core transport](../weighted_core_transport_20261020/core_count.md).
The original formulas are also recorded in
[THM-514, A000568 unit-tail kernel](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md)
and [THM-516, divisor-profile form](../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md).
We do not re-promote the already proved future-gcd quotient as new.

The canonical hostile is equal mass/cycle count with unequal future gcd
responses: histories `(11,5)` and `(9,7)` respond differently to a future
three-cycle. A second inherited warning is that substitution excess is not an
actual permutation cycle. The least-used sidecar here is the exact common
denominator of the **individual histories**, before summing them.

The current
[THM-4741 presentation-count packet](../../../01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md)
retains its correction lineage: rigid stacks and the forbidden Hamilton-path
counts 7 and 21 are inherited, and those presentation counts are not A000568.
This computation does not identify individual tournament classes, strong
components, automorphisms, Hamilton paths, or blackgraph neighbours.

The live board is: raw histories / remaining responses / exact denominators /
signed core excess / physical computational cost. Anchor: faster A000568;
niche: pair-free cores; wildcard: the rational fixed-point-free observer.

**Connection contract.** Source: the inherited finite-horizon response DP.
Target: integer weights and a packed response key. Map: multiply each history
weight by the odd part of `N!`, and encode the response tuple in a declared
binary radix. Preserved: every transition contribution, cancellation, state
merge and final coefficient. Destroyed: no mathematical data; tuple length is
recovered from the retained stage, degree, mode and bound. Needed sidecars:
those four fields and the common denominator. Cheapest decisive tests: exact
row-by-row state/transition identity and all coefficients against the old DP;
too-small radix and excess-as-cycle hostiles reject tempting shortcuts.

## 2. Minimal uniform denominator for the raw histories

Fix `N>=0`. Let `mu` be a multiset of odd cycle lengths with actual mass
`M<=N`, multiplicities `m_l`, and

\[
 z_\mu=\prod_\ell \ell^{m_\ell}m_\ell!,\qquad
 e(\mu)=\sum_\ell m_\ell(\ell-1)/2
       +\sum_\ell\ell\binom{m_\ell}{2}
       +\sum_{\ell<j}m_\ell m_j\gcd(\ell,j).
\]

Write `odd(a)=a/2^{v_2(a)}` and `D_N=odd(N!)`.

**Lemma 1.** Every rational `2^{e(mu)}/z_mu` has odd denominator dividing
`D_N`. Moreover `D_N` is the least positive integer that clears **all these
individual raw weights** of mass at most `N`.

**Proof.** For every `k>=0`, `v_2(k!)<=binom(k,2)`: for `k>=2`, the stronger
`v_2(k!)<=k-1<=binom(k,2)` suffices, and `k=0,1` are immediate. Thus the
within-length contribution to `e` cancels every factor two of `z_mu`.
The centralizer `z_mu` divides `M!`, because `M!/z_mu` counts permutations of
that cycle type. Hence `odd(z_mu)` divides `odd(N!)`. Finally the all-singleton
history at mass `N` has weight `2^{binom(N,2)}/N!`, whose reduced denominator
is exactly `D_N`. QED.

This is minimality for raw history weights, **not** a claim about the smallest
denominator of already aggregated states, and not a lower bound on all
possible algorithms. Removing length one preserves sufficiency but can remove
the witness for minimality.

At fixed stage and committed degree, store the sum
`V=D_N * sum(history weights)`. It is an integer, including signed core
substitution terms. Adding `k` actual cycles of current odd length `ell`
uses the inherited exponent

\[
 E=k(\ell-1)/2+\ell\binom{k}{2}+kR_\ell.
\]

Its exact integer update is

\[
 V' = \frac{V\,2^{E-v_2(k!)}}{\ell^k\operatorname{odd}(k!)}.\tag{1}
\]

The exponent is nonnegative by Lemma 1's local inequality. The division is
exact: each contributing raw history becomes a valid new history of mass at
most `N`, hence has its denominator cleared by `D_N`. Sums inherit this
divisibility. The checker uses `divmod` on every production transition and
rejects a remainder; it never rounds a rational transition.

The final ordinary coefficient is the integer state divided by `D_N`.
The state integer itself is not a class count: at `N=3`, final state weight
6 represents the two unlabelled tournaments.

## 3. Packed responses, with the length retained

Process odd lengths in descending order, starting at 1 in ordinary mode or
3 when fixed points are omitted. At stage `ell` and committed degree `d`, the
live coordinates are the ascending odd lengths

\[
 j\in\{j_0,j_0+2,\ldots\},\qquad j\le\min(\ell,N-d).
\]

For every such coordinate,
`R_j=sum_{a in mu} gcd(a,j)<=sum a=M<=d<=N`.
Choose `b=max(1,bit_length(N))`, `B=2^b>N`, and encode

\[
 R=\sum_i R_{j_0+2i}B^i.\tag{2}
\]

**Lemma 2.** At a declared stage, bound, mode and degree, (2) is injective.
After inserting `k` actual `ell`-cycles, precompute

\[
 G_\ell=\sum_{j<\ell}\gcd(\ell,j)B^{(j-j_0)/2}.
\]

Then the new packed response is `(R+kG_ell) & (B^r-1)`, where `r` is the
number of future coordinates at the new committed degree.

**Proof.** Every retained digit lies in `[0,N]`, so the radix expansion is
unique. In each future coordinate the addition is exactly the inherited
`R_j+k gcd(ell,j)`. Its value is at most the new actual mass, hence at most
`N`; no carry can cross a digit boundary. The mask deletes precisely the
coordinates that no continuation can query. QED.

Only `N,ell,d,j_0` determine the length; it is not encoded in the integer.
For example the empty tuple and `(0)` both encode as zero, which is harmless
only with this sidecar. At `N=2`, an incorrect one-bit radix identifies
`(2,0)` and `(0,1)`; the correct radix separates them.

The tuple implementation recomputes a gcd and a response entry per surviving
coordinate of every transition. This implementation precomputes each gcd row
once, then performs a bigint multiply/add and a prefix mask. Bigint operations
are **not unit cost**: a response may occupy `O(N log N)` bits. The packet
does not claim that the number of states or total bit complexity is polynomial.

## 4. Signed cores and a rational observer

The inherited core specialization is

\[
 Z_C(p)=Z_T\bigl(p_\ell/(1+p_\ell)\bigr),\quad
 (x^\ell/(1+x^\ell))^k
 =\sum_{s\ge0}(-1)^s\binom{k+s-1}{s}x^{\ell(k+s)}.
\]

For `k>0`, multiply (1) by the signed integer binomial and increase committed
degree by `ell(k+s)`. Responses increase by **k**, not `k+s`. For `k=0`, only
`s=0` is allowed. Actual mass remains at most committed degree, so both lemmas
continue to apply. Negative intermediate weights and their cancellations are
retained exactly.

For instance one 3-cycle with excess one has degree six and singleton response
one. Inserting a singleton multiplies its signed weight by 2, not 4.
Treating excess as a second cycle already changes this individual degree-seven
contribution from `-4/3` to `-8/3`.

Omitting length one computes the inherited fixed-point-free automorphism
mass, not the number of a subclass of cores. Its degree-three coefficient is
`2/3`. The API returns exact `Fraction` coefficients in every mode and does
not force the observer into integers. The production recurrence itself uses
integers throughout.

## 5. Reproduction and measured cost

[check.py](check.py) produces [receipt.json](receipt.json), which includes the
full sequences, dependency hashes and per-stage work statistics. Validation
remains active under `python -O`.

```text
python -B 05-knowledge/results/tournament_enumeration_20261021a/check.py
python -B -O 05-knowledge/results/tournament_enumeration_20261021a/check.py --verify
python -B 05-knowledge/results/tournament_enumeration_20261021a/check.py --benchmark 100
python -B 05-knowledge/results/tournament_enumeration_20261021a/check.py --benchmark 200
python -B 05-knowledge/results/tournament_enumeration_20261021a/check.py --benchmark 100 --core
python -B 05-knowledge/results/tournament_enumeration_20261021a/check.py --benchmark 100 --memory
```

Declared mathematical controls:

- All bounds 0–24, each of four ordinary/core and fixed-point modes, compare
  packed integers, packed rationals and the unchanged inherited tuple DP.
  State and transition counts match at every layer.
- All four modes compare unmerged partition/substitution calculations through
  18; opposite elimination order is checked at 24.
- Every odd cycle partition through mass 30 checks the common denominator;
  the all-singleton history checks its claimed minimality.
- Direct invariant-orientation enumeration and a homogeneous-pair test through
  six independently control actual core counts.
- A000568 through 200 agrees with the frozen preceding receipt; all 77 entries
  of the independently stored [OEIS table](../tournament_future_boundary_20261020/a000568-reference.txt)
  through order 76 are checked. Unmerged Burnside through 60 is recomputed.
  There is no independent external count claim for orders 77–200.
- Pair-free classes and the rational no-fixed core observer through 100 are
  recomputed by both implementations and agree with the frozen receipt.
- Literal type, insufficient-radix, missing-length, substituted-cycle and
  integer-versus-observer hostiles are retained.

Single sequential timing observations on this Windows host, Python 3.13.14
AMD64, are recorded verbatim in [benchmarks.json](benchmarks.json):

| A000568 bound | Old tuple / Fraction | Packed / Fraction | Packed / odd integer |
|---:|---:|---:|---:|
| 100 | 1.107 s | 0.980 s | 0.171 s |
| 200 | 50.303 s | 39.992 s | 7.935 s |

The pair-free core run through 100 took 15.958 / 15.461 / 2.009 seconds in the
same three implementations, with 1,697,975 contributions and 5,480 peak states.

At bound 100, a separate `tracemalloc` run measured peak traced Python
allocations of 2,935,995 / 1,672,240 / 1,504,240 bytes respectively. These are
not RSS measurements. Tracing changes timing, so the memory-run times are not
used for the speed ratios. Single runs establish a concrete local improvement,
not a portable speed guarantee or a statistically estimated performance law.

The deterministic order-200 work ledger records 10,767,006 tuple-coordinate
updates replaced by packed updates, against 4,950 gcd entries precomputed for
the new implementation. Peak live tuple-coordinate count is 609,498; peak sum
of the actual packed-key bit lengths is 4,488,501 bits. These logical sizes
are not heap measurements. The common denominator itself has 1,049 bits.

The old and new mathematical state counts are identical. A more ambitious
next target is a closed summation over a final divisor-profile family; it must
preserve actual response distinctions or prove that summation removes their
effect. Merely replacing the tuple by another coordinate system does neither.

## 6. Audit boundary

This packet's all-order proofs concern arithmetic and representation of an
already proved recurrence. Finite reproduction validates the implementation
on its declared universe. No new literature theorem, universal strong-core
classification, symbolic lower bound on every algorithm, or asymptotic
complexity improvement is claimed.

Normal and optimized-mode replay agree at **3,243,955 checks**.
Receipt LF SHA256:
`3a452d88007b1f596e768a1c6c426780e7b739ca76b44a09c8a56816ee48b84c`.
Independent peer review by the `tiling_fold` agent passed the complete proof,
production interfaces and optimized saved-receipt replay at the same count and
hash. It specifically checked raw-history denominator minimality, signed
transition divisibility, packed no-carry updates, retained length and the
core-excess distinction. The benchmark timings were not independently rerun;
they remain the single-host observations stated above.
