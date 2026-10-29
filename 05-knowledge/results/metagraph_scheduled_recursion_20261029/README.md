# Tournament recursion by requested precision and scheduled continuations

**2026-10-29. PROVED scoped / FINITE-EXACT / measured implementation
experiment. Same-agent proof review and alternative exact algorithms; no
independent audit or literature-wide novelty claim.** Base `201fa96604`.

This session improves the *requested-output recursion* rather than claiming
a new Burnside formula or a smaller universal structural carrier. The newest
Session105 integer projective method was already substantially faster than
its Fraction predecessor, but slower at size200 than Session87's richer
whole-prefix output. Both comparisons are inherited, not newly discovered.

## Results

1. **A precision-dependent exact recursion.** To obtain b low binary digits
   after the natural normalization of a tournament counting channel, only
   O(sqrt(b)+1) permutation cycles can contribute, independently of n.
   The implementation retains arbitrary odd cycle lengths, gcd interactions
   and the factorial's binary valuation. For fixed b this gives a polynomial
   bound in the numeric value of n, not in its binary input length.
2. **A sharp first correction at every odd size.** For r tournament layers
   and s binary attributes, put A=r(n-1)/2+s and delta=2r+2s-1. For every
   odd n>=3,

   \[
   v_2\bigl(n F_n(r,s)/2^A-1\bigr)=\delta.
   \]

   In particular the ordinary channel has delta=1 and the aligned-pair
   channel has delta=5. This strengthens Session89's first normalized bit.
   The missing correction comes from exactly the partition(1,1,n-2), with
   the all-one triple handled separately at n=3.
3. **Metagraph growth checks.** The existing root-switch reset Q satisfies
   trace(Q)=n a_n. At odd n>=3 its normalized trace is3 modulo4. At even
   child order N>=4, the number R_N of pointed extension states satisfies
   `v2((N-1)R_N/2^(N/2)-1)=3`. The complete inherited order-ten bank passes,
   including its97,213,472 pointed states. Literal root-reset matrices
   through order5 are independently reconstructed here.
4. **Contract forced-zero stages.** An exact event schedule skips cycle
   lengths that no longer fit the remaining capacity. At bound200 it
   reduces the inherited prefix recurrence from2,994,814 to2,905,899
   transfers. Matched timings show only a small gain; no material or
   universal speedup, smaller response quotient, or memory reduction is
   claimed. Session87 remains the established production baseline.

Proofs: [PRECISION.md](PRECISION.md), [SCHEDULE.md](SCHEDULE.md).
History and accepted connections: [RECOVERY.md](RECOVERY.md).
Validation and limitations: [AUDIT.md](AUDIT.md).

Some exact **partial outputs**, not full counts at these orders, are

\[
\frac{a_{1001}}{2^{500}}=11\pmod{32},\qquad
\frac{F_{1001}(2,1)}{2^{1001}}=57\pmod{64}.
\]

The first uses21,001 partitions with one or three cycles; the second also
agrees with that independent three-cycle evaluation. At n=1,000,001 the
sharp formula gives normalized residues3 modulo4 and33 modulo64 in the
same two channels, without calculating either complete integer. These
statements do not extend the complete A000568 prefix or construct a
million-vertex metagraph.

## Inheritance and portfolio

**Anchor:** exact A000568 recurrence and root-switch checks. **Niche:** the
March enumeration experiments and the role of their admissible permutation
filter. **Wildcard:** precision as a complexity parameter rather than size
alone. The concept board is gcd responses, exact-degree gauge, binary
valuation, root transport, and enabled future operations.

Closest mechanisms are Session89's channel low-bit theorem, Session87's
degree buckets, Session98's projective quotient and Session105's integer
realization. The hostile is the twelve-vertex81/27 pair:67 labelled cyclic
triples agree, but the transporter3 to7 has sizes27 and0. This is replayed.
The corrected near miss is March16's all-partition sum labelled A000568;
the next-day odd-cycle repair is now propagated back as a warning. The
least-used useful sidecar is the binary part of multiplicity factorials:
it determines exactly where the first extra bit appears.

## Reproduce

```text
python -X utf8 -B 05-knowledge/results/metagraph_scheduled_recursion_20261029/verify.py
python -X utf8 -B -O 05-knowledge/results/metagraph_scheduled_recursion_20261029/verify.py
python -X utf8 -B 05-knowledge/results/metagraph_scheduled_recursion_20261029/benchmark.py
python -X utf8 -B 05-knowledge/results/metagraph_scheduled_recursion_20261029/precision_counts.py 1001 5
```

`verify.py` compares3,402 channel/size/precision requests with full counts,
checks10,367 individual summand bounds and independent Fraction sums,
compares complete prefixes through200, and retains literal and inherited
metagraph controls. Normal and optimized execution use explicit checks.
Sources and output hashes are saved in `verification.json`; timing trials
and machine metadata are in `benchmarks.json`.

## Open targets

- **OPEN:** replace the bounded-cycle enumeration by a faster arithmetic
  convolution for three- and five-cycle terms, retaining gcd and repeated
  parts. This could extend useful precision without computing full counts.
- **OPEN:** combine integral projective states with adaptive multi-output
  queries without losing the multiplicity coordinate needed by consecutive
  degrees. A one-degree quotient cannot simply be reused for a prefix.
- **OPEN:** a comparable precision truncation for actual joint attachment
  responses. Burnside's global cycle-energy bound does not automatically
  bound a native port entry or its complement. The81/27 control forbids
  treating the arithmetic recurrence as a structural classifier.
