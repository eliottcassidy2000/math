# LRC n=14 Information Bottleneck, S612

The recent n=14 progress is starting to look less like a hunt for a larger
finite table and more like a sufficiency theorem.

The key distinction is address entropy versus relevance entropy. Address
entropy asks how many bits it takes to name a row in some layer. Relevance
entropy asks how many bits of the proof predicate remain undecided. S610
crushed address entropy: raw `Res_27` subsets need `23.31` bits, the S610
proof ledger needs `3.46` bits, and primitive floor atoms need `1` bit. But
S611 showed that a small quotient can still be lossy if it forgets the wrong
coordinate.

The carry coordinate is the warning label.

For a lift `v=r+27k`, the congruence `27 == -1 (mod 14)` gives

```text
v == r-k (mod 14).
```

So the least-positive `Res_27` shadow is not what the n-clock sees. The
n-clock sees residue plus carry. In the AP/`V*` scalar audit, the carry
indicator has only `0.4138` bits of entropy, but it perfectly separates floor
shadows from strict shadows. That is low Shannon entropy and high proof
relevance.

This reframes the quotient tower:

```text
raw row
  -> sufficient quotient candidate
  -> side-information reattachment
  -> floor predicate.
```

The theorem should not say "compress the row." It should say "forget exactly
the data independent of the floor predicate, and keep every side channel with
nonzero conditional mutual information." In the current n=14 stack, the side
channels are owner labels and the carry cocycle.

From the anti-Poisson frame, this is the rare-event lesson. Floor rows are
rare enough that the Shannon entropy of the pinch status is only `0.0016`
bits, but those rare bits are the proof. An entropy-minimizing quotient that
throws them away is useless. A proof-relevant bottleneck has to preserve the
coimage of the floor predicate, not just minimize the row address.

The clean theorem shape is:

```text
Floor(V) is independent of the raw integer presentation
given (Res_27 proof atom, owner route, carry cocycle, Cprime window).
```

That is the information-theoretic reading of HYP-2166 and HYP-2167. S610
found the bottleneck. S611 found the side-channel leak. S612 says the next
proof should be a sufficiency theorem for the corrected statistic.

The incoming HYP-2168 C2b work strengthens this reading: the multiple-of-n
branch is a positive-measure `p_0>0` problem with `2n-1` tick and small-q clock
probes, so the `Cprime window` belongs inside the statistic rather than after
it.

Artifacts:

- `04-computation/lrc_n14_information_bottleneck_s612.py`
- `05-knowledge/results/lrc_n14_information_bottleneck_s612.out`
- `05-knowledge/hypotheses/HYP-2171-lrc-n14-information-bottleneck.md`
