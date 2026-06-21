# LRC14 Periodmax Certificate Formal Kernel

The complete THM-563 bounded-base periodmax audit changes the role of
single-far work.  It is no longer a search target.  It is an import boundary:

```text
finite exhaustive audit  ->  small arithmetic certificate  ->  proof DAG glue
```

The S77 script `lrc_periodmax_worstrow_certificate_codex_s77.py` is deliberately
not a second 12805-row scanner.  It recomputes the exact plateau margins for
the six per-k worst rows and rechecks k=8/k=9 by brute force because those are
small periods.  The larger periodmax values are imported from the completed
mac-mini audit.  The point is to make the constants that downstream proof work
will cite hard to misquote:

```text
k=8   headroom 303/392                 periodmax=2
k=9   headroom 44585/196196            periodmax=86/49
k=10  headroom 464839/535080           periodmax=3299/1470
k=11  headroom 5417609/2942940         periodmax=6730439/1961960
k=12  headroom 2664689/1177176         periodmax=536399/196196
k=13  headroom 2207057/588588          periodmax=122064/49049
```

The k=8 line matters as a guardrail.  The complete summary's certificate table
prints `pm=1`, but the row ratio `10.8188`, the threshold, and the brute
period scan all require `pm=2`.  This does not weaken the theorem; the
headroom is still `303/392`.  It does mean downstream formal notes should cite
the S77 kernel or the per-k output rather than copying that literal.

The Lean module `LRCPeriodmaxCertificate.lean` mirrors exactly this boundary.
It proves finite arithmetic facts: positive worst-row headrooms, k=9 as the
largest ratio among the six per-k worst rows, the completed count checksum, and
the k=8 normalization.  It intentionally does not claim the full enumeration;
that remains a Python/mac-mini certificate until someone generates a full row
atlas import.

Tournament Analysis also becomes cleaner here.  Useful vertices are proof
carriers:

```text
full_periodmax_scan
worstrow_arithmetic_kernel
endpoint_orbit_dedekind_packet
genuine_wide_survival_currency
q0_boolean_guardrail
Phi_low_boolean_transfer
raw_sumR_or_ETK_abs_bound
runner_vertex_tournament
```

The switch is predicate fidelity after HYP-2793's proof-DAG compression.  The
endpoint-orbit quotient preserves the single-far inequality, while runner
vertices and raw absolute discrepancy forget the signed phase that made the
single-far leg close.  This is a good example of the current LRC14 pattern:
formalize closed finite arithmetic aggressively, but keep genuine-wide
relation/survival data unsquashed until the final cap comparison.
