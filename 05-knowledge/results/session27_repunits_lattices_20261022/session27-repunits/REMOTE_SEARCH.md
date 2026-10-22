# Bounded remote discovery addendum

The frozen REPORT.md and catalogue remain separate from this addendum. On 2026-10-22 the GitHub connector `search` successfully searched repository `eliottcassidy2000/math`, default branch, with `topn=100` and these exact plain-string queries:

- repunit: 100 results, 48 absent locally.
- Mersenne: 100 results, 58 absent locally.
- 1729: 100 results, 36 absent locally.
- Wagstaff: 38 results, 8 absent locally.
- Fibonacci fan: 100 results, 32 absent locally.

The returned URLs pin snapshot 54bf96642934c6cdf01f02e364f3bba00892ae22. All 438 result rows and returned excerpts are preserved in `remote_repunit.json`, `remote_Mersenne.json`, `remote_1729.json`, `remote_Wagstaff.json`, and `remote_Fibonacci_fan.json`. `remote_index.json` combines them with a local-presence flag. They comprise 415 distinct paths, including 174 absent locally and 193 not previously present in the local match catalogue. The latter figure includes new vocabulary/1729 hits, not necessarily new relevant repunit work.

## Search limits

Four queries hit the requested 100-result cap. The connector exposes no cursor/page parameter, total-hit count, or exhaustion flag. The 38 Wagstaff results therefore do not establish exhaustive coverage either. Multiword 'Fibonacci fan' behaved as token matching, with fan matching Fano and many irrelevant results. Numeric 1729 can also match message IDs, not mathematics. This is bounded complementary retrieval, not a full remote-content census or history search. No access denial or API error occurred. No unsupported GitHub-query flags, browser workaround, or bulk blob fetch was used.

## High-value additional recoveries

Three previously absent full files were fetched at the verified commit above; their connector blob SHAs and URLs are in `remote_fetch_manifest.json` and their contents in `remote_fetched/`:

1. `04-computation/experiments/repunits_hyp9357_20261021/audit_S51/REPORT.md`: original independent audit. It confirms reproduced numerical counts while correcting generalized subtraction, pre-registration blindness, low-feedback-arc explanation, and omitted prior art. This is stronger provenance than quoting only the resulting ledger entry.
2. `04-computation/witt_necklace_s90at.py`: March 16 trace-to-Witt/necklace exploration. An early geometric-series/primitive-necklace connection, with explicit exploratory self-corrections. It is historical idea provenance, not a proved identification of Hamiltonian multiplicities with necklace coefficients.
3. `04-computation/mersenne_cyclotomic_ocf.py`: March 14 cyclotomic interpretation of the OCF evaluation point. It predates the October synthesis and explicitly studies Phi_p(2)=2^p-1 and Phi_3(2)=7. Its opening is a proposed connection and finite-order observation, not an all-order spectral completeness proof.

The remote results also recover source routes for `h63_gap_check.py`, `forbidden_values_audit.py`, `universal_zone_s90ar.py`, `formal_group_s90as.py`, and the complete S51 audit scripts. Their returned excerpts are indexed; full-file verification is not claimed.

## 1729 handoff

The 1729 query supplied missing-local routes including `ramanujan_near_miss_kpc1.py` with two verification companions, `modular_tournament_s18e.py`, `moonshine_plus_one_s18g.py`, `ham_excess_spectra_kps_S134.py`, and `ramanujan_1729_bridge_s133.py`. These were sent to the modular/lattice lane to reconcile with its correction history. Avoid treating these exploratory filenames or matching numbers as evidence for a modular-form theorem.

## Reproduction

For each query above call GitHub connector search with repository_name=`eliottcassidy2000/math`, query equal to the exact string, topn=100. Search follows the current default branch, so later reruns may differ; use the saved returned commit URLs and blob SHAs for this snapshot. Full targeted reads used fetch_file with ref=`54bf96642934c6cdf01f02e364f3bba00892ae22`.
