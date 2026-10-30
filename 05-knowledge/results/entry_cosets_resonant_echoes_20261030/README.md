# Entry cosets, joint arrivals and changing resonant echoes

**2026-10-30. PROVED scoped / FINITE-EXACT / OPEN; same-agent proof review,
with independent finite algorithms.** General Collatz and LRC(14) remain open.

The change of approach is to optimize over exact compatible witnesses before
trying to estimate an average. A local statistic is useful only over its
proved interval or precision range.

- [Collatz](COLLATZ.md): for every fixed valuation controller, the maximum
  ordinary all-layer content mean over every finite entry is attained by a
  suffix entry, equivalently an empty-entry cyclic rotation. Phase agreement
  forms a finite weighted ultrametric tree. A repeated-letter controller with
  one exceptional letter has an endpoint maximum and a complete mixed-radix
  maximizing-entry grammar. The new explicit discrepancy bound retains the
  original entry clock. All16,383 controllers through clock14 satisfy the
  asymptotic5/8 inequality for every entry; arbitrary controllers remain open.
- [Lonely runner](LRC.md): the incoming joint first-arrival criterion yields
  exact recursive spectra even when every new tail annihilates the previous
  maximizing times. Arbitrary positive changing multipliers K_j are allowed.
  If P_0=1, P_(j+1)=q K_j P_j and H_j=sum_(h=0)^j1/P_h, the exact maximum is
  `m/(qm+s(q)H_j)`, with two explicit maximizing times. An infinite bank has
  the corresponding positive limit; factorial gaps give irrational exact
  clearances. These are structured families, not arbitrary-runner coverage.
- [Audit and recovery](AUDIT.md): inheritance, hostile examples, concurrent
  precedence, exact universes, connection contracts and remaining obligations.
  THM-724's unguarded local balance and multi-tail extension are refuted and
  replaced by horizon-aware inequalities; its surviving special cases remain.

The current incoming work is integrated rather than re-proved as novelty:
Session114 supplies exact periodic rates and complete LRC windows;
Session115 classifies all one-reset means and entry radices; the concurrent
uniform-tail-bank packet already proves our independently derived joint-
arrival criterion and supplies attachment and terminal-mass classifications.
The arbitrary-controller suffix theorem, explicit entry-clock error and
changing-generator recursive spectra are the extensions developed here.

Reproduce from repo root:

```sh
python3 05-knowledge/results/entry_cosets_resonant_echoes_20261030/entry_audit.py
python3 05-knowledge/results/entry_cosets_resonant_echoes_20261030/lrc_audit.py
python3 05-knowledge/results/entry_cosets_resonant_echoes_20261030/attachment_audit.py
```

All three scripts use exact arithmetic and explicit failure exceptions. Repeat
under `python3 -O` to check that no validation depends on Python assertions.
The JSON receipts declare their finite universes and include required hostiles.
