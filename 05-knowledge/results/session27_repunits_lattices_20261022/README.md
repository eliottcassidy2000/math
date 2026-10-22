# Session27: repunit rarity and typed lattice connections

This first independently audited checkpoint recovers prior repunit work, distinguishes labelled from unlabelled tournament counts, and states precise modular/lattice bridges and obstructions. No theorem IDs are reserved and no canon correction is applied here.

- [Repunit recovery and measure-aware rarity](session27-repunits/REPORT.md): provenance routes, full local-search catalogues, and the automorphism-weighted SCC counting law.
- [Typed modular and lattice synthesis](session27-modular/report.md): theta/eta normalizations, E8 and rank-24 characters, the rational-versus-positive-weight modular obstruction, and carefully separated 1729 mechanisms.
- [Independent audit](session27-repunit-modular-audit/REVIEW.md) and [source-bound manifest](session27-repunit-modular-audit/manifest.json).

The complete plain-text JSON catalogues are included: the initial and expanded passes record their own commits, available worktree files, omissions, line previews and hashes. They are discovery indexes, not exhaustive recovery of every historical Git blob or every use of 31, 63 and 127. Initial fetched sources retain their connector provenance. The later remote-search and E8-holes/1729 addenda are not part of this checkpoint or its audit.

## Reproduction

From session27-repunits run `python3 verify_weighted.py`; from session27-modular run `python3 check.py` and `python3 check_camion.py`. Saved JSON outputs and independent replays are included. To repeat discovery against a checkout, run `python3 discover.py /path/to/math /tmp/catalogue.json`; checkout availability and later commits can change the result.

Run `python3 session27-repunit-modular-audit/independent_checks.py /path/to/math` from this packet directory. The checker needs the referenced Session24 order-eight CSV and catalogue source files in that checkout. It independently enumerates small tournaments, checks the two order-eight H=63 witnesses, root coordinates, oscillator coefficients and the Camion control. Catalogue hash differences from subsequent repository changes are reported separately. Without an argument it locates the containing Git repository.

Only the independent checker's two workspace paths were adapted for distribution: repository root is an explicit argument or Git-discovered, and catalogue is a sibling file. Mathematical code is unchanged. The original audit manifest retains the pre-portability script hash and original absolute source locations; these are historical provenance, not distributed-path requirements. The packaged replay and publication plan record the adapted script bytes. `bind_manifest.py` is the original audit-time provenance utility, not a portable reproduction command.

## Correction proposals and boundaries

`session27-repunits/corrections.patch` and `session27-modular/historical_1729_correction_proposal.md` are proposals only. Their inclusion does not apply them, upgrade historical claims, or change canon status. Read the audit's scope and limits before relying on a claim. Imported identities remain cited results; finite checks corroborate but do not prove all-parameter statements. No result here solves LRC(14) or supplies a universal repunit-rarity law.
