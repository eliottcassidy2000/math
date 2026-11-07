# Reproduce the audited logarithmic refinement

Read RESULT.md for the theorem, audit/REVIEW.md for the independent mathematical review, and audit/RELEASE_REVIEW.md for verification of the correction and this package. Dependencies are pinned and included in dependencies/.

Python 3, standard library only. From this folder:

    python3 check.py > check-rerun.json
    python3 audit/audit.py > audit/portable-rerun.json
    sha256sum -c SHA256SUMS

The scripts also run from other directories when invoked by absolute path. Finite controls do not replace the analytic proof. Frozen original text, the mathematical correction patch and original-input hashes are retained under audit/ for provenance.
