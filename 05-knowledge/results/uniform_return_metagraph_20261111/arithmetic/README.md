# Arithmetic fibres and exact sequence enumeration

**PROVED scoped / FINITE-EXACT; Chan–Pak, Goldbach and Landau OPEN.**

[Proofs and inheritance](PROOF.md) contain the factor-count source bounds,
the complete consecutive-prime maximum family, the exact A002496/A135244
extremal equivalence, and the composition boundary.

- A135244: complete maximum-source prefix2..1,000,000, all999,999 entries
  independently compared. [Producer receipt](maxima_receipt.json),
  [independent audit](maxima_audit.json).
- A002496: complete even-root sieve through1,000,000,54,110 prime values,
  with explicit Lucas and composite-divisor certificates.
  [Receipt](quadratic_receipt.json).
- First-three-core square-plus-one charts: all2047 positive even-square
  targets through2^24-1. Axis-only fails at2501; minimum short flank reaches29.
  [Receipt](near_square_charts.json), [witness table](near_square_charts.tsv).
- Primary OEIS numeric tables agree in full at their pinned ranges.
  [Source record](source_comparison.json). No external submission is made.

Run from this directory (Python standard library; the chart search also
requires a C++17 compiler named g++):

    python -B maxima.py --output maxima_receipt.json
    python -B audit_maxima.py --output maxima_audit.json
    python -B quadratic_primes.py --output quadratic_receipt.json
    python -B compare_sources.py
    python -B near_square_charts.py

To export the complete computed numeric tables, add `--table output.txt`
to `maxima.py` or `quadratic_primes.py`. The stored receipts hash those
complete normalized streams; bulky generated tables and executables are
not committed. The offline source comparison uses the two attributed
numeric b-files saved here; `--refresh` explicitly downloads new versions.

The cited sources and inherited code are listed in the proof and receipts.
The factor-layer and consecutive-prime arguments have separate proofs;
finite controls do not supply their universal quantifiers. The large
prefixes are workspace computations, not claimed world records or new
OEIS terms. Normal and optimized replay agreement is recorded separately
at close-out.
