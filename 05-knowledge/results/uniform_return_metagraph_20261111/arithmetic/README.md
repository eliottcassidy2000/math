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
  [Receipt](quadratic_receipt.json),
  [complete indexed table](A002496_roots_through_1000000.txt).
  The table is independently generated here, with OEIS attribution and the
  pinned10000-term comparison recorded in [the source receipt](source_comparison.json).
  Its last prime is999920001601; no submission or global record is claimed.
- First-three-core square-plus-one charts: all2047 positive even-square
  targets through2^24-1. Axis-only fails at2501; minimum short flank reaches29.
  [Receipt](near_square_charts.json), [witness table](near_square_charts.tsv).
- Primary OEIS numeric tables agree in full at their pinned ranges.
  [Source record](source_comparison.json). No external submission is made.
- [Exact binary chart inversion](INVERSE.md): complete normalized inverses
  in O(sqrt(N) log^2(N)) elementary integer steps for each fixed chart.
  [Controls](chart_inverse.json) compare all90 charts through255 with the
  inherited signed generator. [The +4 ray receipt](ternary_plus_four.json)
  gives target witnesses through k25, without a source-return claim.
- [Two growing flanks](BIFLANK.md): exact product families and a complete
  scoped obstruction for independent fixed-offset growth. The balanced
  growth probe remains explicitly finite. [Receipt](biflank.json).
- [Minimum defect and Fermat candidates](MINIMUM_SQUARES.md): an exact
  square-plus-one intersection of the all-degree extremal root strata,
  plus exclusion of the residual-C3 defect6 sector. [Controls](minimum_squares.json).
- [Complete near-square spectrum at branching defect<=6](TAIL_SQUARES.md):
  every4^m+1,m>=1, and exactly nine additional values. Its prime intersection
  is the Fermat primes of index>=1 and37,101,4357,15877, with the prime
  guard on the infinite family explicit. [Exact certificate](tail_squares.json).

Run from this directory (Python standard library; the chart search also
requires a C++17 compiler named g++):

    python -B maxima.py --output maxima_receipt.json
    python -B audit_maxima.py --output maxima_audit.json
    python -B quadratic_primes.py --output quadratic_receipt.json
    python -B compare_sources.py
    python -B near_square_charts.py
    python -B chart_inverse.py
    python -B ray_inverse.py
    python -B biflank.py
    python -B minimum_squares.py
    python -B audit_defect.py
    python -B audit_cycle_residual.py
    python -B tail_squares.py

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
