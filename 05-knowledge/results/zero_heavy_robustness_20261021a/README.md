# Empty gaps, first shared endpoints, and physical-port robustness

**PROVED scoped / FINITE-EXACT controls; independent audit pending.**

- [Zero-heavy mass layers](zero_heavy.md): at deficit k, the first empty
  gap gives F_k[binom(3k,k)-1] classes for k>=3, with the exact strong
  correction. At deficit k+1, one P3 is the first possible nonmatching
  component; exact all/strong formulas hold for k>=4. A unified binomial formula
  closes every 0<=delta<=2k-3 for k>=3. These are total-
  weight deficits, not heavy-coordinate deficiency or feedback-vertex
  bounds. A minor inherited onset sentence is corrected from k to k+1.
- [Physical-port robustness](palettes.md): for N actual root-only gates,
  the exact repair distance is N-max_c m_c. Complementary palettes attain
  floor(N/2), arbitrarily large, although every whole graph is four-
  colourable. Conditional robustness is not unconditional palette forcing.
- [Inheritance and connection contracts](INHERITANCE.md).
- [Same-agent audit, hostile controls and open boundaries](AUDIT.md).

Reproduce from this directory, or prefix these filenames by the published
packet path from the repository root:

    python zero_heavy_check.py
    python -O zero_heavy_check.py
    python palette_check.py
    python -O palette_check.py
    python verify.py

palette_check.py discovers the repository from its location; alternatively
set MATH_REPO to the pinned checkout. The manifest records all mathematical
source dependencies and authenticates the gate receipt before use. The
standalone zero-heavy checker has no repository imports.

The final verification command validates frozen file/dependency hashes and
reruns both checkers in both modes against the stored JSON receipts. There
are no third-party Python dependencies.
