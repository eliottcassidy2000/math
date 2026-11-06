# Session 208: Rational strong-atom mixtures

The attached paper's rational support-compression mechanism yields an independently proved finite-dictionary LP framework for Session207's inverse Hamiltonian-path order problem. Its central logarithmic spanning-tree theorem is not used.

The concrete construction combines three strong order-eight H=375 tournaments and two strong order-eight H=405 tournaments to realize H=15^11 on forty vertices. Thus the asymptotic minimum-order efficiency for powers of 15 is at most 40/11. Exact rational prime prices prove this rate optimal within the complete strong-order-at-most-eight dictionary, not globally.

- `PROOF.md`: general theorem, explicit example, scope and exact all-dictionary envelope.
- `verify.py`: independent witness checks, complete dictionary extraction from prior census, exact dual and integer DP checks.
- `census_session207.json`: unchanged audited Session207 input census.
- `RESULT.json`, `verify.out`: deterministic verification output.

Run `python3 verify.py` with assertions enabled. No third-party Python dependencies are required. The general LP and support-elimination arguments are standard; the contribution is their precise transfer and the explicit checked atom mixture.

Independent replay: `python3 independent_check.py` should match `independent.out`. `DEPENDENCIES.json` pins the prior census and paper version; `REVIEW.md` records the independent audit. Verify frozen payloads with `sha256sum -c SHA256SUMS`.
