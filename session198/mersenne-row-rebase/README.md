# A balanced rebase for the complete 3*2^b-1 row

For every b>=5, this independently reviewed construction gives an actual balanced C1 tournament with exactly 3*2^b-1 Hamilton paths and b+4 vertices. The intermediate source profiles are primitive. Two linked power-of-two run lengths cancel the intermediate term, covering even and odd exponents, with a separate seed at b=5.

Scalar attainment of 3*2^b-1 was already known from transitive cloning. The contribution here is the balanced C1 realization, linked-run cancellation, and explicit graph order.

This removes an infinite row of the transitive residual arithmetic family. It is a positive all-exponent construction, not a proof of universal or cofinite odd coverage, and it makes no exclusion claim about the three smaller row values.

Run `python3 verify.py` and `python3 -O verify.py`. The frozen receipts agree and include literal source and final graph checks through b=10. The independent audit and unchanged source-grammar/C2 proof inputs are included. Exploratory template searches are not part of this release.
