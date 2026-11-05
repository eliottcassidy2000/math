# Independent review of the rank-one parity corollary

PASS, no correction, bound to `../session197-rank-one-orbits-audit/PARITY.md`, SHA256 `92e483a633710af46b36a4ceb3ad34cfda1a5ed283dde3b5e083caec5aabceb8`.

Every rank-one class has an odd normalization weight. Therefore its unweighted strong-class count has the same parity as the strong weighted series xM'/M from the reviewed orbit correspondence. Over F2, M=1+x+x^2, xM'=x, and x/(1+x+x^2)=(x+x^2)/(1+x^3). Thus the count at positive side order t is even exactly when t is divisible by three. The coefficient at zero is zero, but no empty rank-one object is asserted.

This is a formal-power-series proof at every order, not an extrapolation. An independent exact-integer recurrence d_n=n M_n-sum_(i=1)^n M_i d_(n-i), with d_0=0, checked the parity through n=100. No source file was changed.
