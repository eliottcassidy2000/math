# Independent audit: Gerbicz square-sum reconstruction (Session 124)

## Verdict

PASS. The reconstructed finite certificates plus the stated induction prove
Hamiltonian cycles in the square-sum graph on [1,N] for every N >= 32.
Consequently they prove Hamiltonian paths for every N >= 25 (using the
literal paths at 25 through 31). This is a verification of Gerbicz's classical
construction, not a claim of a new solution to the ordinary conjecture.
No claim of nonexistence at other small orders is certified here.

## Independent method and source boundary

`check_independent.py` is independently written; it does not import or run the
producer's verifier/constructor. It parses only integer initializer data from
three uniquely named arrays with a recursive brace/comma/integer grammar.
Every token and complete terminating initializer is checked. Expected table
sizes are checked before mathematics. The downloaded C is never compiled,
imported, or executed. Its program assertions are not trusted as proof.

Both portable verifiers reject a mismatched source SHA256 before parsing.

Source SHA256:
`602c1919ee51101519a2f19f2c1817407006824f6f09c85a3b772a818a0b736e`.
No original source, archive, or extracted certificate table is included in
this audit directory. Redistribution permission for that source is unknown.
Source provenance and attribution are documented by the producer in
`GERBICZ_RECONSTRUCTION.md`.

## Independently checked finite obligations

* All 3,984 base rows: exact permutation coverage, every edge and closing
  square, endpoint 1 and terminal 3 (odd) / 8 (even).
* All 1,992 adjacent-order pairs at n = 41,...,2032: common labels occupy
  positions of matching parity.
* All 196 universal words: 49 distinct offsets -24,...,24 and 24 distinct
  singleton labels 1,...,24, legal orientations/member choices, exact
  extension cutoff, total length, endpoints, every connector and closure.
* Physical row order agrees with the original dispatch formula
  4(d-24) + 2(n mod 2) + e; no duplicate or missing state can hide.
* All 98 paired transitions preserve singleton position parity and the
  positional-parity shift for each transformed block.
* Literal paths at 15,16,17,23,25,...,31; literal cycles at 32,...,40.
* 196 fully expanded paths exercise all universal states from input pairs
  at 41 and 42. Ten further fully expanded paths exercise recursive pairs
  at 2033,2079,2080,100000,1000000. These are extra integration checks,
  not the justification for universal induction.
* Five negative controls are rejected: missing word, reversed dispatch,
  illegal offset, wrong member, and corrupted singleton.

## Universal reasoning independently reviewed

Write P_0,P_1 for a nice pair of lengths n,n+1. For a local position j,
use T_c(x_j)=49x_j+(-1)^(j+1)c. Every inherited edge has its offsets cancel,
so its sum is 49 times an input square. Reversal does not change this.
For each common x <= n, its local position parity agrees in both members;
therefore all 49 offsets produce exactly [49x-24,49x+24], irrespective of
which member each block uses. These intervals are adjacent and disjoint.
Together they cover [25,49n+24].

The new label n+1 has position parity n+1 in P_1: subtract the common-label
counts in odd and even positions from the total position counts in the two
paths. It is NOT assumed to occupy the last position. Its lifted value is
49(n+1)+(-1)^n c. Thus selecting member 1 precisely when
49+(-1)^n c <= d+e fills exactly [49n+25,49n+d+e]. Adding singletons 1..24
establishes exact interval coverage, with no collisions.

Only the parity of n is needed to evaluate block lengths, endpoints, and
position shifts. For a block after t existing positions, the output parity
shift relative to its original index j is t for a forward block and
 t+L+1 for a reversed block, modulo 2. The checker evaluates these expressions
at one representative of each parity, which is exact symbolic reduction,
not an extrapolation. Matching shifts for every c and matching singleton
parities prove compatibility for every common output label, including
shared extensions of the input pair.

For N >= 2033 define m=floor((N-24)/49) and d=N-49m. Necessarily
24 <= d <= 72 and 41 <= m < N. Nice pairs are supplied for every
41 <= N <= 2032; strong induction therefore supplies every larger pair.
The equality 2033=49*41+24 verifies the lower boundary exactly. There is no
finite gap. Start 1 and terminal 3/8 close because 1+3=4 and 1+8=9.
Small literal cycles extend the scope down to 32.

## Replay

```
python check_independent.py /path/to/squares.c
python -O check_independent.py /path/to/squares.c
```

`audit-normal.json` and `audit-optimized.json` record the same 41,275,807 explicit checks,
206 concrete paths, source hash and rejected negative controls. The checker
uses explicit exceptions; optimized Python cannot remove its checks.
No defect was found in the producer's final proof. The audit prompted an
explicit statement of the extra-label parity-count argument and a stronger
physical-dispatch check, both now incorporated by the producer.
