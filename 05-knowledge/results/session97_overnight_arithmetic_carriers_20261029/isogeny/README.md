# All 66 high-rank snapshot curves lie in distinct Q-isogeny classes

**Status: FINITE-EXACT** (2026-10-29, checked at 07:07 UTC).
This closes a precise gap: distinct rational j-invariants rule out quadratic
(or other geometric) twists, but do **not** rule out Q-isogenies.

## Result

For each of the **2,145 unordered pairs** in the session's 66-curve rank-lower-bound
≥30 snapshot, we exhibit a common good prime p≤71 with unequal Frobenius traces.
Consequently **all 66 curves are pairwise non-Q-isogenous**. There are no
unseparated pairs. This result uses neither BSD nor GRH, and does not depend
on the correctness of the database rank labels or independence of its points.

The last two witnesses in the increasing-prime search are:

| Curve IDs | Prime | Traces |
|---|---:|---|
|398, 5990|71|−15, −16|
|3792, 5990|71|−9, −16|

This is an exact finite separation certificate, not a claim that 71 is the
smallest prime cap under some alternative model, nor does this initial scan alone establish fingerprint minimality; the
last section proves the minimum subset within the stated finite bank. Bad reduction of the supplied integral model is skipped; this can
skip a prime where a different model has good reduction, but cannot create
false separation.

## Why the certificate proves non-isogeny

If two elliptic curves over Q are Q-isogenous, then for every auxiliary prime
ℓ their rational ℓ-adic Tate modules are isomorphic as Galois representations:
an isogeny and its dual compose to a nonzero integer, which is invertible in
Q_ℓ. At a common good prime p≠ℓ, their Frobenius characteristic polynomials
therefore agree. In particular their traces a_p, and hence their point counts
p+1−a_p, agree. A single unequal trace at a common good prime is a rigorous
obstruction to Q-isogeny. This implication does not require a converse theorem.

The producer computes a_p by completing the square in y and using Euler's
criterion. The independent verifier enumerates **every (x,y)∈F_p²**, adds
the point at infinity, and checks smoothness using both partial derivatives.
For a singular Weierstrass cubic the unique singular geometric point is
rational over F_p; the point at infinity is always smooth. Thus this brute
check independently verifies each good/bad flag without importing the
producer's discriminant routine. Both methods agree for all recorded rows.

### Hostile control

The models y²=x³+x²+x and y²=x³−2x²−3x have distinct j-values 2048/3 and
35152/9 but are 2-isogenous. The verifier checks equal point counts at every
common good prime in the bank. Their standard 2-isogeny comes from the
rational 2-torsion point (0,0): for y²=x³+a x²+b x the quotient has
coefficients −2a and a²−4b. This is a concrete control against confusing
distinct j with non-isogeny.

## Information carrier and family interpretation

- **Source:** the 66 named integral Weierstrass models.
- **Target:** good/bad flags and exact Frobenius traces through p=71.
- **Preserved decision:** an unequal-trace witness certifies non-Q-isogeny.
- **Lost:** rational points, Mordell–Weil lattice, family parameters, generic
  ranks, and all local information beyond the tested primes.
- **Sidecar:** ordered source IDs, pinned input hash, and the upper-triangular
  pair-witness prime matrix.
- **Boundary:** equal finite fingerprints never certify isogeny. Distinct
  Q-isogeny classes need not belong to distinct algebraic parameter families.
  Indeed one nonisotrivial family can contain many different isogeny classes.

Thus the 66 examples are not merely 66 isogenous presentations of fewer
rational elliptic curves. But this does not refute an eventual rank≤21
heuristic, establish infinitely many high-rank curves, or identify a theorem
relating rank21 to tournament exclusions {7,21}.

## Provenance and reproduction

Input: `../ranks/phase2/highrank_inputs.json`, SHA256
`8dffce8a9c53ee3a757c18ef4e7449d972d4eaf96139ea4a89bd519ed586903d`.
The sibling input is the compact mathematical snapshot already extracted in
this session from <https://elliptic-rank.icarm.cloud/database.json>.
This research made use of the Elliptic Curve Rank Leaderboard, maintained by
ICARM and supported by NSF Grant DMS 2425401. No live refetch was made for
this computation. Session baseline source pin:
`5d6eb91dd76eee86b91e04d856daa32c32b3d61d`.
The new files do not modify that baseline or earlier checkpoints.

Run:

```sh
python build_certificate.py
python verify_certificate.py
python -O verify_certificate.py
```

`certificate.json` records IDs, the trace bank and every pair's first witness
prime. Both normal and optimized verification pass identically; checks do
not rely on Python assertions.

### Optional full-database integrity check

`check_database_integrity.py` uses the existing local raw database, which is
**not redistributed**. It checks exact rational equation membership for all
**53,498 listed points on 6,127 records**, with no failures, and verifies
6,127 distinct record IDs. `database_integrity.json` pins the raw source hash.
This is only a membership/ID audit: it establishes neither independence,
saturation, rank, dataset completeness, nor uniqueness of isomorphism classes.
The raw-snapshot checker and raw database are intentionally omitted from the
minimal publication packet. The optional local script needs the retained `../ranks/database.json`;
the published non-isogeny certificate does not.

## Minimal decision bank inside this prime window

**FINITE-EXACT:** among the recorded primes ≤71 the unique minimum bank is

`{19,29,31,37,43,47,53,59,61,67,71}` (11 primes).

`minimum_prime_bank.py` independently exhausts subsets in increasing size
(230,949 tested through the first cover). More simply, **each of these 11
primes is forced**: the following pair has no other common-good unequal-trace
witness anywhere in this window.

|Forced prime|Curve IDs|
|---:|---|
|19|273,4031|
|29|273,5981|
|31|3545,3872|
|37|3567,5981|
|43|4032,5984|
|47|4239,5981|
|53|3545,5981|
|59|273,3587|
|61|273,3545|
|67|4239,5984|
|71|398,5990|

Those 11 primes together cover every pair, so the forced witnesses prove
both minimality and uniqueness **within the finite bank**. This gives a
short lower-bound certificate instead of relying on a black-box optimizer.
Prime3 covers no pair. No bad-status mismatch is used as a separator, and
no optimality claim is made allowing primes outside the recorded window.
