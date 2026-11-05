# Signed one-apex responses: exact 2-adic structure (THM-4883, HYP-9476)

kind-pasteur-2026-11-05-S12. Lane L2 of the dispatch "use Session 198's rational-response coverage note as
inspiration": determine S(a, b) = {a n - b g(n) : n >= 3 odd} for every rational slope, where
g(1) = 1, g(2n) = g(n), g(2n - 1) = g(n) + n - (n mod 2). S(a, b) is exactly the set of signed responses aU - bD of
the strong unit-deletion sources (the one-apex sources), so every element has a literal tournament.

**PROVED** (`PROOF.md`), with b = 2^s b', b' odd:

- **T1.** S(a, b) contains a full residue class modulo 2^(s+1) iff 1/2 < a/b < 1 and b is not a power of two; then
  exactly one (C*), and every other class modulo 2^(s+1) meets S(a, b) in a set bounded on one side. For odd b this
  is the full parity class (Session 198 for a even, THM-4882 for a odd; re-proved here by a different route).
- **T2.** C* is the class of a n - b g(n) for every n whose O/Z word ends with the binary digits d_L ... d_1 of a/b,
  L the first position > s with d_L = 1; these n alone already fill C*.
- **T3.** S(a, b) equals C* on a positive half-line iff d_2 = ... = d_(s+1) = 0, on a negative half-line iff
  d_1 = ... = d_(s+1) = 1, with explicit attained thresholds; otherwise on neither. For b = 2 (mod 4) one side is
  always exact: e.g. S(7, 10) n (-15, infinity) = {r = 3 (mod 4)}, S(5, 6) n (-infinity, -3) = {r = 3 (mod 4)}.
- **T4.** For every nontrivial physical left source v and strong unit source q, the C2 offset
  M(v, Jq) - 2 N(v, Oq) is a signed response of q with odd denominator; if E_v < d_v the offsets fill one residue
  class modulo 2 gcd(D_v + 4U_v - 2d_v, d_v - E_v), namely U_v - D_v + E_v + 2d_v (every odd offset when
  d_v - E_v = 1, e.g. v = t, O t), else they are bounded above by an attained constant.
- **T5.** In the unit-deletion class a word moves aU - bD by a source-independent amount iff b A_w = a(2^k - 1); for
  1/2 < a/b < 1 these are exactly the powers of Session 198's clock (none for even b at any slope).

**FINITE-EXACT** (`VALIDATION_NORMAL.json` = `VALIDATION_OPTIMIZED.json`): pre-registered P1-P6 (commit 2c11d95d56).
P2-P5 confirmed (544,115 full-class targets witnessed, 186 half-line cases, six literal C2 pairs with offsets +-1
from five left sources, 59 class densities in (0.097, 0.617)); P1's structure confirmed but its per-class
attainment clause REFUTED (bounds are attained per group of classes, not per class); P6 partly REFUTED (the 8/9
witness lengths exceed the guessed bound 3j + 8 at j = 4 and j = 9).

**AUDIT:** independent, adversarial, own code (`audit/`): PASS WITH CORRECTIONS, 2 majors (a false example in T4;
T5's "iff" stated for every slope), 14 minors, 6 upgrades, all applied (MISTAKE-724).

**OPEN** (HYP-9476): (B) every one-side-bounded class has density strictly between 0 and 1 on its unbounded side;
(C) shortest witnesses have O(log |r|) letters. ((A) records the pre-registered outcomes.)

**Not claimed.** No positive Hamilton-path coverage, no C2/C3 range result, no fixed-parent children (lane L3), no
witness-size bound, no all-odd spectrum consequence.

## Files

- `PREREGISTRATION.md` - predictions P1-P6, pushed before the confirmatory run.
- `PROOF.md` - statements, proofs, outcomes.
- `verify.py` - confirmatory checker (numpy). `dependencies/` - byte-identical copies of Session 198's
  `bilinear.py` and `core.py` (SHA256 in the receipt), used only for the C2 forms and graph constructor.
- `witness_length.c` - exhaustive minimal-witness search (P6); `witness_length_3_5.txt`, `witness_length_8_9.txt`
  its confirmatory outputs; `exploration_witness_length_*.txt` - the same program on 2/3, 4/5, 6/7 (exploratory,
  not pre-registered).
- `VALIDATION_NORMAL.json`, `VALIDATION_OPTIMIZED.json` - receipts.
- `audit/` - the independent adversarial audit: `AUDIT.md`, `independent.py`, `audit_census.c`, `INDEPENDENT.json`
  (`python3 audit/independent.py --exe <compiled audit_census> --out INDEPENDENT.json`, about five minutes).

## Reproduce

    gcc -O2 -o witness_length witness_length.c
    ./witness_length 3 5 34 1024 ZOOZ > witness_length_3_5.txt
    ./witness_length 8 9 34 1024 OOOZZZ > witness_length_8_9.txt
    python3 verify.py --out VALIDATION_NORMAL.json
    python3 -O verify.py --out VALIDATION_OPTIMIZED.json
    cmp VALIDATION_NORMAL.json VALIDATION_OPTIMIZED.json

The C searches take about a minute each; `verify.py` about four minutes and 1.5 GB of memory. Do not commit the
compiled binary. `python3 verify.py --smoke --out SMOKE.json` runs every code path on small ranges in seconds.
