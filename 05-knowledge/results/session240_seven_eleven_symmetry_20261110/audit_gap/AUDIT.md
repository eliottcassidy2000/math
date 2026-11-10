# Independent audit: elevenfold gap and minimum-order 7/11 merger

Status: **PASS / independently VERIFIED exact finite inputs and proof reduction.**
Session240, 2026-11-10. No producer program was imported or executed by this auditor.

## Reviewed claims and dependencies

Reviewed `../gaps/PROOF.md` in full. The all-order reduction is already proved in
THM-5099 (`01-canon/theorems/THM-5099-value-refined-scc-recurrence-and-sharp-seven-symmetry-gap.md`)
and its Session239 `proof_bridges/PROOF.md`, section4. The cyclic profile bound
and equality characterization are THM-5094
(`01-canon/theorems/THM-5094-cyclic-orbit-profile-cost-and-sharp-sevenfold-path-floor.md`),
whose Session238 `cyclic_cost/PROOF.md` was also read. New here is the elevenfold
cover table and specialization, not a new prime-action reduction or a new value
of the inherited minimum mu11=92411. No external priority claim is audited.

For at least two moving11-orbits, the internal OCF product cannot attain equality
in a strong tournament: its equality condition makes the blocks an ordinal sum.
Parity and free action make H and mu11^s odd multiples of11. Consequently the
strict improvement is22, not merely1 or11. This yields8539792943 or larger.
For one moving orbit, fixed vertices see it uniformly, contraction is strong,
and a triangle through its contracted vertex exists. Its lift is an induced
C3[R,1,1]. Deterministic vertex insertion injects old Hamilton paths into new
ones, so the induced witness supplies a non-strict lower bound. The proof
correctly does not infer strict monotonicity or a complete equality classification.

The exact one-run/two-run decomposition gives3H+F2 with **ordered** covers.
The finite input proves the all-order bound1848759 for strong order>11 under an
actual order11 automorphism. The witness D={2,3,4,5,10}, inflated with two
singletons, is strong and supports rotation of the core, proving sharpness.
A transitive singleton padding has H92411, so omitting strongness is false.
Neither divisibility of H by11 alone nor replacing H by H/|Aut| is licensed.

## Independent algorithms and exact results

`audit.py` enumerates all32 skew connection sets. Its first independent method
is the integer walk inclusion–exclusion sieve

    H(A) = sum over nonempty S of (-1)^(n-|S|) 1^T A[S]^(n-1) 1.

Inclusion–exclusion removes walks missing any vertex. A walk with n vertex slots
which uses all n vertices is a Hamilton path. All arithmetic uses Python
arbitrary-precision integers; no floating point or rounded eigenvalues appear.

The second method deletes the first endpoint in subset DP. It obtains every
induced path count, forms the ordered complementary cover sum, and separately
counts paths in each full13-vertex inflation, without using the cover formula
for that inflated count. The independent walk sieve agrees on every core.
Small p3/p5 core and inflated cases additionally agree with exhaustive
permutation enumeration. No classification of arbitrary eleven-tournaments
is assumed. Counts of labelled specified-cycle orientations are:

- H92411, F2=1571526, inflation1848759:10.
- H93027, F2=1574958, inflation1854039:10.
- H93467, F2=1577598, inflation1857999:10.
- H95095, F2=1586310, inflation1871595:2.

Normal Python and `python3 -O` produce byte-identical output. Explicit runtime
exceptions protect gates under optimization. Deliberately wrong gamma expectations
fail in both modes; exit statuses appear in `failure_controls.txt`.

## Additional mixed77 claim

An order77 permutation has cycle lengths1,7,11,77. Below19 vertices it must have
one7-cycle and one11-cycle to have lcm77; these exhaust18 vertices. (Below18 it
is impossible.) The diagonal cyclic action is transitive on their77 cross
pairs because gcd(7,11)=1, so all cross arcs have one direction and strongness
fails. At19 the sole possible profile is(7,11,1). All cross-block arcs are
uniform, and strongness forces the three-vertex quotient to be a directed
triangle. Conversely such an inflation supports an order77 rotation. Thus19
is the exact minimum order of a strong tournament supporting this action.

For C3[A,B,1], after erasing vertex labels the allowable run words are
ABz, BzA, zAB, ABzA, BzAB, ABzAB. These yield respectively
3ab, F2(A)b, aF2(B), and F2(A)F2(B). This is a direct derivation of the stated
specialization of the inherited path-cover kernel.
The audit independently computes all8 seven-circulants and32 eleven-circulants;
the formula minimum is3049061631. Exactly6*10 choices minimize, and the two
cyclic quotient orientations give120 fixed labelled orientations for a specified
(7)(11)(1) permutation. These are not120 unlabelled isomorphism classes.

`audit19.cpp` is a separate direct19-vertex endpoint-deletion DP on
C3[Q7,R11,1]. It does not invoke the substitution formula or an input value for
H7/H11/F2. It returns3049061631 exactly. uint64_t is safe because every state
counts at most19! =121645100408832000, smaller than2^64. Compiler flags are
`-std=c++17 -O2 -Wall -Wextra -Werror`; an intentionally wrong final expectation
also fails. This verifies the minimum-order witness, not a universal all-order77
lower bound. The at19 classification plus finite cyclic inputs proves the
minimum within that order stratum. No all-order equality classification follows.

## Reproduction

From repository root, with P=05-knowledge/results/session240_seven_eleven_symmetry_20261110/audit_gap:

    python3 "$P/audit.py" > /tmp/audit240.json
    python3 -O "$P/audit.py" > /tmp/audit240O.json
    cmp /tmp/audit240.json /tmp/audit240O.json
    cmp /tmp/audit240.json "$P/output.json"
    g++ -std=c++17 -O2 -Wall -Wextra -Werror "$P/audit19.cpp" -o /tmp/audit240-19
    /tmp/audit240-19

Negative controls: both Python commands with `--inject-failure` must exit
nonzero; `/tmp/audit240-19 --inject-failure` must also exit nonzero.
