# Independent audit: exact append-only signed-section targets

**Verdict: PASS in the stated scope.** The candidate's Lemma 1, Theorem 2,
Corollary 3, infinite hostile family, and target-linear cover bound are
mathematically sound relative to the frozen physical signed-section theorem.
No mathematical correction is required. Promotion from “pending independent
audit” is justified for these claims only.

This is a theorem about descendants of a supplied actual source. It does not
prove all-target occupancy, signed-to-physical surjectivity, or the full
positive C1 spectrum. The candidate and README preserve this boundary.

## Audited inputs and preservation

The audit reviewed `PROOF.md`, `README.md`, `verify.py`, both author receipts,
and `dependencies/SIGNED-SECTION.md` in the sibling
`session202-fiber-intersection` directory. Exact SHA256 fingerprints appear in
`VERIFIED-INPUT-HASHES.json`. In particular:

- Candidate proof: `13d883e54eb1b902b7bede93e8536f233637608b2147c2804d3fcd5e306df4cb`.
- Frozen dependency: `92b300d272d529bb21a504a914226b1245a631249d30d22d0e9362625a698182`.
- Author verifier: `d2b22358abcfd374064e6b533ec4384ca22ad8620f2042526fa75825ddd160a7`.

The original proof, README, and dependency fingerprints were recorded before
testing and verified again afterward. No audited file was edited. All audit
artifacts are in this separate directory. The shared repository's existing
unrelated changes were left untouched; this audit made no GitHub writes.

The physical section theorem is an explicitly inherited dependency, not
reproved from tournament graphs here. I checked that the portions imported
by the candidate really occur in its proof, including the forced terminal
letters and the order in which forgotten deletion coordinates are recovered.

## Proof review

### 1. The physical carrier and finite search

I independently transformed the recurrence to nonnegative coordinates
`(a,f,e,p)`. The three operators are

    O(a,f,e,p) = (a+e, a+2f+e+2p, 0, e+p),
    Z(a,f,e,p) = (2a+f+e+p, f, e+p, 0),
    J(a,f,e,p) = (f+p, 0, a+e+p, e).

The initial state is `(0,0,0,1)`. These formulae both verify the stated
behavior of f and supply an implementation independent of the author's raw
`(D,E,U,d)` recurrence. All actual states have `d=e+p>=1`, `D>=d`, and
`U>=D`. Both D and d are nondecreasing. U strictly increases under each
operation, so every height-capped descendant enumeration is finite and
complete when it retains all states reached from the supplied source.
Deduplication by the four-coordinate profile is sound because all future
operations and slice tests depend only on that profile.

### 2. Balanced descendants and the exact classification

A balanced actual state has `U=D`, hence `f=p=0` and `E=d`. For every
reachable v, `v.OJ` is balanced with endpoint count `D_v`. When `f_s=0`,
`s.Z` is balanced with endpoint count `d_s`.

Conversely, fix a presentation of a balanced descendant and inspect its
last appended J, if one exists. That J sets deletion to the D-value of its
actual prefix; later O and Z preserve it. The final balance condition makes
the endpoint equal to this deletion. If no J is appended, an appended O
would make f positive permanently along the remaining Z letters, so it is
impossible. The remaining deletion is exactly `d_s`. This also covers the
empty word and proves Lemma 1 without replacing actual sources by cone
points.

For the section classification, the frozen theorem supplies the exact
predecessor facts needed by the candidate:

- Every exceptional slice point ends in J, with a balanced predecessor
  `(3n,n,3n,n)`.
- Every canonical slice point has three forced terminal J letters, with
  a balanced third predecessor `(m,n,m,n)`.

If the appended word contains all required forced letters, the balanced
predecessor is still a descendant of the supplied source, and Lemma 1
applies. If the prefix boundary falls inside those letters, the entire
appended word is `J^j`, `0<=j<=2`, precisely the set C(s). This resolves the
potential prefix-boundary error: a balanced predecessor must not be claimed
to descend from s after peeling past s. The candidate handles that case
correctly.

For the forward inclusion, W realizes every D-tail value, balanced
endpoint counts are realized by JJJ, and C(s) consists of literal actual
extensions. Both inclusions therefore hold at all heights.

### 3. Sharp minima and short faces

When `f_s>0`, s itself cannot lie in the slice. Nor can s be the immediate
predecessor of a slice point, since both canonical and exceptional immediate
predecessors have f zero. In the two-letter case, the canonical second
predecessor likewise has f zero. The sole remaining case is an exceptional
`s.JJ`. Then `s.J` is balanced, so `E_s=0` from
`U(s.J)-D(s.J)=E_s`; its endpoint, and hence the target, is `D_s`. Thus
`C(s)` adds nothing to the D-tail in this chamber. Monotonicity of D and
W(s) prove the exact minimum `D_s`, not merely a lower bound.

For `f_s=0`, `D_s>=d_s` removes the D-tail from the minimum calculation.
Direct substitution gives exactly the three tests and readouts in (6).
Positivity is asserted only for actual profiles. This restriction matters:
the same formal hyperplanes do contain nonphysical points with nonpositive
readout.

Both extra pieces in the general formula are genuinely needed. Independent
controls found:

- `s=t.OJ=(5,1,5,1)`: target 1 comes from the optional deletion term and is
  absent from both its D-tail and C(s).
- `s=t.ZJ=(5,2,6,3)`: target 1 is a short completion (already in the slice),
  below its deletion 3 and its D-tail minimum 5.

### 4. Infinite hostile family

The recurrence gives `b=t.Z=(3,1,3,1)` and `b.OJ=3b`. Therefore

    t.Z.(OJ)^r.OZ = 3^r (9,1,11,1).

Its f coordinate is `2*3^r`, its signed target is `3*3^r`, and its minimum
reachable section target is `9*3^r`. No appended word can preserve the
signed target while reaching the slice. Independently, `b.(OJ)^r.W`
is an actual section witness for that same target. This is a real
append-only obstruction at already represented targets, not an argument
from an unrepresented target or a finite failed search.

### 5. Complete finite bound

For a target `n>D_s`, a word realizing `D=n` cannot consist solely of O
letters. Strip only terminal O letters from that appended word. The
remaining core is still a descendant of s and ends in Z or J. Its D-value
is unchanged. For a Z output, `2D-U>=3d>=3`. For a J output,

    2D_new-U_new = D_old-2E_old+U_old >= 3.

Indeed `D_old-2E_old=a+p>=0`; if `U_old>=3` this suffices, and the unique
actual state of height 2 is t, where the expression is 3. Hence
`U_core<=2n-3`. Every prefix of that core also satisfies its height cap,
because U strictly increases. This verifies the completeness of the D-tail
decision procedure even if the supplied source has a long O history.

For the section bound:

- If `n=D_s`, W(s) has height `14U_s+6n`.
- If `n>D_s`, applying W to the bounded core gives at most `34n-42`.
- The optional deletion witness `s.Z.JJJ` has exact height
  `14U_s-8d_s`.
- Every short completion has height at most `4U_s`.

Thus `max(34N-42,14U_s+6N)` covers all positive targets at most N in the
append-only spectrum, including N=1. The argument proves existence of a
bounded witness for an already reachable target. It does not assert
reachability of any additional target.

## Independent finite-exact checks

`independent_check.py` was written before reading the author's verifier.
Its primary recurrence and slice predicate use `(a,f,e,p)` coordinates.
It imports no author code and uses explicit exceptions, so optimization
does not disable the tests.

The ordinary and optimized receipts are byte-identical and report PASS:

- Complete geometry bank: 33,691 actual profiles with `U<=2048`.
- 101,073 outgoing operation checks, including strict U growth,
  monotonicity, both core inequalities, and recurrence conjugacy.
- All three short-face equivalences on the complete bank; hit counts are
  488, 1,109, and 2,044 respectively.
- All 555 actual sources with `U<=128`, each checked against every integer
  target `1..63`: 34,965 exact source/target decisions.
- Separate forward slice banks and D-tail banks, with a total of 368,961
  forward profiles examined. Each source's forward cap is
  `B(s,max(63,predicted minimum))`, permitting a direct minimum check too.
- Balanced endpoint sets checked independently on bounded descendant banks.
- Both chamber minima checked for every one of these 555 sources;
  307 sources have positive f.
- Direct construction and height verification for every right-hand-side
  contribution in the declared target windows.
- Twelve source prefixes additionally checked at twice the claimed
  target-window height bound, with no extra low-target section values.
- Exact hostile formulas checked for `r=0..40`; complete bounded descendant
  searches also checked the first three hostile sources.
- Full raw-coordinate formal basis checks for the recurrence conjugacy and
  W, plus a full basis of the f-zero plane for the deletion construction.

These computations are finite-exact evidence. The all-height conclusions
rest on the proof review above, not extrapolation from these banks. In
particular, completeness of the finite decision windows uses the separately
reviewed core and section-height arguments.

## Exact author replay

Both `python3 verify.py` and `python3 -O verify.py` were rerun from the
candidate directory. Their outputs match each other and both frozen author
receipts byte-for-byte, with SHA256
`6a5ceebbd7f66e2d58a21039c566d87efe1a55f140abe90ba4a54b84a326befa`.

That replay includes the complete author bank of 59,723 profiles through
`U<=3000`, ten supplied prefixes with complete target windows `1..151`,
36 formal-basis checks, and hostile formulas `r=0..30`. Source inspection
confirms that the checks survive `python -O` and no cone-only substitution
or unproved height truncation is used.

## Reproduction

From the common workspace directory:

    python3 session202-append-section-audit/independent_check.py
    python3 -O session202-append-section-audit/independent_check.py
    (cd session202-fiber-intersection && python3 verify.py)
    (cd session202-fiber-intersection && python3 -O verify.py)

The four saved receipts and their fingerprints are included in this audit
bundle. No mandatory revision or unresolved mathematical blocker remains.
