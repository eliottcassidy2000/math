# Independent audit: general cyclic certificates and minimum-order row hostile

**Verdict: ACCEPT.** The general positivity and negative multiplier-chamber
criteria are valid sufficient theorems for actual nonempty same-order module
pairs. The strong seven-vertex same-card hostile is genuine, and its minimum
module order is independently established by a complete finite census through
order six. Failure of either sufficient test remains inconclusive.

General arbitrary strong-quotient scalar zero-rigidity remains **OPEN**. This
packet supplies no necessary sign criterion, no universal eventual-sign
algorithm, and no actual counterexample to scalar zero-rigidity.

## Accepted scope and proof audit

Let `d=f_M1-f_M0` for two actual nonempty tournaments of the same order `n`.
The modules need not be strong or share a deletion card. The ordered-cover
profiles have positive coordinates through their order, vanish afterward,
and have zero coordinate at index zero. For any nonempty actual side modules
`A,B`, the outer quotient is the directed triangle with its complete uniform
inter-module incidences. The comparison changes only the first module:

`H(C3[M1,A,B])-H(C3[M0,A,B]) = f_A^T K f_B`.

The kernel on coordinates `1,...,n+1` has diagonal
`d_(r-1)+3d_r+d_(r+1)` and adjacent entries `d_r+d_(r+1)`.
The factor three counts the three possible initial phases of a linear run
word with equal module multiplicities; these are not Hamiltonian cycles.
Other allowed run multiplicities have one initial phase. Ordering the paths
within each module already accounts for all choices, so no extra factorial
belongs in the kernel. Nonempty sides eliminate the zero-index terms.
Truncating side coordinates above `n+1` is exact, including when a side/base
is larger than the matrix.

### Certificate P

Entrywise `K>=0` implies nonnegative response. If `K11>0`, retain its term:

`response >= K11 H(A)H(B)>0`.

The coefficient is sharp at singleton sides. Within this entrywise
nonnegative class, strict positivity for every nonempty actual pair also
requires `K11>0`, because singleton sides isolate this entry. Entrywise
nonnegativity itself is **not** proved necessary over the realizable sides.

### Certificate N

The supplied base must be an actual nonempty tournament `V`, with positive
log-concave profile `v` on `1,...,N`. A pivot `c` lies inside this support.
All rows before it are entrywise nonnegative; all rows after it are entrywise
nonpositive, with zero padding when necessary. The pivot row is unrestricted.
The strict hypothesis is `v^T K v<0`.

The conclusion applies to every *actual* pair satisfying
`P_A=P_V R_A` and `P_B=P_V R_B`, where each `R` is a nonzero polynomial with
nonnegative **monomial** coefficients. It does not assert realization of
arbitrary formal products or apply merely from nonnegative coefficients in
some other basis. Ordinal additional tournament factors are legitimate by
multiplicativity and the nonnegative colored odd-cycle-family expansion.
An eventual-power certificate requires a supplied nonempty actual base
`V=TT_k[F]`; there is no assertion that any successful `k` must exist.

The multiplier proof is correct. In the binomial basis, multiplication by
`m` acts by `(Lv)_r=r(v_r+v_(r-1))`. The three summands in the expanded
adjacent minor for convolution by `(1,1)` are nonnegative by log-concavity;
zero boundaries preserve the inequality. Pointwise multiplication by `r`
also preserves it. The ratios `(Lv)_r/v_r` increase, and iterated ratios
increase on the common original support. Their nonnegative linear
combination gives increasing `w_r/v_r`, strictly positive original support,
and a nonnegative tail. No log-concavity of the final mixture is assumed.

For every nonnegative vector `w`, the row condition makes `Kw` have only a
positive-to-negative crossing, and no positive coordinate beyond the base
support. If `v·h<0`, take the scale `w_c/v_c` at the last positive coordinate
of `h`; positive contributions can only decrease relative to the scaled
base sum and negative contributions can only increase in magnitude. The
tail contributes nonpositively. If there is no positive coordinate, original
support positivity supplies strictness directly. Apply this first to `Kv`,
then to `Kf_A`, using symmetry. This proves strict negativity for both sides
without assuming that either resulting side profile is log-concave.

## Actual hostile and sharper small-side control

The literal shared card is `(32,1,3,7,15,30)`. Adding the new vertex at cuts
`5` and `14` gives the two displayed strong seven-vertex modules. Independently
computed full profiles are

- old: `(51,678,3798,11280,18240,15120,5040)`;
- new: `(39,678,4038,11760,18480,15120,5040)`;
- difference: `(-12,0,240,480,240,0,0)`.

Both anchored-partition DP and literal enumeration of `7!` vertex words per
module agree. For the actual base `K1`, `v=(1)` is log-concave and the base
quadratic is `-36`. But `K22=228>0`; the only allowed pivot is `c=1`, so the
row condition fails. Direct Hamilton-subset counts give:

| Two equal side modules | Old H | New H | Response |
|---|---:|---:|---:|
| `K1` | 831 | 795 | -36 |
| `TT2` | 27279 | 28107 | +828 |
| `C3` | 1408455 | 1476387 | +67932 |

The accepted claimed `C3` hostile has multiplier `m^2+2` from the singleton
base. The additional `TT2` control has the even simpler multiplier `m`.
Thus side order two already witnesses failure if row crossing is dropped;
no change to the minimum module-order argument is needed. Mixed side pairs
were also checked: `(K1,TT2)` gives `-60`, `(K1,C3)` gives `-180`, and
`(TT2,C3)` gives `+5364`, symmetrically. These are finite controls, not a
classification of all side tournaments or a claim that any response is zero.

## Independent minimum-order certificate

`enumerate_profiles.cpp` enumerates all labelled tournaments of orders one
through six. It does **not** use the producer's permutation/descent method.
It extends paths forward by an outgoing edge to count Hamilton paths on
subsets, then takes an unrestricted first-block convolution of subsets to
count ordered path covers directly. It does not use an OCF interpolation,
an isomorphism representative list, or a presumed profile classification.

| Order | Labelled tournaments | Distinct profiles | Unordered profile pairs |
|---:|---:|---:|---:|
| 1 | 1 | 1 | 0 |
| 2 | 2 | 1 | 0 |
| 3 | 8 | 2 | 1 |
| 4 | 64 | 3 | 3 |
| 5 | 1024 | 8 | 28 |
| 6 | 32768 | 24 | 276 |
| Total | **33867** | **39** | **308** |

All 39 order-specific profiles are log-concave. All 308 kernels are
entrywise one-signed. Equal profiles give the zero kernel, so need not be
listed; reversing a pair only negates its kernel. Thus this covers every
ordered pair of actual same-order modules through order six, a stronger
universe than just strong pairs or same-card pairs.

If such a kernel has a negative base quadratic at a positive-support base,
it must be nonpositive and possess a strictly negative entry meeting that
base support. Any admissible multiplier preserves strict positivity on that
support. Hence the subsequent bilinear response stays strictly negative:
no smaller-module negative-to-positive propagation hostile exists for any
actual base or admissible multiplier. This is the precise finite-exact
minimality assertion. It does not claim a census at order seven or a
classification of all minimal witnesses.

The producer's claimed `23,717,429` literal words equals
`sum_(n=1..6) 2^(n choose 2) n!`; our census validates that universe size but
uses DP instead of iterating those words. It must not be described as a
second independent literal-permutation traversal.

## Executed tests and reproducibility

The independent C++ census passes both `-O2` and UBSan; outputs are byte
identical, SHA256
`fa2bc93ac417540b361b67ad1966bbdde07d35b86b49d32a8be3493f72e0fe7e`.
The retained JSON contains every one of the 39 profiles, not just totals.
The Python checker independently rebuilds 14,623 kernel entries for all
308 small profile pairs from literal forced cyclic run words. It also checks
nine direct actual hostile side-pair contexts, four actual certifier controls
(including a base larger than the matrix), and eight invalid-actual-input plus three invalid-formal-input rejections.

For API conformance, 14,400 formal finite cases are checked; 1,642 satisfy
the certificate and 26,272 multiplier-pair controls stay strictly negative.
These vectors are explicitly **formal**, not claimed realizable tournaments.
The producer certifier is imported only as the API under test; expected
profiles, direct H counts, and acceptance predicates are independently
calculated. The low-level interfaces retain their realizability boundary.

Normal and `python -O` runs match, SHA256
`98335b35eae1e6669c7e6aa40f0a2bc6e6cc8e8fe06ec1ca234d2c3d2a60dc9e`.
All four independent scientific replays and both independent C++ builds exit zero with empty stderr.
Timings/commands are retained in `enumeration_receipts.json` and
`python_receipts.json`; O2 census 0.154 s, UBSan 0.731 s, Python normal 1.253 s,
optimized 1.144 s in this environment. These are observations, not promises.

Preserve the sibling directory layout: the checker imports the frozen
independent graph/profile kernel `../independent/audit_profile_dichotomy.py`
and tests `../zero_framework_later/certify.py` with its `profile_tools.py`.
From this audit directory:

```sh
g++ -std=c++17 -Wall -Wextra -Werror -O2 enumerate_profiles.cpp -o /tmp/s138_framework_profiles
/tmp/s138_framework_profiles > profiles_O2.json
g++ -std=c++17 -Wall -Wextra -Werror -O1 -fsanitize=undefined -fno-sanitize-recover=all enumerate_profiles.cpp -o /tmp/s138_framework_profiles_ubsan
/tmp/s138_framework_profiles_ubsan > profiles_UBSan.json
cmp profiles_O2.json profiles_UBSan.json
python3 -B audit_framework.py > normal.json
python3 -B -O audit_framework.py > optimized.json
cmp normal.json optimized.json
```

Use a writable portable copy for replays, leaving frozen artifacts untouched.
Source pins and producer manifest validation are recorded in `SOURCE_PINS.json`.
The failed order-six discovery precursor was correctly retired under the
producer's excluded `exploratory/` path; it is not a valid reproducer of the
accepted order-seven witness. No theorem repair was required in this audit.

## Final freeze and additional source replay

The accepted producer manuscript SHA256 is
`29b85ded2a59a44bc862e33b316d284968c714d9b79406892e2f457ed6e4dfc9`.
Its `public_manifest.json` SHA256 is
`84ea98c7ceb5d7ec18d628bbd2a5251d83b222bce22b44108a4f1b3bc941ff23`;
its `SHA256SUMS` SHA256 is
`10f2c2de720824dfb9f0e86f07504f11680d6133f31ef8617bbbb64f26102a42`.
All 11 curated payload records (27,740 bytes) and all 12 checksum entries
match. Include the two manifests as well as the 11 listed payloads. Final
integer-input guards and the TT2 witness were checked against this freeze.

In addition to the independent tests above, the final producer Python audit
was replayed in normal and optimized modes and its literal C++ census in O2
and UBSan. All four outputs exactly match retained producer outputs, all
exit zero, and all stderr is empty. The two builds also pass. These four
are explicitly **same-source replays**, not four new independent algorithms.
Commands, terminal outcomes, hashes, and timings are in
`producer_replay_receipts.json`. Combining the two replay blocks gives
eight scientific replays and four C++ builds, all passing.

Both earlier frozen producer packets still match their previously recorded
manifest pins and every checksum entry. No frozen producer file, prior
checkpoint audit, Git state, or publication file was changed by this audit.
An unrelated exact Pfaffian hostile is outside this audit's scope and must
not inherit this acceptance.
