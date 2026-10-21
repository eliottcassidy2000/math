# Independent review: route-free residual suffix

Status: **PASS, PROVED conditional suffix lemma; FINITE-EXACT bounded native applications.** The uniform assertion for every N>=4096 remains **OPEN**. This review neither establishes universal entry nor assumes that pending native routes terminate. Repository was read-only; all new review files are in this directory. No Git or remote action was performed.

## Exact scope and inherited carrier

The carrier is the finite source universe S_(64N), the frozen 627-point frontier V, and the original source charge A_(64N)(s). The inherited definition of arrived energy, not the full moving-endpoint Gram energy, is used. Source norms and compatible mod15/mod35 registers are not a substitute for actual native prefixes.

Read dependencies:
- `math/05-knowledge/results/native_cancellation_105_20261018g.md`, sections 1–3 (native charge alphabet, Nyquist transport, fixed-source scale transport, source budget).
- `math/05-knowledge/results/session18_native_annulus/README.md`, Theorem A and its ownership discussion (exact cut ledger and nonarrival convention).
- `math/05-knowledge/results/cross_scale_frontier_credits_20261021/check.py` (reviewed, not imported by the independent checker; extraction conventions).
- Candidate `session22-residual/routefree_probe.py` (reviewed, not executed or imported).

In particular, actual arrival contributes exactly (-1)^t A_(64N)(s). There is no additional rotation along a route at this Nyquist consumer. Rotations in other source/scale lift identities concern changing labels or scales, not hidden transport along the fixed native route. The full integer clock is still authenticated; replacing it by a parity-only certificate would lose an essential native sidecar even though this energy uses only its sign.

The candidate exploratory script calls an eventual-route resolver for every source before finding its cut. Its numbers are correct for independently checked cases, but that script by itself is not a bounded-prefix certificate. The independent `audit.py` replaces this with a genuine cutoff computation.

## Reviewed theorem and proof

Fix a finite cut d0 at which every old source has actually arrived in V. Let W_u be its now-fixed old vector, H_u the new arrived vector at d0, and T_u=W_u+H_u. Include all u in V, with T_u=0 for presently empty terminals. Set

    Delta0 = sum_u [Q(T_u)-Q(W_u)].

Let P be the distinct new sources not yet arrived at d0, with original unsigned charge a_s=A_(64N)(s). Zero charges may be retained or omitted. Define

    D = max_{epsilon_s in {-1,1}} Q(sum_{s in P} epsilon_s a_s),
    L = sum_{s in P} max_{u in V} |beta(T_u,a_s)|.

Then, for every finite d>=d0,

    Delta_N(d) <= Delta0 + D + L.

This holds without knowing any pending source's future terminal, arrival time, sign, or eventual termination. The case d=d0 is included. Nonarrival forever is allowed.

Proof. At d, some subset A of P has arrived. Partition A by its actual terminal. For each terminal let X_u be the sum of its signed charges, where each source has the actual sign (-1)^t_s. The old vectors are constant, so exact polarization gives

    Delta_N(d)-Delta0 = sum_u Q(X_u) + sum_u beta(T_u,X_u).

The linear sum is at most L, by expanding into distinct source contributions, taking absolute values, and then extending the nonnegative sum from A to P.

For the quadratic sum, independently choose fair signs eta_u for the nonempty terminal groups. Independently choose fair signs xi_s for every source in P\A. The vector

    Z = sum_u eta_u X_u + sum_{s in P\A} xi_s a_s

is always an allowed full individual-sign choice in D: each arrived source receives its actual sign times eta_u, and each excluded source receives xi_s. No claim is made that these artificial signs or assignments are dynamically realizable. They merely enlarge the comparison domain, which is lawful for an upper bound. Bilinearity and the independent mean-zero signs give

    E Q(Z) = sum_u Q(X_u) + sum_{s in P\A} Q(a_s) <= D.

Discarding the latter nonnegative term proves sum Q(X_u)<=D. This also explicitly handles never-arriving sources; a terminal or clock need not be invented for them. Combining the bounds proves the theorem.

The same proof works for any positive-semidefinite quadratic form and any finite charge set. The native alphabet is needed only for fast exact D and the coarse corollary.

## Exact twelve-pattern maximum

The nonzero native alphabet in cyclic Euclidean angular order is

    (1,0), (2,1), (1,1), (1,2), (0,1), (-1,1),
    (-1,0), (-2,-1), (-1,-1), (-1,-2), (0,-1), (1,-1).

Write c_j for the multiplicity of vector v_j and S=sum c_j v_j. The convex hull of all individual-sign sums is the zonotope sum_s [-a_s,a_s]. A convex quadratic attains a maximum at one of its vertices. Every vertex is exposed by a generic linear functional, which assigns all parallel copies a common sign. Generic functionals select six consecutive alphabet directions as positive, with the opposite six negative. Thus

    D = max_{i=0,...,11} Q(2 sum_{j=i}^{i+5} c_(j mod12) v_(j mod12) - S).

Missing directions only duplicate candidate vertices; antipodal directions and multiplicities are handled exactly. Zero charges do nothing. No heuristic angular discretization is involved. The independent audit exhaustively enumerates all 2^n individual-sign choices for 220 deterministic examples with n=0,...,10, including repeated directions and zeros, and agrees with this twelve-pattern expression.

If m is the number of nonzero pending charges and C=max_u Q(T_u), the native bound Q(a_s)<=3 gives D<=3m^2 and L<=2m sqrt(3C), hence the coarser sufficient bound Delta0+3m^2+2m sqrt(3C). The exact D+L expression is the stronger result audited here.

## Extraction and nonreuse

At any cut let A be the arrived new-source norm, M=153N-A, and I=Delta-A the original cross interaction. An extracted bank J' must be an actual signed sum of named unordered pair contributions, each counted once; the remaining ledger is R'=I-J'. Then identically

    Delta<=153N  iff  R'<=M-J'.

The suffix theorem therefore implies R'<=M-J' whenever its displayed bound is at most153N. It does not authorize adding a previously extracted credit to the bound a second time. If J'<=0, the stronger sufficient target R'<=M can fail even when the desired original inequality holds; replacing the exact criterion by that stronger target is a loss of sharpness, not a contradiction. Norms and cross credits are different ledger objects.

The suffix proof itself requires no extraction and spends no source collateral twice. Its ownership requirement is essential: a single pending charge (1,0) gives D=1, but illegitimately assigning that same source twice in a terminal gives Q(2,0)=4. The checker includes this hostile and rejects duplicated source rows. Auditing a specific expanded extracted bank would additionally require its pair-address ledger; this review does not certify arbitrary external banks.

## Independent bounded-prefix native computation

`audit.py` imports no producer module and reads no saved producer route table. It computes charges independently by scanning original odd integers and stripping each to its primitive B-ancestor, rather than summing each source's B-chain. It authenticates the frontier hash and its descending rank edges.

For each source it performs at most d0 native steps, or joins a previously certified actual prefix whose combined first-hit depth is at most d0. Failed cutoff searches remain pending; no continuation is attempted for that source beyond the cutoff. Cached prefixes are actual first hits bounded by d0, not eventual assumptions. Caching can reuse a short already-certified tail but cannot accept an arrival whose combined depth exceeds the cutoff. Every old charge is compared with a separately recomputed old-scale charge. Every new source norm is retained in the exact budget, including zero charges. The new norm sum is checked equal to153N.

The audit regenerates all native source IDs at each tested scale, hashes arrived rows including full clocks and both registers, and preserves pending nonzero source IDs and charge multiplicities in `receipt.json`. It independently replays a representative row set and every arrival exactly at the cutoff. The entire census is built from its own native step implementation, not merely checked via the representative subset.

Verified numerical applications:
- N=4096,d0=105: all old sources arrived by67; 386 nonzero new sources still pending. Delta0=421276, D=91876, L=104672. The bound617824 is below153N=626688 by8864. Hence every finite d>=105 passes, irrespective of what the pending sources do.
- N=8192,d0=112: all old sources arrived by86; 615 nonzero new sources pending. Delta0=793287, D=229957, L=188427. The bound1211671 is below153N=1253376 by41705. Hence every finite d>=112 passes.
- N=7,d0=8: Delta0=1084>1071 already fails. The suffix envelope1415 does not claim success.
- N=7,d0=22: no pending nonzero sources remain and Delta0=1097>1071; the exact envelope stays1097. The genuine hostile is not erased.

Previous-cut tests d104 and d111 are also independently regenerated: these envelopes are641801>626688 and1262521>1253376, respectively. This only says these preceding sufficient envelopes do not pass; it does not say the underlying inequality fails there, or establish a globally earliest possible certificate.

Controls: full-clock corruption by +2 is rejected despite preserving the sign; corrupting mod15 by +5 or mod35 by +5 is rejected even though modulo5 compatibility survives; repeated source ownership is rejected. Algebra controls include 220 arbitrary partial arrivals with arbitrary signs and terminal assignments, zero charges, currently empty terminals, omission/nonarrival, and1000 exact extraction-criterion checks.

Reproduce from the workspace root:

    python3 session22-residual-review/audit.py

The lemma is uniform in its explicitly stated hypotheses. Its two successful scale certificates are not a uniform all-N>=4096 theorem. No conclusion is drawn about earlier cuts from the suffix alone, about full moving-state Gram energy, or about universal Collatz entry. Producer claims at N=16385 and32768 were inspected as proposed values but are not independently recertified by this review's default run.

## Addendum: general pending-old gate (independently reviewed)

**PASS, PROVED.** The final producer's extension removes the old-completion hypothesis. Let z_o,u and z_n,u be future old/new contributions from sources pending at the cut. Direct subtraction of the two energies cancels the old quadratic term and yields

    Delta(d)-Delta0 = sum_u [Q(z_n,u)+beta(T_u,z_n,u)
                      +beta(H_u,z_o,u)+beta(z_o,u,z_n,u)].

Let D_o and D_n be the signed-zonotope maxima for the pending old and new sources, L_n the support against T, and L_o the support against H. The previously proved missing-destination lemma applies independently to each source class. Cauchy–Schwarz in the direct sum over terminals therefore gives

    sum_u beta(z_o,u,z_n,u) <= 2 sqrt(D_o D_n),
    Delta(d) <= Delta0+D_n+L_n+L_o+2 sqrt(D_o D_n).

Replacing the final term by ceil_sqrt(4 D_o D_n) is a safe integer upper bound. All old self-energy cancels, so no extra D_o term is needed. Crucially the old linear support is against H, not T or W. For W=0,H=(1,0), one pending old charge (1,0), and no pending new source, Delta rises from1 to3; omitting L_o would be false.

The independent checker now exhausts all2401 possibilities for a fixed two-old/two-new source system, allowing each source absence or either sign in any of three terminals (one initially empty). It directly computes both energies, verifies the expansion, and verifies the general bound. These are algebraic controls, not assertions that arbitrary synthetic assignments occur in native dynamics. Actual native finite applications remain the same two successful old-completed scales documented above.

The latest producer `check.py` was inspected independently. Its source charges reconstruct the three-row formula; its old charge equality and full new source-norm153N checks are correct. Native edge valuation is1 or2 for the primitive-source domain. It may cache verified suffixes longer than h by concatenating a bounded fresh walk with a previously verified suffix; it correctly returns them as pending when their combined depth exceeds h. This is safe: it does not assume any unproved route terminates, and the gate reads only unsigned charges for sources classified pending. The independent checker here is stricter and does not cache such over-cutoff concatenations.

The producer receipt's moment-hostile last row is [149,9629,433053,-8039,626688,425014]. This correctly refutes Y<=105N at N4096 because433053>430080. That receipt and its arithmetic were inspected here; this review's native recomputation does not claim a separate d149 census.

### Optimized replay and pins

All assertions in this independent checker were converted to explicit non-removable `require` calls before the final replay. Normal Python and `python3 -O` generated byte-identical `receipt.json`; normal and optimized logs are retained. This matters because an optimized execution of an assert-based checker would not be a valid independent verification.

Audited SHA-256 pins at this addendum's completion:

- independent audit.py: 7bd6f1f98229667496f460ba060cd03348e097a3cac656995021616af63c5925
- independent receipt.json: 01091e9a6cef159e3355de5962be3c58e507029f6a54776760aa3051ae953aa9
- producer README.md: a60a8eb0dfe475e1907f372f7a503834cde84a085c1511fd39f64585ac228773
- producer check.py: 1d24188af100b3394b9e47b4b27e35e80b15863e9f033c05081868d7ff5a1632
- producer receipt.json: e5ec575906451e4a6ee293a98ed85dc3d90a2a8b18b878a921df151166d9032a
- producer moment-hostile.json: 780baef1f1d4e874ad8c69f6c9a542f84e36f443584ba0eb79348ceb96fedc8a

These pins identify the pending-old Cauchy–Schwarz version. A subsequent exact joint-zonotope strengthening requires a separate source/receipt pin update and must not silently inherit these file-level pins.

## Separate exact joint-zonotope strengthening

**PROVED algebraic strengthening; independently tested.** For pending old/new multisets define

    G = max_{all individual old/new signs} [Q(S_n)+beta(S_o,S_n)].

Then the exact strengthened gate is

    Delta(d) <= Delta0+L_n+L_o+G.

Proof: at the future cut assign one independent fair sign eta_u to each terminal and multiply BOTH its old and new group sums by eta_u. Independently add fair signed charges of every absent old source and every absent new source. Each realization is an allowed pair of full individual signed sums. Its expected objective Q(S_n)+beta(S_o,S_n) equals the desired sum_u[Q(z_n,u)+beta(z_o,u,z_n,u)] plus the nonnegative sum of norms of the absent new sources. Cross-terminal terms vanish, same-terminal old/new mixed terms survive, and all additional mixed terms have expectation zero. Hence the desired quadratic/mixed expression is at most G. The previously established linear estimates finish the proof. Shared group signs between the two classes, rather than independent signs for the old and new groups in the same terminal, are essential to this argument.

Let Z_o,Z_n be the two signed zonotopes. For fixed old sum, the objective is convex in the new sum, so its maximum over Z_n occurs at a new vertex. For that chosen new vertex, it is linear in the old sum, so an old vertex maximizes it. Thus at most12x12 corner pairs suffice, using the exact twelve-pattern vertices already proved. Also

    D_n <= G <= D_n+2 sqrt(D_o D_n).

The lower bound follows by choosing a new maximizer and an old sign choice with nonnegative mixed term (the old domain is centrally symmetric). The upper bound follows from Q(S_n)<=D_n, Q(S_o)<=D_o and Cauchy–Schwarz. When the old pending set is empty, G=D_n exactly, so existing successful native old-completed certificates do not change numerically.

Separate independent `joint_audit.py` compares the144-corner evaluator with exhaustive individual-sign enumeration on144 deterministic multiset pairs, old/new sizes each0..5, with repeated vectors and zeros. It exhausts2401 absence/three-terminal/sign assignments for two old and two new sources. The tested G=17 equals the maximum actual quadratic/mixed expression17, versus the older rounded bound D_n+ceil(2sqrt(D_oD_n))=18. Normal and optimized executions produce identical output, with non-removable checks.

Separate pins:
- joint_audit.py: 9972a63762654a517ad21abd79eb3be0056638f0814c1710ebf2bf795816ebf8
- joint-receipt.json: 339038e30cae1970b57120a2454d71f66eb0a152a31c6b9cc78fdc87d0251842

This section proves and independently tests the formula. It does not pin a subsequently edited producer implementation. The earlier producer pins remain explicitly for its previous Cauchy–Schwarz version until a final producer hash addendum is recorded.

### Final producer mixed-corollary review and immutable replay

The producer kept the base files unchanged and placed the joint refinement in separate `MIXED-COROLLARY.md`, `mixed_gate.py`, and `mixed_receipt.json`. All three were reviewed. The implementation correctly takes the maximum over both vertex lists, including empty-list zonotopes represented by the zero vertex, and its independent-sign exhaustive controls and partial-terminal controls implement the claimed predicates. A reviewer-run optimized immutable replay passed:

    python3 -O session22-residual/mixed_gate.py --verify session22-residual/mixed_receipt.json

This producer replay is additional to, not a replacement for, the independently implemented `joint_audit.py` checks above. The base README/check/receipt hashes were rechecked unchanged. Final additional producer pins:

- MIXED-COROLLARY.md: 0b8184b119ab701b4c5a2376965d45b8623d9323c394511214d1a0cb0f32d7a8
- mixed_gate.py: 1bc34c05094ca2b855d3f50d2fc73a7994f563e6fabfcb75eadbacfc7c5f9e0d
- mixed_receipt.json: 4a65b10d2066e63258eb9be99039a4c12c58d3960e999ab1f3d4668e35d6fbab

**Final verdict: PASS** for the conditional general suffix theorem and separate exact joint-zonotope corollary, with the finite-native independent audit scope and all open universal obligations unchanged.

### Publication portability correction

The independent `audit.py` repository locator was changed, without changing any mathematics or computation, to discover either the staging sibling `math` checkout or the nearest containing checkout with the authenticated frontier file. It therefore runs from a published audit subdirectory as well as this workspace. A relocated optimized N7d8 smoke test passed in a temporary nested checkout layout. The earlier independent audit.py pin is superseded by:

- final audit.py: 9d12d89a3c18dfeed1500408628f3840366478400001ec15f78d2368bf3636f8

All receipt and joint-checker pins are unchanged. Publish REVIEW.md, audit.py, receipt.json, joint_audit.py, and joint-receipt.json; replay logs are optional and duplicate `.normal.json` files need not be published. From a published audit directory, use `python3 audit.py` or `python3 -O audit.py`, and likewise `python3 joint_audit.py` or `python3 -O joint_audit.py`.
