# Scoped review of the new Session259 growth results

Status: REVIEWED elementary arguments; not a canon promotion or an audit of inherited claims.

Scope: the new statements in `multiplicand/RESEARCH.md`, `dichotomies/RESEARCH.md`, and `transport/RESEARCH.md`. The summand manuscript was not present at the initial review. No heavy census was run.

## Findings

No substantive error found in the three inspected manuscripts. One minor convention clarification is recommended: transport initially defines its vertex set as integers at least 5, so its weak variant must explicitly enlarge that vertex set to include 4 before adding the loop 4→4. This does not affect any strict-graph theorem.

The power formulas in dichotomies should explicitly say integer k>=1 next to Q_k; the surrounding phrase “positive power” already supplies that restriction. The recurrence L_0=2 is valid, but the displayed divided-difference positivity and strictness claim should not be read at k=0.

## Load-bearing checks

- **Natural barcode:** the star of 2 is a nested contractible subcomplex. Each outside odd squarefree face m has every proper face in the star, because removal of an odd prime divides its product by at least 3. Hence the relative differential is zero. Inclusion sends the named generator to itself until exactly 2m, then to zero. Naturality of reduced-to-relative homology makes this an actual persistence decomposition, not merely a Betti-count coincidence. The degree-zero case is handled by nonempty contractibility; prime bars and the essential component are separated correctly.
- **Integral representatives:** the boundary on {2} union supp(m) exists at m and has coefficient one on the outside face. It bounds when 2m enters. This also justifies the integral version directly.
- **Strict multiplicand warning:** deleting a→a² restores that relation through a→ad→a² for composite a and a proper divisor d. Prime-square relations cannot be restored. Thus the weak barcode is correctly withheld from the strict reachability poset.
- **Joint chart:** discriminant parity and positivity give the claimed inverse and parent clock. Scaling and powers are injective on positive pairs, but not generally surjective; the manuscript states the image restriction. Input and output filtrations remain distinct.
- **Prime arity:** cyclic rotation has only orbit sizes 1 and ell when ell is prime. The arity-four hostile correctly prevents extending the congruence to arbitrary arity. Symmetric-power division is explicitly confined to characteristic zero.
- **Transport backlog:** the number of strict nonunit sum pairs born by N is floor((N-1)^2/4)-(N-2). The accumulated strict factor count is (D(N)-2N+2-floor(sqrt(N)))/2, with the n=1 correction included. Their difference is legitimate since every strict nonunit product is later than its sum.
- **Transport paths:** the legal a=4 lane at u>=5 gives u→4(u-1)←2(u-1)←u-1. An additional halving step reduces odd u by more than half in four edges. Bases 1,2,3,4 connect using 3→8 and doubling, so the stated excursion and finite-cutoff bounds hold.
- **Conductor bridge:** conditional on the two-seed closure formula proved in the summand lane, shifting the coprime semigroup Frobenius number by a+b gives ab. The restrictions a,b>=2 and, in the gcd version, both reduced seeds>=2 exclude the exceptional unit-generator cases.

No inherited theorem, historical numerical assertion, or literature-priority claim was used as a new reviewed conclusion.

## Summand manuscript follow-up

The summand manuscript subsequently became available and was reviewed. No substantive error found in its new two-seed closure, genus/conductor, cumulative witness deficit, parity-frontier, or exact deficit-coordinate formulas.

- The candidate closure {a,b} union (a+b+H) is closed under distinct-value sums, while its translated monoid can be constructed by repeatedly adding either smaller seed. Numerical coincidences in coefficient vectors do not invalidate this proof.
- Translation adds a+b−1 positive holes; restoring both seeds removes exactly two. Combining with (a−1)(b−1)/2 yields the stated numerator ab+a+b−5. The seed-one and gcd exceptions are correctly separated.
- The cumulative deficit formula includes all weak hole pairs, including diagonal pairs, and requires exactly N>=2F. Its witness certificate counts gates for all outputs, not a minimal circuit for one output.
- The parity-frontier induction is valid. For extra explicitness, the sum of two values from the inferior parity class is at most 2M−2delta−2, hence at most 2M−8 while delta>=3; it cannot displace the three exhibited top even values. The maximal opposite-parity sum is 2M−delta. Once delta=1, the three consecutive top values reproduce themselves under the update.
- The synchronous time T=(b+5)/2 and the much smaller static final hole 2b concern different filtrations. The manuscript expressly avoids claiming the whole conductor-to-frontier interval has filled at T.

## Strict multiplicand extension follow-up

Section 8 subsequently supplied the complete strict-reachability barcode. Its new splitting proof is sound, including the higher corrections [p²m,2p²m) in degree omega(m), and not merely the isolated-component corrections.

The map replacing every prime square p² by p is indeed order preserving in both posets. Prime-square vertices are pairwise incomparable; their strict upper links are copies of the *weak* divisibility order on 2,...,floor(N/p²), because no deleted relation x<x² with x prime can occur among two multiples of p². Thus the relative chain complex is the direct sum of the asserted cone/base complexes. The natural retraction splits the reduced homology exact sequence and consequently the persistence module. Empty links give degree-zero bars [p²,2p²); nonempty links give the suspension of reduced homology, not ordinary H_0. This distinction correctly avoids spurious eternal bars.

The weak equivalence to the common base is natural; the strict retraction is not asserted to be an equivalence. Allowing p to divide m is necessary and correct, including [27,54) for p=m=3. The first higher hostile [12,24) correctly refutes the simpler H_0-only correction. No additional mathematical fix is needed.

## Final synthesis and event-word review

The condensed table, new arc-shadow weld (W), Section 9 strict Euler identity, and event-uniqueness/silence consequences are mathematically sound. No scripts were rerun during this final review; the parent is performing the bounded replays.

For (W), a strict positive additive arc x→z exists exactly when x<z and z≠2x. Among the excluded doubling arcs, strict nonunit multiplicative witnesses exist precisely for x>=3. This proves the disjoint decomposition under the stated conventions, without importing a statement about paired-parent closure.

For Section 9, a correction birth p²m uniquely identifies p as its only repeated prime, with exponent 2 or 3. Weak births are squarefree, so cannot collide with these. Birth 2-adic valuations are 0 or 2 and death valuations are 1 or 3; this proves all event separation and the v_2>=4 silence claim. The empty-link Euler contribution at floor(N/p²)=1 is correctly included as M(1)=1. “Silent” is explicitly homological, not a homotopy-equivalence claim.

One synthesis wording repair was requested: “full symmetric packet or explicit collision data” in Section 6 should read “full elementary-symmetric coordinates or the parent tuple.” Collision partition data do not reconstruct a higher-arity tuple: the synthesis's own (1,6,6) and (2,2,9) share partition 2+1 as well as sum and product. The lane proof itself already uses the correct stronger data.

Optional precision: the radical map retracts to the squarefree subposet, whose order complex is the barycentric subdivision of the prime-product complex; the synthesis can say “naturally homotopy equivalent” instead of describing a literal retraction directly onto that latter complex.

Final disposition: all four clarifications are applied in the delivered copy:
the transport variant explicitly admits vertex4; powers explicitly require
k>=1; the inferior-parity bound is stated in the induction; and the synthesis
uses full elementary-symmetric coordinates/parent tuple and describes the
radical/subdivision map precisely. Six normal/optimized/stored-output replays
passed, as recorded in FINAL_REPLAY.json.
