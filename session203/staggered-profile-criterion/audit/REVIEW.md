# Independent audit of the conditional staggered profile compiler

6 November 2026. **PASS for the conditional implication, finite certificates, and guarded compiler.** The universal hypothesis F(8,137) remains OPEN. This audit neither proves it nor finds a general contradiction to it.

## Frozen inputs and dependency

The reviewed PROOF.md has SHA-256 acff2714fb6fe02a205f09f27add420a893e190a3bdf7fb01386477bfd66f546. The compiler hash is 84482dd4c2abf3293155ba6d3963b5335118fa152f10a1039a2bc89906cc7997; the producer verifier hash is 4f3c54e5ae423086499f49bb6705114e37f64bf5501ea963786c73c4fdbd63a4. All entries of the producer's frozen input manifest passed. The reviewed proof, programs, and exact two data files are retained here.

The physical readout H(C(U))=3B(U)+2s is inherited from the separately audited positive additive construction. Its frozen proof has hash 211623af96e90d516fbce7598450642878f81fc2892446dee5173b55fb41b099, identical to the published session203/free-additive-difference/PROOF.md at commit 150e50fc16d9ec57c596872559b9823b1612f3fe. Its proof preserves the same root-ending two-path-cover convention, including an empty second path. No source-bank completeness is needed here.

## Actual coupled profiles, not free coordinates

For a nonempty rooted tournament U, B(U)=h+L-s is the count of the literal rooted-edge extension, so it is positive and odd. F(J,B0) asks for an actual single rooted graph realizing both B(U)=B-2j and s=r+3j, with j in a fixed finite range and r=1 or 2. The two residues may use different graphs. The condition is not an assertion that separately attainable coordinates can be combined.

For J=8 and B0=137, the shifted base is at least121 and the required endpoint count lies between1 and26. The input graph may have arbitrary order and need not be strong. Its edges and root are copied into the four-core child. The exact identity then gives

    3(B-2j)+2(r+3j)=3B+2r.

The cancellation is an evaluation of one realized graph's count, not subtraction of two graphs or addition to an independently prescribed base.

## Exact induction endpoints

The general implication is sound on the domain where F(J,N) is defined: N is odd, N>=3, and N>2J. The seed interval [N,3N-2] contains every initial odd target. The next odd target is exactly3N, not3N-1.

- If n is divisible by3, n/3 is an odd earlier target at leastN. An ordered triangle join multiplies its count by3.
- For n congruent2 modulo3, r=1 and B=(n-2)/3 is odd. The first such target after the seed interval is3N+2, giving B=N.
- For n congruent1 modulo3, r=2 and B=(n-4)/3 is odd. The first is3N+4, again giving B=N.

Thus the two nonmultiplicative branches never request a profile below the hypothesis's cutoff. At N=137 the exact seed endpoint is409, followed by411=3*137,413=3*137+2, and415=3*137+4. The supplied witnesses23 through135 fill the lower interval after the conditional induction establishes all odd targets at least137.

This proves only F(8,137) implies that every odd count at least23 occurs. F is a sufficient condition, not claimed necessary for the unrestricted spectrum. Failure of this particular profile condition would not refute that spectrum statement. The exclusions7 and21 are not needed for this implication.

## Independent complete replay of the supplied finite data

check_certificates.py imports no producer code. It validates every tournament and uses nonmemoized outgoing path DFS. For profiles, it counts all literal paths on each subset, tracks the actual rooted endpoint counts, and pairs their choices across complementary subsets to compute L. This is independent of the producer's subset dynamic programming recurrence.

The replay confirms exactly:

- 194 seed graphs, one for every odd target23 through409. Their orders are6,7,8, with counts10,58,126 respectively.
- 862 actual profile certificates, both residues for every odd B from137 through997. Their source orders are7 and8, with counts159 and703 respectively.
- 862 literal four-core children, each with the exact required count3B+2r. The largest target is2995.
- 1,918 literal graph calculations and29,609,972 path-prefix visits.

The j frequencies are294 at0,300 at1,140 at2,79 at3,35 at4,10 at5,2 at6,0 at7, and2 at8. This verifies the declared finite bound, not its minimality or an all-order bound. Normal and optimized independent receipts are byte-identical and use explicit exceptions.

Every certificate stands on its own adjacency rows. Neither deletion closure, known isomorphism-class totals, nor completeness of the original bank is used to validate these witnesses. The source-bank statement is provenance only.

## Compiler behavior and hostile tests

check_compiler.py exercises the exact frozen compiler and checks its returned graphs independently. Targets23,409,411,413,415,1239 and2995 all return valid literal witnesses with independently confirmed counts.

The compiler correctly rejects both absent base999 records: target2999 needs residue1 and target3001 needs residue2. These are missing-certificate failures, not nonexistence claims. It also rejects an even target, a target below23, an incorrect seed count, j=9, an invalid root, a false profile equation, and an invalid tournament. Nine rejection controls pass.

A large control compiles413*3^8=2709693 to34 vertices. Instrumentation records a subset-DP count only on its10-vertex profile child. The eight surrounding triangle factors and every cross edge are checked directly, and the base is independently counted as413. The resulting target count is certified by the ordered-join product bijection. No literal enumeration of the34-vertex graph is claimed.

Python integer bitmasks avoid a32-bit vertex ceiling for these joins. Small seed/profile outputs are still subject to the compiler's explicit20-vertex literal-count guard, and ordinary memory, runtime, and recursion resources remain practical limits. The program does not turn F into an implemented oracle or guess a missing profile.

## Open-hypothesis boundary

The review did not establish all-order truth or impossibility of F(8,137). Every supplied source happens to be strong, and 807 of the 862 supplied profiles have gcd(B-2j,s)=1; PROFILE_STRUCTURE.json records these finite observations. They are not an occupancy argument. Proposed SCC normalization deductions are being investigated separately and are not part of this approved release.

Existing fixed-context obstructions do not automatically refute F, since F chooses an arbitrary full input tournament separately for each B and residue. The infinite quantifier remains the substantive unresolved step.

## Reproduction and approved status

    python3 check_certificates.py > NORMAL.json
    python3 -O check_certificates.py > OPTIMIZED.json
    python3 check_compiler.py > COMPILER_NORMAL.json
    python3 -O check_compiler.py > COMPILER_OPTIMIZED.json

Both receipt pairs agree. The independently approved publication status is “conditional implication and finite certificates independently audited,” with F(8,137) still explicitly OPEN. It must not be promoted to an unconditional threshold or an all-order profile theorem.
