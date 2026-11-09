# Independent audit: uniform proper-graft scalar descent

**Verdict: PASS.** The all-order reduction and finite low-target compiler are sound in their stated scopes. No mathematical correction is required. The finite image uses the explicitly **CITED** Busch strong floor; this audit verifies the repository citation and its application, not the primary paper. No all-odd completeness theorem or new unrestricted cutoff follows.

## Objects and proof audit

For finite nonempty rooted A with root a, use D=H(A), E=Hamilton paths ending at a, d=H(A-a), and U=ordered two-path spanning covers whose second path is nonempty and starts at a. Empty first paths count once. Gadget B is any nonempty tournament, with H=H(B) and Q its ordered two-bin spanning cover count, either bin allowed empty. The graft directs B to a and A-a to B. Its scalar output does not depend on a choice of gadget root.

1. **Scalar row and positivity.** Put k=D-E, V=U-d, and f=V-k. The inherited physical inequalities give d>=1, k>=1 for nonsinks, f>=0, Q>=2H, and H>=1. The two displayed forms of N agree: H(D-E+U)+(Q-2H)(U-d)=Hd+Qk+(Q-H)f. Hence N>=H+Q>=3H. Also N-HU=Hk+(Q-2H)V>=H, proving U<=floor(N/H)-1. These are integer inequalities, with no matrix invertibility hypothesis.
2. **Quantitative strict descent.** The inherited endpoint theorem gives k>=(D+1)/2, while U>=D and V>=1. Therefore N>=H(3D+1)/2+(Q-2H), exactly as claimed. The rearranged D bound is strict: since Q>=2H, D<=(2N-H)/(3H)<= (2N-1)/3<N for N>=3. Thus D<N and H<N for every proper graft, not just reachable-root or strong sources.
3. **All gadgets really are bounded.** Each ordered vertex bipartition admits at least one directed path per part. Thus Q>=2^|B|, including empty parts. With H+Q<=N and H>=1 this gives 2^|B|<=N-1. For a fixed target there are finitely many labelled gadget tournaments, with no menu restriction. Grouping them by (H,Q) loses no scalar information. It would lose information about the full output four-state, which is not the claimed decision problem.
4. **Finite inverse quantifiers.** D and U are bounded; d<=D and 0<=E<=d bound the remaining coordinates. The parity and positivity filters are necessary physical conditions, not claimed sufficient. The Session228 representative theorem makes the subsequent full-state membership search genuinely terminating. Every retained witness is constructed literally, proving sufficiency; every proper graft's full input state and gadget pair occurs in the search, proving necessity. At N=1 there are no proper grafts; the stated algorithm starts at odd N>=3. This proof does not assert recursive scalar replacements preserve a full state.
5. **Reachability survives the size bound.** The Session228 compression deletes only singleton components in the inactive initial prefix of A-a. Such a prefix beats the root and every later deletion component, so its vertices are unreachable from the root. When the root reaches all A, this prefix is empty, and the original graph itself already satisfies the representative bound. Consequently filtering the bounded oracle for root reachability remains complete. Conversely, strong output is equivalent to root reachability and |A|>=2 for arbitrary nonempty B, as the inherited Session226 theorem states.
6. **Prime normalization.** SCCs of A preceding the root's SCC form C=>R. A root-starting path cannot enter C; the other bin therefore contains a Hamilton path through all C. All four state coordinates scale by H(C). The root reaches every vertex in R, and if it was a nonsink in A it remains a nonsink in R. The suffix graft count is at least 3. For prime N, factorization N=H(C)H(graft(R,B)) forces H(C)=1, preserving the full state when C is removed. This is normalization of an existing representation, not existence for every prime. The gcd(H,Q) obstruction is also valid by the scalar row; for prime N its divisor cannot equal N because H<N.
7. **Sink boundary.** A sink has state (d,d,d,d), and N=Hd. Singleton A can reproduce H unchanged. Including it would invalidate strict descent and make an unrestricted scalar graft-image decision tautological. Properness excludes exactly this boundary through V>0.

## Complete low-target source universe

For a reachable root, its SCC is the first SCC. It is either singleton or strong. Any later components are strong and ordered. The cited Busch floor H(T)>=25 for strong order at least 7 excludes initial order >=7 whenever total D<=20: component Hamilton counts multiply and are positive. Since appending a component multiplies U by Q>=2, the initial U cannot exceed final U. Therefore enumerating initial rooted strong tournaments through order 6 with D<=20,U<=31 is exhaustive.

Appending any nonempty final SCC with pair (h,q) sends (D,E,U,d) to (Dh,0,Uq,dh): the root cannot be the final Hamilton endpoint, and each path bin independently appends its assigned final component path. Empty paths cause no exception. Since U>=1, every appended component has q<=31; q>=2^order implies its order is at most 4. The only usable strong component pairs found independently are (1,2) and (3,12); the strong four-vertex pair (5,32) is excluded by q. Each append at least doubles U. Thus arbitrary singleton tails are included but cannot make the compiler infinite. Full-state deduplication is safe for this exact continuation law.

For any proper output N<=31, the universal bounds imply D<=20 and U<=30, so the compiler's U<=31 cap is harmlessly generous. Every gadget has order <=4. Singleton starts are correctly included for later tails but excluded as final proper sources. These arguments establish the final compiler's complete universe, rather than an experimentally chosen cap.

## Independent finite verification

`check.py` imports only Python standard-library modules, and no producer or repository implementation. It enumerates all 33,867 labelled tournaments of orders 1 through 6 with root fixed at 0. All labelings with fixed root cover every rooted isomorphism type, so iterating every root again is unnecessary. Fresh endpoint subset DP computes path counts; applying the same recurrence to the reversed tournament computes root-starting counts. Independent literal permutation enumeration checks all four coordinates for the 1,099 tournaments through order 5.

The independent compiler finds 49 initial states and 77 final states. Every retained state is checked against its actual graph and root reachability. All 71 source-state/gadget-pair combinations producing N<=31 are checked by literal graft path DP and graph strongness. The exact image is:

    3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31.

Additional controls check all 5,625 arbitrary source/gadget pairs through order 4, including 900 sink cases. This tests the scalar law beyond reachable-root sources and checks all claimed strict descent inequalities where applicable. The prime-normalization precursor, exact four-state scaling by the root's preceding SCC prefix, is checked in 327 nonempty-prefix cases through order 5.

These finite checks support the implementations and reproduce the claimed image. They do not replace the all-order arguments. The absence of 7 and 21 here is scoped to the reachable-root proper-graft image, not a new proof of their unrestricted impossibility.

## Dependencies, reproduction, and boundaries

The producer PROOF.md and INHERITANCE.md specify the inherited Session225 scalar row, Session228 physical-state representative bound and its independent AUDIT.md, THM-4876 statement 7 endpoint inequality, Session226 strongness theorem, and the October 2 Camion/Busch note section 2. `AUDITED_INPUTS.json` records exact SHA-256 identities of the reviewed inputs. The pending label retained in the historical Session228 producer proof is resolved by that packet's independent audit; no reserved theorem is used.

Run in this directory:

    python3 -I -B check.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O check.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

Both supplied execution modes pass with identical receipts. All correctness gates use explicit exceptions, not assertions. `SHA256SUMS` seals the code, report, receipts, and input manifest. The independent witnesses may differ from the producer's because of root enumeration and representative selection; the complete state count and scalar image agree.

This worker made no Git writes. Repository startup reported existing unrelated dirty paths; they were left untouched. The optional coordination peek could not run because this checkout has no .machine-id; no identity or registration was created. This is not a mathematical blocker.
