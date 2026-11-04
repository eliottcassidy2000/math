# A named, two-branch, all-clock zero-ready construction

2026-11-04. **PROVED all-clock construction with finite-exact suffix witnesses.** This is an infinite family of literal two-birth `3H+2` and `3H+4` returns to the zero-ready class. It is not a closed iterated ternary grammar and does not close the 170807 source.

## 1. Exact affine lifting criterion

For a nonempty tournament A, a strong marked tournament (B,w), and a new physical vertex r, form Z(A;B,w) by orienting A→B, w→r, and r→every other vertex. This is strong. Its canonical zero-ready root is r, with sole in-port w and final SCC B after r is deleted. Write h=H(B), F=E_w(B,2), where E_w counts ordered two-slot spanning path covers, first slot nonempty and ending at w, second possibly empty. Then

    H(Z(A;B,w)) = H(A)h + P_A(2)F.

Here P_A(2) counts ordered two-slot covers with either slot allowed empty. This is the inherited ZERO_READY identity, proved by decomposing covers across A→B. For A=TT_s,

    H(Z_s(B,w)) = h + 2^s F,  s≥1.

**Affine lifting theorem.** Two marked strong suffixes (B,w),(C,v) satisfy

    H(Z_s(C,v)) = 3 H(Z_s(B,w)) + d for every s≥1

if and only if

    H(C)=3H(B)+d,      E_v(C,2)=3E_w(B,2).

The implication from right to left is immediate. For necessity, it suffices to subtract the identities at any two distinct s, then substitute back. In defect coordinates b=F-h, the exact condition is b'=3b-d. In particular b≥d/3 is necessary because strong suffix defects are nonnegative.

For arbitrary common prefix A, the same suffix pair gives offset d H(A), not generally d. Thus the two-number carrier is sufficient for this particular attachment grammar, not for arbitrary complements or port changes.

## 2. Explicit paired branches from one named suffix

Use labels 0,1,2,3 and outward-neighbor bit rows

    B = [8,1,3,6], w=0, a=1.

It is strong and (h,F)=(5,7). Add x=4 with every old vertex except a beating x, and x→a. The intermediate suffix rows are

    X = [24,1,19,22,2].

Then choose one of two births y=5:

* +2 branch: every original suffix vertex beats y and y→x. Result C2 has rows [56,33,51,54,2,16] and (H,E_0(2))=(17,21).
* +4 branch: w→y and y beats all other suffix vertices including x. Result C4 has rows [56,1,19,22,2,30] and (H,E_0(2))=(19,21).

Both outputs are strong. Their restrictions to labels0,1,2,3 are literally B; no re-marking or isomorphism replacement occurs. In the full Z_s attachment, every prefix vertex and r beats each newborn. Therefore the physical root r retains its entire old star and its sole in-port w, and every old name and edge is retained. The new suffix is the last SCC after deletion of r, so both outputs are again zero-ready.

Consequently, simultaneously for every s≥1,

    H(source)=5+7·2^s,
    H(+2 output)=17+21·2^s=3H(source)+2,
    H(+4 output)=19+21·2^s=3H(source)+4.

For example H=19 has outputs59 and61; H=33 has outputs101 and103; H=61 has outputs185 and187; H=117 has outputs353 and355. Source order is s+5, output order s+7. Each literal program has two permanent births, no complements, and no deletions. The actual intermediate full tournament is strong too: its suffix contains the old strong B, x has both an incoming and outgoing old incidence, and the attachment remains unchanged.

This answers the local construction question positively: both small offsets can coexist with an unchanged named zero-ready root, at every transitive-prefix clock. It does not say these particular outputs can be substituted back as sources of the same recipe.

## 3. Discovery census and verification scope

`search.py` enumerates every labeled strong tournament through order6: counts by order1,...,6 are1,0,2,24,544,22320. It retains one actual marked witness for each (h,F), yielding180 distinct pairs. Searching this finite bank gives6 affine +2 pairs and7 affine +4 pairs. The common source pair(5,7) above admits literal induced embeddings for both branches; those embeddings were renamed once in this writeup so every old label is unchanged. The census does not claim completeness above order6 or classify all embeddings of any response pair.

`verify.py` independently counts suffix H and E2 by permutations and one cut position, not subset DP. It verifies all old edges, strongness, and exact two-birth rows. It checks source/output counts and actual zero-readiness at s=1,2,3,4, independently verifying full output counts by permutations at s=1,2. `verification.json` gives exact receipts. The arbitrary-s statement follows from the attachment identity, not finite interpolation.

The whole small census took seconds without installs or external software. The returned row witnesses, not just their response pairs, are the operational objects.

## 4. What remains open and how this relates to 170807

The source suffix at the saved56935 boundary has (h,F)=(18977,18979). A clock-preserving branch would require (56933,56937) for +2 or(56935,56937) for +4. The observed170807 suffix(2195,42153) satisfies neither. For the actual170807 suffix, clock-preserving targets must have (h,F)=(6587,126459) or(6589,126459). These are precise synthesis targets; existence is not established here.

The family above is clock-uniform, not recursively closed. In particular the simple birth recipe, tried with every in-neighbor a of the port on the displayed +4 suffix(19,21), does not satisfy the lifting criterion. Another suffix transformation, another physical mark, or a richer grammar is needed to iterate arbitrary digits.

F49 and F343 coordinates and the81/27 hostile control are observers of fuller physical data. This construction makes no phase-identification claim and does not use equality of phases, counts, or unrooted responses to replace a named source. It is compatible with the earlier warning that two-run summaries do not close arbitrary complements: no complement is used, and all incidence changes and actual names are explicit.

Source identities: pinned repository commit4ad839b20dd280cdcfbadf86cb62b2be3c9615ad, `05-knowledge/results/triangle_return_170807_20261103/ZERO_READY.md`, `RETURN.md`, and `TRIANGLE.md`. The affine criterion and paired all-clock witnesses here are additional results, not restatements of the fixed-clock obstruction there.

## 5. Independent fixed-name audit and the new512423 endpoint

A separate worker's `../older/audit_phases.py` independently applied the two births directly to the full source labels, keeping root index s+4 fixed and appending the newborns. It verified all old induced rows, strong intermediate and output, sole in-port, both offsets, and E_root(2)=H at s=1,...,4, with full permutation controls at s=1,2. It also permutation-checked all180 saved suffix representatives and the six/seven affine pair counts. This avoids confusing attachment-order indexing used for counting with a permitted physical operation. Both `verify.py` and the independent audit use explicit checks that remain enabled under Python `-O`.

The separate marked-return construction now supplies a512423 endpoint. Using its independently computed endpoint responses, `candidate_phases.py` recomputes the inherited field codecs from scratch, validating the implementation against all four saved170807-history states first. The principal physical-port observations are:

* Port39: golden(512423,10021081,10021126), F49/charge(3,3;6), F343 orbit(2,0,3)→(6,0,5)→(1,0,6)→(2,0,3).
* Port50: golden(512423,72196954,72196999), F49/charge(5,0;2), F343 orbit(2,4,0)→(4,1,5)→(3,2,2)→(2,4,0).

Port50 is F49-Frobenius-fixed but still has degree-three F343 codec. The actual current phases at ports6,51,4 are also retained in `candidate_phases.json`. This replays the same distinction as the earlier170807 interface without assuming phases are preserved by the new program. The cubic modulus is X^3+X^2+1 over F7. The saved golden data are obtained from g=(H,(S2−E1)/2,(S2+E1)/2). These are observer codecs, not tournament automorphism groups or the inherited ancestral spectral clock. The five Q4 presentations are not recomputed here.
