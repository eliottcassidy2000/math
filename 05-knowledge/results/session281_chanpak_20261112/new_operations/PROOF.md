# Coprime-deletion central trades beyond the old fifteen channel

**FINITE-EXACT seed classification at order seven; PROVED all-context transport from the inherited fusion law. Chan–Pak OPEN.**

## 1. What was sought, and what was found

The known central +2 trade at H41/43 has deletion count15. We sought physically reusable central trades with deletion count coprime to15. The experiment enumerates all56 isomorphism classes of six-vertex deletions, then every one of their64 root cuts. Every seven-vertex pointed tournament occurs up to pointed isomorphism in this universe, though some pointed classes repeat.

The nontrivial pointed carrier is

    M(T)=[h,U-h,d-h,H-h; 0,Q,0,V-h; 0,0,0,e-h; 0,0,0,h].

Definitions and physical rooted fusion are inherited from [Session280](../../session280_arbitrary_switching_20261112/escape/PROOF.md), citing Session264. Equal(h,Q,U,V,d,e) and count difference4 give carrier difference4E14.

Among1728 presentations with gcd(h,15)=1, exactly four full-carrier fibers have more than one H value (converse pairs):

|h|Q|U|V|d|e|H values|
|---|---|---|---|---|---|---|
|17|224|89|152|11|6|65,69|
|17|224|152|89|6|11|65,69|
|29|320|122|227|20|9|91,95|
|29|320|227|122|9|20|91,95|

No central +2 trade with deletion count coprime to15 occurs at this order. No larger-order absence is inferred. Explicit representatives rooted at vertex6 are

    X65=(96,65,3,7,79,94,12),
    X69=(102,8,66,69,79,30,34),
    X91=(116,81,66,7,12,94,24),
    X95=(114,80,3,7,76,94,12).

The count and all six remaining coordinates of every displayed witness are independently recomputed by whole-tournament endpoint DP, rather than only by the deletion-cut compiler. Converse representatives are in `RESULTS.json`. No global minimality or literature-novelty claim is made.

## 2. Universal physical transport

For any pointed tournaments L,R, including singleton identities, triangular multiplication gives

    M(L) E14 M(R)=h_L h_R E14.

Therefore L*X65*R and L*X69*R have equal full continuing carriers except for count gap4h_Lh_R; likewise X91/X95. Associativity covers arbitrary iterated common fusion contexts. These are actual constructions; subtraction only compares their matrices. Since h_L,h_R are positive odd integers, no common context converts these +4 trades into +2.

## 3. A physical two-axis combination, with the capacity gate retained

Suppose actual B_i, 0<=i<r, have common receiver(h,V,e)=(b,V,e) and H(B_i)=H0+2i. Let X0,X1 be either central pair, deletion a, gap4. Then

    H(X_j*B_i)=N+2a i+4b j.                          (1)

All outputs share receiving coordinates

    h'=ab,
    V'=Q_X V+(V_X-Q_X)b,
    e'=e_X b.

These follow by matrix multiplication; nontrivial X has its third diagonal0. Equality of Q,U,d among the B_i is unnecessary. For the renewable six-bank b15 and the h29 trade, the offsets are58i+60j, 0<=i<6, j∈{0,1}. For i=1,...,5, choices(i,0) and(i-1,1) give consecutive odd counts N+58i and N+58i+2. All five pairs share one receiver and survive common left unit-context transport.

The twelve offsets are

    0,58,60,116,118,174,176,232,234,290,292,350.

They do NOT form an interval. Primitive gcd2 leaves explicit large holes. Frobenius/CRT would still require sufficiently large physically admissible digit ranges at these same interfaces. Equation(1) is a mixed-channel construction, not physical subtraction or a universal target solver.

## 4. Exact universe and reproduction

Canonical representatives are obtained by enumerating all32768 labelled six-vertex tournaments and minimizing adjacency code over permutations within sorted degree cells. Every isomorphism preserves degree, so this is exact canonicalization, not hashing. The root-cut engine splits a path into the subset before the root and its complement; its prefix must end at an incoming neighbor, its suffix start at an outgoing neighbor. U,V,d,e are corresponding full-subset linear endpoint sums.

Run `python3 verify.py` and `python3 -O verify.py`. The deterministic output includes64 literal source/sink/triangle context checks and independent whole-graph profiles for all displayed seeds. Explicit exceptions enforce checks under optimization. The original exploratory two-apex-word route did not reveal reusable constant-gap trades in its bounded tests; no all-word obstruction is claimed from that failure.

This packet was restored after runtime rollback and replayed, rather than presumed byte-identical to lost workspace files. Mathematical statements and literal witnesses match the retained report; see `RECOVERY_RECEIPT.json`.
