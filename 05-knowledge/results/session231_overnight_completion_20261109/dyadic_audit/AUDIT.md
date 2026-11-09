# Independent audit: periodic near-dyadic twins

**Verdict: PASS, in the stated positive-zero-tail actual-source scope.**
No mathematical correction to the theorem is required. The proof below supplies exact-order and asymptotic details not explicitly proved in the producer packet. Universal odd-count coverage and globally minimal tournament order remain open; nothing here proves either.

Audited input: `../dyadic/PROOF.md`, `core.py`, `verify.py`, and `CHECK.json`. Independent checker `check.py` imports none of these. It was run normally and with `-O`; outputs agree byte for byte. All checker predicates use explicit exceptions rather than assertions. The two producer graph cores were independently validated as tournaments, and their cut states recovered using literal permutations of all vertex subsets and ordered complement splits. Full Hamilton counts also use a separately written backward endpoint subset recurrence.

## Independent proof and scope audit

### 1. Physical state and insertion identities

Write (h,Q,A,G,T,S,E) as in the candidate. For a nonempty marked parent L, insert r pointing precisely to the selected vertices. Paths with r not last correspond to ordered complementary path pairs with the second path beginning selected and the first either empty or ending unselected: exactly A−T. Paths with r last contribute h−E. Thus

    F = H(L+r) = A−T+h−E.

The inherited source coordinates therefore satisfy q=A, g=E, and c=q−F−g=T−h. In particular source g must not be confused with seven-coordinate G. For a final unmarked singleton v, ordinal composition gives A'=2A, T'=T, h'=h, E'=0, so H(L⇒v+r)=2A+h−T. More unmarked singletons double A and preserve T,h. These identities follow from actual paths, independently of the historical packet. Integration subsequently verified the same convention directly in `../../square_sum_spectrum_overnight_20261031h/README.md`, sections 1 and 3, at main commit b38a6a7d891eb699c9007f13925b27700b1cfb50. The file was present in Git but absent from the sparse working tree.

For two nonempty ordinal blocks, every path visits the first before the second. Splitting the two covering paths across blocks yields the stated composition law. For T, the possibilities are: the two relevant endpoint constraints are separated (A1 G2); both in the left block (T1 h2); both in the right (h1 T2); and the exceptional concatenating full paths (E1 S2). This gives all four terms without overlap.

Literal enumeration gives exactly

    C: (323,5096,3672,3672,2610,252,254), cut 189;
    D: (323,5096,3672,3672,2610,254,252), cut 190.

For an MSB-first singleton word P, its h=1, A=P, T=f(P), E=P mod 2. Consequently A_L=5096P+3672 and T_L=3672P+323f(P)+2610+(P mod 2)S_core. For odd P, T_D−T_C=2.

The arbitrary binary-context formula inherited from Session229 is also correct. In particular its difference for D minus C is

    H_D−H_C = 2(first suffix bit − last prefix bit),

with absent endpoint bits assigned zero. Both selected boundary bits cancel the difference; neither also gives zero. A selected prefix endpoint and zero suffix give −2. Thus this pair is not a universal pure-T switch. The checker directly tests 12 core/context combinations, including empty, zero, selected, and cancelling boundaries.

For general nonempty ordinal contexts U,V, the corresponding difference is

    H_D−H_C = 2(h_U S_V−E_U h_V),

because A,h are unchanged and ΔT=2(E_U h_V−h_U S_V). This is a useful all-topology explanation; the theorem needs only singleton contexts.

### 2. Exact congruence class and positive prefixes

5096=8·637 and 3672=8·459, with 637=49·13. Modulo 49, 2 has exact order 21: 2^21=1, while 2^7=30 and 2^3=8 are not 1. Modulo 13, its exact order is 12: 2^12=1, while 2^6=12 and 2^4=3 are not 1. Hence ord_637(2)=lcm(21,12)=84.

Since 2^14≡459 (mod 637), the integral solutions to 5096P+3672=2^m, with m≥3, have precisely m≡17 (mod 84). No m<3 works. The first positive solution is m=17, P=25, and all later ones are m=17+84t, t≥0. The quotient (2^(m−3)−459)/637 is odd. Thus the claimed period is exact, not just a sufficient period.

For m≥17, 2^(m−13)≤P<2^(m−12). The lower inequality is weakest at m=17 and follows from 2^m(1−5096/8192)≥3672; the upper follows from 5096>4096. The binary length is exactly m−12.

### 3. Literal graphs, strongness, and all k

Set b=k−m−1≥0. Append b unselected singleton vertices, then add v and r as prescribed. All pairs of vertices receive exactly one arc. The mandatory v supplies an additional unselected singleton, so the actual tail length is b+1≥1, including at the first exponent.

The first prefix bit is selected, so r reaches its initial vertex, which reaches every later ordinal block. Every parent vertex reaches r through v. Hence all vertices lie in one strongly connected component. The order is (m−12)+8+b+2=k−3. The insertion identity yields 2^k−c and 2^k−c−2 for every integer k≥m+1.

The independent literal-graph counts are:

- C, k=18: 15 vertices, 163283 paths.
- D, k=18: 15 vertices, 163281 paths.
- C, k=19: 16 vertices, 425427 paths.
- D, k=19: 16 vertices, 425425 paths.

Independent source counts give (F,q,g,h)=(31957,131072,254,323) and (31957,131072,252,323). Large-parameter graph claims rely on the proved composition identity, not on an alleged direct enumeration of huge graphs.

### 4. Bounds and source optimality, including clock shifts

For a length-ell word, all pair weights in f are nonnegative. The all-one word has sum 2^ell−2 (including the adjacency correction), so 0≤f(P)≤2^ell−2. Substituting ell=m−12 gives exactly the candidate's upper bound; its coefficient 459/637+323/4096 is below 4/5 and its additive constant is negative. The lower bound uses f≥0, and already exceeds 2^(m−1) at m=17. Thus

    2^(m−1)<c<c+2<2^m.

Every tournament admits a Hamilton path (insert a new vertex into an existing path at its first out-neighbor, or at the end), hence every actual source has F≥1. As g≥0, c=q−F−g<q. If q is dyadic, each deficit therefore forces q≥2^m. The constructed actual parent attains q=2^m, proving the claimed global minimum within this source class, without restrictions on the competitor parent's topology or order.

Absorbing d≥0 zero-tail vertices into the parent changes q to 2^d q and leaves c fixed. This cannot lower source cost. A smaller-slope physical parent is ruled out by c<q, even if found by a different topology. Relabeling the exponential clock or permitting formal nonphysical source states would change the definition and does not evade the argument. For a literal equality q2^r=2^k with integer r≥1, q itself is dyadic. Thus arbitrary integer slopes do not produce a smaller dyadic landing with the same fixed deficit through this multiplication alone.

The earliest route entry is m+1 because r≥1. Allowing r=0 is a different route: at zero tail the root count is F, while q−c=F+g. The proof properly does not infer impossibility of unrelated earlier-exponent realizations, zero-tail coincidences, or smaller tournaments.

### 5. Infinite family and actual asymptotics

Let M=2^84, B=459(M−1)/637, R the reversal of the 84-bit word B, and d=459/637. Then

    P_t=(25+d)M^t−d,
    f(P_t)=14+R(25+d)(M^t−1)/(M−1)+t(f(B)+1−Rd).

This follows by summing the candidate's concatenation recurrence. Cross pairs contribute RP_t; the boundary adjacency contributes one since both boundary bits are 1. The values B,R,f(B) reproduce exactly. P and f strictly increase, so c_t strictly increases and the theorem gives infinitely many distinct twin pairs.

In particular c_t=Theta(2^m_t), and more precisely

    lim c_t/2^m_t = 459/637 + 323R/[5096(M−1)]
                 = 57269513676081110771901207/75533314657546670029764568
                 ≈ 0.758202045491184.

Thus “near dyadic” does not mean a vanishing relative deficit at the first admissible exponent. At fixed t and k tending to infinity the fixed deficit is relatively negligible, as usual. The checker verifies parameters t=0 through 8; the displayed identities supply the all-t proof.

### 6. Hostiles and historical-source limitation

Suffix rigidity is valid: if 2^b Z+hR=2^a, Z>0, odd h, 0≤R<2^b, then a≥b since the left side is at least 2^b. Reduction modulo 2^b and invertibility of h imply R=0, then Z is dyadic. No arbitrary binary suffix repair is available.

The named historical `square_sum_spectrum_overnight_20261031h` packet was not present under the supplied worktree or scratch root. This audit therefore does not claim to have inspected that particular file; it independently proves the source convention the candidate explicitly states. Session229's retained dynamic proof was inspected and agrees with the arbitrary-context and composition formulas. No discovery-census completeness or literature-priority claim was used.

## Reproduction

From this directory:

    python3 -B check.py > CHECK.json
    python3 -O -B check.py > CHECK.optimized.json
    cmp CHECK.json CHECK.optimized.json

No producer or repository files were edited, and no commits were made by this audit worker.
