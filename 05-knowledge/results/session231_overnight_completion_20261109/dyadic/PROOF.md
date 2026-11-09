# An infinite optimally affordable near-dyadic twin class

Status: PROVED, with literal finite certificates and an independent all-parameter audit (../dyadic_audit/AUDIT.md). Universal Chan–Pak completeness remains OPEN. No literature-priority claim.

## Inheritance and research contract

The closest constructions are `square_sum_spectrum_overnight_20261031g/README.md`, its 118-deficit window, and packet h's exact source-cost convention. The repair law is Session229 `dynamic/PROOF.md`, independently audited there; endpoint erasure is Session231 `../fibers/PROOF.md`. The canonical hostile is the expensive deficit9 source: existence does not authorize an arbitrary target exponent. The corrected near miss is trying to repair a non-dyadic A using a binary suffix. The least-used sidecar here is the finite power-of-two orbit modulo the cover count Q.

Concept board: pure-T fibers; endpoint erasure across topologies; modular source entry; positive tail length; near-dyadic coverage. Anchor: parameterized physical windows. Niche: sharp entry cost. Wildcard: periodic binary words and their pair-energy cocycle. The successful bridge maps a physical fixed-(h,Q,A,G) endpoint swap into actual dyadic sources via one periodic language of selected prefixes. It preserves the twin difference and actual source slope, but does not preserve a fixed deficit as the prefix parameter changes. The binary prefix and actual graph/cut witnesses are the necessary sidecars. The decisive controls are literal path counts on the two cores and four completed graphs.

## 1. Coordinate convention and the physical pair

Use Session229's seven coordinates (h,Q,A,G,T,S,E). Q is the number of all ordered spanning two-path covers, allowing empty paths. A counts those whose second path starts in the marked set; G those whose first path ends there; T imposes both. S,E count marked starts/ends of full Hamilton paths. For nonempty ordinal blocks,

    h'=h1 h2, Q'=Q1 Q2,
    A'=A1 Q2+h1 A2, G'=G2 Q1+h2 G1,
    T'=A1 G2+T1 h2+h1 T2+E1 S2,
    S'=S1 h2, E'=h1 E2.

A selected/unselected singleton has (1,2,b,b,0,b,b).

The two eight-vertex tournaments are given by outgoing-neighbor bit masks, with vertices numbered0 through7:

    C=[4,65,26,99,43,135,181,31], cut189,
    D=[68,9,50,21,35,139,62,95], cut190.

Their actual states are

    C: (323,5096,3672,3672,2610,252,254),
    D: (323,5096,3672,3672,2610,254,252).

Thus they exchange S and E, while preserving h,Q,A,G,T. The verifier obtains these states separately from endpoint subset DP and literal permutation-and-split enumeration. Neither discovery census completeness nor a large numerical search is needed for the theorem.

For a positive binary integer P, let its bits be the selected marks of a singleton prefix, read in forward order most significant bit first. Put

    f(P)=sum over selected positions i<j of
         [2^(j-i-1)+1_(j=i+1)].

Prepending this prefix to either core gives

    A_L=5096P+3672,
    T_L=3672P+323 f(P)+2610+(P mod2)S_core.

When P is odd, the T values differ by exactly2, although the right endpoint E remains different. A zero suffix ignores this endpoint difference. This is endpoint erasure in the precise contexts being used; it is not asserted to be a universal pure-T switch before that erasure.

## 2. Periodic prefixes and the theorem

For every integer t>=0 define

    m_t=17+84t,
    P_t=(2^m_t-3672)/5096,
    c_t=3672P_t+323 f(P_t)+2539.

**Theorem.** P_t is a positive odd integer. For every integer k>=m_t+1, both

    2^k-c_t,  2^k-c_t-2

are Hamilton-path counts of explicit strong tournaments on exactly k-3 vertices. Both deficits c_t,c_t+2 have minimum dyadic source slope exactly2^m_t, among all actual marked-parent sources in the inherited positive-zero-tail construction. Therefore m_t+1 is their earliest entry exponent in that route. The theorem does not say these counts are unavailable at earlier exponents through other constructions.

**Integrality.** Write5096=8*637, 3672=8*459. The exact congruences are

    2^14=459 (mod637),  2^84=1 (mod637).

For example637=49*13, 2^21=1 modulo49 and 2^12=1 modulo13; both exponents divide84. Hence P_t is integral. Its numerator after division by8 is even minus459, so P_t is odd. P_0=25.

**Actual graphs.** Let L be the ordinal sum of the binary prefix for P_t and one of C,D, retaining all marks. Append b=k-m_t-1 unmarked singleton vertices to L. Add vertices v,r with every L vertex pointing to v, v pointing to r, and r pointing exactly to the marked vertices of L. All remaining arcs from L point to r. This specifies every arc of an actual tournament.

The star-closure identity is H=2A_L+h_L-T_L. Each appended zero doubles A_L and leaves h_L,T_L unchanged. For C the initial T_L-h_L is c_t, and for D it is c_t+2. Since A_L=2^m_t, the stated counts follow. The first prefix bit is1. Thus r reaches the initial prefix vertex, which reaches every later ordinal block; every vertex of L reaches r through v. The output is strong.

The inequalities

    2^(m_t-13) <= P_t < 2^(m_t-12)

hold for m_t>=17: the lower follows from
2^m(1-5096/8192)>=3672, already true at m17; the upper follows from5096>4096. Thus the prefix has length m_t-12. The completed graph has

    (m_t-12)+8+b+2 = k-3

vertices.

## 3. Sharp source cost, with no hidden tail shift

In packet h's source convention, an actual parent L and marked external vertex r have coordinates (F,q,g,h): F=H(L+r), h=H(L), g counts parent paths ending in the marked set, and q counts ordered two-path covers whose second path starts in that set. Thus q=A_L, not the seven-coordinate total-cover count Q_L. Their deficit is

    c=q-F-g=T_L-h.

Appending r>=1 unmarked zero vertices in the inherited orientation produces q*2^r-c. In our displayed graph there are b+1 such vertices, including v. In particular the family is valid at b=0 with genuine positive tail length1; no formula is extended backwards to a forbidden tail length.

For a binary word of length ell, f(P)<=2^ell-2: the full all-one word maximizes the nonnegative pair sum, and summing by pair distance gives equality2^ell-2. Therefore

    c_t+2 <= (459/637+323/4096)2^m_t
              +2541-3672*459/637-646 < 2^m_t.

The coefficient is less than4/5 and the remaining constant is negative. Also

    c_t >= (459/637)2^m_t +2539-3672*459/637
          > 2^(m_t-1),

already at m17 and hence for every larger m. Consequently both deficits lie strictly between2^(m_t-1) and2^m_t.

Every actual source has F>=1,g>=0, so its positive deficit c=q-F-g is strictly less than q. Any dyadic source for either deficit must therefore have q>=2^m_t. Our actual source has precisely that slope. This proves global optimality within the defined dyadic-source class, independently of parent size or topology. This is not a minimal-order theorem for all tournaments and is not a lower bound on arbitrary-slope sources.

At t0 the deficits are98861 and98863. The source certificates are

    C-prefix: (F,q,g,h)=(31957,131072,254,323),
    D-prefix: (F,q,g,h)=(31957,131072,252,323).

At k18 they give strong15-vertex counts163283 and163281. At k19 they give strong16-vertex counts425427 and425425. All four are checked by forward Hamilton subset DP on literal adjacency graphs. The source q is independently recovered as H(first zero extension)-F-g, rather than only by substituting the seven-coordinate expression.

## 4. The binary construction has an exact all-length controller

Put M=2^84 and

    B=(M-1)*459/637 =13937757016090795383111705.

The84-bit word of B is

    101110000111011011110101001001100010110111010000100100111110010010101100000000011001.

Then P_(t+1)=M P_t+B: the prefix is11001 followed by t copies of this fixed word. If R is the integer represented by the reversed84-bit word, concatenation gives

    f(P_(t+1))=f(P_t)+R P_t+f(B)+1,
    R=11485775830610390975213085,
    f(B)=8276249774333076071621348.

The cross-pair contribution is R P_t, with one extra adjacency correction because both boundary bits are1. This also follows directly from the ordinal composition formula. Hence construction and exact deficit evaluation use a finite periodic binary controller, not a separate search at every t. The deficits strictly increase because P_t increases and every energy term in this recurrence is nonnegative. This is genuinely an infinite set of distinct pairs of deficits.

## 5. Hostiles, scope, and replay

A suffix cannot repair a bad congruence. If h is odd, Z>0, 0<=R<2^b, and2^b Z+hR is a power2^a, then a>=b. Reduction modulo2^b forces R=0, so Z is already a power of two. This explains why most incoming fixed-interface pairs/triples do not automatically enter the dyadic source bank. The search instead tested the finite residue orbit2^m modulo Q. The successful pair uses different core topologies and exchanged endpoints.

Repeating independent nontrivial modules still multiplies pure-T differences by their other Hamilton counts. Nothing here evades that theorem: one fixed-size core is used, and only singleton contexts grow. Our window contains two consecutive odd values, not an increasing number of consecutive odds. The118-deficit fixed bank and this infinite cost-optimal twin class do not combine into an all-target cover. In particular, changing t changes both deficit and its required exponent; no arbitrary fixed-deficit existence inference is made.

Run `python3 -B verify.py > CHECK.json`, then `python3 -O -B verify.py > CHECK.optimized.json` and compare. Both have explicit active checks, independent core-state enumeration, literal completed-graph checks, exact large-parameter composition checks for t0..7, cocycle checks, order/cost checks, and a bounded suffix-rigidity hostile. Large-t graphs are certified through composition, not falsely described as directly counted.

`search.cpp` adapts the Session231 endpoint-erasure census to retain fibers whose A lies in the power-of-two residue orbit modulo Q. `SEARCH8.jsonl` has22 such endpoint-erased pair fibers in its declared order-eight input; the search is discovery only. Its external class input is `session231-dense-fibers/classes_8.txt`, generated by the inherited Session230 class generator. The retained mathematical theorem depends only on the two explicit witnesses, so reproducing the census is optional. No file outside this packet was changed.
