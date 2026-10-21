# Native affine phase lifts and sparse unit boundary defects

**Status:** PROVED elementary lift, exact sparse-defect formula, and guarded simultaneous-arrival cancellation. FINITE-EXACT bounded controls. OPEN uniform annulus-increment bounds and universal entry. This is a separate bounded follow-on; the frozen annulus packet is unchanged.

## 1. Inheritance and exact carrier

The source is the same S_N={positive odd s<4N : s not congruent to5 modulo8}, with the Eisenstein charge A_N(s), palettes L=((1,2),(-2,-1),(1,-1)) and U=((0,1),(-1,-1),(1,0)), braid B(s)=4s+1, and Q(a,b)=a^2-ab+b^2.

Closest proved mechanisms:

- `native_cancellation_105_20261018g.md`, section3: A_(64N)(s)=A_N(s) for each named old source, and exact annulus source norm153N.
- THM-4690, `THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md`, statement2: the exact inverse P-preimages and strip clocks. The two lift-chains below are its period-three repackaging, not a new inverse classification.
- `eisenstein_circuit_completion_20261018g.md`: source-disjoint zero packets and the distinction between cancellation and an upper bound on surviving energy.
- `session18_native_annulus/README.md`, equations(A),(D): finite-arrival debt and the old-baseline correction for mixed old/new packets.

Canonical hostiles are the literal band crossing at N=2019,s=1009 and the stopped-frontier triple {1,113,7281}. Board: native inverse fibre / full strip clock / three-phase rotation / exact band boundary / synchronous source packet / remaining fibre interaction.

Write a(s)=2 if s=1 mod8, and a(s)=1 otherwise. If

(3s+1)/2^a=B^r(u),

with u primitive, then P(s)=u and the native edge clock is a+2r. These are primitive P-edge depths and full integer braid clocks, not fully odd U-word lengths or valuation costs. All residue receipts retain the compatible mod15/mod35 registers, jointly encoding mod105.

## 2. Exact native lift and its stopping-set guard

For a in{1,2}, put

c_a=21+7*2^a,   c_1=35, c_2=49,
L_a(s)=64s+c_a,   R(a,b)=(-b,a-b).

R is multiplication by the Eisenstein unit omega, so R^3=I, I+R+R^2=0 and Q(Rv)=Q(v).

**Lift theorem.** For every primitive s, with a=a(s), L_a(s) is primitive of the same branch, and

P(L_a(s))=P(s),   r(L_a(s))=r(s)+3,
edge_clock(L_a(s))=edge_clock(s)+6.                    (1)

Proof: a=1 lifts to residue3 modulo8 and a=2 lifts to residue1. Directly,

(3L_a(s)+1)/2^a=B^3((3s+1)/2^a),
B^3(x)=64x+21.

The three added B-strips are real strips, and the remaining r strips are unchanged. This proves the integer identity, not only clock parity.

For the actual inherited frontier, min(S\V)=1009 and max(V)=63665. Since64*1009+35>63665, an outside-V source and every positive lift remain outside V. Thus, if s reaches V, every lift has the same actual terminal and the same P-depth, with full arrival clock increased by6 per lift. If s never enters V, neither does its lift. Finiteness of an arbitrary V alone would not prove this guard.

For every integer N>=1,

s in S_N if and only if L_a(s) in S_(64N).             (2)

Indeed c_a<64, and the cutoff4N is even. This exact admission equivalence has no fractional endpoint relaxation.

## 3. Exact phase law with only sparse unit defects

For s in S_N, put c=c_(a(s)), p=4^j, b_j=(4^j-1)/3, and m_j=B^j(s). Then

B^j(L_a(s))=64m_j+(c-21)4^j+21.                       (3)

Its residue modulo3 is m_j+c. Both palettes rotate accordingly: L_(r+c)=R^c L_r and U_(r+c)=R^c U_r.

For K in{2,4}, define the exact crossing indicator

e_(K,j)(N,s)=1 if m_j<KN and 64KN<=B^j(L_a(s));
             0 otherwise.

**Phase-boundary theorem.**

A_(64N)(L_a(s))=R^c A_N(s)+epsilon_N(s),              (4)

epsilon_N(s)=-sum_j R^c[U_(m_j mod3)e_(4,j)
                         +(L-U)_(m_j mod3)e_(2,j)].  (5)

Moreover each s has at most one nonzero indicator in(5). Consequently epsilon_N(s) is either zero or exactly one Eisenstein unit: Q(epsilon_N(s)) is0 or1.

Let h_N=floor(log_4(12N)), equivalently the largest h with4^h<=12N. The exceptional source set E_N={s in S_N:epsilon_N(s) is nonzero} satisfies

|E_N|<=2(h_N-1).                                      (6)

Proof of the formula: write each source weight as

U_r*1_(m<4N)+(L_r-U_r)*1_(m<2N).

The lifted normalized point is larger, so only loss of one of these indicators occurs. Subtracting the rotated old weight gives(5).

Proof of sparsity: in source coordinate s, the old threshold is T=(KN-b_j)/4^j. The lifted threshold is T-delta_c(j), where

delta_c(j)=(c-21)/64+21/(64*4^j),   0<delta_c(j)<1.

The two branch intervals are nested inside the c=49 interval. Hence, for each j,K, at most one integer source, and therefore at most one primitive odd source, can cross. At j=0 the old threshold KN is an integer, so there is no crossing. If an old B-chain point lies below4N, then j<h_N. This gives at most2(h_N-1) candidate events.

Proof of the unit assertion: (3) gives

m_j < B^j(L_a(s))/64 <2m_j.

For the upper inequality, L_a(s)<128s and the positive affine B-iterate preserve strict inequality against128B^j(s). Any boundary crossing therefore forces N<m_j<4N. A B-chain has at most one point in that interval, since its next point is4m_j+1>4N. The two cutoffs cannot both cross at that point: the lower crossing requires m_j<2N, whereas the upper crossing requires m_j>2N. Both U_r and(L-U)_r have norm1. Thus every exceptional source has exactly one nonzero unit defect.

The receipt enumerates E_N without scanning S_N: at each j,K take the greatest odd integer s with4^j s+b_j<KN, test primitive membership and the exact crossing inequality. The width proof shows that no other candidate is lost.

### Actual boundary hostile

At N=2019,s=1009, a=2,c=49 and L_a(s)=64625. The unique crossing is j=1,K=2:

B(s)=4037<4038,
B(L_a(s))=258501>=128N=258432.

The old charge is(-1,-2), its expected rotation is(2,1), but the lifted charge is(1,0). The defect is(-1,-1), of norm1. All sources and clocks are native; this is not an arbitrary-weight counterexample.

## 4. Guarded three-source packets with no P-depth boundary lag

Fix a primitive s outside V in S_N. Let s_0=s, s_1=L_a(s), s_2=L_a^2(s), with the same stable branch a. Suppose

epsilon_N(s_0)=0 and epsilon_(64N)(s_1)=0.             (7)

At the single enlarged source scale M=4096N, inherited fixed-source transport and(4) give

A_M(s_0)=A_N(s_0),
A_M(s_1)=R^c A_N(s_0),
A_M(s_2)=R^(2c) A_N(s_0).

Their sum is zero, because c is1 or2 modulo3. The sources are distinct, and by(1) their full clocks differ by6 and12, so all three Nyquist signs agree. Their actual first arrivals have one terminal and one P-depth. A P-depth cut therefore sees none of the packet or the entire packet: this packet has no partial-arrival boundary at any finite P-depth horizon. Removing it preserves every arrived terminal vector and its cross terms with every other source. Its internal cross sum is minus its source norm3Q(A_N(s)).

No universal entry assumption is needed for an earlier fact: the three sources share their immediate P-target, and their edge-clock parities agree. Their transported signed sum is zero after **one P-edge**, whether or not that common target ever enters V. In a stopped system the outside-V guard ensures that all three do take that first edge. At depth0 their distinct physical source coordinates need not have zero full-state energy; the one-edge statement does not erase that initial energy.

This one-edge cancellation is a scoped statement about these packets, not a full moving-endpoint energy bound for the whole source universe. It does not show any individual source reaches V.

### Literal simultaneous packets

At N=256, M=1048576:

| Named sources | Common first terminal | Common P-depth | Full arrival clocks | Charges at M |
|---|---:|---:|---|---|
|1009,64625,4136049|47|1|6,12,18|(-1,-1),(1,0),(0,1)|
|1015,64995,4159715|571|2|4,10,16|(-1,-1),(0,1),(1,0)|

The second packet first merges at the actual moving target1523, outside V, with edge clocks1,7,13. It cancels there at depth1, and all three would first arrive at571 at depth2.

By contrast, using s=1009,N=2019 without the defect guard gives the same three named labels at scale4096N but total charge(0,-1), not zero. At the correct depth cut the residual has norm1 and must remain.

The outside-V guard is independently necessary. At M=4096 the three phase charges of {1,113,7281} sum to zero, but their first terminals are1,113,1, their depths0,0,1, and their actual arrival clocks0,0,14. Different stopped fibres cannot be merged to manufacture a packet.

### Full clock polynomial, not all-frequency deletion

If a good packet has arrival clock t for its first member and unsigned charge A, its complete clock polynomial is

z^t A[1+omega^c z^6+omega^(2c)z^12]
 =z^t A(z^6-1)(omega^(2c)z^6-1).                      (8)

This preserves the full integer clock information. In particular it vanishes at z=-1 and z=1, but is not the zero polynomial unless A=0. Equal P-depths do not mean equal clock cuts; before all clock thresholds are crossed, an integer-clock horizon may expose a partial packet.

## 5. Inherited inverse fibres and a no-reuse compiler

For an actual primitive target u and a in{1,2}, let

r_a(u)=least nonnegative residue of(2-u) modulo3, if a=1;
r_a(u)=least nonnegative residue of(1-u) modulo3, if a=2,
b_a(u)=(2^a B^(r_a(u))(u)-1)/3.

Then every source in P^(-1)(u) belongs to exactly one of the two disjoint chains

{L_a^k(b_a(u)): k>=0},  a=1,2.                        (9)

This is precisely the period-three form of the inherited inverse formula. Integrality forces the displayed r residue. For a=2 the source is1 modulo8. For a=1 and r>0 it is3 modulo8; at r=0, the allowed primitive residues of u give source residue3 or7. Thus no extraneous branch is introduced. Increasing r by3 is exactly the affine lift. Uniqueness follows from the native branch and actual strip exponent.

At a fixed scale M=64N, each admitted adjacent lift pair has its predecessor in S_N by(2). Fixed-source transport changes(4) into the same-scale relation

A_M(L_a(s))=R^c A_M(s)+epsilon_N(s).                  (10)

Remove V members from the chain prefixes, and split the remaining finite chain segments at bad edges s in E_N. Every resulting segment rotates exactly by R^c. Group consecutive members into blocks of three. The blocks are source-disjoint by construction; different chains cannot reuse a source. Each block outside V is a synchronous zero packet as above. At most two ungrouped members remain per segment.

There are at most2(h_N-1) bad edges globally. This is a sparse obstruction to the phase relation, not a small bound on the number of remaining chain segments: there can still be many native targets and many short chains. Their terminal vectors and mutual coherence are untouched obligations. The construction does not assert that a substantial fraction of the annulus budget is cancelled.

## 6. What a lift pays, and the exact remaining bill

For a good lift pair, let v be the old signed charge and w=R^c v the lifted signed charge, with common actual terminal and depth. Since c is1 or2 modulo3,

Q(v+w)=Q(v),   beta(v,w)=-Q(v)=-Q(w).                 (11)

Thus its own new-source diagonal is paid exactly by its internal old/new cross term. If Z is the rest of that terminal vector, the actual increment is still

Q(Z+v+w)-Q(Z+v)=beta(Z,w).                            (12)

The interaction with the rest is not paid by(11). A family of pairs must be source-disjoint, and cross terms between pairs remain unless another exact identity removes them.

A three-source packet can span two scale annuli, or have other old/new membership depending on its position. Its total zero does not entitle one to spend all its norms against one153N annulus. In particular, a mixed old/new zero packet changes the old baseline in the energy difference; the correction in the frozen annulus note, equation(D), remains mandatory.

The bounded result is therefore an exact native cancellation mechanism with sparse, explicit unit boundary errors and synchronized P-depths. The outstanding task is to control the untouched terminal coherence, bad-edge residuals and old/new baseline terms on unbounded scales. Neither this packet nor its immediate-target version proves the full153N annulus bound or universal Collatz entry.

## 7. Exact checks and reproduction

    python3 -B check_phase_lifts.py --write
    python3 -B -O check_phase_lifts.py

The only input is the adjacent `frontier.json`; no repository imports, third-party packages or network are used. Checks remain active under-O.

Declared finite controls:

- every original primitive source for every integer N=1,...,256, including zero charges;
- exact equality of the directly found exceptions and the sparse candidate compiler;
- independently summed original three-real-row weights;
- every primitive source1<=s<65536 reconstructed from its inverse-chain address;
- sparse exact exception controls at N=2019,4096,16384 and10^40+7, without a large orbit census;
- the two literal outside-frontier packets and their actual common moving/arrival targets;
- the native boundary and stopped-frontier hostiles;
- all twelve signed norm1/norm3 types, both rotations, pair payment and polynomial factorization.

Normal and optimized runs match at945113 explicit checks. The receipt hash is d1c9ce40b69d4952dee72b2f10b8b39ba7d960eaa6df8de5ba8e84d18a059de6. The manifest pins the sources and every deliverable. No theorem identifier or priority claim is assigned.
