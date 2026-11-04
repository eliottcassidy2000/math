# Exact reduction to the active strong component

**PROVED:** localized two-port queries factor exactly through the quotient's strong component containing the replaced module. **FINITE-EXACT:** complete interpolation controls on 75 named rooted quotients with two core choices, and a literal eleven-vertex quartic amplification witness. **OPEN:** the strong prime twelve-core hidden-converse residual and the general no-hidden-module converse. No exact full-response collision is claimed.

This is a structural reduction of the ambient search, not a solution of that residual. Its elementary SCC argument is stated explicitly; no literature-priority claim is made. The inherited twin obstruction C2 already closes the homogeneous-pair branch, and the inherited T4 cubic theorem remains limited to whole converse.

## 1. An exact ambient reduction, stronger than residue descent on this scope

Let Q be a tournament with distinguished vertex z. Let K be the strongly connected component (SCC) of Q containing z. For a nonempty tournament M on fixed physical names, write

    T = Q[M],       B = Q[K][M],
    h = product over SCCs D of Q other than K of H(Q[D]).

The distinguished root is replaced by M, including all its names; B retains its literal names as a subset of T. No assumption that M or B is strong is needed. Every h is a positive odd integer.

For an arbitrary named attachment query (a,c) on B, attach x,y with x->y inside the position occupied by K in the SCC order: every old vertex in an earlier SCC beats both ports, and both ports beat every old vertex in a later SCC. Keep the requested attachment bits on B. Denote the resulting actual ambient masks by L(a),L(c).

**Theorem.** For every query, exactly

    F_T(L(a),L(c)) = h F_B(a,c).                       (1)

In particular, for any two same-named cores M,N in this same quotient,

    F_Q[M] = F_Q[N]  implies  F_Q[K][M] = F_Q[K][N].   (2)

These equalities mean the complete named two-port functions. Equation (2) is an implication only; we do not claim its converse. No inference that H,S,E,C alone determine F is used.

**Proof.** The SCCs of a tournament are totally ordered, with every arc from an earlier component to a later one. This follows since opposing directions between components would make a directed cycle in the condensation; a tournament's acyclic condensation is transitive. A Hamilton path in an ordered sum A1->...->As must visit the blocks in that order, in one contiguous run each: after entering a later block it cannot return to an earlier one. Conversely concatenating an arbitrary Hamilton path of each nonempty block gives a unique Hamilton path of the ordered sum. Therefore H(A1->...->As)=product H(Ai).

The port assignment above makes the actual extended ambient tournament exactly the ordered sum of the untouched earlier SCCs, the extended block B+x+y, and the untouched later SCCs. The same block-product argument applies even if B+x+y is not strong. Its other factors multiply to h, giving (1). Apply this formula to both M and N with identical lifted masks. Full ambient equality gives h times every active-query difference equal to zero; h>0 permits exact integer cancellation. QED.

The theorem extends verbatim to any fixed number of ports and any prescribed tournament among them. The old outside SCCs still see every port uniformly, and the same product proof applies.

## 2. Consequences and precise limits

- A hidden internal module converse in any ambient quotient is already hidden in the active quotient SCC. Searching unrelated earlier/later SCCs is unnecessary for existence of a hidden example.
- If z itself is a singleton SCC, then B=M. In this scope ambient invisibility implies *exact* intrinsic full-response equality, not just the inherited universal modulo-four descent. In particular, any intrinsic separating query lifts to every such ambient placement with its difference multiplied by h.
- More generally, if any query on B gives a nonzero integer difference d, the lifted ambient query gives h*d, with the same sign. Ordered exterior padding cannot cancel that witness.
- When the active quotient is strong with at least three vertices, this theorem stops at Q[K][M]. It does **not** cancel that strong context or prove exact descent to M. That is the obstruction left by this reduction.

The theorem allows a prime twelve-core to remain in an arbitrary strong quotient. It does not exclude such a core, prove universal named reconstruction, or manufacture a full-response collision. Existing residue tests and C2 remain necessary additional constraints.

## 3. Quartics persist under exact localization

Use the strong rooted quotient

    R=(8,1,3,6),  root0,

and replace root0 by the pair0->1. Give its other vertices physical names2,3,4. The active tournaments before and after reversing that pair are

    B =(18,16,3,7,12),
    B'=(16,17,3,7,12).

These are the S186 strong five-vertex quartic hostile with a documented change of labels, not a new minimum-order witness. Its complete named difference has

    [a_2 a_4 c_2 c_1](F_B-F_B') = 2.

Now put a directed triangle before B and another after B, keeping B's names0..4, the earlier triangle's names5..7, and the later triangle's names8..10. All arcs point between these three blocks in that order. The resulting rows are

    T =(1810,1808,1795,1799,1804,1887,1951,1855,512,1024,256),
    T'=(1808,1809,1795,1799,1804,1887,1951,1855,512,1024,256).

For lifted masks, set both port attachment bits to1 on names5,6,7 and to0 on names8,9,10. Equation (1) gives

    (F_T-F_T')(L(a),L(c)) = 9(F_B-F_B')(a,c).

Hence the same coefficient in this *specialized* named ambient polynomial is18. It is not presented as a coefficient with all exterior bits zero. The sixteen actual ambient query masks and their integer counts' differences are retained in hostile.json. Literal Hamilton subset DP independently verifies their finite difference18.

The complete B-vs-B' polynomial has degree4. Separately the complete B-vs-B^op polynomial has degree3, as T4 requires. Thus exact SCC localization respects rather than erases the distinction between internal module reversal and whole converse. Interpolation never drops the quartic basis.

## 4. Verification and reproduction

verify.py enumerates all75 labelled tournaments on quotient orders1..4, each with physical root0, with each of the two cores (2,4,1) and (6,4,0). Relabelling the root to0 covers every rooted quotient up to relabelling; the finite packet does not count distinct rooted isomorphism classes. These are150 precisely listed cases and34476 exact query checks, not a reconstruction census.

For every case the checker computes both sides of (1) independently by literal endpoint subset DP on the two actual port-extended tournaments. It evaluates every query with at most two active bits in each port, plus three high-weight controls. This is a complete interpolation basis because the full response has bidegree at most(2,2): a Hamilton word has at most two incidences with each port. All arithmetic is unbounded Python integer arithmetic and checks use exceptions, not removable assertions. Every quotient, component list, active core and factor is retained.

hostile.py independently constructs the eleven-vertex example, evaluates all sixteen finite-difference terms directly, and retains both complete active difference polynomials. The script imports only the local standard-library checker. Normal and optimized runs agree byte for byte.

    python verify.py > receipt.json
    python -O verify.py > receipt.optimized.json
    python hostile.py > hostile.json
    python -O hostile.py > hostile.optimized.json

Inherited sources: S186 converse/PROOF.md (C2 integration and T4 boundary), S186 structure/PROOF.md (the five-vertex hostile), and two_port_twelve_principles_20261103/PROOF.md (general modulo-four descent and exact whole-converse cubic). This packet does not replace their hypotheses or extend their finite residual census.

## 5. Compatibility with the inherited three-state gate

For S185's physical gates G_j(r), the inherited count is

    H(G_j(r)) = 4*3^r + 9*2^r - 3 + 2j,  j=0,1,2.

Putting any fixed nonempty ordered prefix P and suffix S around the whole gate gives

    H(P -> G_j(r) -> S)
      = H(P)H(S)(4*3^r + 9*2^r - 3 + 2j).

More importantly, placing both ports in the gate position gives its *entire named response* multiplied by H(P)H(S), not just this scalar clock. Thus the incoming two-query twin sign, every retained whole-gate converse cubic, and every internal-clone quartic survive this ordered padding without cancellation. The unpadded global-converse cubic is the comparison on the active gate; reversing just the gate inside padding is not the global converse of the padded tournament. No cubic claim is made for its unrestricted ambient attachments. This is an all-order consequence of (1), not an additional gate census or a new digit operation.
