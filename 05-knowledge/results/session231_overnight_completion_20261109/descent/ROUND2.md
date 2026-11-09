# Round2: zero-defect face and the sharp low-growth boundary

Status: PROVED; independent [audit PASS](../descent_audit/ROUND2_AUDIT.md). This is an all-order structural barrier with explicit equality families. It does not establish a coverage sector or Chan–Pak completeness.

## Provenance and exact advance

Write (h,e,u) for a rooted tournament (A,a), and δ=u−h. The inherited bound δ≥r(r−1), where r is the root's outdegree, is THM-4873 (`01-canon/theorems/THM-4873-end-count-lemma-and-c2-digit-heredity-the-strict-ternary-parent-dies-at-the-grandchild.md`, part (ii), item D); the one-exit digit identity is THM-4862 (`01-canon/theorems/THM-4862-one-exit-extension-law-the-readiness-digit-is-the-cycle-count-of-the-entry-attachment-fibonacci-bijection-ready-floors-through-order-10-and-the-full-digit-c2-sweep-realises-every-odd-number-in-307-26551.md`, part b). The nonsink endpoint bound e≤(h−1)/2 and singleton-join floor (3h+1)/2 are THM-4876 (`01-canon/theorems/THM-4876-deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak.md`, statement 7).

THM-4876 audit D reported only FINITE-EXACT through order9 that the floor is attained only at h=3, and that for h≥5 the least count is at least (3h+3)/2. The proof below establishes both statements at EVERY order and classifies equality in the improved floor. The underlying familiar reversed-extreme-arc family is inherited, not a new count family. The sharp face inequality and complete equality selection are the present synthesis.

Board change: the near-extreme interface is far thinner than the numerical endpoint cone suggests. Minimum growth is not an all-target selector: its large equality outputs are all divisible by3.

## 1. Exact sole-exit identities

If a has unique outneighbor w, set B=A−a, d=H(B), q=e_w(B), and R=u_w(B). Directly partitioning paths and covers at a gives

    e=d−q,   h=d−q+R,   u=d+R,   δ=q.                  (1)

For e, append a to any B-path except those ending at w. For h, paths ending at a contribute e; all other paths split immediately before a→w into a possibly empty first B path and a w-starting second B path, contributing R. The u-covers whose rooted path is just a contribute d; the others are exactly the R covers. These identities agree with the inherited one-exit laws.

## 2. Complete zero-defect face

Suppose a is a nonsink and δ=0. The outdegree bound forces r=1, so (1) applies with q=0. A vertex has zero endpoint count iff it is outside the final SCC (Camion provides a spanning path ending at every vertex of a nontrivial final SCC, and the singleton case is immediate).

Therefore B splits uniquely as P→C, where P ends in the SCC containing w and C is a nonempty suffix. Put p=H(P), c=H(C), z=u_w(P), and Q=Q_C(2), the ordered two-bin path-cover count allowing empty bins. Concatenation in the two bins gives R=zQ. Thus

    (h,e,u)=(pc+zQ, pc, pc+zQ),   C nonempty.          (2)

Conversely every such physical P,C, with w in the final SCC of P, constructs a zero-defect nonsink root by adding a→w and having all other vertices beat a. This is a physical characterization, not just a necessary cone inequality.

Since z≥p and Q≥2c,

    e≤h/3.                                           (3)

Equality holds iff z=p and C is a singleton. Indeed Q>2c whenever C has at least two vertices, because at least one cover has both bins nonempty. The sharp boundary can have arbitrary physical P with zero defect; no assumption of all scalar p is made.

A useful further exact slice is e=1. Then p=c=1, so P and C are transitive. (H=1 iff transitive: a nontransitive tournament contains a directed triangle with3 Hamilton paths, and inserting each remaining vertex injectively extends the Hamilton paths.) As w is in the final SCC of P, it is P's sink and z=1. If t=|C|≥1, then

    (h,e,u)=(2^t+1,1,2^t+1).                          (4)

Graphically, after removing a universal transitive prefix, this is the transitive tournament on t+2 vertices with its first-to-last arc reversed, rooted at its former last vertex. Conversely these graphs supply every t≥1. No other scalar occurs on this physical nonsink face with e=1.

## 3. All-order sharpening of the singleton-join floor

Join (A,a) to the singleton by the proper one-arc operation. Its count is

    L=2h+δ−e.

For a nonsink, let b=(h−1)/2−e≥0. Then

    L=(3h+1)/2 + (δ+b).                              (5)

Equality in the old floor forces δ=b=0. By (3), (h−1)/2≤h/3, so h≤3. At h=1 a nonsink has e=0, δ≥1, so it does not attain equality. Thus equality holds iff h=3 and e=1,δ=0. The zero-defect characterization then forces a transitive prefix followed by a directed triangle, rooted in that triangle. In particular the output is5.

For every h≥5,

    L≥(3h+3)/2.                                      (6)

This promotes the finite-only audit D observation to an all-order theorem. It uses actual physicality at the zero-defect face, rather than merely strengthening a numerical cone.

## 4. Complete equality in the improved floor

For h≥5, equality in (6) means δ+b=1. There are exactly two possibilities.

**A. δ=0, b=1.** Then e=(h−3)/2 and (3) gives h≤9. Since h is positive odd, only h=5,7,9 remain; h=7 would require e=2, impossible because zero defect forces e=pc odd. At h=5 the state is (5,1,5), supplied by (4) with t=2. At h=9 it is (9,3,9), and equality in (3) means C is a singleton and P has (H,u)=(3,3); the root w lies in its final SCC. Here P may have a sink root, so it need not be strong or merely a transitive prefix followed by a rooted C3. Retaining the exact physical conditions H(P)=u_w(P)=3 and w in the final SCC avoids that invalid graph simplification. Formula (2) is the complete description of this case.

**B. δ=1, b=0.** Again the outdegree bound forces one exit. In (1), q=1 and e=d−1. The equation e=(h−1)/2 gives h=2d−1, and (1) then forces R=d. Thus B rooted at w has zero defect and endpoint count1. The sink possibility d=1 would give h=1, excluded. Applying (4) gives d=2^t+1 and

    (h,e,u)=(2^(t+1)+1, 2^t, 2^(t+1)+2), t≥1.       (7)

Conversely all these states occur. Their normalized graph is the reversed-extreme-arc transitive tournament on t+3 vertices, rooted at the penultimate vertex. A universal transitive prefix may be added.

Consequently the equality INPUT scalars for h≥5 are exactly

    h=2^k+1, k≥2,

with the extra zero-defect interface at h=9 (and the zero-defect h=5 interface, in addition to (7)). The equality OUTPUTS are exactly

    L=3·2^(k−1)+3, k≥2.                             (8)

All are composite. In particular, for a prime singleton-join target L>5 with source h≥5, equality in (6) is impossible; hence h≤floor((2L−5)/3). More importantly, the sharp lower-growth boundary is lacunary and divisible by3. It cannot provide physical predecessor selection for an arbitrary prime target. This is a structural limitation of near-extreme singleton joins, not of unrestricted two-sided joins or tournaments.

## 5. Transfer to arbitrary two-sided joins

For any second physical triple (k,f,v), the two-sided count satisfies

    N − k(2h+δ−e) = u(v−k)+e(k−f) ≥0.                (9)

Thus every nonsink input with h≥5 gives 2N≥(3h+3)k, strengthening the Round1 bound (3h+1)k. Equality forces the source to lie in the complete list above; since e>0 on that list, (9) also forces v=f=k, so the second port is a sink. Consequently equality targets are composite. If both ports are nonsinks and h≥5, equality is impossible. This is an immediate exact interface transfer, not a new coverage claim.

## 6. Controls and remaining direction

`round2_verify.py` exhausts all33,867 labelled rooted tournaments through order6, checks the sole-exit identities, the complete zero-defect numerical formulas, (3), both floor equality implications, and all predicted equality scalar states in that universe. It also constructs reversed-extreme-arc representatives for a larger parameter range. Exact finite controls validate the implementation; the arguments above supply the all-order quantifiers.

The next positive avenue must leave this extreme face, use defects increasing with the target, or use genuinely two-sided joins with varying physical triples. A finite computation finding small deficits cannot establish an all-target branch. No unrestricted image hole or new infinite arithmetic progression is claimed here.
