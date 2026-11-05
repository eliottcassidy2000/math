# Five C3 scalar forms and a sharp physical source cutoff

Session198. **PROVED, conditional only on the named previously reviewed source/classification dependencies:** a reduction from eight to five C3 scalar forms, nonnegative source coordinates, sharp product bounds, and an improved complete source cutoff. **FINITE-EXACT:** the accompanying support computation. Universal odd-target coverage remains **OPEN**.

## 1. Inheritance and the physical cone

Use the actual rooted start-one sources and exact source grammar of
`dependencies/BILINEAR_C2.md`. A source has profile

    v=(D,E,U,d)=(H(V), e_c(V), R_c(V), H(V-c)),

with root-start count one. Let S be their profile set, o=(1,1,1,1), and S*=S−{o}. Equal profiles may be merged for scalar computations, retaining one physical recipe; they do not identify rooted graphs.

For every source, deleting the final root in an ending-root Hamilton path injects those paths into Hamilton paths of V-c. Hence E<=d. For a nontrivial source, its unique root-start path gives the root an outneighbor. Inserting the root immediately before the first outneighbor on any path of V-c injects those d paths into paths not ending at the root. Thus D−E>=d.

Cutting each Hamilton path immediately before the root injects D paths into the U reserved-root two-covers. There are also d covers whose reserved path is the root singleton. Their intersection with that injection has size E, since concatenation is a Hamilton path precisely when it ends at the root. Therefore U>=D+d−E, and consequently U>=2d.

It follows that the coordinates

    y=(D−E−d, U−D−d+E, E, d−E)

are nonnegative for every nontrivial source. The inverse is

    (D,E,U,d)=(y0+2y2+y3, y2, y0+y1+2y2+2y3, y2+y3).

In particular U=w.y for w=(1,1,2,2), and d=y2+y3>=1. These inequalities are a necessary cone, not a source-membership characterization. Every support computation below uses the physical forward source oracle, not the cone.

Concept board: exact source membership; positive bilinear images; singleton boundaries; scalar versus graph symmetry; local versus exact representability. The inherited hostile is a formal positive tuple without a source recipe. Passing a cone or congruence test never repairs that missing information.

## 2. Nonnegative C3 matrices and product bounds

Use the eight physical C3 forms A,B,C,T,J0,J1,J2,J3 and their domains from `dependencies/C3_TWO_SOURCE.md`. Substituting the inverse coordinates above yields y(v)^t Q_X y(w), where

    Q_A = [[14,10,21,20],[10,7,16,15],[21,16,33,31],[20,15,31,29]],
    Q_B = [[14,10,26,23],[11,8,20,18],[20,15,35,32],[19,14,34,31]],
    Q_C = [[23,18,44,40],[18,14,34,31],[44,34,80,73],[40,31,73,67]],
    Q_T = transpose(Q_B),
    Q_J0 = [[12,12,24,24],[12,12,21,21],[15,12,30,27],[15,12,30,27]],
    Q_J1 = [[12,12,21,21],[12,12,18,18],[18,15,32,30],[18,15,31,29]],
    Q_J2 = transpose(Q_J1),
    Q_J3 = transpose(Q_J0).

Every entry is nonnegative. Entrywise comparison with c*w*w^t proves, for two nontrivial physical sources,

    A,B,T >= 7 U_v U_w,
    C >= 14 U_v U_w,
    J0,J3 >= 6 U_v U_w,
    J1,J2 >= (29/4) U_v U_w.                     (1)

These are exact rational comparisons, not fitted inequalities. The verifier derives every matrix entry independently by evaluating the inherited forms on the four inverse-coordinate basis vectors and checks the comparisons with integer arithmetic.

The constants in (1) cannot be increased uniformly over physical sources:

- Let t_m=(1,0,2^m,1), realized by m universal-sink moves, m>=1. Writing a=2^m,b=2^n, we have A(t_m,t_n)=7ab+a+b−3 and C(t_m,t_n)=14ab+3a+3b−1. Their product ratios tend to 7 and 14 respectively as both m,n grow.
- B(t_1,t_n)=14*2^n+3, whose product ratio tends to 7. T follows by transposition.
- Let z=2^k+1 with k>=1 and p_z=(2z−1,z−1,2z,z). This is the source obtained by M_(2^k), then M_1, from the singleton. J0(p_z,t_n)=12z*2^n+3z. Its product ratio tends to 6 as n grows; J3 follows by transposition.
- J1(t_1,t_1)=J2(t_1,t_1)=29, attaining product ratio 29/4 exactly.

The displayed J0 formula is also a useful independent sharpness check rather than an appeal to an unrealizable cone extremum.

## 3. Singleton boundaries and the sharp common cap

The singleton is not in the nonnegative y cone and is handled separately. The only Type-I singleton-left cases are C and T. In y(w) coordinates,

    C(o,w) = 21 y0 + 16 y1 + 36 y2 + 33 y3 >= 16 U_w,
    T(o,w) = 12 y0 + 9 y1 + 15 y2 + 15 y3 >= (15/2) U_w.

For Type II the only singleton cases are J0(v,o) and J3(o,w), and

    J0(v,o)=T(o,v)=3(D−E)+9U−6d,
    J3(o,w)=T(o,w).

Thus they satisfy the same 15U/2 bound. Two nontrivial sources have U>=2, so the weakest product bound in (1) gives each source U<=H/12. Every source in every admissible C3 representation therefore satisfies

    U <= floor(2H/15).                            (2)

The singleton source itself obeys (2), since the allowed outputs have H>=15. This improves the inherited common cutoff floor(H/3).

The constant 2/15 is sharp even along infinitely many physical representations: for the p_z family above,

    T(o,p_z)=J0(p_z,o)=15z=(15/2)U(p_z).

Here k>=1 is essential: M_1 has an exceptional endpoint subtraction, so the formula is not asserted at k=0. The smallest equality example is separately t_1, giving H=15 and U=2. Sharpness is representation-wise; it does not claim that these targets have no other presentations with smaller sources.

## 4. Eight forms reduce to five scalar images

The matrices prove the all-profile scalar identities

    A(v,w)=A(w,v), C(v,w)=C(w,v),
    T(v,w)=B(w,v), J3(v,w)=J0(w,v), J2(v,w)=J1(w,v).

All domain restrictions must still be respected. For two nontrivial sources the transposes match admissible domains. A T singleton-left boundary equals J0(w,o), which is explicitly admissible. A J3 singleton-left boundary likewise equals J0(w,o). Therefore the complete strong C3 scalar spectrum is exactly

    A(S*,S*) union B(S*,S*) union C(S,S*)
       union J0(S*,S) union J1(S*,S*).             (3)

In C only its left source may be singleton; in J0 only its right source may be singleton. Unordered pairs suffice for A and for the nontrivial portion of C. B,J0,J1 retain ordered pairs. This is a scalar reduction, not a claim that the corresponding graphs or their rooted presentations are isomorphic.

## 5. Certified inverse and coverage computation

To decide a C3 target N, enumerate the exact physical source bank with U<=floor(2N/15), then test the five forms and permitted boundaries in (3), pruning two-source pairs only by (1). Every retained profile has an actual finite predecessor recipe. Equations (1)-(3), together with the inherited complete C3 theorem, prove soundness and completeness of this algorithm. No order cutoff is imposed.

For comparisons, C2 uses the inherited complete two-form oracle and its cutoff floor((N−1)/6). C1 uses the exact scalar readout D−E+U on every nontrivial source with U<=N. Indeed adjoining one vertex whose sole outneighbor is the root closes exactly one cycle and has that path count. Conversely the completed C1 theorem gives an outdegree-one vertex in every C1 graph; deleting it gives just such a source. Since D−E>=1, any source representing a C1 target at most N has U<N. This establishes the C1 cutoff, with `dependencies/C1_COMPLETE.md` as the only additional structural input.

The verifier compares the new C3 image against the inherited eight-form search through 5,001; checks every scalar certificate by a separate recurrence implementation; counts selected actual graphs by subset dynamic programming; and tests the new matrices and sharp families. It writes exact totals and explicit witness recipes. Run `python3 verify.py` and `python3 -O verify.py`; the two runs must agree. No assertion statement is used for verification conditions.

**Scope:** the certified finite support does not prove cofinite coverage. The result reduces exact membership cost and isolates physical bottlenecks, but does not supply an all-target descent or show that the realizable bilinear images cover every remaining odd integer.
