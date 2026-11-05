# Sharp two-thirds ratio for every weak deletion

Session 194, 5 November 2026. **PROVED, independently reviewed; extension kept separate from the first frozen theorem.** Dedicated exact controls and the snapshot-bound independent audit passed; see the companion structural-cut-audit final report.

## 1. Statement

If T is a strong tournament of order at least three and S=T-v is weak, then

    Delta(T,v)=H(T)-H(S) >= 2H(S),
    H(T)>=3H(S),
    Delta(T,v)>=(2/3)H(T).                              (1)

The constant is sharp: T=C3 gives H(T)=3, H(S)=1 and Delta=2. Thus no strictly larger universal ratio is possible. This strengthens the strict one-half inequality in the frozen `PROOF.md`; it does not change the remaining strong-deletion sector or construct the full Hamiltonian-path spectrum.

If S has at least three SCCs, the frozen Theorem 1 already gives (1), because its multiplier 2^(k-2) is at least two. Only two SCCs require a new argument.

## 2. A local cut-response compensation inequality

Let C be any nonempty tournament, let x be a new vertex beating a nonempty set X subset C, and put

    h=H(C), d=H(C+x)-h,
    R_X=sum_(a in X)R_a(C), e_X=sum_(a in X)e_a(C).

Then

    d+R_X+e_X>=2h.                                    (2)

If |X|>=2, every R_a>=h by the path-splitting injection, so R_X>=2h, while d>=0 by OCF and e_X>=0. If X={a}, every root cycle has the form x followed by a nonempty path in C starting at a, then back to x. It must use an even number of C vertices. Conversely every such even path closes to x, since x loses to every C vertex except a. The rooted parity identity therefore gives

    d=R_a-e_a.

Hence d+R_a+e_a=2R_a>=2h. This also treats a singleton C: d=0,R_a=e_a=h=1.

By reversing all arcs, if x loses precisely to a nonempty set Y in a nonempty D, then

    d_D+L_Y+s_Y>=2H(D).                               (3)

No strongness assumption is needed in these two local inequalities. The only monotonicity used is induced deletion monotonicity from OCF, not orientation monotonicity.

## 3. Exact two-component decomposition

Write S=C->D for its two strong components. Set

    X=N^+(v) intersect C,       Y=N^-(v) intersect D.

Both are nonempty because T is strong. Use h=H(C), j=H(D), and

    d_C=H(C+v)-h, d_D=H(D+v)-j,
    R=R_X(C), e=e_X(C), L=L_Y(D), s=s_Y(D).

All these local tournaments retain the actual incidences of v. Then

    Delta(T,v)=j d_C+h d_D+RL+es.                      (4)

Indeed a root cycle either lies entirely in C+v, lies entirely in D+v, or visits C followed by D. The first two classes give j d_C and h d_D, because untouched SCCs independently contribute their complete H. In the crossing class its initial C vertex is a unique member of X, its final D vertex a unique member of Y, and its C/D paths have equal vertex-count parity. Applying the two endpoint parity identities gives exactly RL+es. These are disjoint and exhaustive classes: a directed path in C->D cannot return to C after entering D. All other odd cycles lie within the unused parts of C or D, so their contribution has already been included.

By (2)--(3), d_C>=2h-R-e and d_D>=2j-L-s. Substituting these valid lower bounds, even if one right side is negative, gives

    Delta(T,v)
      >=2hj+(R-h)(L-j)+(h-e)(j-s)
      >=2hj=2H(S).                                   (5)

The last inequality uses R>=h, L>=j, and 0<=e<=h, 0<=s<=j, since e and s count subsets of all Hamiltonian paths classified by their terminal or initial vertex. This proves (1) with all signs explicit.

## 4. Exact equality classification

For two SCCs, put

    sigma_C=d_C+R+e-2h>=0,
    sigma_D=d_D+L+s-2j>=0.

The preceding proof gives the exact nonnegative slack identity

    Delta-2hj=j sigma_C+h sigma_D
                 +(R-h)(L-j)+(h-e)(j-s).                (6)

Equality in (1) holds **if and only if** S has exactly two SCCs C->D and all the following hold:

* X={a} and Y={b} are singleton endpoint sets;
* R_a(C)=H(C), L_b(D)=H(D);
* at least one of C,D is a singleton.

For the two-component necessity, if |X|>=2 then R>=2h and e>0, because every vertex of a strong C is a Hamiltonian-path endpoint. Thus sigma_C>0. Consequently equality forces X singleton; its exact local identity makes sigma_C=2(R-h), so R=h. The argument for Y is the reversed version. The remaining term in (6) is `(h-e_a)(j-s_b)`. On a nontrivial strong tournament every other vertex is also an endpoint, so e_a<h and s_b<j. Equality therefore forces a singleton SCC. Conversely all displayed conditions make every term in (6) zero.

For k>=4 the frozen multiplier gives Delta>=4H(S)>2H(S). For k=3, equality in that bound would require Q(C2)=2H(C2). Every nonempty proper subset contributes positively to Q, so C2 must be a singleton w. Fix the endpoint a,b used in the frozen proof. If v->w, the triangle v,w,b omits C1; if w->v, the triangle v,a,w omits C3. Either is an extra positive marked-cycle contribution disjoint from the selected family, all of whose cycles visit both C1 and C3. Thus the inequality is strict at k=3 as well. This uses no additional cut-capacity identity.

### Sharp equality at every order other than four

The C3 example gives order three. For every m>=1 take D on vertices a,b,z1,...,zm with

    a->b->T_m->a,

where T_m is the transitive z-chain. Then D is strong and

    H(D)=R_a(D)=2^m+1.

To count H, the paths using a->b correspond to the arbitrary subset of chain vertices placed before a, with the rest after b, giving 2^m. The sole path not using a->b is b,z1,...,zm,a. For R_a, the singleton path a has complement count one; every other path starting a is a,b followed by any subset of the chain, whose complement is transitive. This gives 1+2^m.

Take S=D->{d}, and add v with sole out-neighbor a. Then v reaches D via a and reaches d, while every old vertex reaches v. Thus T is strong. The equality conditions hold with C=D and terminal SCC {d}; hence

    |T|=m+4, H(S)=2^m+1, H(T)=3(2^m+1), Delta=2(2^m+1).

This realizes equality at every order n>=5. At order four a weak deletion has three singleton SCCs or a strong order-three SCC; only the former is weak, and k=3 is strict. Hence no order-four equality occurs. These are positive controls for the exact ratio, not a spectrum-filling family.

## 5. Which weak parents can perform a literal ternary step?

**Theorem (exact weak-parent form for ternary births).** Suppose S is weak, h=H(S)>=3, T=S+v is strong and

    H(T)=3h+d,       d in {2,4}.

Then S has exactly two SCCs. The possible presentations are precisely:

1. One SCC is a singleton and the other is a strong C of count h. If S=C->{z}, the root v has a unique out-neighbor a in C, loses to z, and `R_a(C)=h+d/2`. If S={z}->C, reverse this statement: v beats z, has a unique in-neighbor b in C, and `L_b(C)=h+d/2`.
2. The nontrivial component is C3 and the relevant selected endpoint set in that component has size two. Here h=3,d=4,H(T)=13.
3. Both SCCs are C3 and X,Y are singleton endpoint sets. Here h=9,d=4,H(T)=31.

Cases 2 and 3 are exceptional **presentations**, not claims that the resulting scalar H values have no other presentation. Every listed set of conditions constructs a strong tournament of the displayed count.

**Proof of the SCC reduction.** The target is exactly `Delta(T,v)=2h+d<=2h+4`. If k>=4 the frozen bound gives Delta>=4h>2h+4, since h>=3. If k=3 and the middle SCC C2 is nontrivial, its order m is at least three. The identity `sum_a R_a(C2)=Q(C2)-H(C2)` and R_a>=H(C2) give Q(C2)>=(m+1)H(C2)>=4H(C2). The selected root-cycle contribution alone is then at least4h, again impossible.

For k=3 with middle singleton w and both endpoints nontrivial, the extra root triangle used in the equality proof has complement containing the entire opposite nontrivial SCC, whose H is at least three. It therefore contributes at least six beyond the selected family, already at least2h. The only remaining k=3 shape, up to reversal, is S=C->{y}->{z} with C strong nontrivial, and z->v. Put X=N^+(v) intersect C, R=R_X(C), e=e_X(C), and d_C=H(C+v)-h. If y->v, exact root-cycle decomposition gives `Delta=d_C+3R+e`: cycles within C contribute d_C, cycles from C to y contribute R+e, and cycles from C to z, optionally using y, contribute2R. This is at least4h by the local compensation inequality. If v->y, the decomposition instead gives `Delta=d_C+2R+2h>=4h`: the last term is the triangle v,y,z with an arbitrary cycle collection in C. Both are impossible. Thus k=2.

If both SCCs C,D are nontrivial, their counts h_C,h_D are at least three. If |X|>=2, then sigma_C>=e_X>=2, so the slack identity (6) already contributes h_D sigma_C>=6; similarly for Y. Thus slack two or four requires X,Y singleton. If R_a>h_C or L_b>h_D, the corresponding sigma is at least two, again costing at least six. Therefore R_a=h_C,L_b=h_D. The sole remaining slack is `(h_C-e_a)(h_D-s_b)>= (|C|-1)(|D|-1)>=4`. Slack two is impossible. Equality at four forces |C|=|D|=3. Both factors are C3, giving exactly case3.

**Proof of the endpoint-profile classification.** By reversal suppose S=C->{z}. Strongness forces z->v and a nonempty X=N^+(v) intersect C. Formula (4) says

    H(T)-3h=d_C+R_X+e_X-2h.

For X={a}, the exact local identity d_C=R_a-e_a reduces this to `2(R_a-h)`, proving case1 in both directions.

If X=C, then d_C=0 and the slack is Q(C)-2h>= (|C|-1)h>=6, excluded. If X is proper and |X|>=2, then C+v is strong and d_C>=2, while e_X>=|X| and R_X>=|X|h. Slack two is impossible; slack four forces |X|=2, d_C=2, e_X=2 and R_X=2h. OCF shows that d_C=2 forces precisely one triangle through v: a strong root lies in a triangle, and each such cycle contributes at least two. The inherited single-triangle minimum g(n)>=4 at n>=5 rules out |C+v|>=5. Consequently C=C3 and X has size two, which gives case2. Conversely direct C3 values h=3,R_X=6,e_X=2,d_C=2 give H(T)=13. The corresponding cross-side incidences ensure strongness in all cases. QED.

This is an exact classification for a *single induced vertex addition from a weak parent*. It does not restrict reconstruction, old-arc changes, or births from a strong parent. It replaces an arbitrary neutral-padding search by the two small rooted surplus tests `R-H=1,2` (or their duals), with only the two explicitly described exceptional presentations. It does not prove that every parent count has a core with those surpluses.

## 6. Order consequence and scope

The cited Busch minimum gives the stronger order-only bound

    Delta(T,v)>=(2/3)B(n),
    Delta(T,v)>=2 ceil(B(n)/3)

where the second line uses that Delta is even. Both inequalities hold in the weak-deletion sector only. They beat the inherited g(n) at every n>=4, but do not determine its unrestricted minimum.

For the spectrum program, this says a strong birth whose old tournament is weak costs at least twice the old count. It supplies an obstruction and a positive response decomposition, not a tournament for an arbitrary target. In particular it does not prove prime-lane ternary closure.

Reproduce dedicated exact controls with `python3 verify_stronger.py`. The primary `verify.py` remains frozen and unchanged; the extension imports its explicit finite helper definitions, not any external package.
