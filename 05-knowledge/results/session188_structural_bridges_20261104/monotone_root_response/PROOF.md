# One fixed query detects every original-module family at an equal-profile cyclic root

**PROVED:** for any common complete child profile, the single core query a=A1,c=A2 detects every original-module XOR family containing the whole equal-profile cyclic root. There is no child-order bound for the whole-core statement. Whenever the baseline cyclic-quotient separator survives a fixed ambient context, the same query detects every such family there as well. In particular the separately certified small-child arbitrary-ambient results apply. This strengthens, rather than corrects, the earlier sufficient two-query battery.

## 1. Positive replacement lemma for an actual query

Suppose a physically named block B is a module of the old tournament and both ports see it uniformly. Then B is still a module after the two ports are attached. Compressing its maximal runs in an actual Hamilton word gives

    F=H(augmented tournament)=sum_h c_h q_h(B),
    c_h>=0.

The coefficients depend on the fixed exterior, including the actual ports, and not on the inside of B. Therefore replacing B by a same-name tournament B' with q_h(B')>=q_h(B) for every h cannot decrease this actual query value. Several disjoint child blocks may be replaced one at a time. This uses no approximation, signed coefficient extraction, or independence assumption for joint endpoints.

## 2. Root-containing original families

Let T=C3[A1,A2,A3] on three nonempty disjoint named child sets, with A1->A2->A3->A1. Assume every source child has the same complete ordered-cover profile q. Let F be any family of modules of this original T, with repetitions reduced modulo two, and assume the whole root is present. The simultaneous operation toggles each edge once for every selected original module containing its two ends. The proper selected modules may cross or nest.

Every proper module of T lies within one child: a module meeting two children must contain the third, after which any omitted name in either of the first two sees opposite directions to its occupied neighbors. Thus only the selected root flips the quotient. The resulting target is

    U=C3^op[B1,B2,B3]

on the same child sets.

The hereditary original-module monotonicity theorem in `../module_decomposition/PROOF.md`, Theorem1, applies to every induced child. It includes the whole-child reversal contributed by the selected ancestor. Hence

    q_h(Bi)>=q_h(Ai)=q_h     for every i,h.               (1)

Define the intermediate opposite-quotient core

    U0=C3^op[A1,A2,A3],

with the source child interiors unchanged. The unconditional equal-profile identity and positive-square theorem give, for the actual query a=A1,c=A2,

    F_(U0)-F_T=P(q)>0.                                  (2)

Both ports are uniform on every child in U0 and U. Apply the positive replacement lemma successively to replace Ai by Bi, using(1). It gives

    F_U>=F_(U0)>F_T.                                    (3)

This proves that the **one fixed query** detects every root-containing original-module family, even when scalar H has increased. There is no need to first test the all-zero scalar query. It also proves the quantitative bound F_U-F_T>=P(q).

The theorem does not require equal target profiles, a scalar-equality assumption, laminarity, transitive children, or a bound on descendant depth. Target profiles may be larger and unequal. It is detection relative to the stated source, not pairwise injectivity of an arbitrary original-module orbit.

## 3. Ambient consequence with exact scope

Put the source core into a distinguished root r of an arbitrary fixed outer tournament Q. Select only original modules contained in this core, including the whole core. Keep all exterior names and arcs unchanged. Every selected core module is also a module of the full source, because the entire core is a module.

Lift the same query physically: within the core a=A1,c=A2 and x->y; for each exterior z, give both ports the root's incidence, so z->r implies z->x,z->y and r->z implies x->z,y->z. Every individual child remains a module of the entire augmented tournament, even after replacing other children.

If the baseline opposite-quotient comparison satisfies

    F_(Q[U0])-F_(Q[T])>0,

the positive replacement lemma gives

    F_(Q[U])>=F_(Q[U0])>F_(Q[T]).                       (4)

Thus any independently proved baseline ambient separator with this child-uniform query automatically handles all root-containing original families inside that core. In particular:

- `../equal_small_ambient/PROOF.md` proves it for common child order1 through7.
- `../equal_nine_ambient/PROOF.md` extends the independently certified profile domain through9.
- The separately certified6+4 split argument in `../equal_ten_ambient/PROOF.md` provides order10.

Within those proved domains, the baseline bound is F_(Q[U0])-F_(Q[T])>=H(Q)P(q), so the same lower bound holds for every family target U. The proof here is the positive replacement implication; it does not silently assume an ambient sign theorem at larger child orders.

Selections involving exterior names or ancestors of the displayed core are excluded. Arbitrary disjoint exterior changes would need their own justified restoration argument and are not added here. Nor is unrestricted named-response reconstruction asserted.

## 4. Reproducible physical controls

The checker independently enumerates every proper nontrivial module of the9-name source C3[TT3,TT3,TT3], then all512 root-containing families. It retains the one same root-uniform query for **all** targets, including the296 scalar-increasing cases that the older packet detected with a different scalar query. It also checks all8 root-containing families of the pair-free source C3[C3,C3,C3]. Every family is tested both as a whole core and in a fixed cyclic exterior, with actual physical port incidences.

Every source and target query is counted by literal subset Hamilton DP. The complete target rows and values are retained, together with the family, profile domination checks, scalar status and quantitative bound. These finite controls supplement the all-order positive replacement proof, not replace it. No cubic or quartic response terms are dropped to form the actual query counts.
