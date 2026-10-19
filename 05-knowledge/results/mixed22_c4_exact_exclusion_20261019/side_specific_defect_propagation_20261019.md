# Side-specific exterior budgets and propagation

Date: 2026-10-19. Status: proved sufficient refinements, with a separate exact census. This note does not change the frozen 56,630-core analytic cover or the running residual scan.

Use the mixed-(1,2), c=4 core Q, roots x,y, exterior A below both roots and B above both roots, and the one-defect exchange lemma. Every ambient vertex has incomparability degree at most seven. Write s(v)=7-d_Q(v).

## 1. Forced exterior sides

Both roots are above all A and below all B. Every saturated core vertex (d_Q=7) has the same two forced exterior directions. Let T be the roots together with all saturated core vertices, and define reflexive core closures

    F_A = {v in Q: t<=v for some t in T},
    F_B = {v in Q: v<=t for some t in T}.

By transitivity, every A element is below every F_A vertex, and every B element is above every F_B vertex. This includes all upper slaves in F_A because A<x<a_i, and the lower slave in F_B because b<y<B. It can force a side for additional nonsaturated common vertices via a saturated predecessor or successor.

If q_A(v),q_B(v) count the exterior elements incomparable with v on the respective sides, then

    q_A(v)+q_B(v)<=s(v),
    q_A(v)<=alpha(v):=0 if v in F_A, and s(v) otherwise,
    q_B(v)<=beta(v):=0 if v in F_B, and s(v) otherwise.

A core vertex cannot lie below an A element or above a B element: either relation would compare it with both roots, contrary to membership in Q. Thus there are no further exterior comparison directions to account for.

## 2. Sharper directed defects

For incomparable core u,v, let A_v={a in A:a<v} and B_v={b in B:v<b}. The exact decomposition is

    delta_P(u,v)=delta_Q(u,v)+|A_v\A_u|+|B_u\B_v|.

Since A_v\A_u is contained in A\A_u, and B_u\B_v is contained in B\B_v,

    delta_P(u,v)<=delta_Q(u,v)+q_A(u)+q_B(v)
                  <=delta_Q(u,v)+alpha(u)+beta(v).

Consequently two opposite bounds at most one imply an actual ambient balanced pair. Under the separate maximum-span majority-root hypothesis, either bound for x<z or z<y at most one contradicts the required strict >2/3 orientation for a common z. These are sufficient conditions; the separate alpha/beta bounds discard the joint constraint q_A+q_B<=s and other realizability information.

More generally, a bound at most one for u<v forces v before u in the >2/3 majority order of a poset assumed to have no balanced pair. Such forced edges can be combined with the poset relations and x<Z<y. A directed cycle contradicts transitivity of that majority order.

## 3. A maximal upper slave follows y in the majority order

Let a be maximal among the two upper slaves. Its core upset is empty; all actual successors therefore lie in B=U(y). Also A is contained in D(a), while D(y)=A union {b}. Hence

    delta_P(a,y)=|D(y)\D(a)|+|U(a)\U(y)|<=1.

Thus p(a<y)<=2/3. Under the no-balanced-pair assumption this forces p(y<a)>2/3. This orientation is independent of all exterior choices. Similarly, if the lower slave b is below at least one upper slave, then delta_P(x,b)<=1 and the majority order must place b before x.

## 4. Stronger zero-increment rule

The earlier boundary-tree proof needs only all A below a letter for a prefix increment to vanish, and only all B above a letter for a suffix increment to vanish. Therefore every prefix increment whose word intersects F_A is zero, and every suffix increment whose word intersects F_B is zero. As before, all increments after the first such prefix letter, and at or before the last such suffix letter, vanish. This is stronger than deleting words containing a saturated vertex on both sides. The current residual scan uses the older, weaker saturation rule, so its successful certificates remain valid without any change.

## 5. Exact finite observation

The standalone finite count script `count_side_defect_closure_20261019.py`, using the bundled morphology generator, checks all 79,121 degree-valid cores:

- 49,812 satisfy the two-sided balanced-pair criterion;
- 56,637 satisfy the root/common directed maximum-span criterion;
- 58,466 satisfy at least one;
- 20,655 remain under this refinement.

Testing all forced one-sided majority edges gives exactly the same 58,466 cyclic cores in this census; it adds no further cores beyond those two direct tests. This equality is a finite observation, not a general theorem about forced majority graphs. Compared with the already frozen 56,630-core analytic cover, the refinement adds 1,836 cores. It is kept as a separate sufficient refinement and does not replace any recorded scan or proof receipt.
