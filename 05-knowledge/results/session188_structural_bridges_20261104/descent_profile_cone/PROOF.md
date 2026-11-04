# Palindromic descent profiles separate equal cyclic children through order twelve in every fixed ambient

**PROVED by the exact finite coefficient certificate retained here:** for every common child order `1<=n<=12`, reversing a cyclic quotient of three fixed named child tournaments having the same complete ordered-path-cover profile is detected by one physical two-port query in every fixed ambient tournament. Target interiors may change while preserving that common profile. The argument uses the entire cone of nonnegative palindromic descent distributions, not a Hamilton bound or the smaller Stirling-profile cone.

## 1. The exact word-and-cut transform

For a tournament A on n nonempty physical names, let c_d count all permutations `v1,...,vn` having exactly d backward adjacencies, where an adjacency is backward when `v_(j+1)->v_j`. Thus `0<=d<=n-1`, `c_d>=0`, and `sum_d c_d=n!`.

Reversing a word replaces every one of its n-1 adjacency orientations by its opposite. It is an involution on all permutations and gives the exact symmetry

    c_d=c_(n-1-d).                                      (1)

Let q_k(A) count ordered spanning covers by k nonempty directed paths. Concatenating the k path slots and recording the k-1 cut positions is a bijection with permutations and cut sets containing every backward adjacency. A word having d backward adjacencies allows exactly

    binom(n-1-d,k-1-d)

such cut sets. With binomial coefficients interpreted as zero outside their usual range,

    q_k(A)=sum_(d=0)^(n-1) c_d binom(n-1-d,k-1-d).       (2)

No extra factorial occurs: the cover slots are already ordered in the concatenated word. Set q_0=0 and q_k=0 for k>n. The transform is triangular with diagonal one, hence loses no descent information.

For `0<=d<=floor((n-1)/2)`, define a nonnegative finite sequence

    g_(n,d)(k)=sum_(e in {d,n-1-d}) binom(n-1-e,k-1-e),  (3)

where the set has only one member at the central index. Equations (1)--(2) give

    q(A)=sum_d w_d g_(n,d),       w_d=c_d>=0.            (4)

Every tournament profile belongs to this cone. Conversely, its formal generators need not be realized by tournaments, and no such claim is needed. We deliberately discard the mass constraint `sum c_d=n!` and all additional realizability constraints, making the certificate stronger on its stated domain.

## 2. Full ordered-cover contrasts

Use three children with a common complete profile q in the source cyclic order `A1->A2->A3->A1`, and reverse that quotient on the same physical sets in the target. Attach physical ports `x->y`, with root masks `a=1_(A1)` and `c=1_(A2)`, where `a_v=[v->x]`, `c_v=[v->y]`.

For their port-augmented cores, the exact all-cover formula from `../covered_ambient/PROOF.md` is

    D_h(q)=q_h(T+)-q_h(U+)
          =(h/3) L_q[K E0^(h-1)/(1-XYZ)^(h+1)],         (5)

where `L_q(X^i Y^j Z^k)=q_i q_j q_k`, and

    E0=X+Y+Z+XY+XZ+YZ+3XYZ,
    K=X^2Y^2+X^2Z^2+Y^2Z^2
      +X^2Y+X^2Z+XY^2+XZ^2+Y^2Z+YZ^2
      -X^2YZ-XY^2Z-XYZ^2-6XYZ.

The rational kernel is symmetric in its three colours. Formula (5) counts both same-slot and different-slot port placements in the full ordered-cover count. It is not a truncation of the named response.

Substituting (4) into the cubic (5), the coefficient of the weight monomial with unordered index triple `d1,d2,d3` is exactly

    (h/3) [3! / product_j multiplicity(j)!]
       L_(g_d1,g_d2,g_d3)[K E0^(h-1)/(1-XYZ)^(h+1)].     (6)

The symmetric kernel makes all placements equal; the displayed factor accounts for the distinct placements. Repeated generator indices are not treated as distinct slots.

## 3. Exact certificate and separate strictness

For every n from one through twelve, `verify.py` computes every coefficient (6) at every `h=2,...,3n+2`. All are nonpositive. It retains the complete integer lists, including zeros. Powers are truncated only when a variable exponent exceeds n, which is safe because each generator vanishes above n. After the displayed range the kernel's support is empty in that exponent box, so all later coordinates vanish.

It follows for every nonnegative weight vector in (4) that

    D_h(q)<=0 for every h>=2.                           (7)

Strictness is supplied separately by the already proved finite-support positive-square identity in `../covered_positivity/PROOF.md`: for every nonzero finite nonnegative sequence q with q0=0,

    D_1(q)=-P(q)<0.                                    (8)

In particular, no positive Hamilton lower bound is smuggled into the higher-coordinate certificate. An actual nonempty child has a nonzero profile, so (8) applies. Together (7)--(8) prove componentwise domination of the augmented-core ordered covers, strict at the first coordinate.

## 4. One actual query in every fixed ambient

Let Q be any fixed tournament with distinguished root r, replaced by the source or target core. Every exterior vertex and arc remains fixed. Extend the root masks by assigning both ports the original root's incidence at every exterior vertex z: if `z->r`, use `z->x,z->y`; if `r->z`, use `x->z,y->z`. These are the same actual port bits on both sides.

The augmented core is then a common module. Exact module-run substitution supplies nonnegative exterior coefficients b_h with `b_1=H(Q)>0`. Thus the actual physical response difference is

    F_(Q[T])(a,c)-F_(Q[U])(a,c)
      =sum_h b_h D_h(q)
      <=-H(Q)P(q)<0.                                  (9)

The exterior is arbitrary in size and need not be transitive or non-strong. This is the same physical module lift proved in `../equal_small_ambient/PROOF.md`; zero exterior masks are not substituted for the root-incidence assignment. Profile-preserving target replacements do not affect any count because both ports are uniform on every child.

## 5. Boundary, provenance and scope

At n=13 the sufficient coefficientwise test already has a positive coefficient: at h=2, the weight monomial `w0^2 w6` has coefficient `1930071780`. The verifier retains this exact control. A positive mixed coefficient does not imply that the full cubic is positive anywhere in the nonnegative cone, much less on an actual tournament profile. No counterexample or theorem above twelve follows from this boundary.

The word-and-cut bijection and word-reversal symmetry are elementary identities, and no novelty claim is made for them or for the coefficient-sign method. The new certified application is their stated response-theoretic consequence. The imported all-cover formula, scalar strictness and module lift remain separately identified rather than presented as new here.

This is a common-complete-profile cyclic-quotient theorem on fixed named child sets in a fixed exterior. It does not settle unequal-profile substitutions, arbitrary prime quotients, changes involving exterior names, or unrestricted reconstruction from full named responses. Every cubic and quartic interaction is retained in the inherited literal two-port counts.

Run `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare. The verifier uses only the Python standard library, imports no sibling implementation, and keeps its explicit checks active under optimized execution.
