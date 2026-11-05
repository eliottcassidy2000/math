# A three-module identity and a constructive two-prime-power family

Session 197. **PROVED physical identity and explicit inverse construction.** This extends the positive-unit clock from `../positive-clock/PROOF.md`; it does not claim infinitely many prime values or the full `{7,21}` spectrum.

## 1. The physical substitution

Let Q have outgoing masks

    [30,64,10,66,174,335,21,367,95].

This is the rooted-isomorphic presentation of the other lane's positive clock used in `POSITIVE_CLOCK_REVIEW.md`. Replace vertex0 by a nonempty tournament T, vertex2 by a nonempty tournament U, and vertex8 by a nonempty tournament V. Keep every other vertex a singleton. Arcs between different blocks follow Q, and each block retains its internal tournament. Write the resulting actual graph as Q[T,U,V].

For a tournament A, let B_A(d) be its labelled-bin path-cover polynomial, counting d labelled, possibly empty bins each carrying a directed path. Thus B_A(1)=H(A) and B_Tr_m(d)=d^m.

**Theorem.** For all nonempty T,U,V,

    H(Q[T,U,V]) = 18 B_T(3) B_U(2) B_V(2)
                 + H(T) H(U) H(V).                         (1)

This is a scalar Hamilton-count identity. No unproved formula for the entire output bin polynomial or for nonuniform attachments is implicit in (1).

## 2. Uniform run bounds, before interpolation

First put T=Tr_m, U=Tr_s, V=Tr_t, where m,s,t>=1. Consider one of these modules as the variable transitive clone block while fixing the other two actual blocks.

For a rooted tournament, write O for the vertices dominated by its root and I for those dominating it. If a Hamilton path uses k nonempty root runs, its k-1 intervening outside paths go from O to I. Each contains an O-to-I arc. The outside paths are vertex-disjoint, so selecting one such arc from each gives a matching in the bipartite O-to-I graph. Consequently a cross-matching bound nu gives k<=nu+1.

The relevant uniform covers can be read directly from the nine masks:

- At module0, the only O-to-I arcs are1->6,3->6,4->5,4->7. The clones of2 and8 remain isolated in this cross graph. Its matching number is exactly2, uniformly in s,t. Thus there are at most3 runs of module0.
- At module2, O={1,3}. The only O-to-I arcs are1->6 and3->6. Enlarging modules0 and8 changes none of these arcs. The single vertex6 covers this cross graph, so there are at most2 runs of module2, uniformly in m,t.
- At module8, I={5,7}. The only O-to-I arcs are4->5 and4->7. Enlarging modules0 and2 changes none of these arcs. The single vertex4 covers this cross graph, so there are at most2 runs of module8, uniformly in m,s.

For one transitive module of size m and exactly k distinguished run positions, assigning its m labelled, internally ordered vertices onto the k positions gives k! S(m,k). Hence the Hamilton count is a finite linear combination of j^m with1<=j<=k_max. This is the usual run/Stirling expansion, and the preceding covers prove its bounds for all values of the other two variables, not merely for a sampled grid.

Therefore the three-variable count lies in the tensor-product space spanned by

    i^m j^s k^t,    1<=i<=3, 1<=j<=2, 1<=k<=2.              (2)

For completeness, separate bounded-span statements suffice: interpolate in m using its three positive evaluation points; the resulting three coefficient functions are linear combinations of those evaluations and retain the two-variable bounds in s,t. Interpolate those in s and then t. Each positive-base Vandermonde matrix is nonsingular. Thus the 12 values m=1,2,3 and s,t=1,2 determine the count everywhere on positive integer inputs.

## 3. Exact physical grid and arbitrary-module extension

Independent vertex-subset Hamilton DP on the 12 explicitly constructed graphs gives the following rows, indexed by m, with columns(s,t)=(1,1),(1,2),(2,1),(2,2):

    m=1: 217, 433, 433, 865
    m=2: 649, 1297, 1297, 2593
    m=3: 1945, 3889, 3889, 7777.

The function

    18*3^m*2^s*2^t+1                                      (3)

belongs to (2) and has all12 values. The proved uniform degree bounds and Vandermonde uniqueness establish (3) for every positive m,s,t. This is finite exact interpolation with a prior all-order dimension proof, not an extrapolation of a numerical pattern.

For fixed other blocks, whole-module Hamilton substitution is linear on the complete bin polynomial, as proved by the path-fragment kernel in `UNIVERSAL_CONTEXTS.md`. Transitive inputs provide the monomials x^m for every m>=1, which span all polynomials with zero constant term. Every nonempty tournament's B polynomial has zero constant term. Extend (3) successively in T, U, and V by this linearity. This proves (1), including arbitrary nontransitive input graphs. Empty inputs are not needed or asserted by this theorem.

The checker additionally verifies (1) on all64 triples from K1,Tr2,Tr3,C3 by constructing their actual graphs and using endpoint DP, rather than substituting only their predicted scalar readouts.

## 4. An explicit inverse construction for 2^a 3^b+1

For any integers a,b>=3, choose

    T=Tr_(b-2), U=Tr_(a-2), V=K1.

Equation(1) gives

    H(Q[T,U,V])=2^a*3^b+1,
    |Q[T,U,V]|=a+b+3.                                     (4)

Thus a target N with N-1=2^a3^b and a,b>=3 has an explicit physical inverse, with no search or assumption that N is prime. The construction is the block graph in Section1 with precisely these sizes. Factoring N-1 by repeated division by2 and3 decides whether this rule applies.

Q is strong, and substitution of nonempty blocks into a strong quotient preserves strong connectivity. All graphs in (4) are strong.

The source has the anti-involution

    tau=[0,7,8,5,6,3,4,1,2].

It fixes0 and exchanges2 with8. If a is odd, instead choose equal U,V of order(a-1)/2. Reversing each transitive chain and applying tau then reverses every arc; this gives a physically self-converse graph with the same count and order. There is no restriction on b beyond b>=3. For even a, unequal sizes in this construction do not support that specified lifted anti-involution; we make no assertion that no other self-converse realization exists.

### A new explicit 7p member

Take(a,b)=(4,7), so(m,s,t)=(5,2,1). The resulting14-vertex graph has

    H=2^4*3^7+1=34993=7*4999.

The receipt gives every outgoing adjacency mask and independently checks the Hamilton count. Trial division by all primes at most sqrt(4999)<71 certifies4999 prime; all remainders are listed. This member is outside the fixed-amplitude72 family, whose only7p member was217. There is no contradiction: the extra binary module changes the amplitude from72 to144.

The necessary divisibility condition in the two-parameter family is

    7 | (2^a3^b+1)  iff  2a+b=3 modulo6,

because2=3^2 mod7 and3 has order6. This is only a divisibility filter, not a prime-quotient guarantee. For example, the tempting self-converse choice a=b=7 gives6^7+1=7^2*29*197, not7 times a prime. No infinitude claim is made.

## 5. Scope, source and reproduction

The base positive clock is attributed to the orbit lane's frozen proof, SHA256 `b725136569ee84ddd46c5ef4f0b45b40e216eec9c314ffb604609e16a1e96f6f`; its independent rooted isomorphism and physical checks are in `POSITIVE_CLOCK_REVIEW.md`. The general single-module linearity and run framework are inherited, not claimed as new here. This packet adds the uniform three-variable cross covers, the exact trilinear physical identity, and the resulting varying-amplitude inverse rule.

Run `python3 verify_three_module_clock.py` and `python3 -O verify_three_module_clock.py`. The standard-library verifier has explicit exception gates. It checks the cross covers on all12 interpolation graphs, the exact physical grid, all64 nontransitive-capable input triples, the lifted anti-involutions, the actual14-vertex7p witness, and the complete short primality certificate. All-order scope comes from the proofs in Sections2-3, not from testing finitely many exponents.
