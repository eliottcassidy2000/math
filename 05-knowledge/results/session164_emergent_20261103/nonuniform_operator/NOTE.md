# Nonuniform head operators, sharp strong witnesses, and observer limits

**Status: PROVED**, with the separate finite controls described below. The new mechanism is an exact three-shift operator on path-cover polynomials. It gives stronger literal witnesses for the Session 163 interface bound, and an infinite-rank obstruction to compressing its repetitions into a fixed linear observer. It does **not** prove source-closed ternary branches or the tournament Hamilton-path spectrum. No literature-priority claim is made.

## Inheritance and scope

The closest inherited mechanisms are the run-divisible-module and capacity-two matching theorem in [S163 nonuniform/NOTE.md](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session163_ternary_triples_20261103/nonuniform/NOTE.md), the quadratic single-interface and hidden-observer examples in [MARKED_ESCAPE.md](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session163_ternary_triples_20261103/nonuniform/MARKED_ESCAPE.md), and the two physical suffix matrices in [strong_returns/NOTE.md](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session163_ternary_triples_20261103/strong_returns/NOTE.md).

The hostile examples remain important: a protected triangle in an old tail is not unlocked merely by manipulating a new head; gcd one across a family does not imply a single count is a unit; and a count identity is not a physical source-state return. The less-used sidecar here is the *whole* path-cover polynomial, rather than just H or the three-coordinate marked observer.

## 1. A general one-apex tail operator

For a tournament T, let F_T(r) count covers of its vertices by r labeled directed paths, allowing empty paths. Thus F_T(1)=H(T), F_empty(r)=1 including r=0, and F_T(0)=0 when T is nonempty. If C_k(T) counts ordered covers by k nonempty paths, then

    F_T(r) = sum_k C_k(T) binom(r,k).

This is a polynomial. The ordinal-sum identity is

    F_(A ⊕ B)(r) = F_A(r) F_B(r),

because each path has an A portion followed by a B portion. For a directed triangle C and a transitive tournament T_n,

    f(r) := F_C(r) = r^3 + 2r,       F_(T_n)(r) = r^n.

Fix a head tournament W and D contained in V(W). Form K_(W,D)(Z) on W, an arbitrary tail Z, and a new apex z by prescribing

    W → Z → z,
    z → D,        W\D → z.

Edges within W and Z are retained. In particular Z remains an induced tournament.

For r≥1 define these head-only integers:

* E_W(r): covers of W by r labeled path slots with a distinguished left slot nonempty and ending in D.
* B_W(r): covers of W by r+1 labeled slots with a distinguished right slot nonempty and starting in D.
* C_W(r): covers of W by r+1 labeled slots with distinct left and right slots, the left nonempty and ending in D, and the right nonempty and starting in D.

All unspecified slots may be empty. Then the exact operator is

    F_K(Z)(r) = r [ B_W(r) F_Z(r+1)
                  + (F_W(r) − C_W(r)) F_Z(r)
                  − E_W(r) F_Z(r−1) ].                         (1)

**Proof.** Choose the slot containing z, in r ways, and split that path immediately before and after z into left and right slots. Every resulting slot consists of a W path followed by a Z path.

If the right W path is empty, the right Z path must also be empty, since z loses to Z. This initially gives F_W(r)F_Z(r). The only invalid incoming edges occur when the left W path ends in D and its Z path is empty. Subtract E_W(r)F_Z(r−1).

If the right W path is nonempty, it must start in D. All r+1 Z slots are initially available, giving B_W(r)F_Z(r+1). Again the invalid incoming edges have a left W path ending in D and an empty left Z path; subtract C_W(r)F_Z(r). These cases are disjoint and exhaustive. Rejoin the two distinguished slots through z to invert the construction. ∎

This identity remains correct for the empty tail. The F_Z(r−1) boundary term must not be discarded when r=1 and Z is empty.

## 2. A concrete higher-response head

Take W=C⊕C, with D consisting of one vertex in each triangle. Write G(Z)=K_(W,D)(Z). A specified vertex starts or ends the distinguished path in a triangle in

    q(r) := r^2 − r + 1

ways: the distinguished path has length 1, 2, or 3, contributing respectively (r−1)^2, r−1, and 1. An empty distinguished slot contributes f(r−1). Consequently

    B(r) = q(r+1) [f(r+1)+f(r)],
    E(r) = q(r) [f(r)+f(r−1)],
    C(r) = q(r)^2 + q(r+1)^2.

For the joint condition C, the left endpoint and right start cannot be the same selected vertex. If the left endpoint is in the first triangle and the right start in the second, each opposite distinguished slot must be empty, giving q(r)^2. The reverse assignment gives q(r+1)^2.

Thus F_(G(Z)) = L F_Z, with

    (L F)(r) = r [ B(r) F(r+1) + A(r) F(r) − E(r) F(r−1) ],

    B(r) = 2r^5 + 5r^4 + 12r^3 + 13r^2 + 10r + 3,
    A(r) = r^6 + 2r^4 − 2r^2 − 2,
    E(r) = 2r^5 − 5r^4 + 12r^3 − 13r^2 + 10r − 3.       (2)

In particular,

    H(G(Z)) = 45 F_Z(2) − H(Z)                 if Z is nonempty,
    H(G(empty)) = 45 − 1 − 3 = 41.                         (3)

The stronger full-polynomial identity includes

    F_(G(empty))(r) = r^7 + 12r^5 + 24r^3 + 4r,
    F_(G(empty))(2) = 712,
    F_(G(empty))(3) = 5763.

G(Z) is strong for every tail Z, including the empty tail. The apex reaches the first triangle; that triangle reaches the second and every tail vertex. Conversely each old vertex reaches z, directly or via a triangle vertex or tail vertex that beats z. The two triangles themselves are strong, so the selected incidence is enough to connect every vertex.

This construction adds seven vertices, retains Z literally, and can repeat indefinitely. With Z_0 empty and Z_(k+1)=G(Z_k),

    H(Z_1)=41,  H(Z_2)=31999,  H(Z_3)=168693881.

The first two counts and F_(Z_2)(2)=3749464 are independently counted on the literal graphs. The displayed third count is a theorem-backed operator evaluation, not a direct 21-vertex enumeration. For k≥1, H(Z_k) is congruent to (−1)^k modulo 3.

### Positive polynomial calculus

Although (2) has negative terms, L preserves the cone of polynomials with nonnegative integer monomial coefficients. Put

    U=2r^5+12r^3+10r,      V=5r^4+13r^2+3,

so B=U+V and E=U−V. For every n≥0,

    L(r^n) = r^(n+7)+12r^(n+5)+24r^(n+3)+4r^(n+1)
            +2r U sum_(j odd) binom(n,j) r^(n−j)
            +2r V sum_(j even, j≥2) binom(n,j) r^(n−j).    (4)

All coefficients are nonnegative integers. This follows by expanding (r+1)^n±(r−1)^n and absorbing the j=0 terms. It also shows that L flips polynomial parity and preserves monicity while increasing degree by seven. These are algebraic controls on an actual tournament construction, not independent realizability criteria for arbitrary positive polynomials.

## 3. Strong sharpness with no wrapper vertices

S163 proved that t run-3-divisible blocks and s exceptional vertices force

    v_3(H) ≥ max(t−2s,0),

and gave strong witnesses using an additional wrapper. The following removes the wrapper in the range t≥2s, s≥1, while keeping precisely the specified triangle blocks and exceptional vertices.

Let d=t−2s. Let J=G(empty), the seven-vertex 41-path gadget. Set

    Z = J ⊕ ... ⊕ J ⊕ C ⊕ ... ⊕ C,
        s−1 copies of J     d copies of C,

and output G(Z). There are exactly t disjoint triangle blocks, s designated apex exceptions, and no other vertices: its order is 3t+s. Between triangle blocks and all nonexceptional vertices, orientations are uniform. Every apex from a J factor is nonuniform only on its own two triangles; the final apex is nonuniform only on the newly prepended two triangles. The interface graph is s disjoint two-block stars.

If Z is nonempty, equations (2)–(3) and ordinal multiplicativity give

    H(G(Z)) = 45·712^(s−1)·12^d − 41^(s−1)·3^d
            = 3^d [45·712^(s−1)·4^d − 41^(s−1)].           (5)

The bracket is nonzero modulo 3, so v_3(H)=d exactly. The sole empty-tail case s=1,d=0 has H=41 and the same claimed valuation. All outputs are strong by the preceding argument.

Finite literal controls are:

    (t,s,n,H,v_3) = (2,1,7,41,0), (3,1,10,537,1),
                    (4,1,13,6471,2), (4,2,14,31999,0).

The vertex count is minimal *for retaining t disjoint triangles and s distinct designated exceptions*, since those specified sets already require 3t+s vertices. This is not a global minimum-order theorem for tournaments with the same Hamilton-path count. The case s=0 is not covered: ordinally arranging several untouched triangle modules need not be strong.

## 4. Why the old observer is not enough, and no fixed linear replacement suffices

Let v(X,p)=(h,x,y) be S163's observer defined by

    h=H(X),  H(X^a)=h+2x,  H(X^aa)=h+2(x+y).

Every transitive T_n, n≥1, marked at its sink has v=(1,0,1): the first query stays transitive, and the second creates a terminal directed triangle. But

    H(G(T_n)) = 45·2^n − 1.

Therefore the head operation G cannot factor through this three-coordinate observer on all marked sources. The missing F(2) is essential; it cannot be silently reconstructed from (h,x,y).

There is a stronger **infinite linear-rank obstruction**. Define

    a_k(n)=H(G^k(T_n)),       n≥1, k≥0.

Repeated use of (1)–(2) shows that a_k(n) is a linear combination of

    1^n, 2^n, ..., (k+1)^n,

with coefficient

    c_k = product_(j=1)^k j B(j) > 0

on (k+1)^n; c_0=1. Only the upper shift can reach this largest evaluation, which proves the formula by induction. For example,

    a_0(n)=1,
    a_1(n)=45·2^n−1,
    a_2(n)=28350·3^n+7695·2^n−4049.

The functions n↦j^n are linearly independent for distinct positive j, by a Vandermonde determinant. Thus {a_k} has infinite linear rank. If a fixed d-dimensional linear observer v(T_n), a fixed transfer M, and a fixed output column e represented all these repetitions, every a_k(n)=v(T_n)M^k e would lie in the span of the same d coordinate functions of n. This is impossible.

This excludes fixed finite-dimensional **linear** representations for all transitive tails and all repetitions of this specific G; affine representations are also excluded after adjoining a constant coordinate. It does not exclude nonlinear encodings, bounded-depth states, or useful finite states on a specially restricted family. The full polynomial is an exact growing state, not a solved finite ternary return bank. Moreover a full F polynomial cannot literally triple into another tournament's F polynomial: each is monic, and different orders have different degrees.

## 5. Hidden observer content is detected in one suffix

For the inherited physical suffixes q2 and b, let

    P=[[5,3,4],[8,7,8],[6,6,7]],
    B=[[1,2,2],[4,3,4],[4,2,3]],

with row-vector convention v(Xq2)=v(X)P and v(Xb)=v(X)B. Their three observed counts are

    h,       5h+8x+6y,       h+4x+4y.

The observation matrix with these columns has determinant 8. Since H is odd for every nonempty tournament, this yields the exact integer identity

    gcd(H(X),H(Xq2),H(Xb)) = gcd(h,x,y).                    (6)

Indeed the right side divides all three counts, while the adjugate shows their gcd divides 8 times each coordinate; their gcd is odd. Equivalently, for every odd m, all three counts vanish modulo m iff the whole observer vanishes modulo m.

Consequently any hidden failure of full-state divisibility is exposed by at least one of the two one-atom suffixes q2 or b. When X is strong, both suffixes preserve strongness. This is a finite *nonzero-residue* detector, not a guarantee that one of the three counts is a unit modulo an arbitrary composite m.

The full-state content has a second exact rule. Because det(B)=1 and B has an integer inverse,

    gcd(vB)=gcd(v).

Because det(P)=5, if g=gcd(v), then gcd(vP)/g is either 1 or 5. Writing w=v/g, the ratio is 5 exactly when

    w mod 5 = c(1,3,1),      c∈{1,2,3,4}.

This follows from the adjugate and the one-dimensional left kernel of P modulo 5. No suffix word in b,q2 can remove a prime already dividing the entire observer.

For the inherited R4 example, one reversed interface incidence gives v=(1965,1103,1564). Appending q2 yields the literal strong count

    5·1965+8·1103+6·1564=28033,

which is coprime to 15, not merely nonzero modulo 3. Thus the immediate three-edge family's mod-15 obstruction from S163 is overcome after this permitted suffix. This retains the distinction between the restricted immediate family and the enlarged physical family.

## 6. Multi-apex response support and the critical interface layer

Here is a general extension of the quadratic interface criterion. Fix a tournament A and s new exceptional vertices, including all edges between those vertices. For each variable old–exception incidence, use a Boolean variable indicating reversal from a fixed baseline. The Hamilton-path count is an integer multilinear polynomial in these variables.

A monomial can have a nonzero coefficient only if its variable edges form an undirected linear forest: maximum degree at most two and no cycle. In particular its degree at each exceptional vertex is at most two and its total degree is at most 2s. To prove this, sum path indicators over all vertex permutations. Each indicator is a product of x_e or 1−x_e over the variable edges actually used by that path. Every monomial in its expansion is a subset of a Hamilton path's edges. Summation can cancel coefficients but cannot create forbidden supports.

Thus values at assignments with at most 2s reversed incidences determine the whole response by integer Boolean Möbius inversion. If the baseline is zero modulo m, some assignment is nonzero modulo m iff some assignment with at most 2s reversals is nonzero. This is **not** a sparse-unit theorem for composite m. A fixed old physical observer port gives the same degree/support bound for all its queried counts and hence its integral observer coordinates.

Now suppose A contains pairwise disjoint blocks B_i with m_i dividing every ordered nonempty-path cover count C_k(B_i). Assume baseline orientations from every exceptional vertex to each block are uniform, and all remaining external vertices see each block uniformly. For a polynomial coefficient indexed by a reversal set S, let hit(S) be the blocks incident to S. Then

    product_(i not in hit(S)) m_i   divides that coefficient.      (7)

Every Boolean evaluation used in this coefficient's Möbius difference leaves those blocks genuine uniform modules. Their ordered runs factor independently, proving the product divisibility; take the integer difference. Combining (7) with the degree-two condition recovers the S163 capacity-two matching bound, now coefficient by coefficient.

In the critical equal-modulus case t=2s, modulo m the entire count polynomial is therefore homogeneous of degree 2s: every surviving monomial hits every block exactly once, has exactly two incidences at every exception, and corresponds to s disjoint two-block wedges. All other coefficients vanish modulo m. Hence fewer than 2s reversed incidences cannot escape modulo m. The critical layer is an exact response polynomial, although its coefficients may cancel and it need not have a nonzero Boolean value.

For two ordered triangles and one exception, writing k,l for the numbers of reversed incidences into the first and second triangles gives the complete exact polynomial

    H = 9+36k+6l−6[binom(k,2)+binom(l,2)]−10kl.

Modulo 3 only the critical cross term remains: H≡2kl. The choice k=l=1 is the literal strong 41-path gadget used above.

## Controls and remaining boundary

Run `python3 check.py` and `python3 -O check.py`. Failure checks remain active under optimization. The independent endpoint/subset dynamic program and labeled-cover subset convolution verify:

* 324 general one-apex operator identities for all interfaces on four explicit heads, four tails, and r=1,2,3;
* the specialized coefficient polynomials through r=5;
* G(empty), G²(empty), and four sharp strong witnesses through fourteen vertices;
* all 64 interfaces of the two-triangle head and all 64 two-apex interfaces on one triangle, including forbidden monomial supports;
* 285 literal marked-source transfer/content cases, covering every labeled tournament of orders 1 through 4 and every physical port;
* the R4 one-edge state and the literal strong q2 extension with count 28033.

`check.json` and `check_optimized.json` record matching results. The third 21-vertex iterate is clearly identified as formula-only. The general proofs are the slot bijection, coefficient-support argument, module factorization, elementary matrix algebra, and Vandermonde independence; the finite checks do not replace them.

A tail module sees every newly prepended head vertex and apex uniformly, so its existing run-divisibility survives G. This is the relevant stopping boundary: the construction breaks the two *new* triangle blocks using their designated apex, but does not unlock an arbitrary protected old block. The open return problem still requires a physical operation that both reaches the intended ternary count/state and retains a reusable admissible source class.
