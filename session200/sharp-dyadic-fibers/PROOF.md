# Sharp dyadic source fibers and a height-sensitive exact-value bridge

Status: PROVED below from the named physical source identities; finite computations are supplementary. This does not prove universal odd Hamilton-path coverage, an exact-count no-singleton theorem, or multiplicity of unlabelled graph isomorphism classes.

## 1. Source conventions and provenance

An actual nontrivial rooted start-one tournament source has response

    s=(D,E,U,d)=(H(V), e_c(V), R_c(V), H(V-c)).

The complete physical grammar starts at t=(1,0,2,1), with chronological moves

    O(s)=(D,0,2U,d),
    Z(s)=(D-E+U,d,2U-d,d),
    J(s)=(D-E+U,D-E,D+U,D).

O adds a universal sink. Z adds a vertex beating only the current root, retaining the old root. J adds the same vertex and marks the new vertex. The inherited source theorem gives U>=2d, D-E>=d, and U>=D. Induction on this grammar also gives D,d odd and U=E modulo 2.

The two readouts used here are

    F(s)=D-E+U: append J, producing a strong C1 tournament;
    G(s)=D+2U: append OJ, producing a balanced strong C1 tournament.

Both are odd on every physical source. Inherited inputs are Session198/bilinear-descent, Session198/balanced-source-normal-form, Session198/source-digit-blocks, and Session199/physical-adjoint. The first three may equivalently be read through the reviewed frozen dependencies in Session199/source-cut-rearrangement. The substantive new results are Theorems 1-3 below, especially the all-precision sharpness of Theorem 1.

The motivating attachment is *Z is Diophantine in Q*, dated November 5, 2026, pages 75 and 81-83. Its Boolean coefficient-parity argument and its shared-witness height/separation argument are used as ideas, not as assumed theorems. All facts needed here are proved independently. The attachment's arithmetic headline is not a dependency.

## 2. The exact bit polynomial

Fix any actual prefix source s and append a retained-root block of k letters O or Z. Put x_i=1 at Z and 0 at O, with positions i=0,...,k-1. Direct iteration of the source laws gives

    D'=D+U sum_i 2^i x_i-E x_0-d sum_(i<j) w_ij x_i x_j,
    E'=d x_(k-1), U'=2^k U-d sum_i 2^(k-1-i) x_i, d'=d,
    w_ij=2^(j-i-1)+[j=i+1].

The k=0 block is the identity and is treated separately. This is an integer quadratic polynomial on the entire Boolean cube of physical source words. Both F(s.word) and G(s.word) have the same mixed coefficient -d w_ij. Every w_ij is even. Since every value is odd, all nonconstant Boolean-polynomial coefficients are even. Thus

    P(x)=(H(x)-H(0))/2

has integer coefficients and degree at most two, where H is either F or G.
Because d is odd, the coefficient of x_i x_j in P is odd exactly when

    j-i is 1 or 2.

In other words the odd quadratic support is precisely the square of the k-vertex path. This graph is on the word positions, not the tournament vertices.

## 3. A Boolean degree and parity lemma

For an integer polynomial P of Boolean degree at most d, the bit of index j of its integer value has Boolean degree at most d*2^j. This remains true for negative coefficients and negative values, with bits interpreted by reduction modulo powers of two.

Proof. For integer z the bit of index j equals binomial(z,2^j) modulo 2. Expand the coefficient of T^(2^j) in (1+T)^P. For each nonempty squarefree monomial x_A with coefficient c_A, use the identity on Boolean inputs

    (1+T)^(c_A x_A)=1+x_A((1+T)^c_A-1).

The constant term of P contributes only a scalar series. At most 2^j nonconstant factors can contribute to the requested coefficient. Each uses at most d variables. Truncation makes the calculation finite even if c_A is negative. Multilinearization cannot increase degree. QED.

For P quadratic and r>=1, the indicator of a specified residue modulo 2^(r-1) is the product of r-1 bit indicators. Its degree is at most

    L_r=2 sum_(j=0)^(r-2) 2^j = 2^r-2.

For r=1 the product is empty and L_1=0. The sum over {0,1}^k of a Boolean polynomial of degree less than k is zero modulo 2, since every monomial omits a variable.

## 4. Sharp dyadic fiber theorem

**Theorem 1.** Fix any actual nontrivial prefix source, either closure F or G, and an integer r>=1.

1. At block length k=2^r-2, every odd residue modulo 2^r occurs an odd number of times among the 2^k indexed O/Z words.
2. At every length k>2^r-2, every odd residue occurs a positive even number of times.
3. Thus every larger-depth residue fiber has at least two distinct labelled physical tournament constructions. No conclusion about equal exact H values or nonisomorphic unlabelled tournaments follows.

Proof of the even assertion. H=n modulo 2^r is equivalent to P=(n-H(0))/2 modulo 2^(r-1). Apply Section 3: when k>L_r, its indicator has degree less than k, so its fiber cardinality is even.

Proof at the critical length. Write k=2t with t=2^(r-1)-1. Only the product of the top-degree pieces of all r-1 bit polynomials can contribute to the coefficient of x_0...x_(k-1) in the residue indicator; omitting any bit factor decreases the maximal degree. This coefficient is independent of the chosen residue and of all linear terms of P.

The top-degree part of bit j is obtained by selecting 2^j disjoint quadratic monomials whose coefficients are odd. Each selected term must contribute its linear T coefficient. Hence it is the matching polynomial for matchings with 2^j edges in the odd quadratic support graph.

In the product over j=0,...,r-2, each contribution to the full monomial is a perfect matching of that graph, with its t edges partitioned into distinguished groups of sizes 1,2,4,...,2^(r-2). The multiplicity

    t! / product_(j=0)^(r-2) (2^j)!

is odd: the 2-adic valuation of its numerator is t-(r-1), equal to the sum of the denominator valuations. Consequently the full coefficient equals the parity of the number of perfect matchings of the square of a path on 2t vertices.

Let a_t denote that matching count. Then a_0=a_1=1 and a_t=a_(t-1)+a_(t-2). The first vertex is matched either to the second, leaving a_(t-1), or to the third. In the latter case the second must be matched to the fourth, leaving a_(t-2). Thus a_t is Fibonacci F_(t+1). Modulo two Fibonacci repeats 0,1,1, so F_m is odd exactly when 3 does not divide m. Here t+1=2^(r-1), never divisible by 3. The full coefficient is therefore one.

Summing the indicator over the critical cube extracts exactly that full coefficient. Each odd residue fiber is odd, hence nonempty. For r=1 the empty cube has one word, giving the same result.

For k>L_r, fix any first k-L_r letters. They produce another actual source. The critical-length theorem on its last L_r letters supplies every odd residue. This proves nonemptiness, completing the theorem. Distinct indexed words change at least one edge between the retained root and a named newborn; all final closure edges are identical. Therefore the labelled graph constructions are distinct. QED.

This gives fixed-depth physical modular surjectivity without needing the earlier general modular compiler. It also proves the parity threshold is sharp at every r, not only at experimentally tested precisions.

**Example.** From the edge source with F closure, k=2 gives exact counts 9,11,9,9 for OO,OZ,ZO,ZZ. Modulo four the odd fibers have sizes 3 and 1. At k=3 both residue fibers have size 4. At r=3,k=6 all four odd-residue fibers are odd; at every k>=7 they are positive and even.

## 5. Why the residue theorem cannot supply exact small counts

Every retained-root step obeys U_new>=2U-d while d is fixed. Therefore

    U_k >= 2^k(U_0-d)+d >= d(2^k+1).

Since D_k-E_k>=d,

    F_k>=d(2^k+2),   G_k>=2d(2^k+1)+1.

At the critical modular depth k=2^r-2, these lower bounds are exponential in the modulus m=2^r, not in its bit length. Thus the fixed-depth witnesses for a fixed target residue necessarily escape to infinity in ordinary Hamilton-path count as r increases.

This identifies the first failed implication in a tempting Carmichael analogy. A residue fiber can have many physical labelled objects while an exact scalar fiber is empty. The paper's growing address dimension is useful only together with a bound independent of that dimension. Our source extension has no such ordinary-height bound.

## 6. Exactness from one common pair and a growing height bound

For nontrivial physical source responses define the inherited physical C1 pairing

    K(v,w)=D_v d_w+d_v D_w+E_v E_w+U_v U_w
            -E_v d_w-d_v E_w-U_v d_w-d_v U_w+d_v d_w.

Its physical adjoint turns any pair into an actual C1 source word: reverse one word and exchange O with Z, retaining J. It also proves

    0<K(v,w)<=U_v U_w,
    K(v,w)>=U_v U_w/2+d_v d_w.

Every K value is odd. The exact C1 image is K(S,S), and the balanced image other than 3 is K(S,B), where B is the source orbit from the triangle b=(3,1,3,1). A bound on total J cost is preserved by the adjoint.

**Theorem 2 (subcritical modular height criterion).** Fix an odd n>=3 and m>=1, and put q=lcm(m,2). Suppose q>n. If the same physical pair v,w satisfies

    K(v,w)=n modulo m,
    U_v U_w<n+q,

then K(v,w)=n. The same assertion holds in the balanced or bounded-J-cost classes without enlarging their class.

Proof. Oddness strengthens the congruence to modulo q. The positive integer K lies strictly between 0 and n+q. Since q>n, the only positive integer in this interval congruent to n modulo q is n. Physical compilation supplies the graph and preserves any stated budget. QED.

The strict inequality is equivalent to the integral upper bound U_vU_w<=n+q-1. It is enough to have each height at most floor(sqrt(n+q-1)); a nonuniform product bound is stronger. Thus the bound may grow on the square-root-of-modulus scale. Uniform bounded height, as in the earlier profinite repair, is not necessary.

**Exact reformulation of the remaining lemma.** For the physical balanced class of total J cost at most three, the following are equivalent for any fixed odd n>=5:

- n occurs in that class;
- for some m with lcm(m,2)>n there is a physical pair in that class satisfying the two displayed conditions;
- for every sufficiently large m there is such a pair.

The forward implication fixes one exact physical pair and increases m until its product lies below n+q. The reverse implication is Theorem 2. A proof of this condition for every odd n>27781 would close the existing four-block coverage lemma. Neither Theorem 1 nor the generic modular compiler supplies its height inequality.

## 7. The threshold is sharp in a primitive zero-J class

**Theorem 3.** Let h=2^a+1 with a>=1. There are primitive physical sources

    v=w=(h,1,h,1),

where v=t.Z^a and w=b.Z^(a-1). Their total J cost is zero and K(v,w)=h^2. The corresponding balanced C1 graph is literally

    t.Z^a.O^a.J.

For every fixed positive odd n, choose a large enough that m=h^2-n>n. Then m is even,

    K(v,w)=n modulo m,
    U_vU_w=n+m,
    U_v=U_w=sqrt(n+m),

but K(v,w)!=n.

Proof. Iterating Z on the balanced response (3,1,3,1) gives (2^a+1,1,2^a+1,1). The K identity on balanced sources is K(v,w)=D_vD_w. Alternatively, direct source iteration on t.Z^a.O^a.J gives H=h^2 and exactly one Hamilton cycle. Its last O before J makes it balanced; for a=1 this graph is also obtained from the triangle as required. Since d=1, both profiles have content one. The remaining claims follow by substitution. QED.

The smallest useful hostile for n=7 is a=2: h=5, m=18, H=25. Its pair has product 25=n+m, and cannot be accepted by the strict criterion. The independent finite C1 inverse also confirms that 7 is absent using the inherited complete height cutoff floor(2 sqrt(6))=4. No new unrestricted nonexistence theorem is asserted.

Thus replacing '<' by '<=' fails even on primitive actual sources, with no J budget spent, and with perfectly balanced heights. An asymptotic condition U_v,U_w=(1+o(1))*sqrt(m) likewise does not suffice. What would suffice, for example, is a fixed c<1 with U_v,U_w<=c sqrt(m) at arbitrarily large even moduli, but existence of such witnesses is an open constructive obligation.

## 8. Reproduction and interpretation

Run `python3 verify.py` and `python3 -O verify.py`. The outputs must agree. The standalone standard-library verifier checks both closures and several actual prefixes at all critical depths r<=4, checks larger-depth evenness/nonemptiness, checks the sharpness family by direct graph DP, validates small bit polynomials against graph path counts, and confirms the finite C1 hostile at 7 with the inherited exact cutoff. All checks use explicit exceptions and survive optimized Python.

Connection contract: the source mechanism is Boolean parity plus same-witness arithmetic separation in the attachment; the target is the indexed physical source language and its positive C1 pairing; the maps are the exact bit polynomial and the physical adjoint; the preserved predicates are residue-fiber parity and exact readout; the lost data are exact scalar equality under congruence and unlabelled graph identity; the required sidecars are a literal word recipe and an ordinary product-height bound. The proofs above transfer these mechanisms, not the attachment's elliptic or Diophantine claims.
