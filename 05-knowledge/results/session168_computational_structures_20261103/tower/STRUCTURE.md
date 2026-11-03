# Cayley–Dickson transports: exact generator and nilpotency growth

**Status: PROVED conditional only on the explicit sign-kernel theorem stated below; FINITE-EXACT checks through dimension 128.**

This note continues S167's left-multiplication transport model. It concerns permutations of signed basis units, not arbitrary algebra elements, physical tournament classes, or the full left-and-right multiplication group. The group-order theorem is established in the literature: Jenya Kirshtein, *Multiplication groups and inner mapping groups of Cayley–Dickson loops*, [arXiv:1207.4230](https://arxiv.org/abs/1207.4230), Sections 4–5. Here the dimension-sensitive generator count, lower-central filtration, and exact capacity of a fixed base graph are proved from an explicit kernel description. No literature-priority claim is made.

## 1. Setup and inherited kernel

Let A_n have real dimension 2^n, built by the convention

    (a,b)(c,d) = (ac-conj(d)b, da+b conj(c)).

Its signed basis units U_n={±e_x : x in F_2^n} have 2^(n+1) states. Write e_a e_x=(-1)^f_n(a,x)e_(a+x), and let L_a be the corresponding left permutation. Put G_n=<L_a : a in F_2^n>. Each nonzero L_a has square Z, the global sign reversal. In particular, every basis-unit transport remains an invertible orthogonal map in every doubling dimension, although multiplication by general algebra elements does not.

For n>=3, put d=n-1 and write x=(y,t), y in F_2^d. The explicit sign-kernel theorem is

    ker(G_n -> F_2^n) = {D_(h+ct) : sum_y h(y)=0, c in F_2},       (K)
    D_f(s e_x)=(-1)^f(x) s e_x.

The quotient records index translation. Thus |G_n|=2^(2^(n-1)+n). The companion Cayley–Dickson kernel note, `../cayley_dickson/NOTE.md`, Section 2, proves (K) directly; the checker here independently verifies it in n=3,...,7. The low-dimensional exceptions are G_1=C_4 and G_2=Q_8. For n>=1 the global negative transport is already generated.

## 2. A finite-difference description

Regard functions on F_2^d as the group algebra R=F_2[F_2^d], using point indicators delta_y as its basis and convolution as multiplication (not pointwise multiplication). Its augmentation ideal I consists exactly of even-parity functions. A translation acts by shifting y; write

    Delta_a h(y)=h(y+a)+h(y).

For k>=1, I^k is spanned by products of k differences, equivalently indicators of affine k-flats when the k directions are linearly independent. Choose a coordinate basis v_1,...,v_d and put X_i=delta_0+delta_(v_i). Then

    R = F_2[X_1,...,X_d]/(X_1^2,...,X_d^2),
    dim I^k = sum_(j=k)^d binom(d,j),     I^(d+1)=0.

In particular

    I^2={h : sum h=0 and sum_y y h(y)=0},
    dim I^2=2^d-d-1.

The moment condition is a vector equation over F_2. One proof is that every affine parallelogram has both sums zero; modulo their span delta_y equals delta_0 plus the sum of the coordinate edges specified by y, leaving exactly d+1 independent parity/moment coordinates.

For diagonal transports, conjugation and commutation are finite differences:

    L_(a,A) D_(h+ct) L_(a,A)^(-1) = D_(h(y+a)+c(t+A)),
    [L_(a,A),D_(h+ct)] = D_(Delta_a h+cA),

up to the harmless choice of inverse convention for the commutator. Constants lie in I^d and hence in I^2 because d>=2.

## 3. Exact structural theorem

For every n>=3:

1. The derived and Frattini subgroups coincide and consist of repeated functions from I^2:

       [G_n,G_n]=Phi(G_n)={D_h : h in I^2}.

2. The abelianization is (C_2)^(2n). Exactly 2n generators are necessary and sufficient.
3. For 2<=k<=n-1 the lower central series is

       gamma_k(G_n)={D_h : h in I^k}.

   Consequently the nilpotency class is exactly n-1. The center is {1,Z} and the exponent is four.

### Proof of the initial commutator bound

Define q(u)=0 for u=0 and 1 otherwise, and c(u,v)=q(u)+q(v)+q(u+v). The multiplication recurrence is

    f_n((u,A),(y,t))=f_(n-1)(u,y)+t c(u,y)+A q(y)+At.

The sign difference between L_a L_b and L_b L_a is

    p_n(a,b;x)=f_n(b,x)+f_n(a,x+b)+f_n(a,x)+f_n(b,x+a).

Writing a=(u,A), b=(v,B), direct cancellation gives

    p_n = p_(n-1)(u,v;y)
          +(A+B) Delta_u Delta_v q(y)+B q(u)+A q(v).              (1)

Thus p_n is independent of t. The double difference of q=1+delta_0 is an affine parallelogram indicator or zero. The remaining scalar terms are constants. For n=3 the quaternion term p_2 is constant, so all terms belong to I^2. For higher n, p_(n-1) is repeated across the highest coordinate of y, and its parity on the lower half is even. Indeed, formula (1) expresses that lower half as p_(n-2) plus even-parity parallelograms and constants; the full-domain parity of any p_m is zero, because its four sign rows occur in two shifted pairs. Repetition therefore has zero parity and all first moments. This proves p_n in I^2 by induction. An actual group commutator is a translate of p_n, which remains in I^2.

All commutators of left generators consequently lie in D_(I^2), which is normal by translation invariance. Conversely, commutators of D_h, h in I, with translations span I^2. This gives [G_n,G_n]=D_(I^2).

If g=D_(h+ct)L_a, then g^2 is the diagonal difference of h+ct under a, together with L_a^2 (either identity or Z). All of these lie in I^2. Hence G_n^2 <= [G_n,G_n], and the Frattini subgroup equals the derived subgroup. Comparing dimensions gives

    log_2 |G_n/[G_n,G_n]|
      = (2^d+n)-(2^d-d-1)=2n.

This quotient is elementary abelian, proving the lower bound of 2n generators. For a constructive upper bound choose the n coordinate transports L_(2^i), the d sign flips D_(delta_0+delta_(2^j)) on y, and D_t. Conjugating the d edge flips by the coordinate transports gives all edges of the y-cube, whose span is I. These generate the full diagonal kernel, and the n index translations generate the quotient, so these 2n operators generate G_n. Some chosen generators are words in left-unit transports, rather than single L_a.

Further commutation differentiates I^k once, so gamma_(k+1)=D_(I^(k+1)). The product X_1...X_d is the nonzero constant function, proving the class is d. Every square is diagonal and hence has order at most two, while L_a has order four. To see the center, the stabilizer of e_0 contains D_t and every even h with h(0)=0; its common fixed signed units are exactly ±e_0. A permutation commuting with the transitive action is determined by its value at e_0, so the centralizer of G_n on U_n is precisely {1,Z}. Its center is therefore the same pair. QED.

Examples (algebra dimension, signed states, group order, minimum generators, nilpotency class):

- 8: 16 states, 128 transports, 6 generators, class 2.
- 16: 32 states, 4096 transports, 8 generators, class 3.
- 32: 64 states, 2^21 transports, 10 generators, class 4.
- 64: 128 states, 2^38 transports, 12 generators, class 5.
- 128: 256 states, 2^71 transports, 14 generators, class 6.

This separates two growth rates: the number of distinct exact transport summaries is exponential in the algebra dimension D=2^n, requiring D/2+log_2(D) bits, while the number of independent graph-cycle labels needed is only twice the doubling index. The Clifford anticommutation law for all distinct imaginary basis operators does not persist; successive doublings introduce genuinely higher noncentral commutators. A distinguished smaller Clifford family does persist, as described next.

## 4. A smaller regular Clifford transport survives every doubling

Let H_n be generated by the n canonical doubling operators

    L_1,L_2,L_4,...,L_(2^(n-1)).

They obey L_a^2=Z and L_a L_b=Z L_b L_a for distinct canonical indices in every n. For the anticommutator assertion, equation (1) reduces to the previous dimension when both indices are below the newest bit; when one is the newest bit and the other is a lower coordinate index, its value is the constant one. Induction proves the claim.

Every word therefore reduces to a sign times one ordered product of distinct canonical generators. There are at most 2^(n+1) normal forms. All are distinct because the selected subset is exactly the induced index translation, and the two signs differ. Thus

    |H_n|=2^(n+1).

Its orbit of e_0 contains every signed basis unit, so its action is regular. In particular, the n canonical algebra generators reach all sheets, but they do not generate the full left-transport group G_n when n>=3. Their composed permutations must not be identified with nonassociative unit multiplication.

The survival of this smaller Clifford group is compatible with the dimension-16 failure below: L_10 is not a canonical doubling generator.

## 5. Sharp capacity of a fixed graph

Let B be a connected finite graph, keeping loops and parallel edges when present, and beta=|E|-|V|+1. A permutation-voltage cover of B has monodromy generated by beta fundamental-cycle transports. Therefore:

    A full G_n transport cover over B exists iff beta>=2n, n>=3.  (2)

Necessity is the abelianization lower bound. For sufficiency label a spanning tree by identity, label 2n remaining edges by the explicit generators above, and label the rest by identity; use inverse transports on reverse traversals. Transitivity gives a connected cover. This is an existence statement allowing arbitrary G_n operators as edge labels. It does not prove that any prescribed natural edge labeling generates the whole group.

By contrast, a connected cover on all signed-unit states using some G_n-valued labels exists iff beta>=n. Necessity follows because a transitive monodromy subgroup must surject onto the n-dimensional index-translation quotient. Sufficiency follows by labeling n fundamental cycles with the canonical generators of H_n. That construction is a regular 2^(n+1)-sheet cover. Thus connectedness alone must not be confused with full G_n monodromy.

The full-G_n cover in (2) has 2^(n+1) sheets and deck group of order two. Forgetting the sign gives a regular 2^n-sheet translation cover; the sign layer is a regular double cover. Their composite is nonregular for every n>=3. No nonidentity deck translation of the intermediate cover lifts to a deck transformation of the composite.

The smallest simple base graph supporting full G_n has

    v_min=ceil((3+sqrt(1+16n))/2)

vertices, since a complete graph on v vertices has beta=(v-1)(v-2)/2 and attains this maximum. Thus octonion transport first fits on K_5, sedenion transport on K_6.

[S166's verified order-six non-self-converse black quotient Q_1](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/05-knowledge/results/session166_planar_metagraph_20261103/reflection_voltages/REFLECTION_VOLTAGES.md) has 22 vertices and 157 individual edges, hence beta=136. These counts are inherited from S166 and hard-coded here, not independently recomputed by this packet. Consequently it supports this full formal left-unit model exactly for 3<=n<=68. At n=68 the algebra dimension is 2^68 and the formal signed sheet set has 2^69 states; n=69 is impossible on this fixed graph because it needs 138 independent cycles. The same fixed graph supports a merely connected canonical H_n model through n=136, but not n=137. These huge numbers are exact capacity statements, not enumerated constructions. The lifts have 22*2^(n+1) vertices and 157*2^(n+1) edges. No claim identifies these lifted vertices with tournament isomorphism classes or their edges with physical complement lines.

For contrast, the one-skeleton of the three-cube has beta=5 and cannot carry any full G_n with n>=3. The four-cube graph has beta=17 and can support full G_n for 3<=n<=8. Filling the square faces changes the problem: the full cube is simply connected, so arbitrary nontrivial graph monodromy cannot extend across all its faces. Graph and cubical-complex claims must stay separate.

## 6. The norm and Clifford boundary is exactly visible at 16 dimensions

In the displayed basis convention, put

    x=e_1+e_10,    y=e_4-e_15   in A_4.

Exact multiplication gives

    xy=0,        ||x||^2=||y||^2=2.

Thus the Euclidean norm is not multiplicative. Moreover,

    x^2=-2,
    x(x e_4)=-2e_4-2e_15,
    (x^2)e_4=-2e_4,
    (L_1 L_10+L_10 L_1)e_4=-2e_15.

So alternativity and the basis-operator Clifford relation fail, while each separate L_a remains a signed permutation and preserves the norm. The same obstruction embeds unchanged in every A_n with n>=4. [Dimensions 1,2,4,8 retain the familiar composition norm](https://math.ucr.edu/home/baez/octonions/node5.html); the explicit witness establishes failure at every later doubling. A nonzero vector can therefore have positive quadratic norm and still define singular left multiplication. This is not a contradiction to the invertibility of the signed basis-unit transport maps.

## Reproduction and scope

Run `python3 check_structure.py` and `python3 -O check_structure.py`. The script uses integer and binary arithmetic, no third-party packages, and no Python assertions. It checks the multiplication recurrence, all sign-kernel generators, every generator-pair commutator, Frattini ranks, complete lower-central ranks, and the displayed norm/alternativity obstruction. It independently enumerates the group in dimensions 8 and 16; higher dimensions use exact binary subspace closure, not an infeasible full group enumeration.

The all-n structural proofs are not inferred by interpolation from these finite checks. The kernel companion provides the all-n algebraic prerequisite (K). No theorem here solves tournament spectrum classification, proves a CA equivalence, or promotes an observation-state cover into a physical tournament construction.
