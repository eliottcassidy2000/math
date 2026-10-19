# Five and six upper exclusive neighbors: seven-thin exclusion

Date: 2026-10-19. Status: proved by exhaustive exact morphology certificates and independent standard-library replay. No general seven-thin theorem or novelty claim.

## Theorem and family

Under the hypotheses of the four-upper theorem, replace |U(x)\U(y)|=4 by either 5 or 6. These morphologies are also impossible. The order-dual statements hold. Put m=|U(x)\U(y)|, and c=|I(x)∩I(y)|. The root y has incomparability degree 1+m+c, so m=5 allows c=0,1, while m=6 allows only c=0.

Label roots x=0,y=1, common vertices 2,...,c+1 and upper vertices c+2,...,c+m+1. The core consists of arbitrary naturally labeled common and upper orders, relations x<a for every upper a, and any transitively closed common-to-upper relations. There are no other relations. Choosing a linear extension independently on the two nonroot groups proves completeness. Every common predecessor set of an upper vertex is an ideal, increasing along the upper order. The exact catalogue sizes are 357 and 4,824 for (m,c)=(5,0),(5,1), and 4,824 for (6,0). All cores have at most eight vertices, so all satisfy the core degree budget. The independent verifier also generates them by closing all subsets of allowed edges and comparing the entire resulting set.

## 2. Exterior and a sharpened directed-defect screen

Every vertex outside Q is comparable with both roots. It therefore belongs to A below both roots or B above both roots; the other two possibilities would compare x and y. We have A<B. No Q vertex can lie below an A vertex or above a B vertex. Moreover

    A<a for every a∈T,

by A<x<a. This additional forced relation is useful and should not be discarded in a generic slack bound.

For incomparable u,v define the directed defect

    δ(u,v)=|D(v)\D(u)|+|U(u)\U(v)|.

The exchange lemma says δ(u,v)≤1 implies p(u<v)≤2/3. For completeness: when δ=0, swapping u and v injects extensions with u before v into those with v before u. If the only defect is an extra successor t of u, split the former extensions into u<v<t and u<t<v. Swapping u,v works on the first class. On the second, replace the three occupied labels (u,t,v) by (v,u,t). This is valid because D(v)⊆D(u) and every successor of t is a successor of v: any successor of t is a successor of u distinct from its sole extra successor t. Each class injects separately into v<u. Thus the forward class has at most twice the reverse count. A sole predecessor defect follows by order duality.

Let d_Q(v) be core incomparability degree, and s(v)=7−d_Q(v). Define valid upper bounds on missing exterior comparisons by

    α(v)=0 for roots or upper vertices; otherwise α(v)=s(v),
    β(v)=0 for roots; otherwise β(v)=s(v).

Here α bounds the number of A elements incomparable with v, and β bounds the analogous B count. For every incomparable core pair,

    δ_P(u,v) ≤ δ_Q(u,v)+α(u)+β(v) =: H(u,v).

Indeed the exterior differences are |A_v\A_u| and |B_u\B_v|, where A_v is the set of A elements below v and B_v the set of B elements above v. These are bounded by missing A relations of u and missing B relations of v, respectively.

The following two exact screens are therefore sound for every actual exterior:

1. If H(u,v)≤1 and H(v,u)≤1, the pair is balanced and the entire core is excluded.
2. In a proposed majority order M, if u precedes v and H(u,v)≤1, that order is excluded.

Separate slack bounds may overestimate what an actual exterior can realize, which only weakens these sufficient screens.

## 3. Why the projected order test is necessary

With no balanced pair, every incomparable pair has exactly one direction of probability >2/3. Together with comparable pairs, these directions form a transitive tournament: a directed triangle would have three cyclic probabilities summing to more than two, whereas each linear extension satisfies at most two of the cyclic indicators. Its order M extends P.

All exterior A letters precede x and y, and all B letters follow both. Consequently no exterior letter is between the roots, and their span in M restricted to Q equals their full span. Projection cannot increase another core-pair span. Therefore the projected order must be a linear extension of Q with

    rank(y)>rank(x),
    |rank(v)−rank(u)|≤rank(y)−rank(x)

for every incomparable core pair u,v. All such projected orders are enumerated exactly. This is only a necessary condition, so retaining extra unrealizable orders is harmless.

## 4. Exact boundary-generator lift

For a core extension o, let pre(o) precede its first root and suf(o) follow its last root. Define L(u) as the number of interleavings of A with the fixed ordered word u, divided by e(A), and R(v) dually using B. Then the number of ambient extensions restricting to o is

    e(A)e(B)L(pre(o))R(suf(o)).

All A is before the first root, all B after the last. All other core/exterior relations are automatically respected outside these two boundary interleavings, proving this factorization.

Let U_empty=V_empty=1 and, for nonempty words,

    U_u=L(u)−L(u with last letter removed),
    V_v=R(v)−R(v with first letter removed).

These numbers are nonnegative: append the removed prefix letter, or prepend the removed suffix letter, to obtain an injection. Expanding along the two word trees gives a nonnegative linear combination of generators G_(u,v), where G_(u,v) counts the core extensions beginning with u and ending with v and having both roots between these boundary words.

If d_Q(s)=7, all exterior points are comparable with s, hence A<s<B. Every increment whose boundary word contains s vanishes, including increments later along that word-tree branch: after the first saturated letter in a prefix, all A has already appeared, and the normalized prefix count remains constant; the suffix statement is dual. Thus generators whose prefix or suffix includes a saturated vertex may be removed. No independence or realizability of the other generator coefficients is assumed.

For a retained generator g let E_g count its extensions, and N_g(a<b) count those having a before b. An exact certificate consisting of positive integer weights w_ab on forward incomparable pairs proves a contradiction whenever

    sum_ab w_ab [3 N_g(a<b)−2 E_g] ≤ 0

for every retained g. Multiplying by the actual nonnegative boundary coefficients and summing proves the same inequality for ambient extension counts. Yet every certificate pair is forward in M, so each ambient summand is strictly positive. This is impossible.

## 5. Full exact counts and witness location

| m | c | Cores | Analytically balanced cores | Directed-screen orders on remaining cores | Residual orders | Distinct certificates |
|---|---|---:|---:|---:|---:|---:|
| 5 | 0 | 357 | 0 | 15,533 | 473 | 135 |
| 5 | 1 | 4,824 | 602 | 485,760 | 9,468 | 1,255 |
| 6 | 0 | 4,824 | 391 | 527,387 | 7,974 | 1,012 |

All 17,915 residual orders have an exact one-pair certificate a<y for some upper a, with weight one. All 2,402 distinct witnesses satisfy 3N_g(a<y)−2E_g≤0 for every retained boundary generator. No LP call or floating-point certificate is used. This is an existential finite classification over all necessary orders, not a claim that every upper vertex always works.

## 6. Independent replay

The `five_upper` and `six_upper` directories contain the complete compressed certificates, scan summary, independently implemented standard-library verifier, replay receipt and manifest. Run `python -O verify_upper.py` in each directory. The verifier regenerates every core independently, proves the defect-screen inequalities, checks all necessary projected majority orders, builds every retained boundary generator afresh, checks every certificate in integer arithmetic, and verifies coverage digests. No ambient-size bound is used. Explicit guards remain enabled under -O.

## 7. Attribution

The classical setting is G. Brightwell and C. D. Wright, *The 1/3–2/3 Conjecture for 5-Thin Posets* (1992), https://doi.org/10.1137/0405037; see also E. J. Olson and B. E. Sagan, *On the 1/3–2/3 Conjecture*, https://arxiv.org/abs/1706.04985. General majority-order, exchange-injection, counting and nonnegative-certificate methods are not claimed as original here. These external references provide context, not an attribution of this exact receipt to those authors. The necessary proofs are included above, and no literature-priority claim is made.
