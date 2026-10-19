# Four upper exclusive neighbors: complete seven-thin maximum-span exclusion

Date: 2026-10-19. Status: exact finite morphology theorem, with an independent standard-library replay. This is not a proof of the entire seven-thin balance conjecture.

## Theorem

Let P be a finite poset of incomparability degree at most seven. Suppose P has no balanced incomparable pair, where balanced means probability in [1/3,2/3] under the uniform distribution on linear extensions. Let x precede y in its strict >2/3 majority order, and suppose x,y is an incomparable pair of maximum index span in that order. Then it is impossible that

    D(x)=D(y), U(y)⊆U(x), |U(x)\U(y)|=4.

The order-dual exclusion follows as well. The exterior of the displayed neighborhood is arbitrary, with no bound on its size or on its internal order.

Write T=U(x)\U(y), Z=I(x)∩I(y), and c=|Z|. Thus x<a and a∥y for every a∈T. The root y has incomparability degree 5+c, forcing c≤2. The proof below handles all three values c=0,1,2. It does not extend the false four-slave common-neighbor aggregate inequality.

## 1. Core family and completeness

Set Q={x,y}∪Z∪T. Label x=0, y=1, Z={2,...,c+1}, T={c+2,...,c+5}. Choose natural labels separately on Z and T along a linear extension of each induced order. The complete core relation family consists of:

- arbitrary naturally labeled orders on Z and T;
- x<a for every upper vertex a;
- arbitrary relations z<a, closed under both internal orders;
- no other relations involving either root or between Z and T.

The direction a<z is impossible because it would imply x<z. For each a, its common predecessors form an ideal of Z; these ideals increase along the order on T. This gives exactly 40, 357, and 4,824 distinct labeled cores for c=0,1,2. Every core has at most eight vertices, so every core incomparability degree is automatically at most seven. These are complete finite neighborhood families, not a sample of ambient posets.

The independent verifier regenerates this family by a different construction: enumerate all subsets of the allowed forward nonroot relations, compute transitive closure, and deduplicate. Equality of the resulting set with every recorded core is checked explicitly.

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

## 5. Complete finite certificate and witness-location classification

The full computation has these exact totals:

| Common vertices c | Cores | Cores excluded by balanced defect screen | Orders excluded by directed screen on remaining cores | Residual orders certified by generators | Distinct generator certificates |
|---|---:|---:|---:|---:|---:|
| 0 | 40 | 0 | 630 | 9 | 4 |
| 1 | 357 | 40 | 12,320 | 530 | 171 |
| 2 | 4,824 | 813 | 427,489 | 12,251 | 1,787 |

Each of the 12,790 residual orders has a certificate with exactly one pair, weight one. Moreover every one of these pairs is of the form

    a before y, with a∈T.

Thus the surviving orders have a particularly simple witness location: an exclusive upper endpoint and the root y. All 1,962 distinct certificates verify 3N_g(a<y)−2E_g≤0 for every retained generator. No LP certificate, floating-point decision, or multi-pair weighted inequality is used in the result. The search program includes more general fallback code, but the recorded result never reaches it.

This is a complete exact finite classification of witnesses by core and order. It is not an assertion that every a∈T always has p(a<y)≤2/3, nor a claimed short symbolic classification of exactly which a works in every core. The receipt and verifier establish the required existential witness for every remaining order.

## 6. Reproduction and proof status

Files in this directory:

- scan_four_upper.py: deterministic exhaustive search and exact integer checks (NumPy/SciPy are imported, but no solver is used in this successful run).
- certificates.json.gz: full core records, positive comparison witnesses, exact maximum generator margins, and digests of order-to-witness assignments.
- summary.json: scan totals.
- verify_four_upper.py: independent, standard-library-only certificate replay.
- replay.json: independent verification receipt.

Run:

    python -O verify_four_upper.py

The verifier independently regenerates the entire allowed core set, recomputes all ambient-defect upper bounds, enumerates every necessary projected majority order, reconstructs every retained boundary generator, checks every certificate inequality in integer arithmetic, and verifies complete order coverage and assignment digests. Explicit guards remain active under python -O. The finite case check plus the arbitrary-exterior lift proves the theorem; no exterior-size cutoff is involved.

The familiar core x<a1<a2<a3<a4, z<a4, y isolated remains a counterexample to the old common-neighbor aggregate inequality (48/35>4/3). It does not invalidate this result, which uses actual endpoint comparisons and does not assume that aggregate inequality.

## 7. Attribution and priority limits

The underlying problem and thin-poset framework are classical. See G. Brightwell and C. D. Wright, *The 1/3–2/3 Conjecture for 5-Thin Posets*, SIAM Journal on Discrete Mathematics (1992), DOI https://doi.org/10.1137/0405037, and E. J. Olson and B. E. Sagan, *On the 1/3–2/3 Conjecture*, https://arxiv.org/abs/1706.04985. These references provide context, not a claim that they contain this exact finite receipt. Majority-order, exchange-injection, linear-extension counting, and nonnegative-coefficient arguments are general methods, and no originality or priority is claimed for them. All lemmas needed by this certificate are proved here; a comprehensive literature comparison has not been performed. The result is a verified research record of the precisely stated morphology exclusion, not a novelty claim.
