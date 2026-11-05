# The full opposed-binary carrier remains arithmetically sparse

Session 197, 2026-11-05. **All-order scalar-support theorem and exact bounded-image certificate. The full `{7,21}` tournament Hamiltonian-path spectrum remains open.**

Even allowing **every** binary attachment word, two opposed apexes above a transitive chain do not give a cofinite Hamiltonian-path spectrum. Their distinct strong counts through X number

    O((log X)^2 X^(log_3(2))).

After arbitrary ordered joins and transitive padding, the number of distinct counts is still O(X^(2/3)). This is a stronger obstruction than forbidding a particular word pattern. It needs neither Busch's order bound nor the odd-cycle collection formula: its two lower bounds are explicit injections into actual Hamiltonian vertex words. No literature-priority claim is made.

## 1. Physical class and inheritance audit

Fix a->b and a transitive chain x_1->...->x_m. A word w over {1,2} records

    1: a->x_i->b,       2: b->x_i->a.

Write T(w) for this actual tournament and H(w) for its number of directed Hamiltonian vertex words. This is an intrinsic graph class too: T admits two vertices a,b such that deleting them is transitive and every other vertex is beaten by exactly one of them. Relabel the pair so a->b and list the unique transitive order of the remainder to obtain a word.

The empty word and all-1 words are transitive and have H=1. Every word containing a 2 is strong. That 2-vertex gives b->x->a, so the apexes lie in one component. The component reaches the first chain vertex and is reached by the last one. The increasing chain then puts all its vertices in the same component.

The full four-letter two-apex carrier is **not** the class considered here. In that larger carrier the extra letters are `0: x->a,b` and `3: a,b->x`. They are forbidden throughout this theorem, but appear in the explicit escape at the end.

Closest named sources:

- [S191 two-apex transfer](https://github.com/eliottcassidy2000/math/blob/8d04f8557ea50c985521f039d5f58aa335617c9a/05-knowledge/results/session191_tournament_growth_20261104/two_apex/TWO_APEX.md) supplies the inherited exact three-slot path bijection and nonnegative finite transfer.
- [S191 sparsity](https://github.com/eliottcassidy2000/math/blob/8d04f8557ea50c985521f039d5f58aa335617c9a/05-knowledge/results/session191_tournament_growth_20261104/two_apex/SPARSITY.md) handles polynomially many fixed-run templates and their multiplicative closure. The unrestricted binary language here has exponentially many words and unbounded run count.
- [S192 boundary limits](https://github.com/eliottcassidy2000/math/blob/8d04f8557ea50c985521f039d5f58aa335617c9a/05-knowledge/results/session192_surprising_approaches_20261105/boundary_limits/PROOF.md) gives the general deterministic boundary-allocation mechanism and random boundary laws. The middle-allocation idea is inherited. The present sharp alphabet-versus-slot count and resulting full binary scalar sparsity are separate consequences, not claimed in that source.
- [S196 cloning](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/spectral-moves/PROOF.md) includes the inherited cyclic triangle-clone count used below; it is also reproved directly here.

The session's preceding forbidden-pattern packet gave a different odd-cycle lower bound. That result remains valid, but is not an input here. In particular the present theorem improves its Fibonacci-controller conclusion simply because the Fibonacci language is a subset of this full binary language.

## 2. Direct elementary facts

**Induced-path injection.** If U is an induced subtournament of T, then H(U)<=H(T). Fix an order of the omitted named vertices and insert them one at a time. For a new vertex v, insert it immediately before the first current path vertex it beats, if any; otherwise append it. At an internal insertion, the preceding vertex beats v because it was not beaten by v. The result is directed. Deleting the inserted names recovers the old word, so the procedure is injective.

**Triangle-clone count.** In a cyclic triangle, replace one vertex by a transitive chain of length l>=1, retaining all uniform exterior arcs. The count is

    2^l+1.                                          (1)

There are three paths with the clone chain in one run. With two nonempty clone runs, the other two vertices must occur consecutively in their directed order between those runs. Every ordered nonempty partition of the l chain names gives one path, with each part internally increasing; there are 2^l-2 such partitions. No third run is possible. The total is 3+(2^l-2)=2^l+1. This is a literal path count.

## 3. The two outer runs must be short

Take a nonconstant word w of length m. Let a>=1 be the length of its maximal initial constant run, b>=1 the length of its maximal terminal constant run, and c=m-a-b. The two runs are disjoint, so c>=0.

**Lemma 1.**

    H(w)>=2^a+1,       H(w)>=2^b+1.                 (2)

**Proof for the initial run.** If it consists of 2s, retain all its vertices and both apexes: the induced graph is exactly the cyclic triangle-clone graph `a_apex -> b_apex -> run -> a_apex`. If it consists of 1s, its next chain vertex y has letter 2. Retain the whole run, y, and apex a: the incidences are `a -> run -> y -> a`. Apply (1) and induced monotonicity.

**Proof for the terminal run.** If it consists of 2s, again retain both apexes. If it consists of 1s, its preceding chain vertex y has letter 2. With apex b the cycle is `b -> y -> run -> b`. This gives the second bound. The selected exterior vertices are distinct from every run vertex in all cases. QED.

In particular, if H(w)<=X and X>=3, then

    a,b <= L(X)=floor(log_2(X-1)).                  (3)

This controls the boundary runs regardless of how complicated the middle word is. No assumption about cloning the whole input is being made.

## 4. Three independent placements for every middle vertex

Assume now c>=2. Put sigma=w_1 and tau=w_m. Maximality of the two end runs gives the unique decomposition

    w = sigma^a (bar_sigma) z (bar_tau) tau^b,
    bar_1=2, bar_2=1, r=|z|=c-2.                    (4)

Let P be the initial block `sigma^a bar_sigma` and Q the terminal block `bar_tau tau^b`. Both contain a 1 and a 2, and they are disjoint. The arbitrary middle word z lies strictly between them.

**Lemma 2 (literal free-slot injection).**

    H(w)>=3^r.                                      (5)

**Proof.** Use Hamilton words in the apex order a,b, with increasing chain slots

    S_0, a, S_1, b, S_2.

For such a word the only chain/apex conditions are:

- the last vertex of S_0 must have letter 2;
- the first and last vertices of S_1 must have letter 1;
- the first vertex of S_2 must have letter 2.

There is no first boundary on S_0 and no last boundary on S_2.

Choose one 1-vertex p_1 and one 2-vertex p_2 in P, and one 1-vertex q_1 and one 2-vertex q_2 in Q. Fix the following allocation:

    p_1,q_1 -> S_1;     p_2 -> S_2;     q_2 -> S_0.

Put every remaining P-vertex into S_0 and every remaining Q-vertex into S_2. Independently assign each of the r named vertices of z to any one of S_0,S_1,S_2. Sort each slot in its inherited chain order.

Now S_0 ends at q_2, S_1 starts at p_1 and ends at q_1, and S_2 starts at p_2. All four required boundary incidences hold. All slots are nonempty, every chain name appears once, and every within-slot arc is forward. Thus every one of the 3^r assignments gives an actual spanning directed path.

The assignment is recovered from the path by locating each named middle vertex relative to a and b. Different assignments therefore cannot collide. This proves (5), including r=0. QED.

This proof does not confuse path slots with contiguous blocks in the input word: the individual subsets of z can interleave arbitrarily before being placed in their three different path slots.

## 5. Complete presentation count through a Hamilton ceiling

Let B be the set of distinct H-values greater than one occurring in this full binary class. Fix X>=3, and define L as in (3) and

    R(X)=floor(log_3 X),       alpha=log_3(2).        (6)

There are four exhaustive cases.

1. **Constant words.** All-1 words contribute only one and are discarded from B. All-2 words have H=2^m+1 by (1), so at most L nontrivial candidate words remain.
2. **c=0.** There are exactly two constant runs, with opposite letters. At most 2L^2 candidates remain.
3. **c=1.** The first and last letters must agree and the single middle letter is their opposite. Again at most 2L^2 candidates remain. The four-anchor proof is intentionally not applied because its boundary blocks would overlap here.
4. **c>=2.** For each pair (a,b), each pair (sigma,tau), and each r<=R, equation (4) allows exactly 2^r middle words. Lemma 2 forces r<=R whenever H<=X. Thus at most `4L^2 sum_(r=0)^R 2^r` candidates remain.

Every relevant strong tournament presentation lies in this list; the list may contain candidates whose H exceeds X. Summing yields the explicit bound

    # {h in B:h<=X}
      <= L + 8L^2 2^R
      <= L + 8L^2 X^alpha,
    alpha=log(2)/log(3)=0.6309297535714574... .        (7)

The same bound holds for the number of ordinary isomorphism classes of strong binary outputs through that H-ceiling, since every class has at least one such presentation. It is not a count of all A000568 objects, and no uniqueness of apex presentation is assumed.

## 6. Arbitrary ordered products are still sparse

Let <B> be the set of finite products of the distinct integer generator values in B, including one. Ordered joins of any number of binary outputs and arbitrary transitive padding have exactly counts in <B>.

For every real s>alpha, (7) gives convergence of `sum_(b in B)b^(-s)` by dyadic summation. Since every b>=3, the Euler product

    Z_B(s)=product_(b in B)(1-b^(-s))^(-1)

converges. Expanding its finite partial products counts exponent vectors of generator values. Multiple factorizations of an integer only overcount, so

    sum_(n in <B>) n^(-s) <= Z_B(s),
    #(<B> intersect [1,X]) <= Z_B(s) X^s.             (8)

In particular alpha<2/3 because 2^3<3^2, and hence

    #(<B> intersect [1,X]) = O(X^(2/3)).              (9)

The support has natural density zero. Its reciprocal sum is finite, while the prime reciprocal sum diverges, so it omits infinitely many odd primes. This statement concerns the construction class, not whether any particular omitted prime is globally attainable.

The Euler product counts **distinct scalar values**, not ordered graph decompositions. Many nonisomorphic ordered joins can have the same product. No subexponential bound for the number of all unlabelled joined tournaments is being asserted.

## 7. Exact all-order images and a physical escape

The parameter bounds (3),(5) give a terminating, exact bounded-H enumeration without any guessed graph-order cutoff. `verify.py` enumerates the four categories above, computes every candidate's exact H, and then closes the resulting integer generators under multiplication.

For X=19, L=4 and R=2, so the complete candidate list has exactly

    L+8L^2 2^R = 516

words. Its scalar image, including one, is

    {1,3,5,9,11,13,15,17}.

Therefore H=19 is absent from this entire opposed binary class at every order. Since 19 is prime, arbitrary ordered products do not repair the omission.

The extra aligned letter 0 escapes the obstruction. In the usual four-letter convention the actual word `2020` has

    rows=[2,20,57,51,33,3], H=19.                    (10)

This graph is strong. Its scalar count was known before this session; the point here is its role as an explicit escape from the all-order binary restriction, not a claim that H19 is newly discovered.

For X=501 the exact candidate list has 16,392 words. It gives:

- 96 distinct strong generator H-values;
- 97 values after also including one;
- 139 values after arbitrary ordered multiplication, including one.

In particular 67 and173 are also absent from this binary/product class at every order. The separate two-port packet physically realizes them, so the present result is a concrete limitation of the opposed-binary architecture rather than a global exclusion.

`bounded_image.json` records a literal binary word for every positive generator and a generator factorization for every product value. The finite output has all-order completeness only because the boundary parameter bounds prove that every relevant word was included.

## 8. Reproduction and independent controls

Run in this directory:

    python3 verify.py > normal.stdout
    cp receipt.json receipt.normal.json
    cp bounded_image.json bounded_image.normal.json
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout
    cmp receipt.normal.json receipt.json
    cmp bounded_image.normal.json bounded_image.json

The verifier is standalone, uses only the standard library, and keeps all tests active under optimization.

It checks all 511 binary words of lengths zero through eight against independent adjacency subset-DP. It verifies 988 actual induced triangle-clone graphs, and for 396 words explicitly constructs and checks all 8940 distinct free-slot paths specified by Lemma 2. This tests both their arcs and injectivity.

Every one of the 16,392 canonical candidates through501 is counted by two separately implemented exact transfers: the generated general fixed-apex slot matrices and a direct specialized three-slot update. All 97 positive generator witnesses are additionally counted by graph subset-DP. The H19 escape and the exact integer cutoff calculations are independently checked.

The all-order proof relies on the actual path injections, the canonical boundary decomposition, and the convergent Euler product, not on the finite tests. No asymptotic approximation, random boundary law, OCF identity, or external strong-order theorem is required.
