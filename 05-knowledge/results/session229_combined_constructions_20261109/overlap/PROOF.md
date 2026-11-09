# Exact carry fibers, nearest-neighbor defects, and a scalar entropy barrier

**Status: PROVED within the explicitly specified block constructions; finite witnesses are FINITE-EXACT. Chan–Pak completeness remains OPEN.**

## 0. Inheritance and scope

This packet extends, rather than reruns, Session228's equal-alpha construction in `05-knowledge/results/session228_recursive_constructions_20261108/unsaturated_module/PROOF.md`. That packet supplies the physical six-state transfer convention, the `(H,Q)=(19,256)` block, and a separated sparse exponential family. Its statement that singleton separators kill defects is correct; we strengthen it to **every nonempty intervening block**, and determine the interactions when no separator exists. Session228 `mobile_triangle/GENERAL_DIVISIBILITY.md` supplies the existing general SCC divisibility obstruction. Section 5 gives an independent short proof of its repeated-module consequence, without claiming the general theorem anew. Session225's affine-return recovery, Session226's finite-line scope, and Session227's adaptive constructions explain why a large menu is not itself scalar coverage.

Our live concepts are equal-alpha switches, genuine adjacent interactions, sidecar carries, fixed-template entropy, and physical root insertion. Anchor: turn same-sidecar freedom into overlapping scalar choices. Niche: classify exactly how alpha digits can compensate. Wildcard: count all masks, not just the inherited switch menu. Canonical hostiles are the original superincreasing family, the H19 pair with zero interaction, and the H9 carry example below. No unrestricted interval, prime coverage, or completeness claim is made.

All results here concern literal tournaments, not numerically admissible nonphysical states. We make no literature-priority claim.

## 1. Physical definitions

Delete a root r and take an ordered list of nonempty tournaments B_i, all arcs from earlier to later blocks. For B_i let H_i count Hamilton paths, and Q_i count ordered pairs (P,R) of directed paths covering the block, allowing an empty path. Let alpha_i count such covers with R nonempty and starting at a vertex beaten by the root. Consequently

    0 <= alpha_i <= Q_i-H_i.

The excluded H_i covers have R empty. The physical rooted state is (D,E,U,d), where D is the Hamilton-path count, E counts Hamilton paths ending at r, d counts Hamilton paths after deleting r, and U equals d plus the root-started two-path-cover count of the deletion graph. For these blocks,

    d = product_i H_i,
    U = d + sum_i alpha_i (product_(j<i) H_j)(product_(j>i) Q_j).       (1)

Indeed, choose the first block occupied by R. Earlier blocks are entirely P, later blocks carry unrestricted ordered two-path covers. This is the inherited physical sidecar formula, now also used for carries.

For strong witnesses flank the list by singleton p before all B_i and singleton q after all B_i, with r->p and q->r. Then every vertex is reachable from r through p and reaches r through q. E=d=product H_i. For t repetitions of B,

    U = H^t + 2 Q^t + 2 sum_(i=0)^(t-1) alpha_i H^i Q^(t-1-i).      (2)

The flanks change neither equality of the middle alpha expressions nor the scalar interaction formulas below.

## 2. Every nonempty block kills long-range interactions

Use the inherited six states (p,r)=(0,0),(0,1),(1,0),(1,1),(2,0),(2,1), with column-vector convention. Here p records empty P, P ending at a root-inneighbor, or P ending at a root-outneighbor; r records empty R versus R already starting at a root-outneighbor. A local cover category (pp,rr) is in {0,1,2}^2 except (0,0), with rr=1/2 recording a root-outneighbor/inneighbor as its first vertex. Reject rr=2 when old r=0; otherwise update to `(pp if pp nonzero else p, bool(r or rr))`.

Write e_j for coordinate j, and set

    W = span(e1-e5,e3-e5),
    V = W + span(e2-e4),
    rho(z) = z0+z2+z4,  w=e3-e5.

**Theorem 1.** Every nonempty physical block transfer M satisfies

    M(V) subset W,   M|W = H I,   rho M = H rho.                     (3)

On W, all covers with P nonempty cancel because they erase the old p coordinate; covers with P empty retain p and contribute H. On e2-e4, the P-nonempty covers again cancel. The P-empty covers now require a full R starting at a root-outneighbor; if their count is a, the image is a(e3-e5). This proves the first two statements. The last statement counts precisely covers with R empty, of which there are H, regardless of the previous P endpoint. Nonemptiness is essential: an empty block acts as the identity and does not send V into W.

Let F=M'-M compare two modules with the same H,Q,alpha. Let x be the difference in the number of full Hamilton paths starting at a root-outneighbor, and y the difference in the number ending at a root-inneighbor. Their eight category differences necessarily have the form

    (0,1):x, (0,2):-x, (1,0):y, (2,0):-y,
    (1,1):a, (1,2):b, (2,1):-x-a, (2,2):x-b.                     (4)

Equality of alpha gives the seventh entry; equality of Q gives the eighth. Applying the update rule gives

    image F subset V,  F W=0,  F(e2-e4)=x w,
    [coefficient of e2-e4 in Fz] = y rho(z).                        (5)

Therefore for any such defects F,G,

    FG = x_F y_G w rho,   FGH=0                                   (6)

for any third such defect H. Also F C G=0 whenever C is any product containing at least one nonempty physical block transfer. This last claim follows from (3), even if C contains no singleton. This is the precise strengthening over Session228.

## 3. Exact nearest-neighbor scalar compiler

At position i choose any physical module from a finite menu with common H_i,Q_i,alpha_i, and fix baseline M_i. Denote a choice's defect F_i by F_i(c), with zero defect allowed. Let D_0 be the baseline scalar and L_i(c) its one-position scalar difference. Then the scalar of every word is exactly

    D(c_1,...,c_t) = D_0 + sum_i L_i(c_i)
                      + sum_(i=1)^(t-1) J_i(c_i,c_(i+1)),
    J_i(a,b) = y_i(a) x_(i+1)(b) product_(j not in {i,i+1}) H_j.    (7)

This holds with or without the strong singleton flanks. Expand the product of transfers. Every term involving three or more defects vanishes by (3),(5),(6); a term involving two nonadjacent defects has a nonempty baseline block between them and vanishes. The adjacent pair is (6). Every earlier physical block scales rho by H and every later block scales w by H. The initial vector has rho=1 and the final scalar row ell=(1,1,1,1,0,0) has ell(w)=1. These observations give (7). Thus interactions are real but only quadratic and nearest-neighbor: there is no hidden higher-order carry in this fixed-alpha transfer mechanism.

For k choices at each of t sites, minimum and maximum D follow by ordinary chain dynamic programming in O(t k^2) transitions: store the best sum with each last choice and append L_(i+1)(b)+J_i(a,b). Replacing each best value by the set of attainable values gives an exact target-membership algorithm. Set sizes may be exponential, so this is not a polynomial-time target oracle and does not establish interval coverage.

For the inherited block

    B=[40,1,3,6,15,30], H=19,Q=256,
    alpha weights=[21,20,29,35,73,59],
    full-path ends=[1,8,5,3,1,1], starts=[1,1,2,3,7,5],

the complete list of alpha fibers larger than one is

    alpha94: 17,40; alpha108:24,38; alpha114:19,42;
    alpha123:21,44; alpha129:25,39; alpha143:23,46.

All have size two. For 17->40, (x,y)=(0,-2), so F^2=0 even without any separator. For 24->38, (x,y)=(-2,-10), so F^2=20 w rho is nonzero. Two adjacent copies therefore have scalar mixed difference 20. The checker builds all four strong 15-vertex graphs and checks their full Hamilton counts by an independent subset-path dynamic program. Removing separators is not falsely claimed to manufacture dense scalar intervals.

## 4. Exact finite carries and coprime rigidity

**Theorem 2.** Repeat a fixed block t times. Write g=gcd(H,Q), H=g h, Q=g q, so gcd(h,q)=1 and q>h. Two alpha words have the same U if and only if their differences delta_i=alpha'_i-alpha_i satisfy

    delta_i = h c_i - q c_(i-1),  i=0,...,t-1,
    c_(-1)=c_(t-1)=0,  c_i integers.                              (8)

Every such carry arising from allowed digits obeys

    |c_i| < g, hence |c_i| <= g-1.                                (9)

Divide the equal-U equation by g^(t-1). It becomes sum delta_i h^i q^(t-1-i)=0. Reducing modulo h shows h divides delta_0. Write delta_0=h c_0, cancel, and repeat, which proves necessity and integrality. Substitution telescopes, proving sufficiency. Equivalently the forward difference polynomial sum delta_i z^i is divisible by q z-h; (8) fixes the signs without relying on that convention.

For the strict carry bound use the suffix identity

    c_i = -sum_(j=i+1)^(t-1) delta_j h^(j-i-1)/q^(j-i).

Since |delta_j| <= Q-H=g(q-h), its absolute value is at most

    g [1-(h/q)^(t-1-i)] < g.

Thus if g=1 all carries vanish and the entire alpha word is fixed by U. This does **not** imply unique root masks: masks within each equal-alpha fiber remain free. In particular all fixed-U root-mask words on the fixed repeated H19 deletion graph are Cartesian products of singleton or the six displayed binary fibers. Cross-position alpha compensation is impossible there, even if arbitrary root masks are allowed.

A mixed-block sufficient rigidity condition also follows directly from (1): if

    gcd(Q_i, product_(j<i) H_j)=1 for each i,

reduce equal-U expressions modulo Q_last, use |delta_last|<Q_last, then cancel Q_last and induct. All alpha digits are again unique. This condition is sufficient, not necessary.

**Physical hostile showing why noncoprime carries matter.** The five-vertex block

    C=[16,1,3,7,14], H=9,Q=84,
    alpha weights=[9,10,16,25,15]

has alpha words (10,56) and (16,0) with the same weighted sum 1344. Their difference (6,-56) has carry c_0=2, attaining g-1 for g=3. Reversing the comparison gives carry -2. The explicit mask words are (2,28) and (4,0). With strong flanks their scalars are respectively 13149 and 16719, with the same (E,U,d)=(81,16881,81). The checker independently counts both whole graphs. Hence rigidity without its arithmetic hypothesis would be false. The carry automaton has at most 2g-1 states, an exact finite description rather than a heuristic numerical cone.

## 5. What paying for carries can cost

**Corollary (inherited divisibility, specialized).** For t repetitions of a block with g=gcd(H,Q), with arbitrary root masks and optional singleton flanks,

    g^max(t-2,0) divides D.                                      (10)

For completeness, decompose each rooted Hamilton path as P,r,R. If both P,R are nonempty, mark the block containing the last vertex of P and the block containing the first vertex of R. These are at most two exceptional blocks. Every other nonempty block lies wholly in one path (contributing H), or strictly between the marked blocks and may be covered by an arbitrary ordered path pair (contributing Q). Fixing the marked blocks and their endpoint-constrained covers therefore leaves a product with at least t-2 factors divisible by g. Empty P or R only reduces the number of exceptional blocks. Summing proves (10). Singleton flanks contribute H=1,Q=2, and can be exceptional instead of a repeated block, never weakening the bound.

Thus in this homogeneous framework: coprime H,Q eliminate cross-site alpha carries; noncoprime H,Q permit a finite carry alphabet but impose a growing scalar divisor. For prime targets, t>=3 and g>1 leave only the possible prime D=g when g itself is prime; this is already impossible with strong flanks because D>=Q^t>g (Section 6). This dichotomy does not forbid local equal-alpha switches, changing the block template, or arbitrary constructions outside this framework.

## 6. A scalar sparsity barrier encompassing every mask

**Theorem 3.** Fix any nonempty tournament B on n vertices. Form strong rooted graphs using t copies in forward deletion order and the flanks p,q of Section 1, allowing every one of the 2^(nt) root masks on the middle vertices. If Q is its ordered two-path-cover count, then every such graph satisfies

    D >= Q^t.                                                   (11)

For each ordered path cover (L,R) of the whole middle graph, the sequence

    L, q, r, p, R

is a Hamilton path: the flank arcs connect every nonempty segment and q->r->p is fixed. This is an injection, including cases L or R is empty. The middle ordered-cover count is Q^t, by restriction to each forward block. This proves (11) independently of all root-mask details, alpha collisions, and defect formulas.

Consequently the number of distinct attainable D<=X in the union over all t>=0 is at most

    sum_(0<=t<=floor(log_Q X)) 2^(nt)
       <= [2^n/(2^n-1)] X^beta,   beta=n log(2)/log(Q).           (12)

Every tournament has at least one Hamilton path, by vertex insertion. Therefore Q=sum_S H(B[S])H(B[S^c])>=2^n. If B is nontransitive it contains a directed triangle; its three Hamilton paths make at least one summand strictly larger than one, so Q>2^n and beta<1. The entire family then has natural density zero among integers, and also among odd integers. This is stronger in scope than showing one equal-alpha binary fiber is sparse: even all masks, all U, all adjacent interactions, and all lengths fail to cover all sufficiently large odd values.

For H19/Q256 on six vertices, beta=3/4. Thus even the union of **all** 64^t masks at every length has only O(X^(3/4)) distinct scalars <=X. The statement is tied to this fixed deletion template and strong flanks; it is not a universal Chan–Pak obstruction.

There is also a finite-bank version. For allowed fixed blocks B_j with n_j vertices and cover count Q_j, if some 0<s<1 satisfies

    F(s)=sum_j 2^(n_j) Q_j^(-s) < 1,

then all finite words from this bank, with the same flanks and arbitrary masks, have at most X^s/(1-F(s)) scalars <=X. Bound each word by its mask count, weight it by (product Q_j)^(-s), and sum the geometric series over word length. The sufficient condition can fail even though each individual block is sparse. F(1)>=1 is **not** a density or coverage theorem; it merely means this particular sufficient counting bound is unavailable.

## 7. Verification and limitations

Run `python3 verify.py > CHECK.json` and `python3 -O verify.py > CHECK.optimized.json`; the outputs must be identical. The standalone standard-library checker uses no imported repository files and does not rely on Python assertions. It verifies all eight category generators, the universal four-dimensional defect basis, every basis pair and triple, all 64 masks on H19, all binary choices of mixed equal-alpha menus through length nine, exact dynamic-programming extrema, literal whole-graph Hamilton counts for the interacting and carry witnesses, and exhaustive bounded-difference carry equivalence for the listed seven (H,Q,t) universes. Literal permutation cover enumeration and subset-path Hamilton counting are separate computational routes. The proof, not these finite tests, carries all-order quantifiers.

The new mechanisms identify exact barriers and give an exact scalar compiler. They do not produce overlapping intervals of unbounded length, rule out larger multi-choice modules, or settle the missing unrestricted counts. Next escapes must change the template/grammar, use sufficiently rich mixed block banks, or obtain a structural overlap theorem not contradicted by these bounds.
