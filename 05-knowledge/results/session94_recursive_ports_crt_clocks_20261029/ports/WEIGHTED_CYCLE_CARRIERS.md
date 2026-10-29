# From shared named ports to weighted cycle carriers

**Status: PROVED**, with independent literal finite controls. Scope: the explicit cyclic tree actions below; for the ternary case these are the inherited tournament automorphism and reversal-centralizer actions. This is a family-specific extension using elementary cycle-index counting, not a claim of external novelty or a new formula for unrestricted A000568.

Source pin: `88c969273f2b`. Inheritance:
- [`PORT_MOMENTS.md`](../../lucas_prime_reflection_towers_20261028/PORT_MOMENTS.md): exact all-order ordered-port moments and critical ternary law.
- [Session93 port/LCA carrier](../../session93_reversal_lucas_clocks_20261029/ports/PORT_ORBITS_AND_LCA.md): exact folded diagonal action and complete named-tuple invariant.
- [Session79 prime cycle index](../../session79_tournament_enumeration_20261027/PRIME_CYCLE_INDEX.md): aggregate cycle actions, rather than group orders alone, suffice for unrefined composition. Its unrestricted species input is not rederived here.

The closest mechanism is the shared diagonal action; the hostile is the two mirrored blocks treated independently. The corrected near miss is using fixed-point moments to count invariant subsets. The underused sidecar is the whole cycle inventory. Anchor: weighted attachments; niche: unordered port sets; wildcard: composite odd cyclic alphabets, where proper-divisor rotation clocks must survive.

## 1. Full cycle-index theorem

For odd `d>=3`, let `G_(d,0)=1` on a singleton and

    G_(d,k+1) = G_(d,k)^d ⋊ C_d

in its natural action on length-`k+1` words over `Z/d`. Let alpha negate every digit, and let `D_(d,k)` be its centralizer in this group. The case `d=3` is the inherited `A_k,C_k`. For other d the theorem needs only this specified permutation action; identification with the **full** automorphism group of a tournament is a separate obligation.

Write

    Z_G(p) = |G|^(-1) sum_g product_(i>=1) p_i^(c_i(g)),

where `c_i(g)` counts cycles of length i. Put `h=(d-1)/2`. Then

    Z_(G_(d,0)) = p_1,
    Z_(G_(d,k+1)) = (1/d) sum_(l|d) phi(l)
                          Z_(G_(d,k))(p_l,p_(2l),...)^(d/l),       (1)

    Z_(D_(d,k)) = p_1 product_(j=0)^(k-1)
                          Z_(G_(d,j))(p_1²,p_2²,...)^h.         (2)

For d=3, (1) is `(Z_k^3+2 Z_k[p_i -> p_(3i)])/3`.

**Proof of (1).** A rotation of order l has d/l block cycles of length l. Around each block cycle, the product of its independent uniform child automorphisms is uniform in the child group. Products on different block cycles are independent. A child cycle of length i yields a vertex cycle of length li. Exactly phi(l) rotations have order l. Average these conditional cycle laws.

**Proof of (2).** A root rotation centralizing digit negation satisfies `r=-r`; odd d forces r=0. Child0 carries `D_(d,k-1)`. The h pairs of nonzero opposite children each carry one arbitrary `G_(d,k-1)` element and its alpha-conjugate. Thus each child cycle occurs **twice at its original length**, and the h choices are independent. Recurse down the fixed child.

This is a cycle-multiplicity square `p_i -> p_i²`, NOT a cycle-length doubling `p_i -> p_(2i)`. It is also NOT two independent factors `Z_G(p)^2`: the same element governs both paired copies. These distinctions are the exact algebraic form of retained port alignment.

For d=9, root rotations of order3 contribute `(2/9) p_3^3`. An identity/full-9-cycle-only model already fails at the first level. The exact number of three-subsets of a 9-cycle up to rotation is10.

## 2. Weighted palette and attachment composition

A finite color palette has weights `w_a` in a commutative polynomial ring. A coloring has the multiplicative weight `product_v w_(color(v))`. This is invariant under vertex permutations. Burnside with weights gives its orbit inventory by substituting

    p_i = sum_a w_a^i.                                         (3)

Every orbit contributes its weight once; nonnegative integer coefficients count literal colored choices, and signed formal weights give a signed invariant inventory.

For a one-variable palette `Q(t)=sum_a t^(cost(a))`, including repeated colors with the same cost, define

    P_0(t;Q)=Q(t),
    P_(k+1)(t;Q)=(1/d) sum_(l|d) phi(l) P_k(t^l;Q)^(d/l).     (4)

Then the full centralizer inventory is

    B_k(t;Q)=Q(t) product_(j=0)^(k-1) P_j(t;Q(t)^2)^h.         (5)

For signed coefficients in Q, interpret their absolute values as numbers of colors and their signs as color-weight signs. Every cycle length in these odd-order groups is odd, so raising a sign to the cycle length preserves it. This is why `Q(t^l)` is valid for the signed tests too; the same shortcut would not be valid for even cycles. General polynomial color weights should instead use (3) literally.

The seed `Q²` is the palette of **ordered pairs of colors**, not two independent orbit counts. A cycle in the folded child must carry a constant ordered pair, giving `(sum_a w_a^i)^2` exactly.

Different palettes or cost variables on each immutable level/side tag are permitted: replace the pair seed by the product of the two side palettes, separately for each j and each opposite pair. Thus source labels and side costs need not be collapsed. The theorem counts attachments made independently at individual vertices under the specified group action. A response involving edges, Hamilton paths, or joint port guards needs its corresponding equivariant data; palette weights alone do not model it.

Example: binary colors on the nine vertices of W_2 have192 centralizer orbits. The three fixed vertices give factor8. The coupled pair of three-vertex blocks has24 orbits, because `(4^3+2*4)/3=24`. Treating the two blocks as independent binary necklaces gives only `8*4*4=128`; incorrectly replacing every paired cycle by one doubled cycle gives32.

### Exact small-degree implementation

Equation (4) is evaluated directly in `Z[t]/(t^(R+1))`. Every substitution increases degrees and every product has nonnegative degree, so truncation is exact, not an approximation. Formula (5) can be accumulated using one evolving child inventory and one product. For fixed d, naive polynomial convolution takes `O(k R²)` integer arithmetic operations and `O(R)` retained coefficients. For variable d, repeated-squaring gives an upper bound `O(k tau(d) log(d) R²)`; integer bit costs are additional. The full group and the full cycle-index polynomial are never enumerated by this evaluator. A literal group implementation appears separately only for hostile verification.

## 3. Unordered ports: a new polynomial hierarchy

Take `Q=1+t`, and let `U_(d,k,r)=[t^r] B_k`. These count orbits of r-element subsets, not ordered tuples and not multisets.

**Theorem.** For each fixed r>=1, `U_(d,k,r)` is a polynomial in k of degree r, with leading coefficient `(d-1)^r`.

**Proof.** More generally, for a seed Q with constant term1 and linear coefficient q, the coefficient `[t^r]P_k(t;Q)` is a polynomial in k of degree at most r-1, with leading coefficient `h^(r-1) q^r`. Induct on r. In (4), the all-degree-r contribution to `P_k^d/d` reproduces the previous coefficient. The terms using exactly two positive-degree factors have coefficient h and degree r-2. Their combined leading coefficient is

    h sum_(i=1)^(r-1) (h^(i-1)q^i)(h^(r-i-1)q^(r-i))
       = (r-1) h^(r-1)q^r.

Terms with three or more positive factors have lower degree. Terms from nonidentity rotations have l>=3, since d is odd; they only use coefficients through floor(r/l), so also have lower degree. Discrete summation proves the statement, including possible degree drop when q=0.

For (5), use the paired seed `(1+t)^2`, so q=2. In `P_k^h`, the coefficient of degree i has leading term `h^i 2^i k^(i-1)`; all contributions occupying two or more factors have lower degree. Let L_r be the leading coefficient of U. The product recurrence `B_(k+1)=B_k P_k^h` gives

    r L_r = sum_(i=1)^r L_(r-i) (2h)^i,    L_0=1.

Thus `L_r=(2h)^r=(d-1)^r`. Polynomiality follows by the same lower-order induction. The factor Q initializes `B_0`; it does not change the leading coefficients. ∎

For the ternary tournament family, the first exact polynomials are

    U_(3,k,1)=2k+1,
    U_(3,k,2)=k(4k-1),
    U_(3,k,3)=k(2k-1)(4k-3).                                (6)

At k=2 the entire inventory is

    1,5,14,30,46,46,30,14,5,1.

This hierarchy differs from the inherited named-port moments, whose leading coefficient is `2^r r!`. Division by r! recovers only the leading asymptotic coefficient, not the exact count: stabilizers can nontrivially permute a selected subset.

### The first cycle correction is explicit

For a group element in the ternary family, all cycles have odd length. A fixed three-subset is either three fixed points or one 3-cycle. Therefore

    U_(3,k,3) = E binom(Y_k,3) + E c_3(g).

The root recurrence gives `E c_3(A_j)=2j/3`; (2) gives

    E c_3(C_k)=2k(k-1)/3.

Combining this with the inherited first three moments yields (6). At k=2, the fixed-point-only expression is `86/3`, whereas the answer is30. Even the integrality failure warns that an ordered-tuple moment cannot just be reinterpreted as an invariant-subset count. The missing information is exactly the cycle coordinate, not another arbitrary one-bit correction.

**Stronger hostile: every named-port moment can agree.** The regular C9 and C3×C3 actions on nine points have the identical fixed-point distribution: nine fixed points for the identity, zero for each of the other eight elements. Hence every ordered r-port orbit count is `9^(r-1)` for r>=1 in both actions. Their cycle indices are respectively

    (p_1^9+2p_3^3+6p_9)/9,
    (p_1^9+8p_3^3)/9.

Their three-subset orbit counts are10 and12. Thus the full fixed-point-moment hierarchy does not universally recover subset attachment inventories. Within a known one-parameter recursive family, knowing k of course identifies its cycle index by (1); the obstruction concerns transferring moments without that action sidecar.

At depth11 (177147 vertices) the coefficients through12 are

    1,23,473,9471,184690,3534850,66629112,1239471068,
    22791196703,414731863513,7475458522191,
    133567834293993,2367154817599362.

At depth12 (531441 vertices) they are

    1,25,564,12420,267036,5645596,117730932,2426329620,
    49488811394,1000065739530,20039163042624,
    398432373417408,7864960699520088.

These are port-subset orbits for two specified recursive groups, not values of A000568 at11 or12 and not unrestricted tournament class counts.

## 4. Positive depth templates and fast-forwarding the recurrence

The polynomial theorem has a constructive strengthening in the ternary family. For an r-subset, retain its immutable first-nonzero levels. At each such level, jointly fold the suffixes from the two sides as in S93, and retain every branching LCA level among distinct folded suffixes. Call the union the **active depth set**.

**Saturated-template theorem.** Let beta_(r,l) count C_l-orbits of r-subsets whose active depth set is every level0,...,l-1. Then

    U_(3,k,r) = sum_(l=0)^r beta_(r,l) binom(k,l),              (8)

with nonnegative integer beta. In particular the Newton/binomial coefficients are actual compressed attachment templates, not cancellation artifacts.

**Proof.** For each occupied first-nonzero level a, let m_a be the number of distinct folded suffixes among its selected vertices. Their branching trie has at most m_a-1 branching nodes, hence at most that many distinct branching depth levels. Counting the initial level a as well gives at most m_a active levels. Summing yields at most `sum_a m_a<=r`; the zero root needs none.

Delete all inactive digit positions from every selected vertex. The first nonzero position was retained, and every first difference of two folded suffixes in the same tag pair was retained, so no distinct selected vertices merge. The resulting subset at depth l has all l levels active. Its folded pairwise LCA and cyclic orientation data are exactly the original data with depths order-compressed. By the complete S93 invariant, this construction is well-defined on orbits. Conversely, choose l ordered levels from k and insert zero digits at all other positions in a saturated template. Its invariant is the unique original invariant with those expanded depths. Compression and expansion are inverse on orbits, although the expanded representative need not be literally the original vertex set. This proves (8). ∎

The first binomial-coefficient rows, indexed by l starting at0, are

    r=1: 1,2
    r=2: 0,3,8
    r=3: 0,1,28,48
    r=4: 0,0,46,336,384
    r=5: 0,0,46,1312,4672,3840.

For example, the46 nine-vertex four-subset orbits are the two-active-depth templates; the336 new three-active-depth templates appear at depth3. The top coefficient is `beta_(r,r)=2^r r!`, agreeing with the leading-coefficient proof. The construction says exactly which LCA level information is retained; it does not claim that tags alone determine a template.

Equation (8) also removes linear dependence on a huge depth k from fixed-order evaluation. Compute inventories for depths0,...,R incrementally by (4)-(5), then take forward differences in each coefficient column. This costs `O(R^3)` arithmetic operations for fixed d using naive convolution, with `O(R²)` storage. Evaluating the whole inventory through degree R at any specified k costs `O(R²)` further arithmetic operations via binomial coefficients. Integer bit costs depend on log k. For d other than3, the proved degree bound still makes the same Newton evaluator exact; the explicit saturated-template interpretation here is asserted only for the ternary family. The implementation compares fast evaluation with direct depth50 iteration for d=3,5,9, then evaluates at depth1,000,000 without materializing a single vertex.

## 5. Weighted named ports and the exact source law

The inherited fixed-point law also has a fully tag-weighted form. Assign root weight z0 and weights u_j,v_j to the two mirrored tags at suffix depth j. With independent uniform `A_j` elements, put X_j equal to their fixed-point counts. Then the weighted number of fixed vertices is exactly

    S = z0 + sum_(j<k) (u_j+v_j) X_j.                        (7)

For ordered r-port tuples, allowing repeated vertices, `E S^r` is the weighted orbit sum where each orbit receives the product of its port-tag weights. This follows because the weight is constant on each tuple orbit and weighted Burnside applies. Its exponential moment series factors as

    E exp(tS) = exp(t z0) product_(j<k) E exp(t(u_j+v_j)X_j).

For tags on the same pair, fixed-count covariance is `Var(X_j)=2j`; for distinct j it is zero. Thus the two mirrored sides cannot be independently resampled even if their single-port responses agree. The complete S93 folded LCA/orientation invariant is still needed to identify actual tuple classes, construct transporters, or evaluate a native operation not determined by these weights.

For general odd d, the fixed population law is 0 with probability `(d-1)/d`, and the sum of d independent previous populations with probability1/d. It has mean1 and variance `(d-1)k`. There are h independent mirrored pairs at each level, so the centralizer fixed count has

    E Y = 1+(d-1)k,
    Var(Y) = (d-1)^2 k(k-1).

**Source boundary.** All expectations here sample group elements uniformly, or equivalently use Burnside's uniform group average. They do not sample distinct fixed-count values uniformly, and do not supply a law for native Collatz sources, valuation words, seed multiplicities, or clocks. For example, A_2 has fixed-count multiplicities `{0:62,3:12,6:6,9:1}`. Sampling these four support values uniformly gives mean9/2; sampling its81 elements uniformly gives mean1. The script checks the full signed tag-weighted fixed-source distribution, not just its mean. Vertex/color costs are explicit invariant object weights, not a license to change the averaging measure on group elements. An arithmetic transfer must preserve that separate source measure and its clock, or carry their exact sidecars.

## 6. Reproduction and limits

Run from this directory:

    python3 check_weighted_cycles.py
    python3 -O check_weighted_cycles.py

Both runs give identical output and1782 always-active checks. The checks compare full literal cycle laws and weighted palettes for A_0..A_2 and C_0..C_3; d=5 and9 include depth-two centralizers of orders25 and6561. They independently enumerate all512 subsets of the nine-vertex action, test wrong independent/doubled pair models, test the composite9 clock and the same-all-moments C9/C3×C3 control, compare truncations through12 and15, and verify finite differences through subset order8 for d=3,5,9. Actual subset orbits through four ports at depths2 and3 verify the intrinsic active-level counts and the saturated-template identity. Fast-depth checks compare against direct iteration at depth50 and test depth1,000,000. Finite differences are checks of the proved polynomial theorem, not its proof.

The decisive loss ledger is:

| Source | Target | Preserved | Lost | Repair |
|---|---|---|---|---|
| actual group action | fixed-point moments | named-port orbit counts | cycles moving selected ports internally | full cycle index |
| paired child action | independent marginal orbit counts | each single side's inventory | relative alignment | Cartesian-product palette under one action |
| full cycle index | small-degree inventory | weights/costs through R | higher costs and actual representatives | retain action/LCA carrier when future operations need it |
| uniform symmetry samples | arithmetic seed analogy | only a structural recursion resemblance | source masses, carry guards and clock | no transfer until an exact map is supplied |

No assertion here identifies marking loops with black tournament moves, establishes a Collatz drift theorem, solves arbitrary weighted Hamilton-path composition, or reconstructs individual prime tournament classes. The concrete gain is an exact, cheaply truncated attachment inventory and an explicit first obstruction to replacing cycle data by the entire fixed-point moment hierarchy.
