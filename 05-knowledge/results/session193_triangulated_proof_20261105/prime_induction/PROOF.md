# Prime-atom induction, a compulsory reset, and a mixed infinite fan

**PROVED elementary reduction and construction; PROVED finite-certificate obstruction using a CITED order bound.** Session 193, 2026-11-05. The unrestricted tournament Hamilton-path spectrum conjecture with exceptions `{7,21}` remains **OPEN**. The lower arm of the infinite fan is explicitly **CONDITIONAL** on a smaller attainable factor.

## 1. Inheritance and the complementary-method question

Write `H(T)` for directed Hamiltonian vertex words and `M` for the positive odd integers other than 7 and 21. Ordered joins multiply H: all paths in `A -> B` concatenate a path in A with one in B. This gives a physical implementation of multiplication, including disconnected-in-the-strong-sense outputs.

Closest proved mechanism: [Session 192 arithmetic approaches](../../session192_surprising_approaches_20261105/dyadic/APPROACHES.md), section 3, reduces the conjecture to existential ternary closure `h -> 3h+2,3h+4`, with finite seeds. Canonical hostile: [Session 191 adaptive-end obstruction](../../session191_tournament_growth_20261104/two_apex/ADAPTIVE_ENDS.md) blocks fixed-two-apex end growth along `13 -> 41 -> {125,127}`. Corrected near miss: a scalar parent does not supply a common reusable marked interface. Underused sidecar: the multiplicative atom structure of the *desired* support, coupled to rooted two-path-cover coefficients rather than only H.

The concept board is: multiplicative atoms; prime-only descent obligations; feedback-vertex heredity; an old-arc reset; rooted affine wrappers; bounded-template sparsity. The anchor is a genuine mixed-method induction. The niche is the stronger 41-to-127 induced obstruction. The wildcard, which succeeded, is lifting a finite old-arc repair to an exact all-clock affine family.

## 2. The exact multiplicative atoms

Define

    P = {odd primes other than 7},
    Q = {7p : p is an odd prime, p>=5, p!=7},
    A = P union Q union {49,63,343}.

**Theorem 1 (exact atom decomposition).** The multiplicative monoid M has exactly the irreducible elements A. Every member of M is a finite product of elements of A, and no element of A is a product of two members of M greater than one.

**Proof.** Write `n=7^a m`, with `7` not dividing `m`. If `a=0`, factor m into odd primes different from 7. If `a>=2` is even, split off copies of 49. If `a>=3` is odd, split off 343 followed by copies of 49. These leave the prime factors of m.

The remaining case is `a=1`. If m has a prime factor `p>=5`, then `p!=7`; split off `7p`, then the remaining prime factors of m. Otherwise m is `3^b`. Since n is not 7 or 21, `b>=2`, and `n=63*3^(b-2)`. This proves generation.

A prime is irreducible. The proper two-factor splits of `7p` all use 7; those of 49 use 7; those of 63 use 7 or 21; and those of 343 use 7. Thus all displayed generators are irreducible in M. Conversely every generator decomposition of a nongenerator has at least two factors. QED.

**Consequences for the proof program.**

1. If all members of A have tournament witnesses, arbitrary ordered joins realize M. Thus the infinite work can be split into an ordinary-prime lane and a seven-times-prime lane, with three finite exceptional atoms.
2. A least missing member of M must belong to A. Every nonatom is already handled by multiplication of smaller admissible integers.
3. Assuming the inherited exclusions of 7 and 21, a witness for an atom has exactly one nontrivial strong component. Indeed all component counts lie in M, and a nontrivial factorization would contradict irreducibility. It is therefore enough, and necessary after removing neutral singleton components, to find **strong** witnesses for atoms. A nonstrong product does not manufacture a missing prime or 7p atom.
4. This is a target reduction, not a theorem that an arbitrary existing parent tournament can be extended. The witnesses used for different targets, or for different digits at one parent, need not be the same marked tournament.

## 3. A strictly smaller sufficient ternary hypothesis

The following reduction allows unrelated construction methods to handle different residual targets.

**Theorem 2 (prime-lane ternary induction).** Let S be a multiplicatively closed set of positive odd integers containing

    1,3,5,23,49,63,67,343.

Suppose that whenever `h in S`, `h>=3`, `h not in {7,21}`, and `d in {2,4}`, one has

    3h+d in P union Q  ==>  3h+d in S.                 (1)

Then `M subset S`.

**Proof.** Apply strong induction to a member N of M. If N is a nonatom, Theorem 1 supplies a factorization into smaller members of M, so multiplication finishes it. All finite exceptional atoms and the small seed atoms are supplied. Every remaining atom is in P or Q and is nonzero modulo 3. Its unique odd inverse parent is

    h=(N-2)/3 if N=2 mod3,
    h=(N-4)/3 if N=1 mod3.

It is smaller than N. The boundary `h=1` gives 5 or 7, already treated; `h=7` gives 23 or 25, of which the only atom is seeded 23; and `h=21` gives 65 or 67, of which the only atom is seeded 67. All remaining parents are admissible and at least three. Induction gives `h in S`, and (1) gives N. QED.

For actual tournament counts, multiplication is an actual ordered join. The seven nontrivial seeds have these verified presentations:

| H | Word | Apex count |
|---:|---|---:|
| 3 | 10 | 1 |
| 5 | 100 | 1 |
| 23 | 11000 | 1 |
| 49 | 10212 | 2 |
| 63 | 1001100 | 1 |
| 67 | 12301 | 2 |
| 343 | 1012000 | 2 |

One singleton supplies 1. Word conventions are defined in section 5; the verifier generates exact adjacency rows in `receipt.json`. The universal physical implication (1) is still **OPEN**. This theorem does not assume fixed macros, one common witness for both children, return of the output to an input interface, or strong children for composite targets.

A useful multi-method formulation is equivalent: assign any number of direct construction families or conditional operations to the residual targets in `P union Q`; it suffices that their union covers every unseeded target under the smaller-parent induction hypothesis. Methods may cover overlapping sets and may discard every unused state coordinate. What must be proved is coverage of the target lanes and actual tournament existence, not closure of each method by itself.

## 4. Fixed-complexity families necessarily miss both infinite atom lanes

This connects Theorem 1 to the already proved [Session 191 sparsity proof](../../session191_tournament_growth_20261104/two_apex/SPARSITY.md). Suppose a bank B of nontrivial attainable strong counts, all at least two, has

    #{b in B : b<=X} = O((1+log X)^d)

for some fixed d. The cited packet proves, by the convergent Euler product over the generators, that its multiplicative closure C satisfies

    sum_(c in C) c^(-s) < infinity   for every s>0,     (2)
    #{c in C : c<=X} = O_s(X^s).

Finite unions of banks and finitely many extra seeds do not change the conclusion. The hypothesis applies to any fixed feedback-arc budget, a finite set of fixed transitive-clock templates, or a fixed number of constant runs in the two-apex language. Unlimited ordered joins are already included in C.

**Corollary.** Each such construction class misses infinitely many members of P and infinitely many members of Q.

Indeed if it contained all but finitely many members of P, its reciprocal sum would diverge by Euler's elementary theorem `sum_p 1/p=infinity`, contrary to (2) at `s=1`. The same argument applies to Q, whose reciprocal sum is `1/7` times a prime reciprocal sum with only finitely many terms removed. For completeness, prime reciprocal divergence follows because convergence would make `product_p (1-1/p)^(-1)` finite, while its finite partial products dominate the harmonic sum over any fixed initial interval once their primes include that interval's prime factors.

This is a precise boundary on a proof program, not a new global absent Hamilton-path count. The residual target lanes have zero integer density, but neither can be covered by finitely many fixed clock templates plus multiplication. Template or run complexity must become unbounded, or a genuinely different family must enter.

## 5. No induced H41 survives inside any two-apex H127

For a word `w=c_1...c_m` over `{0,1,2,3}`, put apexes a,b with `a->b`, and chain vertices `x_1->...->x_m`. Bit 0 of `c_i` records `a->x_i`, and bit 1 records `b->x_i`. Every tournament with feedback-vertex number at most two has such a presentation, after choosing and ordering the two apexes. This class is hereditary under induced subgraphs. With one apex, use binary letters recording its outgoing arcs to the chain.

**Theorem 3 (all-order induced obstruction).** If T has feedback-vertex number at most two and `H(T)=127`, then T contains no induced subtournament U with `H(U)=41`. T and U may be nonstrong and arbitrarily padded. The apexes may be selected again at any stage; the assertion does not preserve their identities or restrict insertions to chain ends.

**Order cutoff used.** The inherited external theorem is Busch's bound

    H(C) >= 5^((|C|-1)/3)

for strong tournaments C (A. H. Busch, EJC 13 (2006), N3, DOI `10.37236/1141`). The source's full proof is not independently re-proved here. The integer comparisons `5^7>41^3` and `5^10>127^3` imply strong-core orders at most 7 and 10 respectively.

**Complete finite certificate.** The verifier enumerates **every** two-apex presentation with chain length zero through eight, exactly

    1+4+...+4^8 = 87,381 words.

It computes all counts in two independent ways: the inherited 36-state run-assignment transfer, and a new literal-path Boolean-polynomial/zeta calculation described in section 8. Every value, including all negative-filter cases, agrees. There are exactly 40 strong H127 words. This is a presentation count, not an isomorphism-class count.

The complete strong H41 inventory with at most seven vertices consists of these 21 words:

    3200,3220,3320,
    10020,10201,10310,10321,11310,11321,12011,12031,
    13101,13231,20012,20030,20132,21332,30332,31011,
    31031,32331.

There is no smaller such word; hence every strong H41 member of this hereditary class has order six or seven. The verifier separately exhausts all 127 one-apex words through chain length six and finds no H41, which proves that its strong core cannot omit either chosen apex. No uniqueness-of-isomorphism-class claim is used for H41.

For every one of the 40 strong H127 presentations, independently compute H on every six- and seven-vertex subset. Across the inventory this is exactly **3,456** induced-subset checks. None has count 41. The full H127 inventory and each target's subset count histogram are retained in `receipt.json`.

**Why this covers nonstrong padding and all apex choices.** Suppose U were induced in T with H(U)=41. SCC multiplicativity and primeness give a strong component D of U with H(D)=41. D lies inside some strong component C of T. Every induced Hamilton path injects into a Hamilton path of a larger tournament by choosing a legal insertion for each new vertex, so `H(C)>=H(D)=41`. Since H(C) divides 127, it equals 127. Both C and D inherit feedback-vertex number at most two, and the order bounds apply to these strong components, not to arbitrarily padded T or U. C is therefore one of the enumerated presentations up to relabeling; D has six or seven vertices and occurs among the checked subsets. Contradiction.

The induced monotonicity just used is elementary: to insert z into a directed path, put it first if z beats the first vertex, last if the last beats z, or between consecutive vertices where the relation to z changes from beating z to being beaten by z. Fix the first legal slot. Removing z recovers the old word, so this is injective.

Positive controls show this is not an absent target:

    H41:  word3200,   rows [6,12,56,49,35,3];
    H127: word122132, rows [102,216,250,241,225,194,128,1].


Consequently an induced `41 -> 127` construction must leave the feedback-vertex-number-two class. Alternatively, changing an old arc can stay within the general spectrum program. This is stronger than S191's fixed-apex end-growth obstruction. It does not forbid arbitrary reconstruction, old-arc changes, or induced extensions with three or more feedback vertices.

## 6. A literal finite repair and the complementary product

Here is a nonovershooting physical route from H41 to H127. Rows are outgoing bitmasks on names `0,...,n-1`:

| H | Rows | Operation from previous row |
|---:|---|---|
| 41 | `[14,48,26,50,33,5]` | parent |
| 23 | `[14,48,58,50,33,1]` | reverse old `5->2` to `2->5` |
| 37 | `[14,48,122,50,33,1,59]` | add6 with outgoing mask59 |
| 63 | `[14,48,250,50,33,1,187,59]` | add7 with outgoing mask59 |
| 127 | `[14,48,506,306,33,1,443,315,51]` | add8 with outgoing mask51 |

All five tournaments are strong. Each birth preserves every existing name and incidence. Every count is at most 127, but the route deliberately first decreases the count. Full literal permutation enumeration and subset DP independently verify the five counts. The initial old-arc reversal is exactly the operation forbidden by an induced-only program.

The sibling target `125=3*41+2` needs no such ancestry at all: take three ordered-joined copies of the four-vertex one-apex `100` tournament, each with H5. This is an unconditional H125 witness, independently checked by subset DP. Thus the previously blocked pair is resolved by **different actual operations**, not a repaired common end-growth state.

## 7. A nontrivial infinite mixed-method fan

The finite reset suggested searching for matching rooted cover coefficients. For a core A and root r, let F_r(A) count ordered spanning two-path covers `(L,R)`, where L is nonempty and ends at r, and R may be empty. The inherited S190/S191 wrapper W_k(A,r), for k>=1, adjoins a transitive block Z of length k and w, with

    w -> Z -> A,          r -> w -> A-r.

The exact path-cover bijection gives

    H(W_k(A,r)) = H(A)+2^k F_r(A).                    (3)

Paths starting at w traverse Z and then a core Hamilton path. Every other path has predecessor r at w; its two core blocks form a cover counted by F_r, and each Z vertex chooses the block before the first core path or the block immediately after w. This gives `2^k` choices. With a strong core and nonempty Z the wrapper is strong.

Keep physical root r=4 and use these cores:

| Core | Rows | H | F_4 |
|---|---|---:|---:|
| P | `[16,69,121,115,98,67,1]` | 41 | 80 |
| Q | `[18,68,121,115,98,67,1]` | 55 | 105 |
| R | `[146,68,249,243,226,195,129,2]` | 127 | 240 |

All three are strong. The P-to-Q move reverses the single old arc `1->0`; the Q-to-R move adds vertex7 with sole old out-neighbor1. Every other old incidence is unchanged. Under the wrapper, the newborn is also dominated by w and every vertex of Z. Physical root4, all old names, and every old exterior incidence stay unchanged throughout.

Put `x=2^k`, k>=1. Then the wrapped route is

    N_k=41+80x -> 55+105x -> 127+240x=3N_k+4.          (4)

It uses exactly one old-arc reversal followed by one birth. Both initial and intermediate counts are below the final count; for the intermediate the difference is `72+135x>0`. This is an unconditional construction for the parent and upper child for every k.

**The other arm is a different induction step.** Set

    q_k=25+48x < 41+80x=N_k.

Then

    3N_k+2 = 125+240x = 5q_k.                         (5)

The number q_k is odd, positive and never 7 or 21. **Conditional on q_k being attainable**, an ordered join with the H5 seed realizes the lower child. In particular the strong-induction hypothesis that all admissible smaller numbers are attainable discharges this arm. Arithmetic divisibility alone is not being promoted to unconditional existence of q_k. The zero-digit child `3N_k` is unconditionally realized by an ordered join with the H3 seed.

Thus at every parent in this infinite family, the ternary proof obligation splits: one branch has an exact old-arc/birth construction; one is a product of a seed and a strictly smaller admissible target; and one is the standard product with 3. At the atom level, the lower branch is never required because it is composite. The upper target is never divisible by 7: `127+240*2^k = 1+2*2^k mod7`, whose three possible residues are 3,5,2. Only its prime instances are unresolved atom obligations, and all are supplied by (4).

**Not a return theorem.** This fan does not cover arbitrary parent values. Its scalar children do not re-enter the same clock family: parents are `41 mod80`, whereas `3N,3N+2,3N+4` are `43,45,47 mod80`. Even repeating the same literal flip and sole-out-neighbor birth on the unwrapped upper core fails immediately: `127 -> 119 -> 269`, rather than 385. These are exact finite controls. A general covering family or a reusable richer interface is still needed for Theorem 2.

## 8. Reproduction, independent calculation and scope

From this directory run:

    python3 verify.py
    python3 -O verify.py

The stdout and deterministic `receipt.json` bytes agree; no `assert` is used. The complete word count is independently verified as follows. For each chain length, enumerate the two apex orders and every assignment of chain vertices to the three increasing slots around them. Each is a possible literal Hamilton word. Its adjacent apex/chain incidences impose at most four Boolean requirements. Write its indicator as

    product(required-one bits) * product(1-required-zero bits).

Expand into multilinear monomials and add coefficients over all literal words. A subset zeta transform evaluates the resulting polynomial at **every** one of the `4^m` possible two-apex attachment words. This shares neither transfer states nor the transfer recurrence with the 36-state algorithm. Every count agrees throughout the full 87,381-word universe, so the independent algorithm verifies the absence filter as well as positive witnesses.

Additional checks include all one-apex words through length6, all 3,456 relevant induced subsets by subset DP, seven physical seeds, all five finite-reset counts by permutation enumeration, the H125 product, and all three rooted cover profiles by literal permutation/cut enumeration against a subset recurrence. Fifteen actual wrappers (`k=1,...,5` for each of three cores) independently check the affine counts, strongness, literal old-arc reversal and birth. Their finite checks support the displayed formulas; the all-k theorem follows from (3), not extrapolation. The exact atom decomposition and inverse-parent boundary are additionally exhaustively checked through 100,000.

The exploratory `probe_*.py` and `probe_*.cpp` files are provenance, not dependencies of the certificate. The proof uses no new literature-priority claim. Theorem 3 depends on the cited Busch order bound; Theorems 1 and 2 and the displayed constructive fans are elementary. The unrestricted `{7,21}` spectrum, universal prime-lane closure and indefinitely reusable reset remain open.
