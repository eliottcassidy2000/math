# Marked operators, tensor recursion and three different orbit comparisons

**Status:** the three cited OpenAI paper theorems are **CITED / accepted premises for this session**, as requested by the owner. This note does not audit their proofs. The elementary interface identities and obstructions below are **PROVED**; the companion controls are **FINITE-EXACT**. No universal Collatz conclusion, independently measured source floor, or prize submission is claimed.

## 1. Inheritance and the objects that must remain marked

The closest proved mechanisms are [reciprocal four-channel carriers](reciprocal_four_channel_kernel_20261004.md), [the quadratic word-doubling bridge](quadratic_escape_rank_atlas_20261004.md), [recursive four-state tournaments](tournament_recursive_four_state_20261004.md), and [source/child common-future lifts](supplied_state_join_extensions_20261004.md). The recent [cycles, tubes and debts](collatz_cycles_tubes_debt_walk_openai_20261006.md) supplies the separate long-gap arithmetic question. These are different notions of recursion and exponent.

The hostile is a quotient which retains a slope or trace but forgets the carry, named source, or legality guard. The corrected near miss is to confuse a four-dimensional tensor carrier with four possible states. The least-used useful sidecar here is the **marked tensor basis**, together with the order of composition. The concept board is: guarded affine words; marked tensor powers; independent tensor-leg changes; free-group word metrics; adaptive witness precision; graph bipartization.

For a valuation word `w=(a_1,...,a_r)`, put

```
P=3^r, Q=2^(a_1+...+a_r), F_w^+(x)=(P x+B)/Q,
M_w^+ = [[P/Q,B/Q],[0,1]],       lambda=P/Q.
```

The carry starts at zero and a next letter `a` changes `(P,Q,B)` to `(3P,2^a Q,3B+Q)`. An odd integer follows exactly this word iff

```
P x+B = Q (mod 2Q).                                      (1)
```

Indeed final oddness determines the last division exactly, and reducing successively modulo the preceding powers of two forces all earlier divisions. Positivity and first-hit ROOT status are additional predicates. The companion uses a representative `x>2Q`, which prevents an intermediate ROOT: such a hit would imply `3^j x+B_j=2^(a_1+...+a_j)<=Q`.

## 2. The two Collatz orbit comparisons

### 2.1 Two sources with a shared word

If `x` and `y` both obey the same guarded word, then exactly

```
F_w^+(y)-F_w^+(x)=lambda*(y-x).                           (2)
```

Thus the finite-block logarithmic gap multiplier per odd step is `(r log 3-A log 2)/r`, where `A=sum a_i`. This statement stops applying as soon as the valuation words separate. For example,

```
7 -> 11 -> 17,             15 -> 23 -> 35
```

share word `11`; the gap grows from `8` to `18`, by `9/4`. The next valuations are respectively `2` and `1`. There is no continued shared-word exponent without a new guard.

Different words may instead meet. The exact example

```
7 --1123--> 5 <--1-- 3
```

has `(P,Q,B)=(81,128,73)` and `(3,2,1)`. The endpoint equality is `y=(27x+3)/64`, not an equality of slopes. This is the inherited lift `7+256t -> 5+162t <- 3+108t`, with its two separately retained words. A later common suffix cannot recover the discarded initial carry.

### 2.2 Plus and minus dynamics

On all signed odd integers let `U_s(n)=(3n+s)/2^v2(3n+s)`, `s=+1,-1`. Then

```
U_-(-n)=-U_+(n),
M_w^- = R M_w^+ R,          R=diag(-1,1).                 (3)
```

The reflection preserves each actual valuation. It exchanges the positive and negative half-lines. It does not identify the two dynamics restricted to positive inputs. In fact, for every positive odd `n`, exactly one of `v2(3n+1),v2(3n-1)` is `1`, and the other is at least `2`: the two numerators differ by `2` and are both even. Positive minus-dynamics has the cycle `5 -> 7 -> 5`, corresponding to the negative plus-cycle `-5 -> -7 -> -5`.

The fixed word maps have opposite anchors `+B/(Q-P)` and `-B/(Q-P)`, but the same `lambda` and reciprocal trace. Their common multiplier does not retain the sign, carry, positive basin, or legality of repeated words.

## 3. The exact `9/4` bridge, and its information boundary

For word-matrix squaring, `lambda` goes to `lambda^2`. The rational reciprocal trace

```
J=lambda+lambda^(-1)
```

goes to `J^2-2`. For the one-letter word `(1)`, these are

```
lambda=3/2 -> 9/4,       J=13/6 -> 97/36.                (4)
```

This is the inherited exact semiconjugacy between squaring and the third quadratic. It is a statement about repeated word presentations; an actual integer must still satisfy the repeated guards. It does not identify `9/4` with the matrix-multiplication complexity exponent.

The full marked tensor square is stronger than its character. Writing `M=[[lambda,beta],[0,1]]`, the marked entries of `M tensor M` include `lambda^2` and `lambda beta`. Since `lambda>0`, they recover `lambda` and `beta`; all of `M` is recovered. Moreover

```
(M N) tensor (M N)=(M tensor M)(N tensor N).              (5)
```

The four ordered channels `00,01,10,11` therefore carry a lossless recursive operator representation with unbounded exact coefficients. They are not four finite states or intrinsically the four tournament classes. The earlier tournament codec uses a different carrier and its own ordered block addresses.

There is a decisive obstruction to forgetting the marks. Regarding an invertible `2 by 2` matrix as a two-leg tensor, or as a trilinear tensor with a one-dimensional third leg, independent invertible changes on its legs take it to `I_2`. Thus all our word matrices are mutually restricted to one another in that tensor category. Every restriction-monotone tensor character gives the same value on them, even though `(1)` has multiplier `3/2` and `(2)` has multiplier `3/4`. Reciprocal trace is invariant under conjugation, not under arbitrary independent changes on both legs. A transfer from a tensor character to source payment must specify a more restrictive marked category or an additional observable.

## 4. What the Nine-Fourths paper contributes

[Matrix Multiplication Nine-Fourths, October 2, 2026](https://github.com/openai/math/blob/main/preprints/Matrix-Multiplication-Nine-Fourths-October-2-2026/paper.pdf), Theorem 1.1, gives the accepted arithmetic-complexity bound `omega<=9/4` over the complex numbers. Its tensor restriction, direct-sum and tensor-product operations concern bilinear algorithms, not literal time steps or integer bit cost. We use two exact constructions from that paper.

**Shared-label separation (Proposition 3.1).** When `M` sectors share one leg, use `5M` Fourier-labelled copies. The surviving character equation is

```
u-v+2(g-h)=0 mod 5M,      1<=g,h,u,v<=M.
```

Its absolute left side is at most `3(M-1)`, so no nonzero multiple of `5M` survives. Legwise weights `g^2`, `hu-h^2`, `-hv` then sum, on the surviving terms, to `(g-h)^2`. The zero-weight part therefore has `g=h` and `u=v`. It retains `M` complete blocks, each with an `M`-dimensional dot-product index. It is not `M^2` independent blocks. This is a concrete way to preserve a shared label while separating interactions. It does not restore information from a scalar after that label has already been deleted.

**The four-channel representation split (section 4.1).** The multiplication-compatible sequence has kernel generated by `D=uw-vs` and quotient given by diagonal substitution. At its first level, in characteristic zero,

```
V tensor V = Sym^2(V) direct-sum wedge^2(V),             4=3+1.
```

For an `SL_2` operator of trace `X`, the characters are `X^2-1` and `1`. This is precisely the representation-theoretic version of the inherited reciprocal kernel. The antisymmetric channel and its determinant must remain typed; deleting it is not a lossless four-to-three compression of the full operator. Neither construction supplies a source guard or a positive Collatz floor.

## 5. Free-factor isomorphism, orbit metrics and adaptive witness budgets

[An Isomorphism of the Free Group Factors, September 23, 2026](https://github.com/openai/math/blob/main/preprints/An-isomorphism-of-the-free-group-factors-September-23-2026/An-isomorphism-of-the-free-group-factors-September-23-2026.pdf), Theorems 1.1 and 1.2, supplies the accepted trace-preserving factor isomorphisms. Section 5.2 retains earlier word witnesses by choosing each later perturbation budget after its sensitivity is known. The distinction between a limiting distribution and generation of the ambient factor is essential to that construction.

Here is the elementary transferable estimate. If a word polynomial is `F=sum c_v v`, define `L(F)=sum |c_v| length(v)`. For unitary tuples, telescoping products shows

```
||F(A)-F(B)|| <= L(F) max_i ||A_i-B_i||,                 (6)
```

in operator norm or tracial `L_2` norm. Thus, once a witness and a positive margin are available, later refinements can preserve it by charging their cumulative error to this known sensitivity. The same estimate for a finite signed measurement uses its coefficient `ell_1` norm. This transports an already established margin; it does not prove that a missing source margin is positive.

There is also an exact marking obstruction. Suppose a trace-preserving surjective isomorphism `L(F_m)->L(F_n)`, `m!=n`, sent every canonical generator unitary to a canonical group-element unitary. Trace-freeness would make those image elements freely generate a subgroup `H` isomorphic to `F_m`. Their generated von Neumann algebra is `L(H)`. Surjectivity forces `H=F_n`, because a group element outside `H` is orthogonal to `L(H)` in the canonical group basis. Abelianization ranks then give `m=n`, a contradiction. Consequently the accepted unequal-rank isomorphisms do not preserve this entire group marking.

For comparison, spheres in the named free-group word metric have size `2d(2d-1)^(r-1)`. That metric growth is not transported by an unmarked factor isomorphism. Unitary group orbits have constant Hilbert norm; unitary conjugation preserves operator and tracial norms. None of these is the real growth of a Collatz integer. The separate tensor-power arithmetic exponent is a third type of exponent.

## 6. Amenability: a useful group-completion hostile

[Unitarizability Implies Amenability for Countable Groups, September 23, 2026](https://github.com/openai/math/blob/main/preprints/Unitarizability-Implies-Amenability-for-Countable-Groups-September-23-2026/paper.pdf) is the chosen third paper. Its accepted theorem concerns all uniformly bounded Hilbert representations of a group and similarity by a bounded invertible operator. The paper's triangular cocycle representation and its positive rank-one operator weights retain multiplication and collision data. In particular a sum of Gram contributions is nonnegative even when different labelled paths have the same group product; an arbitrary signed scalar carry has no such property.

Our formal inverse affine generators

```
D(x)=2x,       E(x)=(2x-1)/3
```

give a concrete boundary. Their group lies in an abelian-by-abelian affine group (translations and multiplicative slopes), hence is solvable and amenable. Nevertheless it is not the free group on the two labels. With the convention `[D,E]=D E D^-1 E^-1`, the commutator is translation by `-1/3`; its conjugate by `D` is translation by `-2/3`. These commute. The companion records a nonempty freely reduced word representing the identity. Passing to this formal group has forgotten inverse-branch residue guards and positivity.

Conversely the particular two-dimensional matrix representation with generator `F_(1)(x)=(3x+1)/2` has powers with eigenvalue `(3/2)^r`. It is not uniformly bounded and cannot become unitary under one bounded invertible change of norm. The hypotheses of the accepted amenability theorem are therefore not an available source-height bound. A semigroup contraction and a uniformly bounded group representation have different requirements on inverse powers.

## 7. Prize targets: exact questions and one useful family control

These are read-only snapshots of the public pages on October 7, 2026. No submission, contact or payment has been made. The requested Erdős9 URL ended in `erdos-91`; the catalog's accessible task is the `erdos-9` URL below.

| Task | Pinned mathematical target | Displayed bounty |
|---|---|---:|
| [Erdős9](https://conjectures.io/problems/erdos9-erdos-9) | Positive upper density of odd integers not representable as `p+2^a+2^b`, with `p` prime | $6,637 |
| [Erdős11](https://conjectures.io/problems/erdos11-erdos-11) | Every odd integer `n>1` is `k+2^l` for a squarefree natural number `k` | $6,637 |
| [Erdős213](https://conjectures.io/problems/erdos213-erdos-213) | For every `n>=4`, `n` planar points with integral pairwise distances, no three collinear and no four concyclic | $6,637 |
| [Erdős1085, upper d3](https://conjectures.io/problems/erdos1085-erdos-1085-variants-upper-d3) | The maximum number of unit-distance pairs in dimension three is `O(n^(4/3) log log n)` | $6,637 |
| [Borsuk4](https://conjectures.io/problems/borsuk-borsuk-conjecture-four) | The Borsuk conjecture in dimension four: five strictly smaller-diameter pieces suffice | $6,637 |
| [Erdős23](https://conjectures.io/problems/erdos23-erdos-23) | Every triangle-free graph on `5n` vertices becomes bipartite after at most `n^2` edge deletions | $6,639 |

An infinite family is weaker than positive upper density; a construction for a fixed number of points is weaker than all `n`. A planar unit-distance theorem is not the dimension-three bound. The accepted OpenAI Borsuk construction with lines in `R^4` lives in a nine-dimensional projector space, so it does not settle the four-dimensional task. The targeted repo search found no existing proof of the full Erdős23 target.

There is a clean exact Erdős23 boundary family. Blow up the five vertices of `C_5` to independent parts of sizes `a_0,...,a_4`, with complete bipartite adjacency between consecutive parts. Its minimum bipartizing deletion count is

```
min_i a_i a_(i+1).                                      (7)
```

For a cut, let `x_i` be the number of one colour in part `i`. Its uncut-edge count is affine in each `x_i` with the other coordinates fixed. Successively moving each `x_i` to `0` or `a_i` cannot increase it. Thus an optimum has monochromatic parts. The odd five-cycle forces an uncut edge, and any specified one edge can be the only uncut edge. This proves (7). Balanced parts `a_i=n` attain exactly `n^2`, explaining the target constant; it proves no bound for arbitrary triangle-free graphs. The companion exhausts all 243 vectors with `a_i in {1,2,3}` as an independent finite control.

The [platform FAQ](https://conjectures.io/faq) requires a Lean proof of the pinned statement or its negation, the specified mode, validator acceptance and review. It distinguishes successful solutions from attempts and contributions; the displayed dollar values are not guaranteed fiat payouts. Copying an existing Internet proof and exploiting a kernel defect are excluded. These conditions require rechecking before any future authorized submission. No candidate universal prize proof is produced here.

Equation (7) also gives a constructive certificate for every graph admitting
a homomorphism to `C_5`. For class sizes summing to `5n`,
`min_i a_i a_(i+1) <= (prod_i a_i)^(2/5) <= n^2` by AM--GM when all
parts are nonempty. Deleting that entire adjacency class leaves a subgraph
of a path blowup, hence a bipartite graph. An empty class already leaves
a path. Arbitrary triangle-free graphs need not admit this homomorphism,
so this is a certificate class, not the full prize result or a novelty claim.

For Erdős11, a squarefree-values density theorem over a varying polynomial
parameter does not select one of the finitely many `n-2^l` for each fixed
supplied `n`. This is another instance of the missing source-specific
intersection bound. The platform already lists a reformulation contribution;
we have not generated an independently new full proof or submitted a lemma.

## 8. Reproduction and the next precise obligations

[Script](../../04-computation/experiments/collatz_operator_tensor_bridges_20261007.py), [output](collatz_operator_tensor_bridges_20261007.out):

```
python -B 04-computation/experiments/collatz_operator_tensor_bridges_20261007.py
python -B -O 04-computation/experiments/collatz_operator_tensor_bridges_20261007.py
```

The explicit universe is 341 valuation words of lengths zero through four over `{1,2,3,4}`, 441 ordered tensor compositions, 500 positive plus/minus initial controls, 2,275 shared-label tuples for sizes one through six, 243 complete small blowup tests, and the displayed exact/type hostiles. These test our interfaces, not the cited papers' proofs.

The constructive next step is to keep a marked guarded operator while applying a paper's label-separation or adaptive-precision mechanism. To turn this into a new source floor requires an independently justified positive readout and a proved bound on its total future sensitivity. Tensor dimension, an isomorphism of unmarked algebras, or a shared numerical exponent supplies neither premise by itself.
