# A sign-aware augmentation compiler: parity repair, stabilizer automata, and the root-forgetting defect

Status: PROVED below by elementary finite-group arguments, not independently audited or assigned a canon ID. FINITE-EXACT: the supplied standard-library program independently checks all graph types through order five, all relevant automorphisms, local codecs, and a separate tournament census. This is a new instrument in this research packet, not a claim of literature priority and not an unrooted tournament/sign-even-graph bijection.

## 1. Four recovered mechanisms, and what each actually contributes

All paths in this section are relative to the read-only `math/` repository.

1. **Signed cycle-index compression.** `01-canon/theorems/THM-587-per-level-signed-cycle-index-antipodal-euler-number.md` proves the subset-trace formula
   `P_n(x) = (1/n!) sum_g product_pair_cycles C (1+s_C x^{|C|})`.
   Its coefficients are the tournament quotient's Boolean-Fourier level multiplicities; `P_n(1)` counts tournament types and `P_n(-1)` counts self-converse types. This is a symmetric-subset trace, not `det(I+xM_g)`: exterior reordering signs would give the wrong formula. The reusable mechanism is to retain a character while compressing permutation actions by cycle statistics. I import the proved trace identities, not an unspecified chain complex or homology interpretation.
2. **Augmentation plus exact root/card correction.** `01-canon/theorems/THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer.md`, §§2 and 6, proves the root-preserving Burnside formula and the `1/a(T)` root/card incidence weight. Its important warning is that parent attachment orbits are rooted objects; forgetting the new root merges them. The prime-parent classification is not needed here. The inherited mechanism is local orbit generation together with an explicit root-forgetting fibre.
3. **Finite-state certificates with a retained hidden coordinate.** `01-canon/theorems/THM-4250-odd-core-suffix-lift-automaton-and-density-one-compiler.md`, §3, retains the whole carry vector and persistent success flags, obtains an exact transition automaton, and then a rational generating function for each fixed odd core. Here the replacement hidden coordinate is the *actual remaining automorphism subgroup with its sign character*. The graph is frozen just as that theorem freezes the core. No uniform bounded-state claim across all graph orders is inferred.
4. **The universal sequential-realization warning.** `01-canon/theorems/THM-3163-universal-finite-prefix-markov-realization-and-physical-sidecar-boundary.md` proves every finite terminal law has a posterior-prefix Markov realization. Thus saying merely “there is a sequential algorithm” is empty. The instrument below is more specific: its states are physical stabilizers of actual neighborhood bits; all transition multiplicities are nonnegative integers; the exact completion counts come from a signed cycle index; and local ranking never requires a catalogue of all terminal objects.

Immediate target input: `01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md` already proves the relative-even identity for every permutation subgroup, including the rooted form. We specialize and refine its character principle into an augmentation compiler, prove its exact minimal recurrence and coding threshold, and quantify the root-forgetting defect. These are not presented as a new proof of the underlying equinumerosity theorem.

The closest prior explicit map is `05-knowledge/results/tournament_even_bijection_20261011c.md`: SCC word, Lyndon factorization, and ranked connected atoms. This packet does not repeat that construction. Its canonical hostile is the odd complete graph; its corrected near miss is confusing a relatively-even rooted graph with a globally-even unrooted graph. The least-used sidecar is the sign character on the *current* stabilizer, not just an automorphism order or its ordinary cycle index.

## 2. The signed coloring polynomial

Let `G` be a simple graph on `V=[n]`, `A=Aut(G)`, and

`epsilon(g)=(-1)^{#{ij in E(G): i<j, g(i)>g(j)}}`, for `g in A`.

The character is independent of the reference orientation of the edges. “Even” below means automorphism/sign-even, not Eulerian.

For a coloring `c:V->[q]`, let `A_c` be its color-preserving stabilizer. Call the coloring **sign-separating** if `epsilon|A_c=1`. This asks only that remaining symmetries be even; they need not be trivial. Define `D_G(q)` to be the number of `A`-orbits of sign-separating colorings. If `cyc(g)` is the number of vertex cycles, then

**Theorem 1 (signed color inventory).** For every integer `q>=0`,

`D_G(q) = (1/|A|) sum_(g in A) epsilon(g) q^{cyc(g)}.                 (1)`

The convention `0^0=1` covers the empty graph on zero vertices. The expression is a rational-coefficient polynomial, integer-valued and nonnegative at every nonnegative integer. In its binomial expansion

`D_G(q)=sum_(r=0)^n h_r binom(q,r)`,                                (2)

each `h_r` is a nonnegative integer: it counts `A`-orbits of sign-separating surjective colorings into a fixed `r`-element color set.

**Proof.** On each color orbit, the signed stabilizer average

`(1/|A_c|) sum_(g in A_c) epsilon(g)`

is one exactly when the restricted character is trivial and zero otherwise. Summing over all colorings and dividing by `|A|` therefore counts the required orbits. Interchange the sums: a permutation fixes exactly `q^{cyc(g)}` colorings. This proves (1). Group colorings by their set of used colors to get (2). No cancellation assumption is used. QED.

When `epsilon` is nontrivial, put `A^+=ker(epsilon)`. A useful equivalent interpretation is

`D_G(q) = #orbits(A^+ on [q]^V) - #orbits(A on [q]^V).             (3)`

Indeed, an `A`-orbit splits into two `A^+`-orbits exactly when its stabilizer is inside `A^+`; otherwise it stays one orbit. Thus the positive quantity (1), though written as a signed sum, is an **orbit-splitting defect of the orientation double cover**.

## 3. Individually marked roots turn the polynomial into a coding problem

Adjoin `k` individually distinguished vertices `r_1,...,r_k` to `G`. Fix their induced graph once and for all, for example make them independent. The neighborhoods of the roots in `V` are arbitrary. Each old vertex receives the binary word

`c(v)=(1_(v~r_1),...,1_(v~r_k)) in {0,1}^k`.

We quotient only by automorphisms of the old graph, and require every automorphism fixing each distinguished root to be even. This is **relative/rooted evenness**. It does not require the unmarked resulting graph to be globally even.

Order the old vertices before all roots. A root-fixing automorphism reverses no edge between an old vertex and a root and no edge between roots. Its sign is precisely the restriction of `epsilon_G` to `A_c`. Therefore:

**Theorem 2 (marked-root compiler).** The number `b_G(k)` of relatively-even marked-root extension types is

`b_G(k)=D_G(2^k)=(1/|A|) sum_(g in A) epsilon(g) 2^{k cyc(g)}.       (4)`

If the graph between roots is allowed to vary, rather than frozen, multiply this count by `2^{binom(k,2)}`. The recurrence claim below applies to the frozen-root-graph count, not to that superexponential extra factor.

For nonempty `G`, define `d_G=min{q>=1:D_G(q)>0}`. An injective coloring always works, so `1<=d_G<=n`. Monotonicity of existence under enlargement of the color set gives the exact threshold

`minimum number of added marked roots = ceil(log_2 d_G).           (5)`

This threshold is zero exactly when `G` is already sign-even.

### Sharp hostile: complete graphs are injective binary codes

For `G=K_n`, the edge-reversal sign is the ordinary permutation sign. The stabilizer of a coloring is the product of the symmetric groups on its color classes. It is sign-trivial exactly when every color class has size at most one. Consequently

`D_(K_n)(q)=binom(q,n)`,     `b_(K_n)(k)=binom(2^k,n)`,

and the least number of roots is exactly `ceil(log_2 n)`. Thus there is **no fixed bound on the number of marked roots that can repair arbitrary graph parents** in this model. This is a bounded-root obstruction, not a general computational-memory lower bound.

For the empty graph `I_n`, every permutation is even in the edge-sign character, so

`D_(I_n)(q)=binom(q+n-1,n)`.

The graphs `K_n` and `I_n` have the identical permutation group `S_n` on vertices, with the identical ordinary cycle index. Their root-repair counts differ because the sign character differs. At `n=3`, for `k=0,1,2,3,4`:

- `K_3`: `0,0,4,56,560`;
- `I_3`: `1,4,20,120,816`.

Automorphism order, vertex-orbit partition, and ordinary cycle index alone are insufficient state summaries.

## 4. A nonnegative stabilizer automaton and its exact observable compression

After some root neighborhoods have been fixed, let `H<=A` be their common stabilizer. The next root chooses a subset `S` of old vertices, up to the action of `H`, and the new state is `H_S`. Keep an actual subgroup of `A`; alternatively one may identify states conjugate inside `A`, retaining the vertex action and sign character. Abstract group isomorphism alone is not enough.

Choose one representative of every `H`-orbit on `2^V`. Set

`M_(H,K) = #{chosen S: H_S=K}`,                                  (6)

or sum by the chosen conjugacy class if using conjugacy-reduced states. Accept precisely those states with `epsilon|H=1`. This is a finite nonnegative integer transition system. The trivial group is reachable by giving vertices distinct bit words.

**Theorem 3 (positive transfer and minimal signed observer).** The number of accepted length-`k` paths from state `H` is

`F_H(k)=(1/|H|) sum_(g in H) epsilon(g) 2^{k cyc(g)}`.             (7)

In particular,

`F_H(k+1)=sum_K M_(H,K)F_K(k)`,

`F_H(0)=1_(epsilon|H=1)`,     `b_G(k)=F_A(k)`.

For the initial state write

`a_j=sum_(g in A:cyc(g)=j) epsilon(g)`, `J={j:a_j!=0}`.

Then the generating function and exact minimal rational constant-coefficient recurrence are

`sum_(k>=0)b_G(k)z^k=(1/|A|)sum_(j in J) a_j/(1-2^j z)`,         (8)

`product_(j in J)(E-2^j)b_G=0`,                                  (9)

where `E` shifts `k` to `k+1`. No lower-degree nonzero rational recurrence annihilates the sequence, even eventually. For `n>=1`, the observable recurrence order is at most `n`.

**Proof.** A tuple of root neighborhoods modulo `H` is first classified by the `H`-orbit of its first subset, then by the orbit of the remaining tuple under its stabilizer. This proves the transfer recursion without signed weights. Formula (7) is Theorem 1 with the restricted character and `2^k` colors. Grouping its sum by cycle count gives (8) and (9). Distinct exponentials `(2^j)^k` are linearly independent by a Vandermonde determinant. Any annihilator `p(E)` must satisfy `p(2^j)=0` whenever `a_j!=0`, proving minimality. QED.

The positive automaton and the compressed observer are different objects: the latter may have negative coefficients `a_j`. We do not claim a nonnegative automaton with at most `n` states. This is precisely the quotient-loss boundary that the carry-automaton and posterior-Markov source theorems warn about.

For `K_n`, conjugacy-reduced stabilizer states can be recorded by integer partitions of `n`, namely the sizes of the current binary-code classes. One new bit splits each block of size `s` into sizes `u,s-u`, `0<=u<=s`. Different assignments can reach the same partition and must be retained as transition multiplicities. The accepting partition is `1+...+1`. This is an explicit finite partition-refinement automaton for injective binary codes. Its count nevertheless has the `n`-term exponential compression `binom(2^k,n)`.

## 5. Exact local rank and unrank, without a terminal atlas

Order subset masks increasingly. At state `H` with `r` roots remaining, list the minimal subset representative `S` of every `H`-orbit. Its branch has exactly

`F_(H_S)(r-1)`

accepted completions. Skip zero branches. The branch intervals give the following finite inverse pair.

- **Unrank:** choose the interval containing the input integer; subtract its left endpoint; append its `S`; continue from `H_S`.
- **Rank:** transport the first subset to its minimal `H`-representative, transport all remaining subsets by the same group element, add the preceding branch sizes, then recurse under `H_S`.

Different transporters taking the first subset to the same representative differ by an element of `H_S`; induction makes the final rank independent of that choice. Each input orbit is reached exactly once by the orbit/stabilizer decomposition in Theorem 3. Thus the algorithms are mutually inverse and invariant under the original `A`-action. For an unlabeled input parent, first use a fixed canonical representative of its graph type, carrying its transporter if vertex identities are needed.

This is local enumerative coding on genuinely legal neighborhood moves. It uses ranks, but not a rank pairing of two complete terminal catalogues. It does **not** by itself identify the accepted paths with tournament classes.

## 6. Why even-only hereditary generation and root forgetting fail

### The first attachment-fibre obstruction is already at order two

There is one tournament type on two vertices; its automorphism group is trivial. It therefore has four one-root extension types. There is one sign-even graph type on two vertices, `I_2`, and

`b_(I_2)(1)=binom(3,2)=3`.

The missing fourth relatively-even graph branch comes from the **odd** parent `K_2`, for which `b_(K_2)(1)=1`. Consequently no pairing of the even graph parents with the tournament parents can preserve these one-root extension fibres. Global equality requires contributions from odd graph parents. This does not rule out arbitrary unrooted bijections or constructions with richer intermediate states.

Evenness is also not hereditary: deleting a leaf of the sign-even path `P_3` gives the odd graph `K_2`.

### Exact root-defect identity

For any graph `X` let `A_X=Aut(X)`, `epsilon_X` its sign, and

`r(X)=#orbits(A_X on V(X))`,

`s(X)=#{vertex orbits O: epsilon_X is trivial on the stabilizer of v in O}`.

A relative-even root contributes one precisely on those orbits. For even `X`, `s(X)=r(X)`. For odd `X`, write `A_X^+=ker(epsilon_X)`; then

`s(X)=#orbits(A_X^+ on V(X))-r(X)`.                              (10)

Every contributing vertex orbit in an odd graph has even cardinality: its stabilizer is inside an index-two subgroup and its `A_X`-orbit splits into two equal `A_X^+`-orbits.

**Corollary 4 (root-forgetting defect).** Let `t_n^bullet` count tournament types with a distinguished vertex. Then

`t_n^bullet = sum_(X even) r(X) + sum_(X odd) s(X)`.              (11)

**Proof.** THM-4721's relative-even theorem for the root-fixing subgroup says rooted tournament types equal relatively-even rooted graph types. Partition the latter by their unrooted graph and then by the vertex orbit. Equation (10) is the same index-two orbit-splitting argument as (3). QED.

The full-sign-even rooted graph count is just the first term of (11); omitting the second is incorrect. Exact values at `n=2,3,4,5` are respectively

`2=1+1`, `4=3+1`, `12=7+5`, `48=33+15`.

For actual unrooted even-graph enumeration by attachment, the safe root/card correction is

`e_(n+1) = sum_[G] sum_([S] in Aut(G)\2^V, epsilon_G|Aut(G)_S=1)
             1_(X=G+root(S) is globally even)/r(X)`.             (12)

Every globally-even child supplies exactly `r(X)` rooted classes across all parent types, so its total weight is one. The extra global-evenness test cannot be replaced by the root-stabilizer test. One can check it with full child automorphisms; alternatively, root-moving transporters are the missing sidecar.

## 7. Exact tests and reproduction

Files in this packet:

- `signed_root_augmentation.py`: self-contained standard-library verifier and local codec;
- `signed_root_augmentation.json`: output;
- `signed_root_augmentation_optimized.json`: optimized-Python replay, byte-identical output.

Commands:

`python3 signed_root_augmentation.py > signed_root_augmentation.json`

`python3 -O signed_root_augmentation.py > signed_root_augmentation_optimized.json`

`cmp signed_root_augmentation.json signed_root_augmentation_optimized.json`

The program uses explicit checks, not removable `assert` statements. Exact universe: every unlabeled simple graph at `n=0,...,5` (53 total), every automorphism, and all coloring orbits for `q=1,2,3,4` (212 direct comparisons). Independent sign implementations compare inversion parity with self-paired edge cycles. It constructs 247 reachable actual-subgroup states, checks transfer recursion and recurrence identities, checks every local rank/unrank codeword for two added roots (9,725 objects total), and tests rank invariance under every automorphism. A separate enumeration of tournament orbits checks (11). Total: **1,324,371 checks, PASS**, with normal and optimized outputs identical.

For parents of order `n=0,...,5`, the total numbers of two-root relatively-even extension types are

`1,4,16,88,704,8912`.

For `n=1,...,5`, the histograms of `d_G` across all graph types are

- `n=1`: `d=1: 1`;
- `n=2`: `d=1: 1, d=2: 1`;
- `n=3`: `d=1: 2, d=2: 1, d=3: 1`;
- `n=4`: `d=1: 4, d=2: 5, d=3: 1, d=4: 1`;
- `n=5`: `d=1: 12, d=2: 15, d=3: 5, d=4: 1, d=5: 1`.

These are finite data, not a claimed all-order classification of the histograms. The complete-graph threshold is proved for all `n`; the program also checks its arithmetic boundary through 16.

## 8. What changes for the bijection project

**Recovered instrument:** exact sign-aware canonical augmentation with a physical stabilizer state, positive transition multiplicities, a minimal signed cycle-count observer, and a tested local rank/unrank codec.

**New obstruction:** an even-only one-vertex compiler cannot preserve tournament extension fibres, already at parent order two; and boundedly many added roots cannot repair all odd parent graphs.

**New cross-field bridge:** graph automorphism parity becomes a color-stabilizer separation problem; added roots become binary code coordinates; orbit splitting supplies positivity; finite-state transfer compresses to a signed exponential observer.

**Preserved:** the fixed parent graph, actual neighborhoods of ordered roots, remaining automorphism action, and triviality of the restricted sign character.

**Lost unless retained:** the sign character under ordinary cycle-index compression; positivity under minimal linear compression; root-moving automorphisms under the marked quotient; graph labels under type canonicalization.

**Remaining frontier:** combine this richer signed augmentation category with the tournament side so that the odd-parent branches and root-moving transporters cancel or match constructively. Equations (10)--(12) say exactly what must be accounted for. A global set bijection is not supplied here; counting equality and a local codec are not promoted to one.
