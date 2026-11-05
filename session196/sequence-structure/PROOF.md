# Certified variance expansions and deterministic planted-spine witnesses

Session 196, 2026-11-05. **INHERITED exact moment/even-run identities; PROVED uniform remainder bounds, higher-order expansion, and constructive consequence; FINITE-EXACT verification is separate.** The unrestricted raw Hamiltonian-path spectrum with exceptions `{7,21}` remains **OPEN**. No literature-priority claim is made.

## 1. Definition and inheritance check

[OEIS A000568](https://oeis.org/A000568), checked directly, counts ordinary isomorphism classes of tournaments on n vertices, with offset zero. Its initial terms are `1,1,1,2,4,12,56,456,6880`. It does not count possible Hamiltonian-path values.

The [Session 195 guide](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/00-navigation/S195.md) routes the current constructive interfaces and explicitly leaves full coverage open. The inherited [fixed-count packet](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/05-knowledge/results/spectrum_basins_20261020e/spectrum.md) already proves fixed-H rationality, logarithmic strong cores, exact eventual-polynomial degrees and leading terms. Those results are **not** presented as new here. Its [later inventory](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md) supplies further finite strong-count data.

The exact second-moment formula, even-run cancellation, and `2/n` concentration rate were already proved in [THM-589, even-run slug](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/01-canon/theorems/THM-589-cvH-even-run-closed-form-and-2-over-n-rate.md). The odd-composition/formal-series identity is also inherited from [THM-589, variance slug](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/01-canon/theorems/THM-589-metagraph-H-variance-closed-form-nud-succession.md), which credits the owner's form. The [2026-10-14 run-kernel atlas](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/05-knowledge/results/oeis_transfer_atlas_20261014h.md) already independently proves a general formal-series compiler and a seven-step recurrence for W_n. None of those identities, the sequence, or the concentration rate is a new S196 result.

Relative to these inspected sources, this packet adds an explicit uniform remainder bound supporting expansion to every inverse-power order, higher coefficients (including the vanishing `n^-2` coefficient), and a planted-spine conditional-expectation algorithm producing actual high-count graphs. The orbit-mass formulas are elementary corollaries. The inherited proof is included to make the extensions and their implementation independently checkable; this comparison is not a repository-wide or literature-wide novelty certificate.

A terminology control matters. The “even graphs” in [Royle, Praeger, Glasby, Freedman and Devillers, arXiv:2204.01947](https://arxiv.org/abs/2204.01947), linked from A000568, are defined by the parity of edge reversals under automorphisms. They are **not** simply graphs with all degrees even. No Eulerian-graph interpretation of that equinumerosity is used below.

## 2. Inherited exact pair theorem, restated with proof

Let T be uniform among the `2^(n(n-1)/2)` labelled tournaments on `[n]`, n>=1. Let H(T) count directed Hamiltonian vertex words, and put

    mu_n = n!/2^(n-1).

Define W_0=1 and, for n>=1,

    W_n = n! + sum_(1<=j<=floor((n-1)/2))
                 sum_(1<=r<=min(j,n-2j))
                 2^r binom(j-1,r-1) binom(n-2j,r) (n-2j)!.       (1)

**Theorem 1 (INHERITED: the two THM-589 sources above).** For every n>=1,

    E[H] = mu_n,
    E[H^2] = n! W_n / 4^(n-1),
    Var(H)/mu_n^2 = W_n/n! - 1.                                (2)

Equivalently, as a formal ordinary generating series,

    sum_(n>=0) W_n z^n
      = sum_(m>=0) m! [z(1+z^2)/(1-z^2)]^m.                   (3)

The initial W_n, n=0,...,9, are

    1, 1, 2, 8, 32, 158, 928, 6350, 49752, 439670.

**Proof.** Each of the n! possible vertex words requires n-1 independent arc orientations, giving the first moment.

For the second moment, fix one word as `1,2,...,n`, and let pi be the second permutation. Let f(pi) count its consecutive arcs `i -> i+1`, and let b(pi) count its consecutive arcs `i+1 -> i`. If b(pi)>0 the two words are incompatible in a tournament. Otherwise the union of their required arcs has `2n-2-f(pi)` distinct arcs. Consequently

    E[H^2] = n! 2^(-2n+2) W_n,
    W_n = sum_(pi in S_n, b(pi)=0) 2^f(pi).                   (4)

It remains to evaluate this permutation sum without cancellation. For each spine edge `{i,i+1}`, let F_i and B_i be the indicators that pi uses it in the forward and backward direction. Since F_i B_i=0,

    1_(b=0) 2^f = product_(i=1,...,n-1) (1+F_i-B_i).          (5)

Expand (5). A selected subset of the spine edges is a disjoint union of runs. A summand can occur in a permutation only when the direction is constant along each run: a direction change at a shared vertex would demand two predecessors or two successors. A run of length ell can therefore be all forward, of weight +1, or all backward, of weight (-1)^ell. If s edges in all are selected, each compatible orientation of the runs contracts to n-s ordered blocks, giving exactly `(n-s)!` permutations. The combined weight of a run is `1+(-1)^ell`. Thus odd-length runs cancel and every even-length run contributes 2.

Write s=2j and let r be the number of runs. Their positive half-lengths form a composition of j, giving `binom(j-1,r-1)` choices. Placement of the r runs among n-1 edges, with a nonempty gap between successive runs, gives `binom(n-2j,r)` choices. Multiplying by `2^r (n-2j)!` proves (1). The no-edge subset contributes n!.

For (3), a contracted block is either a singleton, of weight z, or an odd block of size at least three, of weight `2z^3/(1-z^2)`. An ordered sequence of m such blocks, multiplied by m!, gives precisely the preceding expansion. This is an identity of formal series; no analytic convergence of the factorial series is asserted. QED.

## 3. The refinement of A000568

Let `T_n` be the set of ordinary isomorphism classes counted by A000568(n). Write g(T)=|Aut(T)|. Every automorphism fixing a Hamiltonian word fixes every vertex, so its action on words, and on ordered pairs of words, is free. In particular, g(T) divides H(T).

**Corollary 2 (elementary orbit-mass restatement).** For every n>=1,

    sum_(T in T_n) H(T)/g(T)   = 2^((n-1)(n-2)/2),
    sum_(T in T_n) H(T)^2/g(T) = 2^((n-1)(n-4)/2) W_n.       (6)

For n=2,3 the power of two in the second identity has a negative exponent; the product displayed is still an integer. Equivalently the identity can be stated without negative exponents as

    4^(n-1) sum_T H(T)^2/g(T) = 2^(n(n-1)/2) W_n.

**Proof.** A class T has n!/g(T) labelled realizations. Divide the sums of H and H^2 over all labelled tournaments by n!, then apply Theorem 1. The free-action observation also identifies the sums as the numbers of unlabelled tournaments carrying a distinguished word or an ordered pair of words. QED.

Thus A000568 is the unpointed class count, while (6) gives its exact first and second Hamiltonian-pointed masses. Neither H nor H/g should be substituted for the other in the raw `{7,21}` problem.

## 4. Controlled expansion to every asymptotic order

Put R_n=`Var(H)/mu_n^2`. A summand of R_n in (1) is

    t_(j,r)(n) = 2^r binom(j-1,r-1) binom(n-2j,r)/(n)_(2j),   (7)

where `(n)_k=n(n-1)...(n-k+1)`. All terms are nonnegative.

**Lemma 3 (uniform tail).** For n>=12 and integer J>=0,

    sum_(j>J) sum_r t_(j,r)(n)
      <= (2/3) (6/n)^(J+1)/(1-6/n),                         (8)

where only the valid indices of (1) are included.

**Proof.** Let m=n-2j. The last r factors in `(n)_(2j)` are m+1,...,m+r, so

    binom(m,r)/(n)_(2j)
      <= 1/[r! (n)_(2j-r)]
      <= 1/[r! (n)_j],

because r<=j and all omitted positive integer factors are at least one. Therefore

    sum_r t_(j,r)(n)
      <= [sum_(r=1)^j 2^r binom(j-1,r-1)]/(n)_j
      = 2*3^(j-1)/(n)_j
      <= (2/3)(6/n)^j.

The last inequality follows from j<=floor((n-1)/2), so every factor of `(n)_j` is at least n/2. Summing the geometric tail proves (8). QED.

**Theorem 4 (full expansion).** For each fixed nonnegative integer J, an asymptotic expansion of R_n through n^(-J) is obtained by expanding the finitely many rational terms (7) with j<=J. The error is O(n^(-J-1)). In particular,

    R_n = 2/n - 14/(3n^3) - 58/(3n^4) - 916/(15n^5)
                 - 6986/(45n^6) + O(n^-7).                 (9)

The coefficient of n^-2 is exactly zero. The first elementary estimate is

    2(n-2)/(n(n-1)) <= R_n
      <= 2(n-2)/(n(n-1)) + 48/n^2    (n>=12).             (10)

**Proof.** Formula (8) makes the discarded j>J terms O(n^(-J-1)). Each remaining rational term has a convergent expansion in 1/n near zero. Its coefficients are extracted from

    [2^r binom(j-1,r-1)/r!] x^(2j-r)
      * product_(i=0)^(r-1) [1-(2j+i)x]
      / product_(i=0)^(2j-1) [1-ix].                       (11)

Adding the finitely many coefficients proves the general claim and gives (9). The verifier performs exact rational coefficient extraction. The j=1 term is `2(n-2)/(n(n-1))`; retaining it and using (8) with J=1 gives (10). QED.

Chebyshev now gives, for every epsilon>0,

    Pr(|H/mu_n-1|>=epsilon) <= R_n/epsilon^2 -> 0.          (12)

This concerns a uniform **labelled** random tournament, not a uniform member of the unlabelled A000568 classes. No transfer between those probability models is assumed.

## 5. An all-order constructive consequence

**Theorem 5 (planted-spine conditional-expectation algorithm).** For every n>=1 there is a deterministic algorithm constructing a tournament T on n vertices with

    H(T) >= W_n / 2^(n-1)
          = mu_n [1+2/n+O(n^-3)].                         (13)

One implementation uses O(n^4 2^n) exact rational arithmetic operations and O(n 2^n) storage. This is an arithmetic-operation bound, not a claim of unit-cost bit complexity.

**Proof.** Fix `1->2->...->n` and leave every other pair independent and uniformly oriented. The mean H in this planted-spine distribution is

    E[H | the standard spine is a path]
      = sum_pi 1_(b(pi)=0) 2^(-n+1+f(pi))
      = W_n/2^(n-1).                                    (14)

At any partially fixed orientation, give a required directed edge weight 1 if fixed present, 0 if fixed absent, and 1/2 if unfixed. The weighted Hamiltonian sum is exactly the conditional expectation of H, because a word uses any unordered pair at most once. Compute it by the usual subset-endpoint dynamic program

    D({v},v)=1,
    D(S,v)=sum_(u in S-v) D(S-v,u) w(u,v).

For an unfixed pair, the current mean is the average of the two means after fixing its direction. Choose a direction with no smaller mean, with a specified deterministic tie rule. Continue through all pairs. The final mean is the actual integer H(T), never smaller than (14). There are O(n^2) choices and each weighted DP takes O(n^2 2^n) arithmetic operations, proving the bounds. QED.

The verifier uses integer weights 0,1,2 in place of 0,1/2,1, so its weighted sum is `2^(n-1)` times the expectation. It records every decision and confirms the mean is nondecreasing, then recomputes final H by ordinary subset DP and, at small orders, literal permutation enumeration.

This supplies actual graphs at every order and realizes the inherited second-moment improvement over the elementary first-moment bound by a simple explicit deterministic algorithm. It is **not claimed to be a best-known extremal bound or a new asymptotic for the maximum H**. It does **not** tune H to an arbitrary prescribed odd integer; it is a high-count construction, not a solution of the support problem. This algorithm and its stated guarantee were not present in the three compared moment/run-kernel sources; no broader priority claim follows from that limited comparison.

For a verified scope comparison, the inspected primary abstracts of [Wormald, arXiv:2610.16245](https://arxiv.org/abs/2610.16245) and [Divoux and Wang, arXiv:2610.23755](https://arxiv.org/abs/2610.23755) state substantially stronger asymptotic extremal results, including Hamiltonian paths. Their proofs were not audited here and are not dependencies of this packet. The point of Theorem 5 is its simple explicit weighted-DP construction, not a claim to improve those extremal results.

## 6. Scope and reproducibility

Run `python3 verify.py` and `python3 -O verify.py` from this directory. The output receipt records:

- exact enumeration of all labelled tournaments through order six, including the H and H^2 moments;
- an independent signed-adjacency permutation evaluation of W_n through order nine;
- a separate formal-series evaluation of (3), exact through order thirty;
- orbit-weighted moment checks and A000568 class-count controls through order six, using actual relabelling orbits;
- exact asymptotic coefficients from rational formal series, plus finite checks of the uniform tail bound;
- actual deterministic witness graphs and conditional-expectation traces through order ten.

No removable `assert` statements are used. Exhaustive checks validate implementations on their explicitly finite universes; the all-order statements follow from the proofs above. The inherited all-order omissions 7 and 21 are unchanged, and no other value is declared absent or universally attained.

For background, the first-moment formula is classical and appears in [Moon, *Topics on Tournaments*, Section 10](https://www.gutenberg.org/files/42833/42833-pdf.pdf). The present packet makes no claim that the pair formula or its consequences are absent from the literature.
