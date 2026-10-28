# Consecutive tournaments encode a square bit matrix—but not its usual quotient

**Status:** self-contained proofs below; finite checks are separately scoped in
`check_encoding.out.json`. These are elementary consequences of a standard
orientation torsor and Burnside enumeration, not claims of literature novelty.

**Scope:** every `n >= 0`; all `2^(n^2)` binary square matrices, with no symmetry,
regularity, switching, transitivity, or extension condition imposed.

## 0. Inheritance and research contract

- Closest proved mechanism: **THM-479 +
  `01-canon/theorems/THM-479-level-holonomy-and-branch-split.md`**. Its setup
  already gives the affine `S_n` action on tournament orientation vectors;
  its antipodal-pair obstruction explains the even-cycle failure below.
  The file read here agrees byte-for-byte with commit
  `03408f84ac43fcdcfb0887e5d6bc1fc8792415cf`
  (SHA-256 `bd591649d996130d4e7f0b0cfaa8575ab0e90ddabb9e8c0e805328ac65cc16e0`).
- Canonical hostile: the transposition of the two vertices of a single edge.
  It fixes equal matrix entries but fixes no orientation of that edge.
- Corrected near miss: `n^2 = binom(n,2)+binom(n+1,2)` gives a labelled
  bijection, not a congruence for whichever quotient one later chooses.
- Least-used sidecar here: the actual orientation torsor, retained after the
  two copies' difference becomes an ordinary graph.
- Applied method card: `00-navigation/META-PATTERNS.md`, “Controlled forgetting
  and unlabeled quotients require a sidecar.”
- Concept board: fixed-frame bijection / group action; root / attachment mask;
  discrepancy graph / orientation sidecar; transpose / core exchange;
  independent component classes / relative alignment.

The anchor is exact lossless encoding; the niche is the affine-action audit;
core exchange and its square-permutation fixed-point law form the wildcard.
No result below settles a major open problem or proves an optimal runtime.

## 1. The exact labelled bijection

Let `V=[n]`, with its usual order, and let `*` be a distinguished new vertex.
For each matrix `M in {0,1}^{V x V}`, and each `i<j`, put

```
x_{ij}=M_{ij},       y_{ij}=M_{ji},       d_i=M_{ii}.
```

Define `T` on `V` and `U` on `V union {*}` by

- `i -> j` in `T` iff `x_{ij}=1`;
- `i -> j` in `U` iff `y_{ij}=1`, for `i<j` in `V`;
- `* -> i` in `U` iff `d_i=1`.

Every unordered pair receives exactly one direction; zero means the reverse,
not “no edge” or a tie. The inverse reads these same three kinds of arcs.
Thus this is a bijection

```
Phi : {0,1}^{n x n} <-> Tourn(V) x Tourn(V union {*}).
```

The cardinality identity is exact:

```
2^(n^2) = 2^binom(n,2) * 2^binom(n+1,2).
```

The second tournament is rooted at `*`, and its other vertices are **aligned**
with the vertices of the first tournament. These are data, not optional
interpretive choices. An equivalent triple is `(T,S,d)`, where `S=U-*` and
`d in {0,1}^V` is the root-attachment mask.

Crucially, consecutive *sizes* do not mean a consecutive *extension*: `T`
and `S` are independent. Imposing `T=S` retains only
`2^(binom(n,2)+n)=2^(n(n+1)/2)` matrices, exactly the symmetric matrices in
this convention. It drops `binom(n,2)` bits. The empty case `n=0` is the
unique empty matrix paired with the empty and singleton tournaments.

## 2. Three relabeling problems are different

### 2.1 Ordinary simultaneous row/column relabeling

For `pi in S_n`, let

```
(R_pi M)_{ab} = M_{pi^{-1}(a), pi^{-1}(b)}.
```

This is the standard conjugation action by permutation matrices. It regards
`M` as a directed binary relation with optional loops, not as two independent
orientations on each unordered pair.

### 2.2 Natural aligned-pair relabeling

Relabel `T`, `S`, and `d` by the same permutation, keeping `*` fixed. This
is the natural isomorphism notion for an aligned rooted pair. Let `alpha_pi`
be its action transported back through `Phi` to the matrix coordinates.

For `E=binom(V,2)`, let `P_pi` permute edge coordinates, and let `Q_pi`
permute vertex coordinates. At the output pair `e={a<b}`, set

```
c_pi(e) = [pi^{-1}(a) > pi^{-1}(b)].
```

All additions below are modulo two. The transported action is

```
alpha_pi(x,y,d) = (P_pi x + c_pi, P_pi y + c_pi, Q_pi d).       (1)
```

The correction is necessary because transporting a directed edge can reverse
its relation to the fixed numerical order. It satisfies the cocycle law

```
c_{pi sigma} = c_pi + P_pi c_sigma.                           (2)
```

Consequently `(1)` really is an affine action. Formula `(2)` also follows by
transporting the reference transitive orientation in two stages.

In contrast, at an output edge whose input order is reversed, `R_pi` sends
its two entries `(x,y)` to `(y,x)`; `alpha_pi` sends them to
`(1-x,1-y)`. If the input order is not reversed, both preserve the two
slots (and merely permute their pair address).

### 2.3 Independent row and column relabeling

The action of `(pi,rho) in S_n x S_n` is

```
M_{ab} -> M_{pi^{-1}(a),rho^{-1}(b)}.
```

Its native carrier is a bipartite graph with distinguished row and column
parts. It can send diagonal entries to off-diagonal positions and does not
preserve the present split into core arcs and root attachments. It is not
natural paired-tournament relabeling either.

For reference, its orbit count is

```
B_n = sum_{lambda,mu partitions of n}
      2^(sum_{i,j} gcd(lambda_i,mu_j)) / (z_lambda z_mu),
```

because one row cycle and one column cycle have `gcd` position orbits.
At `n=2`, this count is `7`, whereas ordinary simultaneous conjugacy gives
`10` and natural aligned-pair isomorphism gives `8`.

## 3. No alternative bijection can repair the ordinary equivariance

**Proposition.** For `n>=2`, the `S_n`-set of ordinary binary matrices and
the `S_n`-set of natural aligned rooted tournament pairs are not isomorphic.

**Proof.** A transposition fixes the zero matrix under ordinary conjugation.
It cannot fix any tournament on `V`, since it reverses the edge between its
two exchanged vertices. Thus it fixes no aligned pair. Equivariant bijections
preserve fixed-point counts for each group element. Contradiction. More
locally, exchanging the two off-diagonal matrix slots fixes `00` and `11`,
while complementing two tournament orientation bits fixes no pattern. QED.

The same obstruction shows that the orientation cocycle cannot be removed
by translating the origin in the fixed-frame edge vector space: such a
translation would provide an `S_n`-invariant tournament.

This is stronger than the failure of our particular encoding. A common
cardinality is insufficient even if one is allowed to choose a different
bijection.

### A faithful repair with an explicit frame

Include a total order `<` as part of both objects, and define the encoding
relative to that order. Then

```
(M,<) <-> (T,U,*) with the same order on the nonroot vertices
```

is equivariant when the order is transported together with the vertices.
An inverted pair relative to the old numerical order is not inverted relative
to its transported frame. This restores a genuine equivariant *bundle*.
It does not remove the frame: it enlarges both labelled sets by `n!`, and
forgetting it reintroduces precisely the problem above.

## 4. The discrepancy graph cancels the cocycle, but not the orientation

Put `r=x+y` and retain `t=x`. Then `r` is the graph of edges on which `T`
and `S` disagree. This is another exact coordinate system since `y=t+r`.
For natural pair relabeling,

```
(t,r,d) -> (P_pi t+c_pi, P_pi r, Q_pi d).                      (3)
```

The discrepancy graph is ordinary linear data; the first tournament remains
an affine orientation torsor. For ordinary matrix relabeling the exact rule is

```
(t,r,d) -> (P_pi t+c_pi*(P_pi r), P_pi r, Q_pi d),             (4)
```

where `*` is coordinatewise multiplication. Thus even though `r` transforms
in the same way in both problems, it is not enough to decide which action
acts on the retained orientation.

Two sharp boundaries explain the distinction:

1. `r=0`: the matrix is symmetric and the core tournaments coincide.
   Ordinary relabeling is the graph action on `t`; natural relabeling still
   has the tournament cocycle. This is **not** an equivariant locus.
2. `r=1` on every edge: `M_{ij}+M_{ji}=1` for `i!=j`, and the two core
   tournaments are converses. Equations `(3)` and `(4)` now agree. Ordinary
   tournament adjacency matrices, with arbitrary diagonal bits, form a genuine
   equivariant subfamily of size `2^(n(n+1)/2)`. Indeed this is the exact
   agreement locus: if `alpha_pi(M)=R_pi(M)` for every `pi`, exchange the
   endpoints of any selected edge; equations `(3)` and `(4)` force `r=1`
   on that edge. Since the edge was arbitrary, all opposite entries must
   be complementary.

Neither operation discards a bit if `(t,r,d)` is retained. Keeping only
`r,d` loses exactly `binom(n,2)` labelled bits.

## 5. Exact orbit enumeration: the missing even partitions

For a partition `lambda=(lambda_1,...,lambda_k)` of `n`, put

```
ell(lambda)=k,
z_lambda = product_a a^(m_a) m_a!,
D(lambda) = sum_{i,j} gcd(lambda_i,lambda_j).
```

For an all-odd partition also put

```
e(lambda) = sum_i (lambda_i-1)/2 + sum_{i<j} gcd(lambda_i,lambda_j).
```

Then `D(lambda)=2e(lambda)+ell(lambda)` for all-odd `lambda`.

**Theorem.** The number `P_n` of natural aligned rooted pair classes is

```
P_n = sum_{lambda partitions n, all parts odd}
      2^(2e(lambda)+ell(lambda)) / z_lambda
    = sum_{lambda partitions n, all parts odd} 2^D(lambda)/z_lambda.    (5)
```

The number `C_n` of ordinary matrix conjugacy classes is

```
C_n = sum_{lambda partitions n} 2^D(lambda)/z_lambda.                  (6)
```

**Proof.** If a permutation has an even cycle of length `2h`, its `h`-th
power exchanges a pair of antipodal vertices in that cycle. It cannot fix a
tournament, so it fixes no pair. If every cycle is odd, it has `e(lambda)`
edge orbits with no reversed return. Each core tournament has one free
orientation per edge orbit. Root attachments have one bit per vertex cycle.
The fixed-pair count is therefore `2^(2e+ell)`.

For a matrix, all ordered position orbits may be assigned a bit, including
loops. A cycle of length `a` and a cycle of length `b` contribute `gcd(a,b)`
ordered-pair orbits, giving `2^D` for every partition. There are `n!/z_lambda`
permutations of a given cycle type. Burnside's lemma proves both formulas.
QED.

The classical single-tournament fixed-point calculation is the lemma and
Theorem 41 in Moon, *Topics on Tournaments*, Section 29, printed pages
113–114; we use it for two independent core tournaments and a vertex mask.
[Primary text](https://www.gutenberg.org/files/42833/42833-pdf.pdf).

The exact positive gap is

```
C_n-P_n = sum_{lambda having an even part} 2^D(lambda)/z_lambda > 0
```

for every `n>=2`. This is not a rounding artifact or a missing root factor.

| n | natural aligned pairs P_n | ordinary matrix conjugacy C_n |
|---:|---:|---:|
| 0 | 1 | 1 |
| 1 | 2 | 2 |
| 2 | 8 | 10 |
| 3 | 88 | 104 |
| 4 | 2752 | 3044 |
| 5 | 279968 | 291968 |
| 6 | 95458560 | 96928992 |
| 7 | 111698291584 | 112282908928 |

The `n<=4` rows have independent exhaustive orbit checks. Larger displayed
rows are exact partition computations, with per-cycle-type Gaussian
fixed-point checks through `n=12`; they are not exhaustive object censuses.

For `n>=3`, the exact contribution in `(5)` from a single 3-cycle is
`2^(n^2-4n+6)/(3(n-3)!)`. Ordinary matrices instead also have the
transposition contribution `2^(n^2-2n+2)/(2(n-2)!)`. These are individual
terms, not a claimed complete asymptotic expansion or runtime bound.

## 6. Transposition is core exchange, not arc reversal or literal pair swap

Matrix transposition acts by

```
tau(x,y,d)=(y,x,d).
```

On `(T,U)`, this means: take `S=U-*` as the first tournament and extend the
old `T` by the unchanged root mask as the new second tournament. The two
full tournaments have different sizes, so a literal factor swap is not an
endomorphism of this carrier. Core exchange is the precisely typed swap.

It commutes with natural relabeling and ordinary matrix conjugation. Its
fixed labelled objects are `x=y`, so there are `2^(n(n+1)/2)` of them.
It is not tournament converse, root switching, or relabeling a vertex.
Quotienting by it requires a separate choice of equivalence relation.

**Square-permutation law.** For a permutation `pi`, write `lambda^(2)`
for the cycle type of `pi^2`: an odd part stays unchanged, and an even part
`2a` becomes two parts `a`. Then

```
|Fix(alpha_pi tau)| =
  0,                                      if some lambda_i is divisible by 4;
  2^(e(lambda^(2))+ell(lambda)),           otherwise.               (7)
```

**Proof.** The fixed equations are `T=pi S`, `S=pi T`, and `d=pi d`.
Choose a tournament fixed by `pi^2`; then the other core is determined,
and choose a mask fixed by `pi`. A tournament fixed by `pi^2` exists iff
all cycles of `pi^2` are odd, equivalently no original cycle is divisible
by four. The counts in `(7)` follow. QED.

Consequently the number of classes after also allowing core exchange is

```
(P_n + sum_{lambda: no part divisible by 4}
       2^(e(lambda^(2))+ell(lambda))/z_lambda) / 2.                 (8)
```

It gives `1,2,7,58,1480` for `n=0,...,4`, all independently enumerated.
The analogous ordinary-matrix-plus-transpose quotient is different: at
`n=2` it has `9` classes, not `7`.

## 7. What an isomorphism class and a root do—and do not—encode

For an unlabelled tournament representative `A` on `n` vertices, its
labelled realizations are the cosets `S_n/Aut(A)` and number
`n!/|Aut(A)|`. Choosing a single root gives rooted isomorphism types
`V(A)/Aut(A)`, not all those labelings. For the transitive tournament,
`Aut(A)` is trivial: there are `n!` labelings but only `n` possible roots.
Already at `n=3`, class-plus-one-node cannot specify all labelled versions.
A root is also not a license to reverse the tournament's edges.

For a fixed first representative `A`, the correct residual enumeration of
aligned pairs is

```
(Tourn([n]) x {0,1}^[n]) / Aut(A).
```

Its exact count is

```
(1/|Aut(A)|) sum_{g in Aut(A)} 2^(e(g)+ell(g)),
```

because every automorphism of a tournament has odd order. This follows
from the same even-cycle obstruction and the elementary fact that a finite
group of even order contains an involution. Summing over representatives
`A` gives `(5)` again. This is a safe bridge to stabilizer and double-coset
algorithms: it does not require exploring all `n!` relabelings of the first
component repeatedly.

If the two component isomorphism classes are stored independently, their
relative alignment is lost. Restoring it requires explicit identifications,
or an appropriate double coset under their automorphism groups, followed by
the residual orbit of the root mask. If the root itself is forgotten,
vertex orbits of the larger tournament must also be accounted for.

A fixed-frame labelled encoding contains exactly `n^2` independent bits.
No worst-case lossless fixed-length encoding of every matrix uses fewer
than `n^2` bits. Class-based enumeration can reduce repeated work when the
target is explicitly an isomorphism quotient; it cannot silently discard
labeling, alignment, or attachment data.

## 8. Sources, status, and reproducibility

The affine orientation mechanism is inherited, not introduced here:

- **THM-479 + the full slug above**, its setup and pair-holonomy proof.
- Babai–Cameron, *Automorphisms and enumeration of switching classes of
  tournaments*, EJC 7 (2000), R38. Lemma 3.1 and Section 7 distinguish odd-order
  tournament automorphisms from switching-class symmetries. That paper treats
  switching classes, not the present full paired-matrix universe.
  [Primary paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v7i1r38/pdf/).
- Moon's primary text, Section 29, cited above, supplies the classical
  single-tournament counting background. The additional-product formulas,
  action comparison, and transpose law here are proved directly.

The parent session owns repository publication, integration, and global
ledgers. This lane made no Git writes. The bounded startup ran read-only;
processor message lookup reported no `.machine-id`, so no registration or
shared message state was changed.

Run from this directory:

```bash
python3 check_encoding.py > check_encoding.out.json
```

The standard-library checker contains three distinct verification paths:

1. actual adjacency-matrix encoding and direct vertex relabeling;
2. the explicit affine formula and Gaussian fixed-point equations;
3. partition formulas, independently compared with exhaustive orbit sets.

It verifies all labelled round trips through `n=4`; all objects and all
permutations in direct-vs-affine comparison through `n=3`, plus every affine
basis vector and permutation through `n=4`; all composition laws through
`n=4` on an affine basis; all 272 partitions from `n=0` through `12`, for
natural, ordinary, and natural-plus-transpose fixed counts; all orbit sets
through `n=4`; and the full complementary-core equivariance locus through
`n=4`. It uses no sampling, external packages, or unrecorded data. All controls use explicit
`need(...)` calls: normal and Python-optimized (`-O`) runs retain the same
active checks and must produce identical receipts.

**Open scope:** benchmark practical relative-alignment enumeration at larger
orders; determine which quotient the user's application actually needs;
retain the frame if ordinary matrix isomorphism is the target. None of the
counts alone proves an optimal enumeration algorithm or mathematical novelty.

## 9. Constructive continuation

`ODD_SELF_CONVERSE_LIFT.md` supplies a class bijection to odd self-converse
tournaments, upgrading the inherited Session22 numerical channel identity.
It retains an anti-involution/centralizer sidecar and does not identify full
automorphism groups. That continuation also gives an explicit tournament-
isomorphism reduction, without claiming a new complexity equivalence.
