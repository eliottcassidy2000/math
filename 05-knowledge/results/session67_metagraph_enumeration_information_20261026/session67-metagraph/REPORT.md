# Complete weighted families: recovered closure, exact core-size census, and its limit

**Status:** inherited PROVED classification/counting; derived elementary refinement
and bounds; independent FINITE-EXACT checks. No new closed formula for A000568,
novelty claim, polynomial-time enumeration claim, or full metagraph compression.

## Useful result first

The repository already **does** contain an all-order, complete weighted-family
recurrence. It was recovered from sparse Git objects, not merely rediscovered:

- `weighted_core_transport_20261020/core_count.md`, §§1, 6–7: unique weighted
  pair-free cores, full automorphism actions, alternating augmentation, and the
  exact root-forgetting symmetry weight.
- `defect_boundary_20261020/boundary.md`, §§3–4: signed ordered-list inversion
  and the entire odd-cycle operator system, including its fixed-point-free
  boundary.
- `tournament_future_boundary_20261020/README.md`: the continuation-response
  quotient computes ordinary counts through 200 and core counts through 100,
  without constructing individual cores. This remains Burnside arithmetic.

The executable contribution here is a **self-contained independent census of
ordinary classes jointly by total order and canonical core order**. It verifies
actual core families, their automorphism cycle indices, and duplicate freedom
through order seven, and computes the full integer triangle through order forty.
The new quantitative consequence is an exact all-order core-deficit bound and a
uniform large-core asymptotic. A finite family bank is not where almost all
unlabelled classes live: almost every class has no clone contraction at all.

All inherited source bytes are pinned in `source_manifest.json` to
`b4234f728c0516b5fc521c001d69e26f3d2508d5`. Exact local source snapshots are included.
Use the source path plus this commit, not a historical theorem number alone.

## 1. Inheritance, hostiles, and what counts as a family

The closest proved mechanism is THM-4463,
`01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md`,
§2: the graph joining homogeneous vertex pairs is a disjoint union of consistently
oriented paths. Its connected components are exactly the maximal transitive
modules. A homogeneous pair has identical incidence to every *other* vertex;
its two full adjacency rows are not equal.

The canonical hostile is the transitive triple: pairs 01 and 12 are homogeneous,
but 02 is not. Pair homogeneity is **not** an equivalence relation. Its connected
components, rather than its pair relation, are the canonical chains.

The corrected near miss is arbitrary weighted quotients. For the transitive
two-vertex quotient, weights (1,2) and (2,1) both give TT3, despite the quotient
having trivial automorphism group. Quotient twins must be merged first.

The least-used relevant sidecar is the **permutation action** of the core
automorphism group. Group order alone loses weight-orbit counts. Two inherited
order-six pair-free cores both have automorphism group of order three:

```
Q_A = (0,33,11,19,7,29),  nonidentity cycle type (1,1,1,3)
Q_B = (10,20,33,38,13,19), nonidentity cycle type (3,3).
```

Their numbers of positive weight orbits of total seven are respectively 4 and 2.
The code rechecks both examples directly.

The board is **binary relation / canonical chains / core automorphism action /
positive weights / pointed path allocation / core-size distribution**. The
Anchor is exact ordinary class enumeration; the Niche is its core-deficit
triangle; the Wildcard is uniform asymptotic exhaustion by the full-core layer.

## 2. Recovered exact, duplicate-free family decomposition

Call a tournament **pair-free** if it has no homogeneous pair. This is not
modular indecomposability, strong connectivity, or gcd-one inflation scale.
Let C_q contain one representative of every pair-free class on q vertices.
For Q in C_q and a positive vector w, replace vertex i by an ordered transitive
chain of length w_i, with all cross-chain arcs directed by Q.

Every nonempty tournament is uniquely Q[w], up to a core isomorphism carrying
weights. For fixed Q the classes are exactly positive weight vectors modulo
G=Aut(Q), and

```
Aut(Q[w]) = Stab_G(w).
```

Indeed every isomorphism preserves canonical maximal chains, hence induces the
stated weighted core isomorphism. Each chain's order has exactly one possible
order-preserving map. Conversely, a weighted core isomorphism lifts uniquely.
Inflating a pair-free quotient creates no new homogeneous pair across chains;
a cross-chain pair would project to a homogeneous core pair. This establishes
both coverage and duplicate freedom, including all mergers.

For g in G with vertex cycle lengths d_1,...,d_s, a fixed weight is constant
on each cycle. Thus the family's exact ordinary generating function is

```
F_Q(x) = x^q/|G| * sum_(g in G) product_(cycles c of g) 1/(1-x^|c|).
```

This is local Burnside on Aut(Q), not a new Burnside-free formula. The complete
positive family identity is

```
a_n = sum_(q=1..n) sum_(Q in C_q) [x^n] F_Q(x),   a_0=1.
```

The recovered alternating augmentation closes the bank C_q itself: delete one
vertex from a pair-free child, contract the parent's canonical chains, then
restore the child using bits that alternate along each parent chain. There are
2^q phase candidates on a q-chain parent; reject any new homogeneous pair and
identify isomorphic children. This is an exact well-founded **set recurrence**,
not a recurrence on the scalar a_n and not output-polynomial compression.

## 3. Exact refinement by core order

Let A[n,q] count ordinary n-vertex classes whose canonical core has q vertices.
Let c_lambda be the number of pair-free labelled tournaments fixed by one
permutation with odd cycle type lambda, and z_lambda its centralizer order.
The inherited cycle-index inversion supplies c_lambda exactly. Then

```
A[n,q] = sum_(odd lambda partition q) c_lambda/z_lambda
             * [x^(n-q)] product_(d in lambda) 1/(1-x^d).                 (1)
```

All terms in (1) are nonnegative rational numbers; the total is an integer.
The coefficient automatically retains repeated cycle lengths. It is the sum
of the complete positive family functions over cores of size q, not an
assertion that those individual cores have identical structure.

For an equivalent all-order bivariate form, let Z_T and Z_C be the full cycle
indices of tournaments and pair-free cores, including empty constant one.
The recovered identities are

```
Z_T(p) = Z_C(p_d/(1-p_d)),   Z_C(p) = Z_T(p_d/(1+p_d)).
```

A core d-cycle contributes d to the number of core vertices, so marking that
number by u gives

```
A(x,u) = sum_(n,q>=0) A[n,q] x^n u^q
       = Z_C(u^d x^d/(1-x^d) : d>=1)
       = Z_T(u^d x^d/(1-x^d+u^d x^d) : d>=1).                            (2)
```

These are formal series. Setting u=1 recovers A000568. The diagonal q=n
recovers pair-free cores. No analytic convergence is assumed.

For explicit coefficient calculation, write m_d for multiplicities and

```
L(m,k)=m!/k! * binom(m-1,k-1),  1<=k<=m; L(0,0)=1.
e(lambda)=sum_i (lambda_i-1)/2 + sum_(i<j) gcd(lambda_i,lambda_j).
c_m = sum_(k<=m) 2^e(k) product_d (-d)^(m_d-k_d) L(m_d,k_d).               (3)
```

All other L values are zero. Equation (3) is inherited signed ordered-list
inversion, not a newly claimed transform. Formula (1) followed by summation
and formula (2) expanded directly are independently implemented here.

Selected rows (unlisted core sizes have zero count):

| n | Contributions q: A[n,q] | a_n |
|---:|---|---:|
| 5 | 1:1, 3:2, 4:4, 5:5 | 12 |
| 6 | 1:1, 3:4, 4:8, 5:17, 6:26 | 56 |
| 7 | 1:1, 3:5, 4:16, 5:47, 6:132, 7:255 | 456 |
| 8 | 1:1, 3:7, 4:26, 5:107, 6:444, 7:1693, 8:4602 | 6880 |
| 9 | 1:1, 3:10, 4:40, 5:210, 6:1168, 7:6664, 8:36348, 9:147095 | 191536 |
| 10 | 1:1, 3:12, 4:60, 5:374, 6:2600, 7:19860, 8:162648, 9:1319141, 10:8228360 | 9733056 |

## 4. A quantitative all-order limit on family compression

**Proposition (core-deficit tail).** For 1<=d<n,

```
sum_(q<=n-d) A[n,q] <= binom(n-1,d) a_(n-d).                             (4)
```

**Proof.** Every tournament on the left has at least d surplus chain vertices.
Remove d of them while leaving each chain nonempty, obtaining an (n-d)-vertex
parent. The child can be reconstructed by inflating the parent's individual
vertices into positive chains whose lengths sum to n. There are exactly
binom(n-1,d) such compositions for a fixed representative parent. This is a
surjection, not an injection; its duplications only improve the upper bound.
Conversely these inflations have core order at most n-d. QED.

In particular, for r_n=A[n,n],

```
a_n-r_n <= (n-1) a_(n-1).                                                (5)
```

A fixed finite bank of cores of order at most b contributes at most
sum_(q<=b) r_q binom(n-1,q-1)=O_b(n^(b-1)). But a_n>=2^binom(n,2)/n!, so
that bank's proportion tends to zero. This does not exclude an unbounded
symbolic grammar with a small description, such as the tribonacci family.

For completeness, an elementary quantitative Burnside estimate proves the
stronger statement about **ordinary unlabelled classes**, without importing
an unverified random-tournament heuristic. Put L_n=2^binom(n,2)/n!. Then

```
a_n=L_n(1+epsilon_n),
0<=epsilon_n <= E_n := sum_(m=3..n) (n!/(n-m)!)
                               * 2^(-ceil(m(2n-m-1)/3)).                (6)
```

Only all-odd permutations fix tournaments. A nonidentity such permutation
moves m>=3 vertices. Every unordered pair with a moving endpoint belongs to
an orbit of length at least three (a fixed unordered pair would be a fixed
pair of points or a transposition). Of m(2n-m-1)/2 affected pairs, at most
one third as many pair orbits remain. Relative to the identity term, its
fixed-tournament count therefore loses at least ceil(m(2n-m-1)/3) binary
choices. There are at most n!/(n-m)! permutations with that support size.
Summing proves (6).

Writing rho_n=n*2^(-(n-1)/3), for rho_n<1 the latter sum is at most
rho_n^3/(1-rho_n), hence epsilon_n tends to zero. Substituting in (5) gives

```
0 <= 1-r_n/a_n <= n(n-1) 2^(-(n-1)) (1+epsilon_(n-1)),
so r_n/a_n -> 1.                                                        (7)
```

Thus almost every ordinary class already has core order n. The old labelled
union bound alone did not justify this ordinary-class conclusion.

There is also a useful uniform version of the entire column law. Let R_q
count rigid pair-free classes. Every nonrigid tournament has odd automorphism
order at least three, so

```
R_q >= a_q-(q-1)a_(q-1) -(3/2)(a_q-L_q).
```

For every n>=q, rigid cores have all binom(n-1,q-1) distinct weightings;
nonrigid cores have at most this many. Consequently

```
R_q binom(n-1,q-1) <= A[n,q] <= a_q binom(n-1,q-1).                       (8)
```

Equations (6)–(7) give R_q/a_q->1. Therefore, **uniformly for all n>=q** as
q tends to infinity,

```
A[n,q] ~ a_q binom(n-1,q-1)
       ~ 2^binom(q,2)/q! * binom(n-1,q-1).                              (9)
```

This is an asymptotic for a growing core order; it must not be substituted
for a fixed-q asymptotic, where the individual automorphism weights matter.
The proof of (8) needs no complexity assumption and remains uniform even
when the surplus n-q grows arbitrarily quickly.

## 5. Ports and metagraph scope

The weighted core determines the unmarked tournament exactly, but does not
name a Hamiltonian-path presentation. The full line metagraph needs that
extra coordinate. The recovered allocation carrier stores the ordered list
of chain visits and the allocation of **ordered clone ranks** to those visits,
modulo the simultaneous core action. A mere block word loses information:
D2(C3) has ABBCCA presentations whose complements have different degree
multisets, even after converse merging.

Ordinary tile complement can also leave every bounded-core family: the
unique TT_n presentation has line partner with pair-free core order n for
all n>=5. This exact inherited hostile is rechecked through n=30. The
classification and (1) count metagraph **nodes**, not colored line incidences.
They must be combined with marked-path fibres/relocations for those edges.

The full-cycle-index growth hostile is stronger still: the recovered rigid
order-six cores `(0,5,25,35,11,23)` and `(8,17,3,38,13,23)` both have index
p1^6, but the automorphism-order profiles of their 52 pair-free extension
instances are respectively `(50 of order1, 2 of order3)` and
`(51 of order1, 1 of order3)`. The source contains the exhaustive check.
We do not claim to rerun this particular 52-instance profile here; it is
inherited evidence for why aggregate counting indices are not individual
continuation states.

## 6. Verification and precise remaining work

Run from this directory, with Python 3.10+ and no third-party packages:

```
python3 verify_core_families.py --max-order 40 --class-order 7 --out receipt.json
python3 -O verify_core_families.py --max-order 40 --class-order 7 --out receipt-optimized.json
cmp receipt.json receipt-optimized.json
```

The arithmetic universe is every odd cycle partition needed through order40;
no sampled or assumed core list enters that calculation. The separate structural
universe is every tournament class through order7, generated by arc-flip BFS.
The verifier checks full-permutation minima through6, labelled orbit masses,
all actual automorphism cycle-index coefficients, homogeneous-chain round trips,
Aut(Q[w])=weight stabilizer, all explicit positive-weight orbits, no duplicate
family classes, and exact equality of weighted families with the independent
BFS census. Direct bivariate expansion audits (1) through14. All inequalities
(4) and (8) are checked against the exact triangle, without floats.

Remaining open: a substantially smaller all-order structural family grammar,
its closure under every requested line/extension operation, and an asymptotic
complexity advantage. Neither this refinement nor the already-recovered exact
recurrences supplies those claims. No repository edits, commits, or pushes were
performed by this subagent; all deliverables are in its owned directory.
