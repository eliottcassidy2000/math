# No-three-in-line after exact row/column conditioning

**Date: 2026-09-06. Status: PROVED elementary identities + FINITE-EXACT
conditioned census, independently audited.** This is a new repository synthesis and instrument;
no external priority claim or asymptotic solution is asserted.

The [independent proof review and permutation-pair replay](overnight_20260906_lrc_audit_no3line.md)
accepts the matching identities, sampling weights, and the complete `n=4`
mean/variance/zero-event table by a separate no-import implementation.

## Inheritance and source correction

The user supplied [Guy--Kelly (1968), The No-Three-In-Line Problem](https://doi.org/10.4153/CMB-1968-062-3).
Its proved input is the number of grid triples
`t3(n)=(3/pi^2)n^4 log n+O(n^4)`. The paper then assumes independence
to motivate an asymptotic conjecture; that assumption is not a theorem.
[Voutier (2026), arXiv:2603.00215v2](https://arxiv.org/abs/2603.00215v2)
explains the classical correction from the original cube-root constant to
`pi/sqrt(3)`: after replacing `2n` by `kn`, the entropy term must also
change from `2n log n` to `kn log n`. The corrected claim remains conjectural.
Both primary PDFs were read, including the original counting and heuristic
pages; PDF text extraction loses some inequality and exponent symbols.

[Prellberg (2026), arXiv:2602.07751v1](https://arxiv.org/abs/2602.07751v1)
already imposes exact row counts in a constraint model and notes that
rotational symmetry exposes the weakness of the independence heuristic.
The conditioning itself is therefore inherited, not a discovery here.
The new calculation below isolates what that conditioning retains and what
the first scalar moment still loses.

Closest repository mechanism:
[THM-3166 / falling-factorial matching response](../../01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md)
retains a whole rook polynomial rather than one count. Its tournament target
does not transfer here; the shared exact operation is counting matchings.
The hostile is two cycle types with identical triple means but different
zero-event probabilities. The corrected near miss is the Guy--Kelly entropy
coefficient, separately from the unproved independence premise. The least-used
sidecar in this application is the cycle decomposition of the saturated
row/column incidence graph.

Board: grid line hypergraph; primitive direction and span; degree-two
bipartite incidence; matching polynomial jet; cycle-dependent tail events;
owner-colored LRC midpoint closure.

## 1. The faithful saturated carrier and its sampling measure

Let `G_n={0,...,n-1}^2`, `n>=2`. A `2n`-point no-three-in-line set must
have exactly two points in each row and each column. Represent a point
`(i,j)` by an edge from row vertex `i` to column vertex `j`. The result is
a simple degree-two bipartite graph, a disjoint union of even cycles of
length at least four. Write `c_(2r)` for the number of cycles of length `2r`.

Every such graph is a union of two disjoint permutation matrices. On each
cycle its alternating edges have two choices of which permutation comes
first; hence a graph with `c` components has exactly `2^c` ordered
decompositions. Uniform ordered permutation pairs do not sample uniform
uncolored boards. Dividing each pair's weight by `2^c` repairs that measure.

In particular the number `B_n` of labelled saturated boards is exactly

```text
B_n = n! * sum_(rho a derangement of n) 2^(-cycles(rho))
    = (n!)^2 [z^n] exp(-z/2) (1-z)^(-1/2).
```

The last identity follows by expanding independent permutation cycles:
cycles of length `r>=2` contribute `z^r/(2r)` to the exponential generating
function. This is a classical labelled-cycle calculation, proved here to
make the sampling correction explicit. Row and column labels remain essential:
arbitrary relabelling preserves the skeleton but not Euclidean collinearity.

## 2. All cycle types have the same first collinear-triple moment

Fix **any** simple degree-two bipartite graph on `n+n` vertices, `n>=3`,
and independently label the two shores uniformly by `0,...,n-1`. Let `X3`
be the number of collinear triples among its grid edges.

**PROVED.** If `T3(n)=t3(n)-2n*binom(n,3)` counts non-axis grid triples,
then, regardless of the graph's cycle type,

```text
E[X3] = T3(n) * 4(2n-5)/(n(n-1)^2(n-2)).             (1)
```

Proof: three non-axis collinear points have distinct rows and columns,
so form a matching. There are `binom(n,3)^2*3!` labelled three-edge
matchings in the full grid. A degree-two graph on `2n` edges has exactly

```text
m3 = binom(2n,3) - 2n(2n-2) + 2n
   = 2n(n-2)(2n-5)/3.                                (2)
```

Indeed there are `2n` adjacent edge pairs; each extends to a triple in
`2n-2` ways. A triple with two adjacent pairs is a three-edge path, of which
there are `2n`. There are no triangles. Inclusion-exclusion gives `(2)`.
Uniform shore labels make each fixed matching equiprobable; multiplying
its inclusion probability by `T3(n)` proves `(1)`.

The same expectation holds for any mixture of skeleton types whose labels
are conditionally uniform, including both uniform boards and uniform
disjoint ordered permutation pairs. It need not hold after imposing
geometric symmetries on the labels.

Using only the cited triple asymptotic, `(1)` gives
`E[X3]=(24/pi^2)n log n+O(n)`. Exact row/column conditioning leaves this
leading coefficient unchanged. This is an expectation theorem, with no
independence or zero-probability conclusion.

## 3. The first missing cycle coordinate appears at four matching edges

Let `R_G(t)=sum m_k(G)t^k` be the matching polynomial. Put, as formal series,

```text
lambda_+=(1+sqrt(1+4t))/2, lambda_-=(1-sqrt(1+4t))/2,
q=lambda_-/lambda_+.
```

The elementary path recurrence for matchings yields
`R_(C_L)(t)=lambda_+^L+lambda_-^L`, `L>=3`. Consequently

```text
R_G(t)=lambda_+^(2n) * product_(r>=2)(1+q^(2r))^c_(2r). (3)
```

Since `q=-t+O(t^2)`, cycle length `2r` first changes the matching jet
at degree `2r`, with coefficient exactly `c_(2r)`. Thus the `k`-jet
depends on the total size and cycle counts of lengths at most `k`.
Those short cycle counts can also be recovered successively from the jet.

For `n>=3`, the first case is

```text
m4(G)=n(n-3)(2n-7)(2n-5)/6 + c4(G).                  (4)
```

To check the universal term, a single cycle of length `2n` has
`m_k=(2n)/(2n-k)*binom(2n-k,k)` when `k<=n`; its coefficient at four
is the polynomial in `(4)`. For `n=3` both sides at degree four vanish.
Each four-cycle contributes one additional `t^4` in `(3)`.

If `T4(n)` counts non-axis collinear grid quadruples and `n>=4`, the same
matching-incidence argument gives

```text
E[X4] = T4(n) * m4(G)/(binom(n,4)^2*4!).              (5)
```

This is a precise route from a forgotten graph coordinate to a geometric
clumping statistic. It does not recover the entire variance of `X3`: pairs
of triples may share vertices or use other row/column incidence types.

## 4. Exact hostile: equal means, different tails and sampling bias

On the four-by-four grid the two possible cycle types give:

| Incidence graph | Labelled boards | E[X3] | Var[X3] | E[X4] | P(X3=0) |
|---|---:|---:|---:|---:|---:|
| Two four-cycles | 18 | 2 | 56/9 | 1/3 | 1/2 |
| One eight-cycle | 72 | 2 | 25/18 | 1/6 | 1/36 |

These rational numbers come from the exhaustive independent-label universe,
not Monte Carlo. The zero-event probabilities differ by a factor eighteen
despite identical means. Uniform boards have zero probability `11/90`;
uniform disjoint ordered permutation pairs give `5/27`, since the former
cycle type has four decompositions and the latter two.

This refutes a finite-cycle-blind tail model, not the asymptotic conjecture.
Nor is more short-cycle structure uniformly beneficial: at `n=5`, all
600 boards of cycle type `(4,6)` fail, while 32 of the 1,440 ten-cycle
boards pass. Geometry and full event dependence still matter.

## 5. Exact reproduction and controls

Run `python3 -B 04-computation/overnight_20260906_no3line.py` and its `-O`
variant. The matching output is saved beside this report. The universe is
every binary row/column-degree-two board for `2<=n<=6`, with no symmetry,
cycle-type, or collinearity prefilter. Counts are `1,6,90,2040,67950`.
The recursion prunes only exhausted or unfillable column capacities.

Two independent geometric paths are used: canonical pair-line multiplicity
and literal three-point determinants (the latter exhaustively through `n=5`).
Whole-grid triple/quadruple counts through `n=12` are independently obtained
from primitive direction and endpoint span:

```text
sum_(primitive oriented v) sum_(s>=k-1)
    binom(s-1,k-2) (n-s|v_x|)_+ (n-s|v_y|)_+.
```

Each endpoint pair spans `s` primitive steps, with `k-2` interior choices.
Axes are included or omitted explicitly. Direct matching enumeration through
`n=5` checks `(2)` and `(4)` for every board. All cycle-class means through
`n=6` check `(1)` and `(5)`. Square-symmetry orbit counts reproduce the
original small-grid controls `1,1,4,5,11` for `n=2,...,6`.

## 6. Connection contract and next decisive questions

| Map | Preserved predicate | Lost information / next test |
|---|---|---|
| Saturated grid points -> bipartite edges | Row/column capacities and each labelled point | Abstract cycle type loses geometry; retain both shore labels |
| Two permutations -> their union | Saturated occupancy | Decomposition multiplicity; weight by `2^(-components)` |
| Fixed skeleton -> random shore labels | Matching incidence probabilities | Mean loses joint line events; the n4 hostile is decisive |
| Matching jet -> short-cycle counts | Exact triangular identity `(3)` | Long cycles appear only at their own degree |
| Primitive endpoints -> grid line tuples | Euclidean collinearity and scale | Reducing direction modulo a prime loses actual area zero |

**Completed continuation:** the [joint-incidence theorem](overnight_20260906_moments_pairprofile_theorem.md)
preserves entire unions and proves, for every `n>=6`,
`Var(X3)=A_n+B_n*c4+D_n*c6`. More generally, the kth factorial moment
sees only cycles of length at most `3k`, with the same weighted degree
bound. The [finite compiler](../../04-computation/overnight_20260906_no3line_pairprofiles.py)
checks every skeleton cycle profile through `n=8`; both proof and finite
normalizations passed independent review. Crossing-line events are included.

Next: test whether the third factorial moment actually has nonzero `c8`
or `c4^2` coefficients, and explain the geometric sign changes in the
second moment. These moments alone do not prove a positive zero-event
probability. Independently, the [LRC full-cap theorem](overnight_20260906_lrc_cap_carriers.md)
uses owner-colored convex midpoint closure to prove its bound; that is a
separate arithmetic mechanism with a full-dictionary hypothesis.
