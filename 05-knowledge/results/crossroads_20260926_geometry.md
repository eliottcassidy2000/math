# Collatz crossroads: thin tails, arithmetic shells, paired valuation rank

**Status: PROVED elementary consequences of THM-4476 after an independent
proof audit; CITED unit-equation bound; FINITE-EXACT experiments; two OPEN
strengthenings. No claim to prove Collatz.** Geometry lane, 2026-09-26.

Reproduction: `python3 04-computation/experiments/crossroads_20260926_geometry.py`.
Matching output: [crossroads_20260926_geometry.out](crossroads_20260926_geometry.out).
The script prints its own SHA256. All substantive experimental comparisons
use integers or `Fraction`, not floating-point thresholds.

## 1. Inheritance and concept board

Anchor: Collatz integer height together with residue chronology. Niche:
Newton/valuation matrices and multiplicative rank. Wildcard: logarithmic
shell occupancy instead of a global Lyapunov function.

* Closest proved mechanism: THM-4476,
  `01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md`:
  Terras cylinders, an injective landing pigeonhole, and a scale bootstrap.
* Canonical hostile: the positive `3n-1` cycle `5,7`; infinite harmonic
  mass certifies eventual periodicity, not the identity of the cycle.
* Corrected near miss: THM-2991,
  `01-canon/theorems/THM-2991-pf-infinity-arbitrarily-delayed-newton-ratio-return.md`.
  Even strong local positivity can postpone a different global extremum
  arbitrarily far. A finite observer is not a global height certificate.
* Underused sidecar: THM-3152,
  `01-canon/theorems/THM-3152-multi-place-newton-degree-barcode-and-euclidean-flag-census.md`.
  Keep the joint prime-coordinate capacities, including the coordinate
  root. Separate prime summaries forget compatibility.
* Flatness source: THM-3743,
  `01-canon/theorems/THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction.md`.
  Its actual input is a lattice-point-free body. Each finite Collatz word
  has an infinite arithmetic progression of positive integer sources, so
  flatness has no input until a real height interval is imposed.
* Log-curvature source: THM-3023,
  `01-canon/theorems/THM-3023-newton-ratio-transform-dynamics.md`.
  Its Newton transform kills geometric factors, hence also kills precisely
  the affine log-drift that a Collatz argument must control.

These comparisons retained six coordinates: actual integer source, exponent
word, affine carry, minimum future height, paired prime valuations, and
first-descent/periodic target. No tournament structure is intrinsic here.

## 2. Independent audit of thin divergence

The core proof of THM-4476 passes a line-by-line independent audit. Its
parity bijection is valid for fixed odd b; the carry cutoff, binomial tail,
landing multiplicity, and strong induction all have the stated scope.
It also applies uniformly to finite injective trajectory segments, with
the final k points treated as boundary exceptions.

Two corrections were communicated to the parent for repair in the sources.
Incoming commit `777e4e137` independently repairs the harmonic wording and
mentions the bounded-product repair; the parent retains that lineage and
tightens the displayed inequality and symmetric threshold accordingly.

1. The statement that Collatz is equivalent to harmonic divergence on all
   positive orbits omits exclusion of other positive cycles. Harmonic
   divergence is equivalent to **no divergent positive orbit**. Add cycle
   exclusion to obtain full Collatz. For the half-step convention the
   allowed cycle is `1->2->1`.
2. Corollary 9's plus-sheet estimate cannot upper-bound the partial sum of
   `2^Delta_i` using only a lower bound on Delta. Its conclusion survives:
   earlier corollaries give bounded normalized products, and this directly
   supplies the needed upper bound on the orbit height.

Specifically put alpha=log_2(3), S_j=sum of the first j odd-step valuations,
Delta_j=S_j-j alpha, and

```
m_j = 2^(-Delta_j) c_j,
c_j = n product_(i<j) (1+1/(3m_i)).
```

For an injective positive orbit, thinness gives `c_j -> c<infinity`.
Therefore `Delta_j >= -a log_2 j-O(1)` implies `m_j<=C j^a`.
For every `a<1/h*`, with `h*=h(log_3 2)`, this contradicts thinness.
This is a one-sided drift obstruction, stronger than the earlier symmetric
plus-sheet constant 0.02634. It does not exclude exponential growth.

## 3. A quantitative real interval at each future minimum

**PROVED, conditional only on the audited thinness theorem.** Fix
`h*<gamma<1`, and its uniform count constant C_gamma. For any finite or
infinite injective positive orbit segment all of whose entries are at
least M,

```
sum_j 1/m_j <= C_gamma/(1-gamma) M^(gamma-1).             (1)
```

Proof: its counting function vanishes below M and is bounded by
`C_gamma X^gamma` thereafter. Integrate against `X^-2`. This includes the
atom at M by taking the lower endpoint from below.

On the plus sheet the accumulated normalized product satisfies

```
1 <= product_j(1+1/(3m_j))
  <= exp(C_gamma/[3(1-gamma)] M^(gamma-1)).               (2)
```

If the first term itself is M, all finite normalized affine prefixes and
the infinite limit, when it exists, lie in

```
M <= m_j 2^S_j/3^j <= M+O_gamma(M^gamma).                (3)
```

Thus real normalized centers occupy a uniform sublinear interval while
their exact dyadic itinerary is retained. This is the concrete transfer
from counting geometry to a real/2-adic carrier. It still allows an
unbounded interval; rounding to a single integer is unjustified.

## 4. Why changing only the window length does not improve h*

Take a forward window of length `k=floor(c log_2 X)`, `0<c<=1`.
The same parity entropy count has limiting exponent

```
1-c+c h(log_3 2) = 1-c(1-h*),                          (4)
```

as the dip threshold parameter tends to zero. Its minimum occurs at c=1.
Improving the landing multiplicity from log X to a constant changes only
subpower factors. For c>1 the parity modulus exceeds X: bijectivity alone
provides no equidistribution of the relevant residues in this short
integer interval. A two-place argument likewise needs a joint residue
height bound, not multiplication of marginal savings.

## 5. OPEN hypothesis G1: logarithmic shell occupancy

For every injective positive `3n+1` trajectory segment A, ask whether

```
#(A intersect [M,2M)) <= C log(2M)                     (G1)
```

holds with an absolute C. Test the same question on `3n-1` as a hostile
sheet. This is stronger than current thinness and does not assert descent.
If true, summing dyadic shells gives, for a segment staying above M,

```
sum 1/m_j = O(log(2M)/M),
c-M = O(log(2M)) at a future-tail minimum M.             (5)
```

This is a useful target: it would compress (3) from power width to
logarithmic width without pretending it is already a discrete rank.

**A constant shell bound is REFUTED structurally.** For arbitrary r choose
`S_j=floor(j log_2 3)` and exponent letters `S_(j+1)-S_j` in {1,2}.
Their exact source cylinder is

```
3^r n+C_r = 2^S_r mod 2^(S_r+1),
C_r=sum_(i<r) 3^(r-1-i)2^S_i.
```

Choose its representative plus `2^(S_r+1)(r+1)4^r`. Then

```
n <= m_j < 2(n+j/3) < 3n,       0<=j<=r.
```

There are no repetitions: any repeated subword would give a positive
cycle value at most its carry, bounded above by `r4^r`, whereas every
term exceeds that. The interval [n,3n] meets at most three standard dyadic
shells, so one shell contains at least `(r+1)/3` distinct points.
Here log n=O(r), making logarithmic occupation unavoidable in scale.

**FINITE-EXACT:** these lifts were verified for r=8,16,...,512. On every
positive odd source below 100001, full distinct paths had maximal dyadic
shell counts 24 (plus source53529, shell32768) and25 (minus source9351,
shell8192). No 2000-step cap was hit. The separate 5n+1 hostile universe
has 444 capped paths among500 and supplies no termination conclusion.

## 6. A different carrier: paired prime valuations of transitions

Every plus transition provides the exact unit equation

```
x_j=-3m_j,       y_j=2^k_j m_(j+1),       x_j+y_j=1.     (6)
```

Let r_N be the torsion-free rank of the multiplicative subgroup of
`(Q*)^2` generated by the first N transition pairs. It is the rational
rank of a matrix with rows indexed by `(coordinate slot, prime)` and
columns indexed by transitions. Both slots must be kept: prime support
alone forgets their coupling.

**CITED + PROVED application.** Beukers--Schlickewei, Theorem1.1, bounds
the number of solutions of x+y=1 in a rank-r subgroup (indeed its
Q-closure) by `2^(8r+8)`. Distinct orbit sources give distinct pairs, so

```
r_N >= (log_2 N)/8-1.                                 (7)
```

If s_N counts all primes in the odd nodes through time N, then
`r_N<=2s_N+2`, by adjoining the two fixed coordinate generators for3 and2.
Thus `s_N >= (log_2 N)/16-3/2`. In particular an injective infinite orbit
cannot use a fixed finite prime bank. This observation does not establish
descent: new primes are perfectly compatible with divergence.

Primary source: F. Beukers and H.P. Schlickewei, *The equation x+y=1 in
finitely generated groups*, Acta Arithmetica78(1996),189--199,
[publisher](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/78/2/109078/the-equation-x-y-1-in-finitely-generated-groups),
[indexed primary first page](https://matwbn.icm.edu.pl/ksiazki/aa/aa78/aa7826.pdf).
The primary search extraction supplies the full statement of Theorem1.1;
direct PDF requests returned403. The theorem is cited, not independently
proved here, and no priority is claimed for application (7).

**OPEN hypothesis G2:** transition pairs along every finite injective
positive plus orbit segment are multiplicatively independent (`r_N=N`).
This stronger statement is deliberately separated from the cited theorem.
Exact rational rank probes on1000 source paths through1999 and nine
selected long paths are recorded in the output. A separate modular rank
computation at the prime1000000007 certifies FULL rank on all50000
complete distinct plus paths from odd starts through99999. Full rank
modulo a prime implies full rational rank, since the matrix entries are
integers and a nonzero modular minor is a nonzero integer minor. There
were no deficient modular cases requiring a second-prime replay.
Rank129 occurs on the129-edge path starting77031. Minus cycles1,5/7,
and the seven-cycle through17 retain full rank when each closing edge is
included exactly once. Repeated laps would duplicate columns and are
deliberately excluded.

**Exact hostile: local edge injectivity is insufficient.** Put

```
P={11,43,65,253},   N={13,23,121,215}.
product_(n in P) (3n)   = product_(n in N) (3n),
product_(n in P) (3n+1) = product_(n in N) (3n+1).        (H)
```

The eight nodes in sorted order have distinct successors
`17,5,35,65,49,91,323,95`. Thus even a finite injective edge set can carry
a nontrivial paired multiplicative relation. They are not all comparable
along one chronological orbit: their later branches merge. G2 must keep
full chronology, not merely distinct source and target labels.

**Fresh-prime shortcut REFUTED.** The literal rising run `55->83->125`
has valuations1,1, but125=5^3 introduces no new prime:5 already divides55.
In the general rise formula `m_j=3^j 2^(H-j)u-1`, standard Zsigmondy for
fixed-base differences does not automatically apply. The prefactor
`2^H u` depends on the run height; the needed uniformity is absent.
See the recovered primitive-divisor thread
`collatz_mod6_20260917_zsigmondy_triad.md` for the correctly scoped
fixed-base theorem and its exceptions. The exact rank carrier can retain
multiplicity and coordinate position even when prime support does not grow.

Even G2 would not prove Collatz; its useful consequence would be a linear
lower bound on fresh prime-coordinate complexity. Cheapest next tests:
larger paths, rational cycles, and solving for short integer dependencies
with paired prime valuations retained.

## 7. Independent audit of the parent's sharp-price proof

**PASS, algebraic and counting proof audit.** This is a finite-horizon
pairing-family statement, independent of Collatz convergence.

Retain source parity words of length L with every prefix slope
`w_k=3^e_k/2^k` between1 andK. Their affine maps are

```
T^k(n)=w_k(n+h_k),  h_k=sum_(i<k,bit_i=1)1/(3w_i),
0<=h_k<=k/3.
```

For a fixed actual hub v and fixed(k,e), every actual ancestor lies in
`[v/w-k/3,v/w]`. There are at most `floor(k/3)+1` integer ancestors.
There are at most `floor(log_3 K)+1` exponents e. Consequently a vertex
serves at most

```
M=(floor(log_3 K)+1) sum_(k<L)(floor(k/3)+1)
```

retained sources, and a flipped pair serves at most2M. Deterministic
chronology prevents different residue completions from counting the same
ancestor twice. Every L-step-descent strategy must encounter a flipped
pair before step L along each retained Collatz source path.

Sources n<=X hit only vertices at most `K(X+L/3)`, so their hit pair indices
are at most `Y=(K(X+L/3)+1)/2`. If F(Y) counts flipped pairs and rho is the
retained source density, then

```
rho X+O_L(1) <= 2M F(Y),
liminf_(Y->infinity) F(Y)/Y >= rho/(KM).                (8)
```

Both pair endpoints and the factor1/2 in Y are accounted for. Natural
density need not exist; lower density is enough.

For the entropy construction put b=floor(sqrt L),
e=ceil(b log_3 2), t=floor(L/b), r=L-tb. Among all binomial(b,e) words,
at least binomial(b,e)/b have nonnegative log-slope prefixes: rotate a
word to a minimum prefix; its total log-slope is nonnegative and less
than log3. Concatenate t such blocks and append r ones. Every resulting
word stays in the tube with `K=3^(t+b+r)`, and

```
rho >= 2^(-(1-h*)L-O(sqrt L log L)),
KM = 2^O(sqrt L) poly(L).
```

Equation(8) supplies the matching lower exponent1-h*. Together with
THM-4475's existing upper bound, this closes the stated sharp-price
exponent target. The rotation, exact integer ancestor counts, carry
normalization, first-hit implication, pair normalization, and asymptotic
limits were all audited independently of the parent's implementation.

The accompanying script independently enumerates cyclic-minimum words
through b14 and actual source/hub incidence through source5000 at
L4,8,12 and K2,4,16. These are controls, not substitutes for (8).

## 8. Independent audit of arbitrary vertex edits

The same proof extends to any deterministic positive-integer map G that
agrees with T outside an edit set E and gives every n>=2 descent withinL.
Each retained source must first encounter E before stepL. A hub serves at
mostM sources and is at mostK(X+L/3), so

```
rho X+O_L(1) <= M #E(K(X+L/3)),
lower_density(E) >= rho/(KM).
```

For the upper bound, edit exactly the actual no-descent set
`Bad_L={n>=2: T^j(n)>=n for 1<=j<=L}` by putting G(n)=1 there.
An originally bad source descends immediately. An originally good source
either follows T to its old descent or first hits an edited vertex and
then reaches1<its source inside the horizon. Thus G has the required
property. The actual set differs by finitely many integers from the
slope-undecided residue union: each fixed word with some slope<1 has
only finitely many small carry exceptions, while words with all slopes>1
never descend. Equality of a nonempty slope with1 is impossible because
3^e!=2^j forj>0. Therefore its density is rho_L.

The arbitrary-edit optimum has the same sharp exponent1-h*, with upper
boundrho_L, whereas the pairing construction has upper bound2rho_L.
This audit finds no additional analytic-density or first-hit assumption.
