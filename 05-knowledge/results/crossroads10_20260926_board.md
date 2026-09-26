# Ten vertices, 223 and 233: finite objects that retain the arithmetic

**Current synthesis, 2026-09-26. PROVED / FINITE-EXACT / OPEN, separately
labelled below. Collatz and G2 remain OPEN.** This continues the
[233 session](crossroads233_20260926_board.md), rather than restarting its
rank-height, Bernoulli-boundary or residue-quotient experiments.

The prompt's difference233-223=10 led to an actual two-edge operation.
The strongest mathematical progress came from two different finite objects:
a weighted prime-incidence graph that preserves multiplicative relations,
and the legal pairing tree whose scale averages have a complete certificate
calculus. Neither is licensed by a numerical coincidence alone.

## 1. The exact tournament connection

**FINITE-EXACT, exhaustive neighborhood.** The inherited eight-vertex
staircase has233 Hamiltonian paths. None of its28 one-edge flips has223;
exactly two of its378 two-edge flips do. For one such pair the response is

| Reversals | none | first | second | both |
|---|---:|---:|---:|---:|
| Hamiltonian paths |233|291|123|223|

The interaction is223-291-123+233=42. **PROVED:** when the two original
arcs form u->v->w, this mixed difference counts paths containing the
contiguous forward or reverse triple. Here it is31+11=42. Contracting
that triple gives the exact smaller object; it retains the two boundary
orientations and can cease to be a tournament.

The full odd-cycle formula explains the net difference10:
the number of directed triangles and disjoint triangle pairs is unchanged,
five-cycles fall by3, seven-cycles fall by4, and disjoint3+5 pairs rise by1.
Thus the change is2(-3-4)+4=-10. This is a compatibility effect involving
longer cycles, invisible to the triangle count alone.

**PROVED ten-vertex examples:** adjoining a source and a sink transports
these counts bijectively to tournaments on ten vertices. Their component
sizes are1,8,1. This does not give a strong ten-vertex example or a unique
role for ten; any larger order admits the same padding. The literal
ten-vertex staircase instead has2489 paths, or2265 after the same flips.

[Exact graph note, arrays, proofs and controls](crossroads10_20260926_graph.md).
Its odd-cycle formula is inherited from THM-002 and the
[Grinberg--Stanley paper](https://arxiv.org/abs/2307.05569); the local
contraction mechanism is inherited from THM-082. No priority claim is made.

## 2. A graph genuinely carries number-theoretic relations

**PROVED, THM-4497.** For any finite family of positive integers greater
than1, gcd refinement constructs pairwise coprime atoms. Their exponent
matrix has exactly the same relation kernel as full prime factorization.
Rows supported at one node ground that node; rows supported at two give
weighted edges; the remaining rows form a residual matrix after graph
elimination. This preserves every multiplicative relation, not just rank.

An edge equation u c_i+v c_j=0 forces opposite coefficient signs. An odd
cycle therefore kills its whole connected component. An even cycle can
also kill it when the endpoint exponent ratios fail to balance. The
unweighted graph by itself is insufficient: the same ten-cycle can have
rank9 or10 after changing just one exponent.

The actual odd Collatz orbit of27 contains

    91=7*13,     175=5^2*7,     325=5^2*13.

These form a prime-incidence triangle. They are independent despite having
no private primes and failing strict gcd dominance at every vertex. This
is a direct improvement over the simple sufficient tests in THM-4493.
Seven other nodes of the same orbit give an actual ten-node example.
The consecutive ten-node prefix from199 contains another such triangle.

**FINITE-EXACT universe:** all9446 distinct-node ten-odd-node prefixes
with odd source3..20000. First-slot independence is certified in6487 by
private primes,8481 by one graph pass,9061 by iterative removal, and9350
by the exact weighted graph plus residual matrix. Retaining BOTH paired
coordinates certifies rank10 in all9446. This is not an all-source theorem.

The common factor3 is essential. One actual orbit contains507 and13;
these are independent, but their first slots1521 and39 satisfy1521=39^2.
Their full pairs remain independent. Graph incidence, exponent weights,
the common prime, and the two separate coordinates all have work to do.

[THM-4497: precise statement](../../01-canon/theorems/THM-4497-coprime-graph-rank-compression.md),
[proof, standard gcd-free-basis attribution and exact census](crossroads10_20260926_arithmetic.md).

## 3. A complete family of finite scale certificates

**PROVED, THM-4496.** In the repo's globally legal two-step pairing model,
normalize a finite positive cutoff kernel to a probability vector p on
scales r^k, r=2/3. Its limiting common-assignment minimum B(p) is concave,
shift invariant, and Lipschitz in total variation. Hence

    B(p*q)>=max(B(p),B(q)).

If u_m is uniform on the first m consecutive scales, then

    lim_m B(u_m)=B_*=sup_(all finite kernels p) B(p).

Thus the entire finite positive-kernel method is exhausted by one explicit
family. Every global legal pairing has lower logarithmic density at least

    B_* >=821510388809/2677850419968 =0.30677978974599224... .

This uses the inherited exact rational certificate; the new result extends
its consequence to lower logarithmic density and explains the kernel family.

The strongest new result identifies and attains the global cost:

    lim_(X->infinity) H(X)/log X = B_*,
    one global legal pairing has logarithmic density exactly B_*.

Here H(X) is the minimum finite harmonic cost sum_(i in F,i<=X)1/i.
Prescribing any legal prefix throughY costs at most log Y+3 extra.
More precisely, splitting the tree into Q logarithmic phase cells aboveY
costs at most18Q/(Y-1) to repair, uniformly in the final cutoff. Averaging
shifted phase weights and then assembling long shells proves the equality.
The small tree controls the infinite construction because both boundary
costs are proved. This is an exact variational principle for this model.

A compact family of actual integer subsets with a shared logarithmic phase
still has every grid-kernel value sqrt(6)-2 but logarithmic density1/2.
It refutes a generic compactness argument; it lacks the pairing tree's
quantitative phase repair. The equality above required that extra structure.
The value of B_* beyond the certified lower bound and natural-density
attainment remain OPEN here.

[THM-4496](../../01-canon/theorems/THM-4496-pairing-kernel-convolution-and-log-density.md),
[full proof, harmonic extension and phase analysis](crossroads10_20260926_flow.md).

## 4. Recovered connections and the active research board

Anchor: the global pairing cost and its phase dependence. Niche: exact
weighted prime-incidence reduction. Wildcard: the223/233 tournament move.
The following six concepts were compared after the main experiments.

| Source -> target and map | Preserved predicate | Lost data / next decisive test |
|---|---|---|
| Tournament edge flips -> boundary-contracted paths | Exact mixed path response | Scalar H loses which compatible blocks contribute; retain block ends |
| Integer family -> coprime exponent rows -> weighted graph + R | Entire multiplicative relation kernel | Unweighted graph loses exponents; deleting R loses shared primes |
| THM-2521 LRC signless potential operator -> equal-exponent prime edges | Same linear equations c_i+c_j=0 | LRC's landing predicate is not transported; no prize implication |
| THM-3387 cyclic-sheet gcd graph -> arithmetic incidence lens | Intrinsic symmetric relation | Different target, sheet blocking versus rank; do not force a tournament |
| Pairing cutoff costs -> kernels -> uniform windows -> repaired phase cells | Exact minimum logarithmic cost | Scalar limits alone lose phase; finite repair and phase averaging restore it |
| Harmonic finite tree -> one infinite legal pairing | Actual logarithmic density at the minimum B_* | Natural density and full convergence require additional control |

The Sun reflection carrier in THM-4246 is bipartite; the arithmetic graph
here has actual odd cycles. This distinguishes two uses of graph parity
instead of treating every appearance of an odd cycle as interchangeable.
The LRC signless operator supplies an exact algebraic analogy, not an LRC
proof. The ten-vertex centrality census THM-4137 concerns certified floors
on strong tournaments; padded1+8+1 examples lie outside that hypothesis.

The Bernoulli B_1 boundary idea survives as a prompt to measure, rather
than ignore, the finite head. Here the head's harmonic extension cost is
bounded explicitly, while a shared logarithmic phase can persist in the
bulk. This is a precise boundary-versus-bulk distinction, not a new identity
involving Bernoulli numbers or a proof from their odd-index vanishing.

For223 and233 themselves, ord_p(2) is37 versus29 and ord_p(3) is222
versus232. The quadratic-residue relation is antisymmetric at223 and
symmetric at233. Their multiplicative2-subgroups have indices6 and8.
These are exact local distinctions; subtracting the primes does not
supply an arithmetic map to a ten-vertex tournament.

## 5. Strongest next obligations

1. **OPEN, sharp value and stronger density.** The phase-sensitive
   equality h=B_* is now proved. Compute useful certified upper bounds
   for B_* alongside the lower bound0.3067797897, and determine whether
   some logarithmically optimal pairing also has natural density B_*.
   Phase homogeneity was not needed: averaging shifted weights after
   quantitative repair closed the comparison. The next numerical work
   should exploit this complete one-parameter family and its boundary
   bounds, rather than search arbitrary unrelated coefficient lists.
2. **OPEN, carry-constrained residual rank.** After graph elimination,
   test the remaining matrix on chronological paired families. Shared
   primes divide the actual subword affine carries (THM-4493), so arbitrary
   synthetic graphs need not be realizable. Use those divisibility
   constraints to exclude a surviving relation. Appending a node changes
   prime supports; finite ten-node certificates alone do not iterate.
3. **OPEN, higher-horizon compatibility.** The exact two-step tree no
   longer describes four-step legality. Retain the genuine longer-horizon
   clauses and their boundary states before transporting any density or
   phase theorem. The known four-step counterexample remains a hostile
   control. None of these targets yet supplies a contradiction for a
   hypothetical divergent Collatz orbit.

## 6. Audit and reproducibility

All three lane programs use exact integers/rationals and explicit checks
that remain active under Python optimization. The root
[reproduction program](../../04-computation/experiments/crossroads10_20260926_audit.py)
compares normal and optimized runs with retained outputs, independently
enumerates the four8! tournament path counts, evaluates the arithmetic
triangle determinant, and exhausts the1024 assignments of the ten-bit
pairing carrier. Its [output](crossroads10_20260926_audit.out) records the
verdicts. A companion manifest records committed-content hashes.

Root proof audit accepted the kernel limits, harmonic extension, phase
cuts and averaging/assembly equality, gcd refinement, weighted residual
matrix and sign torsion. Geometry and
arithmetic independently audited each other's core statements; geometry
also audited the kernel calculus. All claims stay within their declared
universes. Reused META-PATTERNS cards: retain local incidence through
counting; type every analogy; find the missing second coordinate; retain
support separately from bounded-height coverage. No unearned new general
research rule was promoted.

Incoming work was checked: THM-4494's binomial monotonicity correction
removes an earlier numerical lead. During this session another author
promoted THM-4495/4498 with refined polynomial factors for no-descent and
dip-spectrum counts. These concern the finite-word population and do not
alter the pairing clauses, phase repair, or prime-incidence kernel; they
are not imported as dependencies here. Its second reservation collision
was repaired upstream by renumbering its theorem4496 to4498. Our pushed
reservation4496 remains the pairing theorem. The checkout interrupted during creation
was repaired from its own HEAD, preserving all authored files and the
unrelated shared checkout. Only this session's explicit paths are staged.
