# SOP trace torsors, the Rule 30 odometer, and the moving-observer boundary

**Research reflection -- 2026-08-15.**  Exact statements are in
THM-3456/3458/3459/3463; the external model-theory theorem and Rule 30 prizes remain
cited/open as marked below.

## Outcome first

The productive merge was not an attempted identification of model-theoretic
SOP, cellular-automaton complexity, and the three mathematical frontiers.  It
was a three-layer separation:

```text
free-input language       distinguished orbit        target observer
trace-prefix torsor   +   fixed seed/boundary    +   center k=t
```

Rule 30 is completely transparent at the first layer: every binary trace is
realized, uniformly, and the trace plus the initial right half is a coordinate
system.  Its right edge is also completely transparent at a second, different
layer: it is a 2-adic odometer on the seed-orbit closure.  The prize sequence
evades both closures because the seed removes the free trace coordinates and
the center reads a coordinate whose address grows with time.

This session now routes through four proved repository objects:

- [THM-3456](../01-canon/theorems/THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary.md):
  the finite-alphabet trace/right-half homeomorphism, direct enriched
  `SOP`/`SOP2`/`SOP3` witnesses, and exact single-seed boundary word;
- [THM-3458](../01-canon/theorems/THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary.md):
  the packed right-edge recurrence, compatible seed periods, 2-adic odometer,
  fixed-offset rational sequences, a width-six moving-observer hostile, and de Bruijn
  C-finite inverse-count compiler;
- [THM-3459](../01-canon/theorems/THM-3459-rule30-ternary-intersection-factorial-truth-lift-and-keller-boundaries.md):
  the mask-intersection compiler, factorial moment formulas, and the
  characteristic/representative Keller boundary.
- [THM-3463](../01-canon/theorems/THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary.md):
  the center-section compiler, dual nonlinear-current charts, period-lift
  restriction, spatial Walsh atlas, and robust-input complexity boundary.

None settles any Rule 30 prize, LRC(14), FC(3), or JC(2).

## External inputs, accurately typed

[Chernikov's preprint](https://arxiv.org/abs/2608.13291) proves
`SOP2 => SOP3`, hence equality of the two theory classes (and, with the cited
recent work, `SOP1=SOP2=SOP3`).  The proof begins with a treetop-indiscernible
SOP2 witness, splits on a mixed partial type, and retains witness--parameter
pairs in both branches.  The reusable methodological point is typing: a bare
tree vertex is not the proof object.

The [2019 Rule 30 announcement](https://writings.stephenwolfram.com/2019/10/announcing-the-rule-30-prizes/)
asks about the isolated-single-seed center column.  The
[current official page](https://rule30prize.org/) still actively lists all
three prizes and accepts submissions; on that dated evidence the repo treats
the questions as open.  The third question has to be handled with special care: its prose,
big-O wording, and displayed predicate should not be blended without fixing a
machine and bit-cost model.  Nothing here lower-bounds the fixed-seed Prize-3 task.

## Inheritance pass

| Lane | Closest proved mechanism | Canonical hostile | Corrected near miss | Least-used sidecar |
|---|---|---|---|---|
| LRC anchor | THM-3395 typed owner/coset/star cochain; THM-2050 global phase exit | `q=6`: `(2,8,14)` and `(2,8,10)` share the same three masks but only one has a star | mask partition is necessary-looking but loses affine gaps/common time | labelled star cochain and physical row |
| FC niche | THM-2810 growing factorial carrier; THM-3362 odd-profile moment detector | Boolean `x^2=x` versus factorial square-zero at `p=2` | Rule 30 truth lift has `L(g)=0` but `L(g^2)=6` | polynomial representative, characteristic, all powers |
| JC niche | THM-3367 separates integrability, determinant, and coefficient image | same Boolean rule has a tame ANF lift and a faithful non-Keller lift | formal `Jac=I` in characteristic two survives noninjectivity | boundary condition, characteristic, representative gauge |
| sequence wildcard | THM-3359 fixed finite state gives periodic support/density | `R_7=R_15 mod 64` but center bits differ | every fixed edge column is rational, but the center is their diagonal | observer address/scale |
| SOP wildcard | THM-3395 already imported typed witness pairs as motivation | Rule 90 has the same prefix-cylinder tree but trivial seed center | causal address tree collapses `LR,CC,RL` to one cell | full address and chosen language |

The closest mechanisms were therefore already warning against quotienting.
Rule 30 supplied an unusually clean laboratory in which the missing coordinate
could be named exactly.

## Concept board

The session kept six objects live:

1. the left-permutive trace-prefix tree;
2. the fixed single-seed section through that tree;
3. the packed edge permutation tower;
4. the moving center observer `k=t`;
5. labelled mask/cochain witnesses for LRC;
6. polynomial representatives of one Boolean truth function.

Every new identity was compared against all six.  This prevented three
tempting but false promotions:

- free-input trace fairness is not deterministic seed balance;
- an odometer with periodic continuous observables does not make a moving
  diagonal observable periodic;
- a truth-table identity does not select a characteristic-zero Keller lift or
  a factorial carrier.

## Pull 1 -- the exact trace torsor

For any radius-one rule

```text
f(l,c,r)=l xor h(c,r),
```

the center at time `t` has the form

```text
(F^t x)_0=x_(-t) xor H_t(x_(-t+1),...,x_t).
```

Thus the center trace and positive initial half-line reconstruct the negative
half-line successively.  At horizon `T`, every trace word of length `T+1` has
exactly `2^T` cone preimages.  Under iid fair initial data, the vertical trace
is exactly iid fair.

This is stronger than a heuristic mixing statement and weaker than the prize
in precisely the right way.  The prize seed fixes the positive half-line; the
trace is then no longer a free coordinate.  Rule 90 and Rule 150 are decisive
controls: the same torsor supports respectively an eventually-zero and a
constant-one single-seed center.

### SOP port

In the explicitly expanded two-sorted structure

```text
configuration x --E--> finite trace prefix w,
```

prefix cylinders form an SOP2 tree: every branch is realized and incomparable
prefixes are disjoint.  Chernikov's theorem gives SOP3 for the complete theory
of that expansion.

The port is exact but almost tautological.  It depends on naming the trace-
prefix incidence relation.  It applies to Rule 90 and therefore carries no
single-seed complexity.

Separately, any future **causal-path-indexed** family faces an endpoint hostile:
`LR`, `CC`, and `RL` are distinct paths with the same cell `(2,0)`.  This is not
a factorization claim about the binary `E`-witness above.  It is a clean
instance of the paper's typed-pair lesson: retain path plus parameter when a
construction actually uses paths, rather than a bare endpoint.

## Pull 2 -- the right-edge odometer

Reading the Rule 30 row inward from the right edge gives

```text
R_(t+1)=R_t xor ((2R_t) or (4R_t)),       R_0=1.
```

Modulo `2^w` this is a unitriangular Boolean permutation.  If `P_w` is the
seed period, then

```text
P_(w+1)/P_w in {1,2},
P_(w+1)=2^(epsilon_w)P_w,
epsilon_w=bit_w(R_(P_w)).
```

The periods are unbounded powers of two, so the inverse limit of the seed
cycles is `Z_2`, and evolution is addition by one.  This is a genuine exact
closed subsystem, not merely a fitted recurrence.

Every fixed offset `bit_k(R_t)` is periodic.  At offset five the word is

```text
(00010101)^infinity,
```

with rational OGF `(z^3+z^5+z^7)/(1-z^8)` and density `3/8`.

The prize center is instead

```text
c_t=bit_t(R_t).
```

The width-six hostile

```text
R_7=R_15 mod 64,     c_7=0, c_15=1
```

proves that this fixed quotient state does not determine the center.  The
center is the diagonal of an array whose columns are individually periodic.
This names the missing observer coordinate; it does **not** prove failure for
every width, which would already bear on the open nonperiodicity prize.

### Restoring the moving observer rather than discarding it

The edge map has an exact three-state, low-bit-first Mealy presentation:

```text
A=(A,B),       B=(C,B)sigma,       C=(A,B)sigma.
```

Here `A` is the packed Rule 30 map and `R_n=1 B^n(0^infinity)`.  Taking the
zero-section `n-1` times restores the center exactly.  With `S` denoting
strict suffix XOR, the resulting compiler is

```text
p_0=p_1=1^n,
p_(k+1)=S p_(k-1) or S p_k,
c_n=xor p_(n-1).
```

This costs `O(n^2)` bit operations.  It does not accelerate the known
computation by itself, but it localizes the obstruction: replacing OR by XOR
is a Lucas/Frobenius transport, and all remaining information lies in the
intersection masks `(S p_(k-1)) and (S p_k)`.  Scalar parity and low Hasse
marginals already fail at width three because they omit precisely that overlap.

There is a dual spacetime chart.  With `L=1+x+x^2`, the packed row polynomial
satisfies

```text
B_(t+1)=L B_t+D_t,
```

where `D_t` marks adjacent `11` collisions.  Rule 150 supplies a constant-one
center; every departure from that baseline is the Green-weighted parity of those
collisions.  A single Green weight has a two-carry-state binary digit compiler.
At dyadic times the sideways linear contribution vanishes entirely, so the
center is pure intersection current in both charts.

This changes the complexity question.  The linear transport is already
cheap; the live object is the labelled nonlinear current.  A shortcut must
compress it, while a lower bound must specify a machine model in which that
compression is impossible.

### Two exact restrictions short of the prizes

The edge-period lift word can never contain `111`.  Two consecutive lifts make
the two newest bit columns exactly uniform over a four-block phase cube, so
the next cocycle parity vanishes.  Hence

```text
P_w <= 2^(ceil(2w/3)-1),
```

improving the previous `2^(w-1)` ceiling for `w>=3`.  Runs ending at a horizontal step
parse into `0,10,110`, a faithful ternary run-length address with multipliers
`1,2,4`.  This resembles the prior three-spine work only at the representation
level; no Berggren arithmetic or branch symmetry is present.

At innovation depths, arrival-normalized traces across spatial phase form an
exact Walsh cube: every joint bit pattern occurs once.  That proves perfect
spatial balance and independence, yet the prize center is the single marked
phase `h=0`.  This is the same warning seen in LRC owner quotients: an exactly
uniform atlas does not control a distinguished address inside it.  It also
omits noninnovation depths; THM-3458's depth-five word supplies the `3/8`
unbalanced control after a harmless phase translation.

Finally, on the `n+1` checkerboard perturbations

```text
{n,n-1,n-3,...,1-n}
```

around the seed, the time-`n` center has full Boolean degree `n+1`.  Its exact
decision-tree complexity is `n+1`, and circuits with XOR/NOT free of charge
need at least `n` binary AND gates.  This is genuine robust-input hardness.
The one unperturbed seed point can still be hardwired, so a reduction from the
binary time index is the missing bridge to Prize 3.

### Sequence compiler that remains lawful

For a fixed cyclic output block `w`, the number of periodic Rule 30
predecessors of `w^k` is

```text
tr(M_w^k)
```

for an explicit `4x4` de Bruijn matrix `M_w`.  It is therefore C-finite of
order at most four and has an eigenvalue closed form.  Arbitrary one-step
inverse counts take `O(N)` fixed-dimension matrix multiplications in the
arithmetic-operation model; no linear bit-complexity claim is made.

This genuinely answers the user's efficient-closed-form theme, but for the
right sequence: repeated-block spatial inverse counts.  It is not an `n`th
center-bit algorithm.

## Pull 3 -- a lawful LRC prefilter and its stopping reason

For masks `A,B,C`, Rule 30 gives

```text
W_A=A xor(B union C),
P=A xor B xor C,
W_A xor P=B intersect C.
```

Three cyclic rotations recover all pair intersections.  Ternary owner codes
of depth `ceil(log_3 r)` therefore test pairwise disjointness of `r` masks in
three gates per digit; equality of the union with the ground set certifies a
partition.

This is a useful early sieve for typed cover search.  It is not a tournament:
overlap is a symmetric edge colour.  It is not an LRC exit either.  The two
`q=6` controls share masks `{0,3},{1,4},{2,5}` while disagreeing on star-
cochain existence.  The compiler must be followed by the affine common-time
equations.  Procedurally:

```text
mask-disjointness prefilter -> labelled star/cochain solver -> physical clock
```

is lawful; replacing the middle or final arrow by Rule 30 is not.

Periodic boundary conditions are an independent obstruction.  The all-one
word on an `N`-cycle has a Rule 30 predecessor exactly when `3|N`.  Hence the
13-slot all-safe word is unreachable in a naive one-step cyclic encoding.
Open fibres do not glue automatically.

## Pull 4 -- factorial moments without a factorial transfer

The unique real multilinear truth lift is

```text
g=l+c+r-cr-2lc-2lr+2lcr.
```

For the factorial functional,

```text
L(g)=0,       L(g^2)=6.
```

The first zero is real structure, not a counterexample.  The second power is
the cheapest decisive hostile.  More usefully, all moments admit a finite
derangement compiler:

```text
M_k=sum_j (-1)^j binom(k,j)(!j)^2,

L(g^m)=sum_q binom(m,q)q!
          sum_s binom(q,s)(-2)^s M_(m-q+s).
```

This avoids expanding `g^m`.  The ANF lift has a separate formal OGF

```text
sum u_n z^n=(1-z)^(-3) sum_d d![z/(1-z)^2]^d,
u_n=L(p^n)/n!,
```

and an order-four polynomial-coefficient recurrence.

The carrier conflict is exact.  Boolean truth wants `x^2=x`; the mod-two
factorial quotient wants `x^2=0`.  Identifying them collapses `x`, and `p=2`
is outside the good-prime live frontier.  The useful output is therefore a
new exact test family and sequence compiler, not movement on FC(3).

## Pull 5 -- complementary polynomial lifts and the JC boundary

The open-edge ANF lift is triangular over characteristic zero, has Jacobian
one, and is tamely invertible.  It is not faithful to real Boolean output: one
local value is `3` instead of `1`.

The unique faithful real multilinear lift has diagonal factors

```text
1-2(x_(k-1) or x_(k-2)),
```

so its determinant is nonconstant.  Closing the edge into a periodic ring
also destroys triangularity; the periodic ANF determinant takes different
values at zero and at all minus one.

Finally, in characteristic two the Frobenius-modified representative has
formal Jacobian `I` while remaining noninjective on the Boolean cube.  The
same truth function therefore supports:

| Boundary / lift | Truth fidelity | Jacobian status | Consequence |
|---|---:|---:|---|
| open edge, ANF over `C` | no | tame Keller | positive control only |
| open edge, multilinear over `R` | yes | nonconstant | naive faithful lift stops |
| periodic ring, ANF over `C` | no | nonconstant/vanishing | cyclic closure stops |
| periodic, Frobenius gauge over `F_2` | yes | formal `I`, noninjective | characteristic boundary |

No row is a JC(2) reduction.  The table is nevertheless useful because it
isolates three independent choices that a future bridge must type: boundary,
characteristic, and representative.

## Tournament check

There is no honest tournament hidden here.  The local Boolean influences of
Rule 30 distinguish the left input and tie the center/right inputs.  Mask
defects are symmetric intersections.  Forcing the tie into an orientation
would destroy the exact symmetry used by the two-boundary inverse.  The right
objects are a rooted preorder, a ternary gate, and a trace-cylinder tree.

## New procedural questions

The session leaves six precise next probes rather than a vague analogy.

1. **Unique Green channels.**  Find an infinite family in which one collision
   event reaches the center with no cancelling partner.  This is a concrete
   route toward nonperiodicity, not a randomness heuristic.
2. **All-depth marked-phase discrepancy.**  First compile the noninnovation
   readouts on the innovation cube, then compare `h=0` with their spatial
   averages.  Both steps are needed before the atlas can bear on density.
3. **Intersection-current complexity.**  Either compress the Mealy/spacetime
   overlap ledgers, or prove they resist compression in a uniform binary-index
   model with advice, preprocessing and bit cost stated explicitly.
4. **LRC two-stage pruning.**  Benchmark ternary mask gates as a fast owner-
   disjointness sieve, but retain the full labelled star equations and phase
   witness for every survivor.
5. **FC moment geometry.**  Use the exact `L(g^m)` compiler as a hostile family
   for proposed finite detectors; do not treat its first zero as persistent.
6. **Truth-lift classification.**  Classify polynomial representatives modulo
   the Boolean ideal that are Keller in characteristic zero.  The two
   canonical lifts are separated, but exotic representatives are not ruled
   out and would still need a planar target map.

## Reusable meta-pattern

The cross-domain lesson is:

```text
closed state dynamics + moving/typed observer does not imply a closed target sequence.
```

Trigger it whenever every fixed-depth or fixed-modulus observer is periodic,
rational, C-finite, or otherwise solved, while the target samples a depth,
owner, phase, or boundary that changes with the index.  The cheapest test is
to find two times with the same proposed compressed state but different
observer addresses and target values.  The Rule 30 pair `(7,15)` is now the
canonical small specimen for width six, not an all-width theorem.

Counterindication: if the observer is fixed and continuous on the closed
state, the odometer/finite-state analysis is lawful and may yield the exact
sequence.  That is why every fixed edge offset and every repeated-block
inverse count does admit a closed form.
