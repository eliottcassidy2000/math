# Three clocks, two cubic determinants, and the observer a bouquet forgets

**Status:** research synthesis and next-operation map.  The proof sources are
[THM-3357](../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md),
[THM-3358](../01-canon/theorems/THM-3358-admissible-composite-parabolic-compiler-and-hensel-normal-offset-atlas.md),
[THM-3359](../01-canon/theorems/THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar.md),
and [THM-3365](../01-canon/theorems/THM-3365-factorial-vandermonde-all-power-detection-benchmark.md).
The bouquet preservation facts are cited from Petrovic--Thoma--Vladoiu,
arXiv:1507.02740; the graph-ideal boundary is cited from
Garcia-Marco--Marquez-Corbella--Tatakis, arXiv:2504.06216.  This reflection is
not an additional proof source.  `LRC(14)`, `FC(3)`, and the proposed
moment-decorated bouquet constructions below remain open.

## Outcome first

The session separated three objects that had been travelling under the vague
name “finite-state structure”:

```text
word clock       state --letter--> state
index clock      n --+1--> n+1, observed through a finite residue state
value carrier    accepted word or sequence term --evaluation--> integer
```

Only the index clock is intrinsically unary.  A regular-language level count
becomes unary only after aggregating the letter transitions into one level
matrix; reading arbitrary words remains a semigroup action.  This is why a
modular C-finite index support has a rational harmonic coefficient, while a
regular ternary address-value set or the square values can have unbounded gaps
and finite reciprocal mass.

At the same time, two different nonlinear three-variable transfers produced
two exact factorial benchmarks:

```text
Berggren transfer determinant
P=(x-y-z)(x+y-z)(x+y+z)

Vieta coefficient-map Jacobian
J=(x-y)(x-z)(y-z).
```

Both lose every odd factorial moment and are detected at the second power,
but for different reasons.  `P` uses an ordered-triangle-to-square
factorization into two odd one-dimensional moments; `J` is the alternating
`A_2` collision arrangement and a transposition kills its odd moments.  The
hostile `xyz` proves that “a product of three linear forms” is not the
mechanism.

Finally, the toric audit identified exactly what algebraic compression does
and does not see.  The exponent toric ideal of `J` remembers a primitive
integer relation, but its bouquets are all singletons.  It cannot recover the
alternating coefficient character or the factorial functional.  The useful
next object is therefore not a bare toric ideal: it must be a relation carrier
decorated by its observer.

## Anchor / Niche / Wildcard portfolio

| lane | exact gain | missing coordinate |
|---|---|---|
| **Anchor -- LRC(14)** | THM-3357 gives a sufficient outer-to-middle gate rule in one fixed ordered rank-two deck; THM-3358 compiles every admissible exact-grade U-spine normal residue | owner, phase, physical Haar mass, seed/exit, and a global row |
| **Niche -- finite-state harmonic mass** | THM-3359 identifies the logarithmic coefficient with accepted cycle density | value support, address order, orientation, and collision multiplicity |
| **Wildcard -- cubic factorial/toric** | THM-3357 and THM-3365 give two all-order determinant moment laws | a classification of parity-killed cubic arrangements and an observer-preserving compression |

The wildcard overtook the original analogy because it produced both a theorem
and a stopping rule.  Bare bouquet algebra is not a shortcut to factorial
positivity; it tells us what extra data a lawful shortcut must carry.

## Concept board

| object | representation | invariant | native operation | quotient/loss |
|---|---|---|---|---|
| C-finite sequence | companion-matrix orbit mod `m` | eventual cycle and accepted positions | unary time shift | forgets integer values and recurrence lift |
| regular ternary language | DFA plus length transfer matrix | C-finite level count | append one letter / advance one level | level count forgets accepted addresses |
| Hensel-normal U-spine shell | `R_(N^2) x (Z/NZ)^x` | `2^k phi(N)` exact residues | prime toggle and normal-unit action | conjugation halves labels, not pulled-back index mass |
| Berggren triple transfer | `xT_L+yT_M+zT_R` | determinant `P` | weighted branch superposition | determinant forgets word and owner |
| Vieta root transfer | `(x,y,z)->(e1,e2,e3)` | Jacobian `J`, discriminant `J^2` | root permutation | coefficient triple forgets root order/sign |
| exponent toric carrier | six columns of `A_J` | kernel, circuits, Graver basis | conformal lattice addition | forgets polynomial coefficients and factorial observer |

In the U-spine row, `N=product_j p_j^(e_j)` is admissible odd with every
`p_j=1 mod 4`, and `k=omega(N)`.  Every new connection below is tested against
all six objects rather than only
against the syntax it shares with one of them.

## Pull 1: cycle density belongs to indices, not values

For a C-finite sequence modulo `m`, the length-`r` state vector lives in a
finite set.  Its unary orbit is eventually periodic.  If `h` of `p` eventual
states are accepted, then

```text
#(H intersect [1,X])=(h/p)X+O(1),
sum_(n<=X,n in H) 1/n=(h/p)log X+C+O(1/X).
```

This mechanism applies unchanged to tournament walk-count sequences and to
modular predicates on level counts of a regular ternary language.  It does not
say that the integers represented by the accepted words are periodic.

The same-language control is decisive.  For `L=1{0,1}^*`, the length-`n`
count is `2^(n-1)`, whose residue-one support modulo three is the odd lengths
and has coefficient `1/2`.  Yet the represented base-three value set has only
`2^(n-1)` points in `[3^(n-1),3^n)`, so its block reciprocal mass is at most
`(2/3)^(n-1)`.  It is harmonically summable and has unbounded gaps.

The square hostile is even sharper.  The sequence `a_n=n^2` is C-finite, but
its **value** support has reciprocal sum `zeta(2)`.  The theorem instead sees
modular **index** supports such as `{n:n^2=0 mod 5}=5N`, with coefficient
`1/5`.  “C-finite” modifies a time evolution; it does not make the image of
that evolution ultimately periodic.

## Pull 2: the normal atlas turns state count into harmonic mass, and nothing more

For the admissible odd grade `N` and `k` just defined, THM-3358 identifies the
exact-grade U-spine normal residues modulo `N^2` with

```text
R_(N^2) x (Z/NZ)^x,
```

of cardinality `2^k phi(N)`.  Consequently any subset `B` of labels has
positive-index harmonic coefficient

```text
delta_B=|B|/N^2.
```

This is the cleanest bridge between the prime-toggle atlas and THM-3359: a
Dirichlet residue becomes normalized finite-state count.  It is also a sharp
loss invoice.  Two label subsets of equal size have the same residue even if
their prime allocations, normal units, Gaussian lifts, and Berggren words are
disjoint.  Passing to conjugation-parent labels halves the label set, but the
integer-time pullback contains both antipodal residues, so its mass is not
halved again.

Changing the observable exposes the boundary.  Along the same ray,
`C_t` grows quadratically, and `sum 1/C_t` converges.  Periodic index density
therefore supplies no physical LRC mass and no owner.  The exact compiler is a
language theorem; the harmonic residue is a clock statistic; neither is a
positive-Haar LRC row.

## Pull 3: two cubic determinants, two parity mechanisms

For the Berggren triple transfer, THM-3357 proves

```text
calL(P^r)=0                              (r odd),
calL(P^r)=(3r+2)!/[2(r+1)^2]            (r even).
```

The three factors of `P` are linearly independent.  Radial/simplex
coordinates reduce the normalized integrand to `(uv)^r` on the ordered half
of `[-1,1]^2`.  Swap symmetry makes this half the full-square integral; that
integral factors, and either odd one-dimensional moment vanishes.

For the Vieta map, THM-3365 proves

```text
calL(J^(2m+1))=0,
calL(J^(2m))=((2m)!)^2(3m)!/m!.
```

The three factors of `J` have rank two and share the diagonal translation
direction.  Their square is the cubic discriminant.  Here a root
transposition preserves the exponential orthant and reverses `J`; an ordered
chamber then reduces the even integral to one beta integral.

| feature | `P` | `J` |
|---|---|---|
| origin | determinant of three branch-transfer matrices | Jacobian of roots-to-coefficients Vieta map |
| factor-form rank | three | two |
| odd cancellation | ordered half-square factorization | `S3` alternation |
| first detector | `calL(P^2)=2240` | `calL(J^2)=24` |
| geometric zero set | three independent planes | three root-collision planes sharing the diagonal |
| FC consequence | fails premise at power two | fails premise at power two |

There is already an exact positive family around `P`.  Off the measure-zero
switch wall `x=z`, the two unimodular branches

```text
T(x,y,z)=(x+y,z-x,x)       when x<z,
T(x,y,z)=(z,x-z,y+z)       when x>z                          (B)
```

are inverse, preserve the positive orthant and `x+y+z`, and exchange the two
open cones.  Give the switch wall any fixed measurable tie convention; it has
zero Gamma measure and is not used in the moment identity.  Almost everywhere
they satisfy `P(Tx)=-P(x)`.  More generally,

```text
F_(b,c,k)=k(x+y+z)[b y+c(z-x)][b(z-x)-c y]
```

obeys `F(Tx)=-F(x)` and

```text
L(F_(b,c,k)^2)=224 k^2(b^4+8b^2c^2+c^4)>0              (C)
```

for real `k!=0` and `(b,c)!=(0,0)`.  The original `P` is `(b,c,k)=(-1,-1,1)`.
Equations (B)--(C) were independently symbolically reduced to zero in this
session; the proposed completeness of this family for the piecewise mechanism
still deserves a separate proof audit before canonization.

This pair suggests a bounded classification problem: classify real products
of three homogeneous linear forms whose exponential-orthant moments vanish
for every odd power.  Only coordinate permutations preserve the fixed
factorial measure among linear orthant automorphisms; a common scalar on the
polynomial merely rescales moments, while an independent positive diagonal
change of variables must transport the Gamma measure as a sidecar.  Indeed
`x->2x` changes `L(P)` from zero to `36`.  The cheapest other hostiles are
`xyz` (no odd vanishing), a complex scalar (no ordered positivity), and a
generic factor perturbation (breaks cancellation).

## Pull 4: the same periodic tail can have a different proof mechanism

For every fixed modulus `M`, both closed even-moment formulas are eventually
zero modulo `M`.  More precisely, for

```text
U_m=L(P^(2m)),        V_m=L(J^(2m)),
```

the prime valuations are

```text
v_p(U_m)=v_p((6m+2)!)-[p=2]-2v_p(2m+1),
v_p(V_m)=2v_p((2m)!)+v_p((3m)!)-v_p(m!),               (A)
```

and both tend to infinity.  Their nonzero modular supports are therefore
finite and have harmonic coefficient zero.  This tail is itself ultimately
periodic; what differs from THM-3359 is the certificate: factorial valuation
extinction, not a fixed companion-state orbit.  The even subsequences are
hypergeometric, the parity-interleaved full sequences are P-recursive, and
their superexponential growth rules out C-finiteness.

For an odd prime `p`, the exact all-later extinction thresholds on the even
index `m` begin at `(p+1)/2` for `U_m` and `ceil(p/3)` for `V_m`.  These clocks
do not give a uniform ranking at prime powers: carry corrections can reverse
their order.  This is a precise state-plus-valuation refinement, not a new
C-finite companion.

This gives a useful diagnostic of four proof/growth mechanisms:

```text
eventual periodicity from a finite state orbit
eventual modular extinction from growing valuations
sparse divergent mass from slowly growing gaps
sparse convergent mass from exponential/quadratic value growth.
```

For a general residue predicate modulo `M`, after extinction the support is
cofinite with coefficient one if zero is accepted and finite with coefficient
zero otherwise.  A harmonic coefficient alone cannot identify which proof
mechanism produced its periodic tail.

## Pull 5: two toric models of one sign polynomial do different jobs

There are at least two natural relation matrices around a Vandermonde
polynomial.

1. The **monomial-exponent matrix** used in THM-3365 has columns
   `210,201,120,102,021,012`.  Its alternating sign vector is Graver but not a
   circuit, and its Gale rows give six singleton bouquets.  It records exponent
   balance and gives no compression.
2. A **signed tournament score-incidence matrix**, derived here from
   THM-1805's score map, records edge choices and equal-score moves.  At `K3`,
   the cyclic difference is one mixed graphic relation.  It retains a Boolean
   tournament fibre only if that fibre and the parity character are attached
   separately.  This directed-incidence toric ideal is not the unsigned
   simple-graph toric ideal of the 2025 paper.
3. The Berggren determinant `P` has the full ten-monomial cubic Veronese
   support.  Its coefficient vector is not even an exponent-kernel relation:
   `A_P c_P=(-1,-7,-1)`.  A rank-seven Gale basis has no zero or proportional
   rows, so these bouquets are again singleton.  Replacing its coefficients
   by their absolute values keeps the support/toric ideal but changes the
   first factorial moment from zero to `32`.

The same polynomial therefore need not have a canonical “toric model.”  One
matrix sees exponent balance; the other sees score-preserving edge moves.
Neither alone sees the factorial functional.  Replacing the alternating
coefficients by all plus signs preserves the first toric ideal but changes
`calL` from zero to twelve at the first moment.  Replacing the coefficient
character on the tournament fibre similarly turns cancellation into counting.

The cited bouquet theorem is exact but narrower than the hoped-for transfer:
integer kernels, circuits, and Graver bases always correspond under bouquet
encoding; Markov and Grobner data require the stable/non-mixed scope.  Neither
statement transports an external nonlinear observer.  The 2025 graph-toric
MG criterion is additionally a graph-incidence theorem, with its
equal-cycle-length implication in the bipartite scope.  It is not a theorem
about the exponent matrix above or about chronological automata.

## Pull 6: a lawful automaton toric ideal must pay for time

The tempting construction sends every automaton transition `q->r` to the
edge `q_src--r_tgt` of a bipartite graph.  Its toric kernel preserves source
and target margins, but a DFA may have parallel letter-labelled transitions;
the 2025 theorem concerns ordinary simple graphs.  Collapsing to simple
support bills the lost letter identity and multiplicity.  In either version,
chronological composability is lost: `L->M->R->L` becomes three disjoint
matching edges, whose toric kernel is exactly zero.

Two partial encodings expose what must be paid:

- signed incidence with affine right side `Bx=e_t-e_s` retains boundary flow,
  but still admits path-plus-cycle or disconnected flows and forgets order;
- a time-layered network can retain a length-`h` path only after adding
  start/end right sides, nonnegative unit-flow/capacity conditions, and the
  transition labels; its matrix grows with the horizon.

A noncommutative path/category encoding is a third possibility.  None of
these conclusions follows from the bare toric kernel alone.

The cheapest decisive test for any proposed finite-state bouquet compiler is
therefore not its Hilbert basis.  It is two equal-margin multisets with
different word order, or one composable path versus one disconnected or
noncomposable multiset.  Terminal-only acceptance differs only when the
start/end or acceptance sidecar has also been omitted.

## New project board

### A. Classify parity-killed cubic arrangements

Search products `Q=l_1l_2l_3` under the symmetry group of the exponential
orthant/simplex.  Record factor rank, oriented matroid, and a measure-compatible
cancellation witness (global involution or separable half-square factorization), and the
hypergeometric ratio `calL(Q^(2m+2))/calL(Q^(2m))`.  Positive controls are
`P` and `J`; `xyz` is the first hostile.  The target is a classification of
the cancellation mechanism, not a search for an FC counterexample.

### B. Build a moment-decorated toric carrier

Start from the exponent configuration but attach:

```text
coefficient character + factorial/Gamma observer + power index.
```

Test whether a stable bouquet map intertwines the finite Macaulay moment
matrices after these decorations.  The all-plus/Vandermonde pair is the
minimal negative control.  A positive result should be stated as an
observer-intertwining theorem, never as a property of the toric ideal alone.

### C. Compare the two determinant sequences modulo primes

Compute minimal extinction thresholds, valuation profiles, and generating
recurrences for the `P` and `J` moment sequences.  Compare these with C-finite
cycle bounds to separate “finite because periodic” from “finite because
factorial.”  A useful invariant may be the smallest state-plus-valuation
automaton rather than a companion matrix alone.  A zero at one index need not
be absorbing for `U_m` (already modulo five), so monotone extinction thresholds
must be proved rather than inferred from a first zero.

### D. Time-layer the prime-toggle compiler

For one small composite grade, unfold the exact Berggren compiler into a
layered transition-margin matrix and audit its bouquets.  Ask whether the
  normal-unit coordinate becomes a stable bouquet or remains a coefficient
  sidecar.  The required hostile is two equal-margin transition multisets with
  different word order, one of them noncomposable.

### E. Keep the LRC consequence honest

The only current LRC gain is pruning: THM-3357's sufficient Horn rule in one
fixed ordered deck and THM-3358's exact U-spine compiler reduce candidate
structure.  Harmonic index density, factorial
positivity, Graver minimality, and bouquet stability do not supply owner,
phase, physical mass, or a global exit.  Any proposed return to LRC must first
exhibit a map from one of these carriers to a labelled determinant deck and
then pass an owner-switch hostile.

## Stopping certificates

- **No syntax transfer:** “C-finite,” “regular,” and “finite automaton” do not
  make a represented value set ultimately periodic.
- **No factor-count transfer:** three linear factors do not force odd moment
  vanishing.
- **No toric-observer transfer:** the same exponent toric ideal can support
  factorial moments zero and twelve.
- **No graph-scope transfer:** a bipartite graph-ideal theorem does not apply
  to a non-graph exponent configuration or a chronological transition system.
- **No harmonic-to-physical transfer:** normalized residue count is not LRC
  Haar mass.

These failures are constructive.  Each names the exact missing coordinate:
value map, measure-compatible cancellation witness, coefficient/Gamma
observer, chronological path witness, or owner/phase height.  That is the
useful frontier left by the session.
