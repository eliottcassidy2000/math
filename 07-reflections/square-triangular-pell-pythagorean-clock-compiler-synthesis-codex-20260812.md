# The parabolic spine, its Pell selector, and the clocks it cannot remember

**Status:** research synthesis and next-operation map.  The proof source is
[THM-3335](../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md),
with the ambient Gaussian carrier in
[THM-3333](../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md),
the parabolic Berggren operation in
[THM-3334](../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md),
the second Pell selector and Gaussian-square transplant in
[THM-3341](../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md),
and the golden/Cassini contrast in
[THM-3339](../01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md).
This reflection is not an additional proof source.  `LRC(14)`, `JC(2)`,
`FC(3)`, and existence at the new skew-EW candidate orders remain **OPEN**.

## Outcome first

The session found a useful compiler architecture rather than another scalar
coincidence:

```text
full unit-gap Berggren spine
  Phi(n+1,n)=(2n+1,4T_n,4T_n+1)
                 |
                 | select T_n=q^2
                 v
Pell state (x,q), x^2-8q^2=1
  |              |                 |
  | half-Hadamard|                 | consumer maps
  v              v                 v
fixed-2 Markov   square-leg triple tournament/LRC/EW addresses
(a,b)            (x,(2q)^2,(2q)^2+1)
```

The selector is lossless through the Pell, consecutive-spinor,
square-triangular, and ordered fixed-two Markov faces.  It becomes lossy only
when a downstream consumer forgets labels: a tournament arc count forgets
orientation, a Gaussian norm quotient forgets represented angle, the
compressed tuple `(Phi,C,D,owner,F)` forgets labelled clock placement, and
even-power decimation forgets the Cassini parity used by the planar-JC
metallic lane.

This architecture produced five concrete gains:

1. a proved two-state `O(log k)` operation-count compiler for the square-
   triangular indices and their Pythagorean/Markov/EW coordinates;
2. a global separation theorem: nontrivial square universal tournament arc
   count is incompatible with skew-EW attainment;
3. seven exact LRC benchmark families with maxima `1/7,...,1/13` but identical
   compressed Gaussian/Kelvin scalar state on six labelled planes; and
4. a global identification of the `70` cannonball event as the unique positive
   square-pyramidal intersection of the selector, after importing Bennett's
   classification; and
5. a second, transverse Pell view in which every non-root square-triangular
   row joins consecutive square-hypotenuse nodes of the same Berggren spine;
   the root row joins its algebraic boundary to the first physical node.

## Anchor / Niche / Wildcard portfolio

| Lane | Inherited object | What moved | Honest boundary |
|---|---|---|---|
| **Anchor — LRC(14)** | THM-2056 determinant/Kelvin gate and THM-3333 consecutive spinors | one spinor ray now feeds a seven-rung exact clock ladder and a six-plane quotient hostile | every row is already clock-safe; no open ledger row closes |
| **Niche — primitive triples and sequences** | Euclid/Gaussian square, square-triangular Pell equation, THM-1880 evaluator | one lossless Pell/Markov/Pythagorean compiler, complete orbit descent, rational GFs, and a sparse Berggren selector | square leg without unit gap is a much larger two-sheet locus |
| **Wildcard — sums of squares / tournament orders** | cannonball `70`, even/odd square split, THM-476 EW gate, universal arc count | the exceptional row is globally unique; two superficially similar order maps are proved distinct | arithmetic EW candidates are not designs; arc count has no orientation data |

The niche supplied the anchor's best hostile control rather than a proof of
the anchor.  That is still progress: it identifies exactly which coordinate a
Gaussian/Kelvin replay must retain before it can predict a phase maximum.

## Inheritance pass

| Required item | Closest controlled source |
|---|---|
| proved mechanism | THM-3333: `Phi(m,n)` and `<Phi(u),Phi(v)>_L=2 det(u,v)^2` |
| native operation | THM-3334: the parabolic Berggren `U` step on the complete `C-B=1` spine |
| canonical hostiles | `(5,12,13)` has unit gap but no square leg; `(63,16,65)` has square leg but gap `49`; norm `65` has inequivalent representations |
| corrected near miss | THM-1880's swapped coupled recurrence, repaired by MISTAKE-367; THM-2142 had the correct crossed formula |
| least-used sidecar | HYP-3075's fixed-two Markov/cannonball row `(29,70,169)`, previously address-only |
| literature upgrade | Bennett's primary-source application of Theorem 2.1: the only positive square-pyramidal squares are `1` and `4900` |

THM-900 is inherited only in its corrected Guy-identity scope.  Its historical
both-clean search conclusion was superseded and is not a dependency.

## Live concept board

| Object | Representation | Invariant | Native operation | Destroyed coordinate / next test |
|---|---|---|---|---|
| parabolic spine node | spinor `(n+1,n)` / triple `Phi(n+1,n)` | primitive, `C-B=1` | Berggren `U` | square-leg predicate; apply Pell selector |
| selector state | `(x_k,q_k)` | `x^2-8q^2=1` | matrix `[[3,8],[1,3]]` | none internally; retain depth under consumers |
| square-hypotenuse state | `C_t=m^2`, `X^2-2m^2=-1` | negative-Pell orbit | multiply by `3+2sqrt(2)` | scalar square loses Berggren ancestry fibre |
| fixed-two branch | `(a,b)=(x-2q,x+2q)` | `a^2+b^2+4=6ab`, `ab=4q^2+1` | Vieta mutation `(a,b)->(b,6b-a)` | unordered product loses mutation address |
| square-sum row | `s=2q`, `s^2=sum r^2` at `s=70` | even/odd split and cannonball height | elliptic/quartic classification | height `N`; scalar `s` alone is not an orbit law |
| tournament order | `R=n+1` or `W=4T_n+2` | square arc count or EW square gate | change order / construct sign matrix | all orientations and Gram/PAF data |
| LRC row | labelled plane `(j,n)` | exact `M=1/j` | move tail placement `j` | `(C,D,F,owner)` forgets the clock label |
| metallic lane | full Pell index versus even subsequence | Cassini sign | multiply by `1+sqrt(2)` or its square | even decimation kills alternating parity |

The board's most important distinction is **selector versus consumer**.  A
selector can be perfectly invertible while every consumer map downstream is
many-to-one.

## Pull 1: square-triangular numbers select square legs on one native spine

For every `n>=1`,

```text
Phi(n+1,n)=(2n+1,2n(n+1),2n^2+2n+1)
           =(2n+1,4T_n,4T_n+1).                         (1)
```

The identity `C-B=(m-n)^2` makes this the complete positive primitive family
with ordered even leg and `C-B=1`; it is not an arbitrary subfamily of
Pythagorean triples.  THM-3334 supplies the native operation

```text
U Phi(n+1,n)=Phi(n+2,n+1).                               (2)
```

The extra predicate `T_n=q^2` is equivalent to

```text
x=2n+1, s=2q,       x^2-2s^2=1,
Phi(n+1,n)=(x,s^2,s^2+1).                               (3)
```

Hence square-triangular solutions are a sparse selector on the parabolic
Berggren spine:

```text
n=1,8,49,288,1681,...,
U^(n-1)(3,4,5)=(3,4,5),(17,144,145),(99,4900,4901),... . (4)
```

This clarifies two graph scales.  The full spine advances by one Berggren
step and lies in the Farey fan about `(1,1)`.  Successive selector nodes skip
`7,41,239,...` spine edges.  They remain distance two through the common Farey
cusp, but their determinants and Lorentz pairings grow without bound.  An
unweighted Farey distance is therefore a deliberately lossy visualization,
not the selector metric.

The ambient square-even-leg locus is larger.  For primitive Euclid parameters,

```text
2mn is a square
iff (m,n)=(u^2,2v^2) or (2u^2,v^2).                     (5)
```

The unit-gap condition intersects these two valuation sheets in the two Pell
forms.  This explains why neither “Pythagorean” nor “square leg” alone is a
useful enough search predicate.

## Pull 2: the Markov and sequence faces are genuinely computational

The half-Hadamard coordinates

```text
a=x-s,              b=x+s                               (6)
```

give

```text
a^2+b^2+4=6ab,      ab=s^2+1.                           (7)
```

Mutation `(a,b)->(b,6b-a)` is conjugate to multiplication by the Pell unit:

```text
[x_(k+1)]   [3 8][x_k]
[q_(k+1)] = [1 3][q_k].                                 (8)
```

The inverse-unit descent `(x,q)->(3x-8q,3q-x)` strictly lowers `q`, so (8)
enumerates every positive solution rather than merely generating examples.
Binary matrix powering computes row `k` in `O(log k)` exact multiplications;
the output integers have `Theta(k)` bits.

The same state compiles

```text
n_k=(x_k-1)/2,
B_k=4q_k^2,
C_k=B_k+1,
W_k=B_k+2,                                               (9)
```

with rational generating functions and scalar recurrences of characteristic
constants `6` and `34`.  This is the useful meaning of “closed form” here:
the formula retains a small state, exact operation, inverse descent, and typed
outputs.  A scalar list without its state transition would be much less useful.

## Pull 3: even and odd square sums locate one globally exceptional row

At selector depth `k=3`,

```text
(a,b)=(29,169),
(x,s^2,C)=(99,4900,4901),
29*169=4901=70^2+1.                                    (10)
```

The square has the exact parity split

```text
70^2=sum_(r=1)^12 (2r)^2 + sum_(r=1)^12 (2r-1)^2
    =2600+2300.                                         (11)
```

Equivalently,

```text
99^2 + (sum_(r=1)^24 r^2)^2 = 4901^2.                  (12)
```

Bennett's cited classification of the cannonball equation proves that
`sum_(r<=N)r^2` is a square only at `(N,s)=(1,1),(24,70)`.  Because positive
selector `s_k` are even, (10) is the unique positive intersection of the
selector with the square-pyramidal locus.  The old scalar coincidence is now
a globally isolated address in an infinite exact compiler.

A different intersection remains open in this session: `q_k` itself being
triangular.  The exact `k<=30` scan sees only `q=1=T_1` and `q=6=T_3`.  It
should be treated as an exponential/elliptic intersection, not inferred from
the cannonball classification.

## Pull 4: every non-root selector row stitches two physical square-hypotenuse depths

THM-3341 reveals a transverse use of the same Markov coordinates.  On the
full unit-gap spine put

```text
C_t=2t^2+2t+1.
```

The square values `C_t=m^2` form the negative-Pell orbit

```text
t=0,3,20,119,696,...,         m=1,5,29,169,985,... .   (12a)
```

Here `t=0` is the algebraic boundary; the physical positive-parent U-spine
starts at `t=1`.

If `T_N=R^2` is a square-triangular row, define

```text
M_-=2N+1-2R,       M_+=2N+1+2R,
t_-=2R-N-1,        t_+=2R+N.                           (12b)
```

Then `C_(t_-)=M_-^2`, `C_(t_+)=M_+^2`, and

```text
M_- M_+-(2R)^2=1.                                      (12c)
```

These are consecutive negative-Pell roots and precisely the adjacent
fixed-two Markov coordinates.  The root row `T_1=1` joins the boundary
`t_-=0` to the first physical square node `t_+=3`; every later selector row
joins two physical U-spine nodes.  At `T_49=35^2`, the edge is

```text
(t_-,t_+)=(20,119),       (M_-,M_+)=(29,169),
29*169-70^2=1.                                           (12d)
```

So the cannonball value `70` is not merely attached to one Pythagorean node:
it is the Cassini carry across an edge between two square-hypotenuse nodes.
Gaussian squaring maps the positive/middle Berggren ray onto these sparse
nodes, but with variable depth drift; it is a state-dependent transplant, not
a branch homomorphism.  Moreover the equal-hypotenuse Boolean fibre ranks stay
unbounded inside this selector, so retaining the square value still loses
ancestry.

The companion scalar `Q_t=2C_t+1` is never square and is triangular on two
norm-`17` Pell orbits.  Those triangular values provide efficiently generated
tournament-size labels, but no orientation or design.

## Pull 5: the same selector feeds two different tournament order maps

Every tournament on `R` vertices has `T_(R-1)` arcs.  Thus

```text
R_k=n_k+1=2,9,50,289,1682,...                           (13)
```

are exactly the orders at which the universal arc count is a square.  This
classification retains no orientation or isomorphism data.

The skew-Ehlich--Wojtas arithmetic address is instead

```text
W_k=s_k^2+2=2,6,146,4902,166466,...,
2W_k-3=x_k^2.                                           (14)
```

It satisfies the proved necessary square gate from THM-476; it supplies no
design.  The two order maps intersect at a structural no-go: if an order has
both square universal arc count and skew-EW attainment, coprime factorization
forces the order to be `2`.  This is global, not a bounded search.

The first larger selector search address beyond `W=6` is `W=146`, with header
`(17,144,145)` and block size `73`.  A useful next computation would
search lawful skew-EW sign/PAF structures over `Z/73`, explicitly recording
which multiplier ansatz is tested.  Failure of a restricted ansatz must not be
reported as nonexistence.

Tournament Analysis is appropriate only after an intrinsic pairwise relation
appears.  The universal arc count is not such a relation; it is an order-level
scalar.  THM-3339 supplies a legitimate comparison tournament on six weighted
edge products, but that is a different object and must not be identified with
the orders in (13) or (14).

## Pull 6: one spinor ray exposes a seven-clock LRC quotient loss

For `7<=j<=13`, `n>=1`, and `a=n+1`, put

```text
S_(j,n)=a({1,...,13}\{j}) union {(j+1)a-1}.              (15)
```

These are seven labelled one-tail planes driven by the same spinor
`d_n=(n+1,n)`.  Exact phases prove

```text
M(S_(j,n))=1/j.                                         (16)
```

The boundary is sharp: at `j=6,n=1`, the surviving core multiple `12a`
destroys the six-clock extremizer and the exact maximum is `2/23`, not `1/6`.

For `j=7,...,12`, all rows at fixed `n` share the compressed tuple

```text
Phi(d_n), C(d_n), D(d_n)=13n, owner=-c_13,
F(d_n)=2n^2-1181n+1,                                   (17)
```

yet their true maxima are `1/7,...,1/12`.  Their full labelled column
configurations and polar polygons differ.  Therefore `(C,D,F,owner)` is an
excellent replay state inside one fixed plane but is provably incomplete
across planes.  The missing coordinate is the labelled tail/clock placement
`j`.

The Pell selector meets the Kelvin residual at precisely

```text
n=1,8,49,288,                                            (18)
```

on every rung, but all four are closed by the clocks in (16).  This is a
positive/hostile benchmark for future LRC state compressors, not progress on
an open physical residual.

There is also a cheap stopping theorem: each of the nine positive streams

```text
n_k,q_k,x_k,s_k,B_k,C_k,W_k,a_k,b_k                    (18a)
```

is at least `7/3`-lacunary, so any thirteen terms from one named stream are
already LRC-safe by THM-928(A').  This excludes neither the undecimated Pell
sequence nor arbitrary mixed-coordinate packets.  Hard selector use must mix
coordinates or attach the address to an independently distinguished physical
path; a one-stream packet from (18a) cannot be a counterexample family.

## Pull 7: a decimation no-go for the planar-JC metallic route

THM-3339's golden/Fibonacci carrier keeps the full Cassini alternation and
turns its sign into an oriented channel-order sidecar.  The square-triangular
compiler instead takes even powers of the silver unit:

```text
1+sqrt(2)           has norm -1,
(1+sqrt(2))^2=3+2sqrt(2) has norm +1.                   (19)
```

Consequently the full Pell Cassini sign `(-1)^r` becomes the constant

```text
q_(k-1)q_(k+1)-q_k^2=-1                                (20)
```

on the selector.  Every defined Newton ratio is then greater than one; the
maximal alternation used by the planar-JC metallic circuit stratification has
been erased.  The missing sidecar is the parity of the original, undecimated
index.  This is a precise no-go: square-triangular projection is useful for
sequence compilation but is the wrong quotient for that JC mechanism.

The same lesson appears in the factorial lane.  THM-3300's torus-invariant
factorial algebra retains `C=|z|^2` and kills the weight-two represented-sum
coordinate `z^2=A+iB`.  The norm collision

```text
8^2+1^2=7^2+4^2=65                                    (21)
```

has determinant shells `1` and `4` against the same endpoint.  A direct
Pythagorean-to-`FC(3)` transfer therefore needs an `(A,B)` sidecar outside the
factorial quotient.  This direct THM-3333 represented-sum/Farey port through
THM-3300's norm quotient cannot move `FC(3)` without it.

## Typed transfer ledger

| Source | Target / map | Preserved | Destroyed | Required sidecar or cheapest control |
|---|---|---|---|---|
| unit-gap Berggren node | Pell selector | primitive triple, gap one, square leg at selected depths | intervening spine nodes | retain `n_k` / Berggren depth |
| Pell state | fixed-two Markov pair | conic, mutation, product/hypotenuse | branch order if pair unordered | seed plus mutation word |
| selector state | exact sequences | recurrences, GFs, random access | consumer meaning if outputs untyped | keep state/output schema |
| selector root | cannonball equation | the `70` row | general orbit structure | Bennett classification plus height `24` |
| square-triangular row | adjacent square `C_t` roots | Markov pair and Cassini carry | intervening U-spine and ancestry fibres | negative-Pell depth plus Berggren word |
| middle Berggren ray | square-hypotenuse U-spine nodes | primitive triple under Gaussian squaring | constant branch length / homomorphism | state-dependent target depth |
| selector | EW order `W` | necessary square gate | sign matrix, Gram, PAF | exact construction/ansatz ledger |
| tournament order `R` | universal arc count | `T_(R-1)` | every orientation | never infer structure from count alone |
| labelled LRC plane | compressed Kelvin tuple | norm, determinant maximum, owner, defect | tail label and phase | six-plane `1/7,...,1/12` hostile |
| full silver orbit | even selector orbit | magnitude, norm `+1` recurrence | original parity/Cassini alternation | retain full index parity for JC |
| represented Gaussian spinor | factorial norm quotient | radial norm | angle, represented sum, determinant | weight-two `(A,B)` sidecar |

## What directly moved, what only narrowed

- **PROVED:** the selector equivalences, orbit exhaustion, valuation sheets,
  recurrences/GFs, Berggren sparse-selector law, global square-arc/skew-EW
  incompatibility, seven-clock all-`n` theorem, adjacent negative-Pell-root
  law, and state-dependent Gaussian-square transplant.
- **CITED + DEDUCED:** Bennett's cannonball classification makes the selector
  intersection globally unique at `(k,N,s)=(3,24,70)`.
- **VERIFIED-EXACT:** the stored companion independently checks direct
  square-triangular hits through `100000`, all primitive Euclid parameters
  through `m=500`, both exact LRC max engines, and the stated positive and
  hostile controls.
- **NARROWED:** one-stream selector packets are outside the hard LRC frontier;
  even decimation is outside the JC metallic-alternation route; this direct
  Pythagorean norm-only factorial transfer is impossible.
- **OPEN:** an EW construction at `146`; the all-range `q_k=T_j`
  intersection; any lawful mixed-coordinate LRC import; any reconstruction of
  a factorial weight-two sidecar; `LRC(14)`, `JC(2)`, and `FC(3)` themselves.

## Best next operations

1. **EW address `146`:** compile a transparent hierarchy of two-circulant,
   multiplier, and unrestricted sign/Gram constraints over `Z/73`.  Each empty
   layer should state exactly what larger construction space remains.
2. **LRC sidecar compiler:** on actual THM-2056 owner cusps, record
   `(Pell depth, physical horocyclic index, labelled tail, owner/tie, F, clock
   histogram, endpoint word)`.  Require an interval certificate before
   extrapolating away from selected rays.
3. **Mixed-coordinate packets:** combine coordinates from different selector
   streams or splice a selector address into a nonlacunary physical core.  Pure
   streams are already closed, so this is the first nonvacuous LRC experiment.
4. **Triangular-root intersection:** treat `q_k=T_j` as an elliptic/exponential
   Diophantine problem.  The finite hits `1,6` are controls, not a conjectural
   proof.
5. **Consumer-aware sequence API:** expose the matrix state, inverse descent,
   rational GFs, and typed output maps together.  This is the useful reusable
   closed-form object; publishing only scalar prefixes would discard the
   operations that made the synthesis work.
6. **Cross-selector ancestry atlas:** compile each square-triangular edge
   together with its two negative-Pell depths, Markov carry, Gaussian-square
   source, and full fixed-hypotenuse Boolean fibre.  The first question is
   whether fibre-grade changes admit a finite-state sidecar; scalar square or
   triangular labels alone demonstrably do not.

The main procedural lesson is sharp: search first for a native orbit, then for
a sparse selector, and only then attach consumers.  Every consumer must carry
its own loss ledger.  That ordering turned several old coincidences into one
proved compiler while preventing the compiler from being overclaimed as an
LRC, tournament, factorial, or Jacobian theorem.
