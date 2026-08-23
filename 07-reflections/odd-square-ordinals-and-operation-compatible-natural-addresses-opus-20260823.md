# Odd-square ordinals and operation-compatible natural addresses

**Session:** opus / 2026-08-23  
**Anchor:** odd-square ordinals on primitive Pythagorean triples  
**Niche:** triangular addresses and centered differences on the live LRC
profile/detector interfaces  
**Wildcard:** operation-compatible ranks for tournaments, ternary walks,
quadratic graph heights, and Poisson exponent lattices

> **Truth boundary.**  This is a research reflection and task generator, not a
> proof source.  At the time of this checkpoint, THM-3756 is a
> `RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT`; use its final
> frontmatter rather than this chronology.  It was subsequently promoted to
> `PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED` after a separate
> census found zero converse failures among all `792` primitive ordered
> triples with `C<=5000` and after both reproduction defects found by the
> audit were repaired.  The exact cross-frontier probe is
> `04-computation/ordinal_encoding_operation_probe_opus_20260823.py` with
> semantic digest
> `7a9920d83f94b98cd962d8418fb390ca37a4c001604f2887153b65f3fe20b60e`,
> raw-LF script SHA256
> `7b6eda61e0352eebbcd392f7d02a21c10c55f24d4039f2e5efde75fd11f13909`,
> and raw-LF output SHA256
> `80eb36e9e58683576438b94cfa5acc16b29b1d0f59de26b2be6ecbc635f63051`.

## Outcome first

The motivating outer rank is mathematically real and is now proved in
THM-3756, `odd-square-ordinal-berggren-affine-descent`:

```text
B+C=(2r-1)^2  ->  call the selected odd square r.        (1)
```

It is not a tree address.  The hidden second coordinate is of exactly the
same type:

```text
C-B=(2s-1)^2.                                           (2)
```

Together `(r,s)` give a lossless chart of primitive Pythagorean triples,
turn all three Berggren branches into affine maps, and admit a strict inverse
descent.  Projecting to `r` has exact fibre size `phi(2r-1)/2`; the first
collision is already the two outer children of `(3,4,5)`.

The more general lesson is sharper than “replace a selected value by its
index.”  An ordinal replacement is structurally useful when at least one of
the following holds:

1. the consumer is genuinely order/shell-only; or
2. the native operation closes on the selected set and pulls back to a cheap
   law on ordinals.

If neither holds, the ordinal is only a scheduler.  It usually destroys gaps,
divisibility, weights, reciprocal mass, coefficient cancellation, or semantic
ownership.

## Inheritance pass

| role | inherited object |
|---|---|
| closest proved mechanism | THM-3333, `gaussian-square-farey-pythagorean-triangular-light-cone`, gives the integer triangular refinement and ordered Pythagorean carrier |
| tree mechanism | THM-3357, `berggren-three-branch-walsh-level-collapse-and-parent-circuit`, fixes the `L/M/R` parameter convention and warns that aggregate branch data loses words |
| canonical hostile | THM-3382, `fibonacci-ray-dual-index-harmonic-bifurcation-and-ternary-heap-addresses`, proves that time, depth, and heap injections of the same objects have different harmonic class |
| corrected near miss | MISTAKE-418: parity/content normalization must precede any primitive Berggren slope transport |
| least-used sidecar | the complementary odd square `C-B`, rather than another statistic of `B+C` |

MISTAKE-222 blocks turning shared triangular syntax into an LRC/JC bridge.
MISTAKE-209 blocks treating indexed multiplicity, support, and reciprocal
weight as one object.  Those failures became design constraints rather than
reasons to avoid ordinal coordinates.

## Active concept board

| lane | object / representation | predicate | operation | lost coordinate | cheapest test |
|---|---|---|---|---|---|
| Anchor | primitive triple as odd-root pair `(2r-1,2s-1)` | primitivity and ancestry | Berggren `L/M/R` | projection to `r` loses `s` | rank-three sibling collision and inverse-cone replay |
| Niche | THM-2349 profile `(1,b,c)` as triangular shortlex rank | lawful row address | only proved row transports, not ordinal successor | `b`, owner, word, root | test `4:(4,5) -> 5:(1,6)` |
| Niche | THM-3713 defect profile on balanced offsets | positivity/detection | centered triangular difference | origin, orientation, augmentation, colours | `5 delta_1-3 delta_2` |
| Wildcard | odd tournament invariant `H=2r-1` | ordered-join composition | multiplication | SCC factors and path-colour fibre | equal-`H` masks with different responses |
| Wildcard | exponent pair `(a,b)` as diagonal triangular rank | monomial landing | product / Poisson bracket | signed convolution owners | `{x+y,(x+y)^2}=0` |
| Wildcard | linear-width layer `(k,j)` | lossless enumeration | child append | canonical within-layer order | compare ternary endpoint cells with a genuine ternary tree |

Every subsequent pull was compared with this board.  The common mechanism is
not “triangular numbers occur”; it is “prefix size plus a within-fibre ordinal
transports a native operation, with a named loss ledger.”

## Anchor synthesis: two odd-square ordinals

With `B` the ordered even leg, define

```text
q=sqrt(B+C)=2r-1,                 d=sqrt(C-B)=2s-1.      (3)
```

Then

```text
(A,B,C)=(qd,(q^2-d^2)/2,(q^2+d^2)/2),                  (4)
Omega={(r,s):r>s>=1, gcd(2r-1,2s-1)=1}.                (5)
```

The map `Omega -> primitive triples` is bijective.  The primitive converse
does not need Euclid's parametrization: `C+B` and `C-B` are coprime positive
odd factors whose product is `A^2`, so each is a square.

The user's centered identity is the decoder:

```text
T(z+h)-T(z-h)=h(2z+1),
[T(r+1)-T(r-3)]/2=2r-1.                                (6)
```

In THM-3357 coordinates `0<m<n`,

```text
m=r-s,                         n=r+s-1,                 (7)
L(r,s)=(r+2s-1,s),
M(r,s)=(2r+s-1,r),
R(r,s)=(2r-s,r).                                        (8)
```

The inverse cones are

```text
r>=3s                 -> L-parent (r-2s+1,s),
2s<=r<=3s-2           -> M-parent (s,r-2s+1),
s<r<=2s-1             -> R-parent (s,2s-r).              (9)
```

The missing wall `r=3s-1` is the component-root wall.  If
`g=gcd(2r-1,2s-1)`, the decoded triple has content `g^2`, the branch maps
preserve `g`, and descent terminates at

```text
root_g=((3g+1)/2,(g+1)/2)  ->  g^2(3,4,5).              (10)
```

Thus the full triangular cone is a forest of odd-square-scaled Berggren
trees; the primitive tree is the `g=1` component.  The nonprimitive slots are
structured objects, not rejected noise.

At fixed outer rank,

```text
|{s:(r,s) in Omega}|=phi(2r-1)/2.                       (11)
```

The proof pairs reduced residues `d` and `q-d`; exactly one is odd.  Prime
`2r-1` gives `r-1` overlaps, so coarse fibres are unbounded.

### Four useful address gauges

| gauge | formula | use | cost |
|---|---|---|---|
| coarse shell | `r` | queries depending only on `B+C` | unbounded collision fibre |
| lossless pair | `(r,s)` | affine branches and descent | not one scalar |
| ambient triangular | `T(1-r)+s=T(r-2)+s` | bijects the whole content forest with positive naturals; holes expose nonprimitivity when restricted | decode before acting |
| selected-shell | `sum_(k<r)phi(2k-1)/2 + j_r(s)` | bijects primitive nodes with positive naturals | totative selection hides affine branch law |
| ternary heap | recursive word address | branch append / ancestry | destroys odd-square shell order and changes analytic weights |

The correct gauge is chosen by the consumer.  There is no encoding-invariant
notion of density or harmonic mass.

## Niche pull 1: a triangular rank for all 165 LRC profiles

THM-2349, `first-depth-one-delayed-shallow-restart`, has the exact valuation
profile universe

```text
(1,b,c),             5<=c<=19,       1<=b<c.             (12)
```

The formula

```text
rho(b,c)=T(c-2)-6+b                                  (13)
```

is a bijection onto `{1,...,165}`.  Its inverse chooses the unique `c` with

```text
T(c-2)-5<=rho<=T(c-1)-6.                               (14)
```

This is a lawful **address** of an already proved finite universe.  It is not
a row operation.  The first shell transition is the decisive hostile:

```text
rho(4,5)=4,                    rho(1,6)=5.               (15)
```

No theorem transports owner/current/word data from the first row to the
second merely because their ordinals are consecutive.

The centered difference makes the limitation visible:

```text
rho(b,c+2)-rho(b,c-2)=4c-6.                             (16)
```

It is independent of `b`, so it cannot distinguish repeated-first from
strict behaviour, much less detect the semantic owner.

The abstract repeated-first lane `b=1`, indexed intrinsically by `n=c-4`,
has ambient addresses

```text
rho_n=T(n+2)-5=(n^2+5n-4)/2.                            (17)
```

Hence `sum 1/rho_n` converges while `sum 1/n` diverges.  The live ledger ends
at `c=19`, so (17) is an encoding warning, not an LRC asymptotic theorem.

## Niche pull 2: the triangular flux reduces to an old detector

For a THM-3713 defect profile `d_u` on balanced lifts `u in [-6,6]`, define

```text
Q_T(d)=sum_u [T(u+2)-T(u-2)]d_u.                        (18)
```

Since the bracket is `4u+2`,

```text
Q_T=4L+2A,
L=sum_u u d_u,                         A=sum_u d_u.      (19)
```

On THM-3713's typed five-colour control, positive coefficients occur at
`1,2` and negative coefficients at `-1,-2,-3`; every term of (18) is then
positive.  But this supplies no new complete detector:

```text
d=5 delta_1-3 delta_2 !=0,              Q_T(d)=0.        (20)
```

The profile is rational and anchored at `d_0=0`, so THM-3713's cyclotomic
mechanism still makes every nonzero deep colour live.  The centered
triangular functional sees only the already present first moment plus
augmentation.  It also depends on the physical origin and orientation, so it
does not descend through arbitrary rotation.

This is a genuine negative result: the triangular potential is a cone
separator with sidecars, not LRC(14) progress.

## Wildcard pull 1: integer quadratics are linear-width prefix potentials

Every integer-valued quadratic on `Z` is uniquely

```text
Q(z)=A binom(z,2)+Bz+C,                 A,B,C in Z.      (21)
```

It obeys

```text
Q(z+h)-Q(z-h)=h[A(2z-1)+2B],
Q(z+h)+Q(z-h)-2Q(z)=Ah^2.                              (22)
```

The forward difference

```text
w(k)=Q(k+1)-Q(k)=Ak+B                                  (23)
```

is affine.  Therefore `Q(k)` is exactly the prefix size of a layered family
whose layer `k` has width `Ak+B`.  If `j` is the within-layer ordinal,

```text
iota(k,j)=Q(k)+j,
iota(k+1,j+delta)-iota(k,j)=Ak+B+delta.                 (24)
```

This explains when triangular and square addresses are structural: the
underlying layer width is linear.  It also gives a stopping reason.  A
genuine `q`-ary tree has exponential level size `q^k`, so its lossless
breadth-first address is radix/heap, not quadratic.  A quadratic can be a
height potential there, but not a collision-free enumeration of nodes.

Centered linearization alone does not characterize quadratics.  The hostile

```text
Q_tilde(z)=Q(z)+(-1)^z                                  (25)
```

has exactly the same centered first differences for every integer `h`.
A fixed step has an even larger periodic kernel.  Constant forward second
difference, or an explicit parity sidecar, repairs the converse.

### Graph-height corollary

For an integer height `h(v)`, write `delta_w=h(w)-h(v)` and use the graph
Laplacian convention `Lf(v)=sum_(w~v)(f(w)-f(v))`.  Then

```text
L(Q o h)(v)
 =(A h(v)+B)sum_w delta_w+A sum_w binom(delta_w,2).      (26)
```

On a `+/-1` graded graph with `c` children and `p` parents,

```text
L(Q o h)=(Ah+B)(c-p)+Ap.                                (27)
```

Height alone destroys the local incidence counts.  Equal-height vertices
with different `(c,p)` are the canonical hostile; `(c,p)` is the minimal
sidecar.

## Wildcard pull 2: ternary endpoint shells end at squares

The endpoint quotient of length-`k` words in `{-1,0,1}` has states
`(k,s)` with `-k<=s<=k`, hence linear width `2k+1`.  Formula (24) becomes

```text
iota(k,s)=k^2+k+s in [k^2,(k+1)^2-1],                  (28)
iota(k+1,s+epsilon)-iota(k,s)=2k+2+epsilon.             (29)
```

Pushing word multiplicities onto the endpoint fibre gives

```text
W(k+1,s)=W(k,s-1)+W(k,s)+W(k,s+1),                      (30)
```

the central-trinomial walk underlying THM-2438's Newton halves.  The square
boundary is the cumulative odd width, not a numerical coincidence.

The endpoint quotient loses word order and prefix barriers:

```text
(+1,-1) and (-1,+1) both land at (2,0),                 (31)
```

but only the first keeps nonnegative partial sums.  Multiplicity `W`, or the
full word, is the sidecar required by word-sensitive consumers.

## Wildcard pull 3: odd tournament values inherit a join law

Let `H(T)` be the odd number of Hamiltonian paths of a tournament and define

```text
r_H(T)=(H(T)+1)/2.                                      (32)
```

For the ordered tournament join, Hamiltonian paths concatenate, so
`H(T -> U)=H(T)H(U)`.  Pulling multiplication through (32) gives

```text
r star s=2rs-r-s+1,
2(r star s)-1=(2r-1)(2s-1).                             (33)
```

The law is commutative and associative with identity `1`.  Cyclic triples
are additive under ordered join because every cross-component triple is
transitive.  Thus the pair `(r_H,c_3)` composes as `(star,+)`.

This is a positive example: replacing the `r`th odd value by `r` preserves a
native operation through a cheap polynomial law.  It still loses the ordered
strong-component factorization and every invariant inside an `H`-fibre.
THM-1370, `h-spectrum-omits-7-21-all-n`, makes the exact holes `H=7,21`
become ordinal holes `4,11`; its coverage through `H=609` becomes coverage
through rank `305` except `4,11`.  Global completeness remains conjectural.

## Wildcard pull 4: diagonal exponent ranks transport landing, not cancellation

The diagonal shell address

```text
rho(a,b)=T(a+b)+b                                       (34)
```

bijects nonnegative exponent pairs with `N_0`.  If `u=(a,b)`, `v=(c,d)`,
`h=|u|_1`, and `h'=|v|_1`, then monomial multiplication obeys

```text
rho(u+v)=rho(u)+rho(v)+hh'.                             (35)
```

For a nonzero monomial Poisson bracket, whose determinant gate is
`ad-bc !=0`, the landing exponent is `u+v-(1,1)` and

```text
rho(u+v-(1,1))
 =rho(u)+rho(v)+hh'-2(h+h').                            (36)
```

This gives a useful exact address of landing shells.  At polynomial support
level it forgets coefficients and contributing-pair owners.  The minimal
hostile is

```text
{x+y,(x+y)^2}=0,                                        (37)
```

where four nonzero raw bracket channels land in two ranks and cancel pairwise.
A signed labelled convolution fibre is the required sidecar.

## How to decide whether “the selected value is just n” is safe

| question | keep |
|---|---|
| only which selected shell? | coarse ordinal |
| reconstruct the object? | ordinal plus complete fibre selector |
| apply a native operation? | prove and use the pulled-back operation law |
| enumerate every object densely? | selected-shell/shortlex address with round-trip decoder |
| retain ancestry? | heap/word address or parent sidecar |
| study gaps, divisibility, density, reciprocal mass, or analytic weights? | native selected value, gap/Jacobian, and collision multiplicity |
| use a quotient in LRC/JC? | owner, phase/sign, root/word, coefficient fibre, and target predicate as applicable |

Examples separate the boundary cleanly:

- squares preserve multiplication after rank extraction;
- odd values preserve multiplication through `star`;
- Pythagorean odd-square shells need the complementary ordinal `s` for tree
  operations;
- primes are not closed under multiplication, so prime index is only a
  scheduler there;
- ranking a generic support destroys precisely the analytic information in
  THM-2000/THM-2438/THM-3382.

## Incoming-signal integration: depth and target value are also quotients

Two proved, independently audited JC near-misses arrived on `origin/main`
during this session.  They do not identify JC objects with Pythagorean ones;
they test the same quotient discipline on different native operations.

| incoming proof | proposed natural coordinate | preserved predicate / operation | destroyed information | required sidecar / cheapest hostile |
|---|---|---|---|---|
| THM-3757, `pell-chebyshev-three-charge-hyperelliptic-obstruction-tower` | Chebyshev construction depth `n` | order in the tower and the recurrence-driven Pell/transport construction | the actual profile divisor, generic-fibre degree/genus, and the class of `dz/Y` | retain `(psi_n,Delta_n,div(dz/Y))`; compare the two nonzero residues at `n=1` with the nonzero holomorphic class at `n>=2` |
| THM-3758, `quadratic-radial-carrier-rational-exact-split-fibre-nonentry` | target value `L=Q`, or an ordinal assigned to target fibres | the fixed-`L` Hamiltonian derivative and complete rational primitive torsor | the two components of `Q=0` and the sign of the primitive's principal part on each | retain the component-labelled principal-part vector; the exact hostile is `(-c a0/(2a1),+c a0/(2a1))`, while a target-only correction adds the same entry to both |

The first row says that `n` is a lawful construction scheduler, but not an
obstruction invariant: the obstruction changes type between depth one and
the positive-genus tail.  The second is a sharper quotient failure.  Generic
fibre integration survives scalarization by `L`, yet global polynomial
regularity lives in a signed fibre over one value.  Merely ranking those two
components `1,2` still loses the sign; the component ordinal must carry a
divisor coefficient.

## Procedural task generator

For any discrete family `X`, selected value set `S={e(n)}`, and candidate
rank `rho`, generate tasks in this order:

1. **Fibre task:** compute `rho^(-1)(n)` exactly; find the first collision and
   decide whether multiplicities are bounded.
2. **Complement task:** search for another signed/dual value of the same type
   that reconstructs the object.
3. **Operation task:** pull every native generator through `rho`; classify it
   as affine, polynomial, finite-state, nonlocal decoding, or not closed.
4. **Descent task:** partition image chambers, derive inverse maps, and seek a
   strictly decreasing potential.
5. **Rejected-slot task:** restore filtered states and classify their content,
   scale, or obstruction components.
6. **Address task:** compare ambient shell, selected-shell, heap, and canonical
   orbit addresses; require round-trip decoding.
7. **Analytic task:** recompute density and harmonic/Dirichlet profiles for
   each injection rather than transferring them.
8. **Consumer task:** name the exact predicate preserved, destroyed
   information, minimal sidecar, and one hostile before using the rank.

Each positive or negative result feeds its missing coordinate back into the
board and generates the next sidecar task.

## Ranked follow-up portfolio

### Anchor

1. **COMPLETED in THM-3756: formalize the odd-square ordinal forest.**  The
   theorem now contains (3)--(10), the `phi/2` fibre, ambient/selected
   addresses, exact ordinary/optimized companions, and an independent
   hostile audit.  Its completion promotes tasks 2--4 to the live anchor.
2. **Primitive-address density and Dirichlet profile.**  Prove the density of
   coprime odd pairs in triangular shells and compare native shell versus heap
   weights.  Do not infer it from a finite census.
3. **Depth distribution inside one totient fibre.**  Minimum, maximum, and
   mean Berggren descent length at fixed `r`; test prime and highly composite
   `2r-1` as hostiles.
4. **Heap/shell transducer.**  Determine whether bounded residue sidecars make
   the map from `(r,s)` shell order to ternary word regular; the U/R boundary
   rays and unbounded L-descent are first controls.

### Niche: LRC

5. **Lawful operation graph on the 165 profiles.**  Add an edge only when a
   proved scalar-row operation transports the relevant current.  Test whether
   `rho`, depth, or another potential is monotone.  Consecutive ordinals and
   the `b`-blind difference (16) are mandatory negatives.
6. **Adaptive selector ordinal atlas.**  For each THM-3718 witness, record the
   complete atom bits, terminal word, full mask/necklace, transporter,
   defect coordinate, target, owner/root provenance, and shortlex task code.
   Search for a fibre with a common semantic interpretation.
7. **Cover-derived triangular-flux stop test.**  Evaluate (18) on the actual
   selectors.  Stop on a (20)-type cancellation; a positive signal matters
   only after proving a cover-specific sign cone and origin/orientation.

### Wildcard

8. **Poisson signed rank-collision atlas.**  Apply (34)--(36) to the live
   W004--W006 monomial atlases; compare raw landing-rank multisets with signed
   aggregated support and record first cancellation-fibre dimensions.
9. **Tournament `star` atom bank.**  Factor attained odd ranks through ordered
   strong components, compute bounded `star` closure, and seek a conductor
   argument beyond holes `4,11` without assuming the open spectrum
   completeness claim.
10. **Ternary endpoint-square quotient.**  Compute the endpoint-positive and
    prefix-positive languages separately; (31) is the minimal lost-barrier
    hostile.  Do not invoke THM-3499 before proving regularity.
11. **Rule30 commuting-parent unit scout.**  Distinguish suffix deletion from
    highest-set-bit deletion, then test their mixed unit rather than another
    valuation.  An exact finite scout found the valuation shell-flat, so the
    valuation-only lane has a stopping reason; a finite-state unit law remains
    open.
12. **Rewrite-normalized Cohn shortlex queue.**  Apply proved reductions before
    ranking words; every closure theorem deletes a task, and the least
    surviving decoded word becomes the next concrete JC obligation.  Retain
    polynomial parameters and seam holonomy.
13. **Signed special-fibre component atlas.**  For every smooth JC candidate
    with a rational primitive, rank the irreducible components of each
    exceptional fibre only after recording multiplicity, intersections, and
    the full principal-part vector.  Test whether available correction
    channels span that vector modulo the diagonal target-only subspace.
    THM-3758 is the two-component hostile and makes this a value-5/cost-2
    exact linear-algebra task.
14. **Pell-depth obstruction-state queue.**  Rank Pell/Chebyshev construction
    words by depth, but immediately compute the supported resultant,
    squarefree generic degree, genus, and divisor class of the forced time
    form.  Delete a task as soon as residues or a holomorphic class survive.
    THM-3757 supplies positive smoothness controls and both obstruction types;
    the live question is whether a non-Chebyshev word pays both debts.

## Meta-pattern use and candidate promotion

Cards used:

- **Type every analogy and every implication.**  Rank, matrix rank, tree depth,
  and address were separated throughout.
- **Find the hidden second coordinate.**  `C-B` repaired the coarse `B+C`
  quotient.
- **Controlled forgetting requires a sidecar.**  The fibre selector, word,
  origin, coefficient convolution, or incidence counts were named explicitly.
- **Fill operation columns.**  Branch, join, append, multiplication, bracket,
  and Laplacian laws produced more structure than scalar matches.
- **Classify response-state growth before naming a closed form.**  Quadratic
  prefix, radix heap, and selected-shell addresses were kept distinct.

Candidate new card: **rank a selected set only after pulling back its native
operations**.  Evidence now comes from distinct threads: affine Berggren laws,
tournament join multiplication, layered append, Poisson landing, and the LRC
profile/flux hostiles.  Its counterindication is equally clear: if the
operation exits the selected set or the consumer uses gaps/divisibility/
weights, retain the native value and its analytic sidecar.

## Honest remaining frontier

- THM-3756 has passed independent promotion; the remaining anchor questions
  concern shell/heap analytic profiles, fixed-shell depth, and transducers,
  not the correctness of the affine forest itself.
- No LRC row was excluded; LRC(14) remains open.
- The triangular flux is a typed cone separator, not a complete detector.
- The tournament spectrum beyond the proved finite coverage and permanent
  holes remains open.
- Poisson/exponent ranking has not yet been run on the live W004--W006 atlas.
- None of the ordinal maps transports JC coefficients or cancellations
  without the signed convolution sidecar.
