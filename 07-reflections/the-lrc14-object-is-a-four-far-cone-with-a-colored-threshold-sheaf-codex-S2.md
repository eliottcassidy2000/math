# The LRC14 Object Is A Four-Far Cone With A Colored Threshold Sheaf

*codex-2026-07-14-S2. A creative information-preservation session following
HYP-6780, the completed Lean THM-755 chain, HYP-6785, and the exact endpoint
sidecar audit.*

## 1. First challenge the word "four-dimensional"

There is a literal four-dimensional object at the current frontier, but it is
not the whole LRC14 moduli space.

THM-738 closes the rows with at most three speeds above 14. The first
unclosed far-count stratum therefore has a fixed nine-speed core
`P subset {1,...,14}` and four ordered far speeds. In gap coordinates those
four speeds are controlled by

```text
(a,d1,d2,d3) in Z^4,
w=(a,a+d1,a+d1+d2,a+d1+d2+d3).
```

This is the honest four-dimensional chart. There are `binom(14,9)=2002`
possible small cores, which explains the natural 2002-body rung behind
THM-741.

The global primitive 13-speed space modulo dilation has twelve ratio degrees
of freedom, and the `f>=5` strata have larger far cones. Any claim that the
entire problem is four-dimensional already assumes a descent theorem that the
repository does not possess. The first reframe is therefore negative and
useful:

> Do not reduce the dimension in prose. Exhibit the functor that performs the
> reduction and list the fiber information it must retain.

## 2. The outer object

For a fixed core `P`, the admissible far-speed set `X_P` is an arithmetic
lattice cone. Covering depends only on divisibility by `2,...,14`, and hence
on residues modulo `360360`. Primitivity is also periodic on this base. Thus
`X_P` is a finite union of residue classes intersected with the positive
ordered cone.

That observation does not make LRC periodic. Inside one residue class,
changing height changes:

- the rational endpoint phases;
- the cyclic event order and simultaneous tie blocks;
- the safe component lengths;
- the pair-sum moduli used by peak witnesses;
- the THM-755 ratio and Bernoulli discrepancy;
- the embedded integer relations and their coefficient heights.

The point of the cone is not finite enumeration by itself. It supplies the
base over which the predicate-bearing information varies.

For `g in X_P`, the speed vector gives a cocharacter

```text
phi_g : T -> T^13,
t |-> (s*t mod 1)_{s in S(P,g)}.
```

Let `K=[1/14,13/14]^13`. The universal safe incidence

```text
Z_P={(g,t): phi_g(t) in K}
```

is the cleanest exact formulation of the four-far problem. LRC on this chart
is surjectivity of `Z_P -> X_P`.

This is a toric statement in a precise sense: speeds are character exponents,
the runner motion is a one-parameter subtorus, and the target is a fixed
coordinate box. It is not yet a toric proof.

## 3. The inner object

Many hard infinite families are affine rays

```text
V(c)=cP+R.
```

Put `u={ct}`. Then exactly

```text
(c*p_i+r_i)t = p_i*u+r_i*t mod 1.
```

One two-torus field

```text
Phi(u,t)=min_i ||p_i*u+r_i*t||
```

contains every member of the family as the slope-`c` slice `u=ct`. Adding
the integer slope `c` and clearance level `lambda` produces the mixed
four-coordinate suspension `(u,t,c,lambda)`.

This recovers several old views as literal coordinate projections:

- `R=0`: a cylinder, hence dilation invariance;
- bounded `R`: transverse shear, not negligible error;
- threshold `lambda`: the persistence direction;
- fixed `c`: the ordinary exact circle arrangement;
- varying `c`: arithmetic tomography by closed geodesic slopes.

THM-742 already proved this strand-versus-area idea for an additive cluster
chart. HYP-6815 is the arbitrary owned-offset version, with the threshold
coordinate retained.

## 4. The dual object

For every coefficient vector `z in Z^13`, project the relation to

```text
(m,n)=(z dot P,z dot R) in Z^2.
```

On the two-torus its character is `exp(2*pi*i*(m*u+n*t))`. On slope `c` it
becomes `exp(2*pi*i*(c*m+n)*t)`. Therefore the relation survives exactly on
the dual line

```text
c*m+n=0.
```

This is the exact bridge among four repository languages:

```text
normalized shape        -> m
owned residue/detuning  -> n
scale slope              -> line c*m+n=0
Fourier/relation mass    -> lattice points on that line.
```

The shape alone invents relations. In the executable audit, the AP shape has
36 support-three relations `e_i+e_j-e_k`; one owned `+1` shear destroys 11 if
attached to owner 1, 9 if attached to owner 7, and 6 if attached to owner 13.
The missing information is not "a residue" in the abstract. It is the
pairing between the residue coordinate, its owner, and the relation
coefficient.

The relation lattice must also remain embedded and marked. A code weight
enumerator, lattice rank, additive energy, or successive minimum forgets how
the relation meets the thirteen coordinate facets of the safe cube. That
incidence is where the observer lives.

## 5. There is no universal sufficient quotient

The phrase "what information must be preserved?" is incomplete until the
next question is named.

### To decide LRC truth at one fixed threshold

Start at `t=0`, where every runner is dangerous. Record the cyclic endpoint
events. Each simultaneous block retains entering owners and exiting owners.
Integrating these flips reconstructs the cubical path

```text
beta(t)_i=1_{||s_i*t||<1/14} in {0,1}^13.
```

The row is lonely exactly when an open cell or exact tie state reaches the
zero vertex. For this Boolean question, metric gap lengths can be forgotten.

Owner colors and tie blocks cannot. If one runner exits when another enters,
their scalar currents cancel. The aggregate current would erase a boundary
state that may be the only equality witness.

### To preserve exact clearance

The fixed-threshold event word is insufficient. Retain the whole threshold
filtration or an exact peak witness and value. A peak certificate determines
`M`, but it does not determine safe measure or component topology.

### To preserve measure and autocorrelation

Retain rational endpoint phases and gap lengths. The signed owner endpoint
divisor reconstructs the Bernoulli `B2` formula. The unweighted endpoint
tournament sees interlacing but not metric discrepancy.

### To preserve covering

Retain the divisor mask. Neither the runner tournament nor endpoint
tournament determines it.

### To preserve the capped-envelope decision

Retain the projective ratio `v*|G|/r`, or its exact comparison bit. A single
lift in the endpoint audit changed the cap status with zero endpoint edge
flips.

### To preserve deletion, peel, or scale transport

Retain owner labels and the action on the scale/residue fiber. A quotient can
be truth-safe for the current row and illegal for the next observer.

### To preserve a proof

Retain the next operation, certificate availability, discharge mode, and
named residual debt. This is proof state, not extra mathematical truth, but
without it local certificates cannot be composed honestly.

The hierarchy is:

```text
truth   = initial state + signed owner event blocks
metric  = truth + rational phases or threshold filtration
functor = metric + arithmetic action + marked relation embedding
proof   = functor + next observer + certificate/discharge
```

## 6. The exact canaries

The companion 552-row audit gives a compact impossibility theorem for several
tempting compressions.

Raw runner tournaments have mixed fibers for covering, `M`, cap status, and
discrepancy. Raw endpoint tournaments also mix all four. Adding the divisor
mask fixes covering. Adding the cap sign fixes the THM-755 bit. Signed
endpoint phases fix `B2`. Exact `M` still mixes until a peak witness is
attached.

The cleanest pair is

```text
{1,...,13}       M=1/14
{1,...,12,26}    M=2/27.
```

They share the endpoint tournament, divisor mask, and cap sign. A quotient
which keeps order, covering, and the terminal cap route still loses the
clearance.

The affine audit supplies a different canary. For the same shape, same owned
offset, and the same `c mod 14`, scales `c=2` and `c=16` both have `M=2/15`,
but their safe measures are `5/84` and `115/1904`, with 4 and 30 components.
Even exact `M` plus residue data is not a metric-topology carrier.

## 7. What unrelated threads contributed

### Tournament switching and tilings

The same cube with a different group folding has different invariants. The
transferable datum is the action and a gauge section, not the underlying
bijection. For LRC the analogue is scale/residue action plus the marked
observer.

### Coding theory and matroids

Identical weight enumerators can hide decomposable and indecomposable support
graphs. Relation support counts are not relation realization. LRC needs
support incidence, coefficient height, and the embedding against safe-cube
facets.

### Ising and Lee-Yang

Total odd-cycle mass and compatible cycle packings separate at higher
interaction order. A scalar moment or `p0` loses the organization of
obstructions. The analogue here is the full threshold or miss-count
filtration. Root surfaces can classify phase changes, but do not replace an
endpoint certificate.

### Resolvent folds

A symmetric inner block becomes exact only after center coupling,
antisymmetric leakage, and boundary sectors are retained. Projective shape
is the inner block; owned shear and endpoint boundary data are the leakage
sidecars.

### CRT and p-adic trees

Residue skeleton and valuation height are independent channels. The affine
audit gives an archimedean version of the same lesson: equal residue and owner
can retain `M` while safe length and component count move.

### Cut/cycle topology and anti-local testing

Local views can satisfy every small consistency check and still fail global
realizability. The endpoint loop and endogenous blocker complex are two
global presentations that must be glued, not two scalar scores to average.

### Observer-cut ledgers and proof circuits

A quotient is legal only relative to its next observer. Every changed
predicate needs reconstruction, exact/coboundary cancellation, dual
annihilation, family descent, a boundary stop, or named residual debt.

## 8. A different style of mathematical thought

The conventional question asks for an invariant of a speed set. The more
productive question asks for a contract between representations:

```text
What operation comes next?
Which fibers does this quotient merge?
Which theorem-facing predicates vary inside those fibers?
What action survives on the forgotten coordinate?
Which sidecar repairs the first failure?
```

This turns failed invariants into useful output. A failure names the next
coordinate rather than merely disqualifying an analogy.

It also suggests replacing one linear proof ledger by a bicomplex:

```text
d_geom  = cross an endpoint/contact wall
d_arith = change residue, valuation, scale, or detuning depth.
```

If the two moves commute, a quotient may descend. If they do not, their
commutator is curvature: an owner, holonomy, coefficient-height, or proof-route
sidecar that was erased. This is an exact research program once both maps are
implemented on endpoint packets; the word "curvature" alone proves nothing.

## 9. Recursive view of infinity

HYP-6780 showed why raw maximum speed is the wrong induction variable. The
four-cone suggests a replacement.

At infinity, keep a flag:

```text
leading projective shape P
first nonzero owned detuning R
residue/valuation address
cluster tree or relation circuit
terminal certificate or smaller descendant.
```

Pure dilation is the zero-detuning cylinder. Additive clusters, multiplicative
packs, and hierarchical clocks are boundary faces with known dispatches. A
genuinely incoherent interior point should either acquire a uniform metric
floor or emit a bounded-height relation/blocker circuit. The desired theorem
is not "the box ends at 500". It is:

> Every infinite path in the normalized cone either enters a known coherent
> face, gains a uniform lonely margin, or strictly descends in a well-founded
> packet order.

That statement would make the finite computations terminal leaves of a proof
rather than samples from an infinite box.

## 10. Concrete next lemmas

1. Formalize the semilinear four-cone and colored endpoint-loop criterion.
2. Prove cubical zero-state reachability equivalent to an empty protected edge
   in the endogenous pair-sum blocker complex.
3. Test a finite-jet-at-infinity theorem on `cP+R`: leading shape plus first
   detuning, finite residue address, and certificate should eventually decide
   every ray, or expose the first counterexample.
4. Prove a four-circuit localization lemma: outside known coherent faces, a
   blocker-complete point forces a bounded-height marked relation involving all
   four far coordinates.
5. Build an observability matrix over the actual `f=4` cone, not a curated row
   bank, and find the smallest sidecar portfolio separating every
   truth/metric/peel-changing fiber pair.
6. Use the existing 2002-core runner as a base case only after its output and
   coverage protocol are present and resumable.

## 11. What this session did not prove

The affine suspension is an exact reparameterization, not a compactness
theorem. The semilinear residue base does not make LRC residue-periodic. The
threshold sheaf has not been globally trivialized. The four-circuit
localization and cone descent remain open. Toric, sheaf, persistence, and
curvature language is retained only where an exact object or falsifiable lemma
has been stated.

The durable conclusion is narrower:

> The predicate-bearing object is a four-far arithmetic cone emitting a
> marked toric loop. Its primal data is an owner-colored threshold event
> path; its dual data is an embedded relation/blocker complex. Any quotient
> must preserve the action on its forgotten fiber, and the amount of retained
> information must be indexed by the next mathematical question.
