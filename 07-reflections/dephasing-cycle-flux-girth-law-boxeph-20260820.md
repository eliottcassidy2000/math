# The cycle-flux girth law for the raw dephasing Schur walk

**Status: PROVED as a finite-dimensional algebraic statement; C3--C6 and
hostile controls are FINITE-EXACT.**  This note concerns raw Schur-complement
words.  It does not identify all convention-dependent terms in a time-local
normal form.

## Statement

Let `G` be a finite simple graph, let `H` be a Hermitian hopping matrix
supported on its edges, and allow arbitrary real onsite energies.  Assume
uniform off-diagonal dephasing `-Gamma Qcal`.  Let `Pcal` project matrices
onto their diagonal, put `Qcal=1-Pcal`, and use

```text
K=-i[H,-],  B=Pcal K Qcal,  C=Qcal K Pcal,  D=Qcal K Qcal.
```

For `m>=2`, define the raw length-`m` Schur word

```text
A_m = B D^(m-2) C.                                      (1)
```

If `G` has girth `g`, every `A_m` with `m<g` is independent of all hopping
phases.  The first raw word that can contain a nontrivial Wilson monomial is
`A_g`.  Formally independent phase variables show that every girth cycle
does contribute a nonzero monomial to `A_g`; special numerical fluxes can
make its evaluated coefficient vanish.

For `Gamma>||D||` (or formally as a series in `Gamma^(-1)`), the
zero-frequency Schur expansion is

```text
B (Gamma-D)^(-1) C
  =sum_(m>=2) Gamma^(-(m-1)) A_m.                        (2)
```

Thus a girth-`g` phase first becomes available at physical-generator order
`Gamma^(-(g-1))`, or at order `Gamma^(-(g-2))` after the leading-rate time
normalization.

## Why girth is the exact threshold

Every matrix element of `(1)` is a sum over paired walks.  The two paths
begin and end on diagonal matrix units, so their union is a closed contour
with `m` Hamiltonian factors and at most `m` edge traversals; onsite factors
are stationary steps.  Its phase is the oriented edge product of that
contour.  A contour with fewer than `g` traversals has no reduced cycle:
after deleting backtracks its support is a tree, and every oriented edge
exponent cancels.  It is therefore phase blind.

The inserted `Qcal` projections do not change the first phase-sensitive
part.  A term removed by one of them returns to the diagonal at an
intermediate time and hence factors its closed contour into two proper
shorter closed walks.  Each factor has fewer than `g` Hamiltonian factors,
hence fewer than `g` edge traversals, and is phase blind.  Consequently the
phase-dependent part of `A_g` equals that of `Pcal K^g Pcal`.

At Hamiltonian-factor length exactly `g`, a phase-sensitive contour must use
all `g` factors as edge traversals; its nontrivial reduction is an orientation
of a simple girth cycle.  Extra edges and onsite factors cannot occur in its
Wilson monomial because all `g` hopping factors have already been used.
This also shows, with independent formal edge phases, that different girth
cycles cannot erase the statement by an accidental symbolic cancellation.

## Exact formula on the unit cycle

Index `C_g` by `Z/gZ`.  Give every consecutive edge unit hopping, put the
whole Laurent Wilson variable on

```text
H_(g-1,0)=z,  H_(0,g-1)=z^(-1),
```

and set the onsite energies to zero.  Let `S` be the cyclic shift
`(S p)_j=p_(j+1)` and define

```text
R_g=(I-S)^g.                                             (3)
```

Then

```text
A_g(z)-A_g(1)
 = i(-1)^((g-1)/2)(z-z^(-1)) R_g,        g odd,          (4a)
 =   (-1)^(g/2)(z+z^(-1)-2) R_g,         g even.         (4b)
```

For `z=exp(i Phi)`, these factors are respectively

```text
2(-1)^((g+1)/2) sin(Phi),
2(-1)^(g/2)(cos(Phi)-1).                                (5)
```

Thus odd cycles first return as an oriented, time-reversal-odd circulation;
even cycles first return as a reciprocal, time-reversal-even correction.
Explicitly, if `d=j-i mod g`,

```text
(R_g)_(ij)=(-1)^d binom(g,d),                    d!=0,
(R_g)_(ii)=1+(-1)^g.                                    (6)
```

It follows at once that

```text
R_g^T=(-1)^g R_g,  R_g 1=0,  rank(R_g)=g-1.             (7)
```

The same formula for a weighted chordless cycle is multiplied by the
product of its hopping magnitudes.  Formula `(4)` follows by expanding

```text
K^g=(-i)^g sum_(r=0)^g (-1)^r binom(g,r)
                    L_H^(g-r) R_H^r
```

and retaining the two complementary simple arcs from the input vertex to
the output vertex.  The binomial coefficients in `(6)` record the possible
interlacings of the left and right arc steps.

This local formula also gives the whole graph-wide first phase term.  Let
`H_abs` have the same onsite entries and hopping magnitudes as `H`, with all
edge phases set to zero.  For every unoriented girth cycle `C`, choose a
cyclic order, let `J_C` embed its `g` population coordinates into the full
vertex space, put `R_C` equal to the product of its hopping magnitudes, and
let `z_C` be its Wilson phase.  If `f_g(z)` denotes the scalar factor in
`(4)`, then

```text
A_g(H)-A_g(H_abs)
 = sum_(girth cycles C) R_C J_C [f_g(z_C)(I-S)^g] J_C^T. (7a)
```

Reversing the chosen orientation changes both factors by a sign when `g` is
odd and changes neither when `g` is even, so `(7a)` is orientation
independent.  Overlapping girth cycles simply add: at total length `g` there
is no room for an extra edge or for a mixed longer contour.

The first four controls are

```text
C3: -i(z-z^-1) (I-S)^3,
C4:    (z+z^-1-2) (I-S)^4,
C5:  i(z-z^-1) (I-S)^5,
C6:   -(z+z^-1-2) (I-S)^6.
```

There is also a sharp spectral size formula for this first raw signal.  Since
`S` is unitary,

```text
||R_g||_2 = 2^g,                              g even,
          = (2 cos(pi/(2g)))^g,               g odd.     (8)
```

On the leading-rate slow scale, the contribution is
`Delta A_g/(2 Gamma^(g-2))`.  On one chordless cycle, if all its hoppings
have magnitude `h`, the elementary bounds on the sine/cosine factors give

```text
||Delta G_slow||_2 <= 4h^2(2h/Gamma)^(g-2),   g odd,
||Delta G_slow||_2 <= 8h^2(2h/Gamma)^(g-2),   g even.    (9)
```

Thus the displayed raw coefficient bound decays exponentially in `g` once
`Gamma>2h`; the binomial growth of the cyclic finite difference does not
defeat strong dephasing.  This coefficient bound is not by itself a
convergence criterion for the full Schur series, which still requires the
usual resolvent separation from the spectrum of `D`.

### The weighted `K_2,2` determinant readout

There is a literal determinant bridge, not just a cycle metaphor.  Given

```text
A=[a b; c d],
```

make a Hermitian hopping graph on row vertices `r1,r2` and column vertices
`c1,c2` by setting `H_(r_i,c_j)=A_ij`.  In the cyclic order
`r1,c1,r2,c2`, the Wilson product is

```text
W=a conjugate(c) d conjugate(b)=ad conjugate(bc),
R=|W|=|abcd|.                                           (10)
```

The leading dephased rates give the four conductances `|a|^2`, `|b|^2`,
`|c|^2`, and `|d|^2`.  Comparing the raw fourth word with the phase-free hopping
matrix having the same magnitudes gives

```text
Delta A_4=(W+conjugate(W)-2R)(I-S)^4.                   (11)
```

In particular,

```text
(Delta A_4)_(00)=4(Re(W)-R),
Re(W)=R+(Delta A_4)_(00)/4,                             (12)
|det A|^2=|ad|^2+|bc|^2-2R-(Delta A_4)_(00)/2.          (13)
```

Thus the leading resistor network supplies the two square terms and the
first `C4` flux correction supplies exactly their missing interference term.
For a Jacobian matrix this reconstructs `|det DF|^2` pointwise.  The bridge
is analytic rather than algebraic: it uses complex conjugation and does not
recover the phase of `det DF`, so it is a diagnostic sidecar, not a new
Jacobian-conjecture implication.  If an entry vanishes, the four-cycle term
vanishes and `(13)` reduces to the conductance data alone.

For a Keller matrix this becomes a useful squeeze rather than merely a
reconstruction formula.  Put

```text
r=|ad|,  s=|bc|,  kappa=det A.
```

Then the scalar in `(11)` is exactly

```text
q=2Re(W)-2rs=(r-s)^2-|kappa|^2,                         (14)
Delta A_4=q(I-S)^4.                                    (15)
```

The reverse triangle inequality gives

```text
-|kappa|^2 <= q <= 0.                                  (16)
```

After normalizing `kappa=1`, the entire raw fourth-order phase deficit has
`||Delta A_4||_2<=16` and its displayed diagonal readout lies in `[-2,0]`,
regardless of how large the four derivatives become.  At the same time,
when `rs>0`,

```text
1-cos(Phi)=[1-(r-s)^2]/(2rs),
|sin(Phi/2)| <= 1/(2 sqrt(rs)).                         (17)
```

Consequently any high-gradient Keller regime with `rs` large forces the
determinant Wilson phase into a very narrow window around constructive
alignment.  The absolute raw deficit stays bounded, but it is tiny relative
to the phase-free four-cycle sector of size `rs`; after imposing dephasing strong
enough for the largest hopping, its slow-scale signal is smaller still by
`Gamma^(-2)`.  This is the dynamical form of
[THM-3548](../01-canon/theorems/THM-3548-planar-keller-conductance-shadow-gates.md)'s
dark-plaquette branch: there the fixed determinant is hidden inside two
increasingly large, increasingly aligned channels.  Its channel-degenerate
or triangular shadow, where `r+s` stays bounded while opposite entries
separate in size, remains a distinct regime and is not squeezed by `(17)`.

The loss of orientation is forced, not an artifact of stopping at `A4`.
On any zero-onsite bipartite hopping graph, the sublattice sign gauge sends
`H` to `-H`.  Complex conjugation sends the dynamics of `H` to that of
`-conjugate(H)`, while fixing diagonal preparations, diagonal observations,
and local dephasing.  Combining the two symmetries shows that the entire
population channel for `H` equals that for `conjugate(H)`.  Hence pure
`K_2,2` population spectroscopy can recover `Re(W)` but never choose between
`W` and `conjugate(W)`.  An odd-cycle reference could expose the missing sine
coordinate, at the price of an extra calibrated hopping phase.  Determinant
modulus is special because it needs only `Re(W)`.

The exact companion verifies every lower word is phase blind, the formulas,
transpose parity, conservation, and rank for `C3` through `C6`; it also
checks `(11)` with four independent symbolic edge magnitudes:

- [script](../04-computation/dephasing_cycle_flux_girth_control.py)
- [stored output](../05-knowledge/results/dephasing_cycle_flux_girth_control.out)

Reproduction commands:

```bash
python3 04-computation/dephasing_cycle_flux_girth_control.py
python3 -O 04-computation/dephasing_cycle_flux_girth_control.py
```

Both executions reproduce the stored output byte for byte.  The current
SHA-256 hashes are

```text
f873637c59da1db34084efd714c45bfc17c5e1465a72b441bf8815c37e463ae9  script
c0630fce8a5c0790d3d53e26d7b5f14ce105dad95db1491c4b8240fb358e4137  output
```

## Hostile boundaries

1. **Trees.**  Every edge phase is removable by a vertex gauge, so the full
   node-population dynamics, not merely its Schur series, is phase blind.
   The companion checks a formal phase on `P6` through `A8` as a finite
   hostile control.
2. **Backtracking.**  Backtracks may contribute large phase-blind terms, but
   their two orientations cancel in the Laurent exponent.  They do not
   lower the girth threshold.
3. **Chords.**  The relevant length for a marked edge phase is the shortest
   cycle containing that edge, not the length of a visually preferred outer
   circuit.  On `C5` plus chord `0-2`, the phase on edge `4-0` is invisible
   in `A3` but already visible in `A4`.  A long chorded-loop holonomy is also
   a product of shorter cycle-basis holonomies, so it has no basis-independent
   claim to a later first-detection order.
4. **Parallel coherent channels.**  They form an effective two-edge cycle:
   amplitudes add before squaring, so their relative phase can already alter
   `A2`.  The simple-graph girth statement deliberately excludes this case.
5. **Special fluxes.**  At `Phi=0` every difference in `(4)` vanishes.  For
   odd `g`, it also vanishes at `Phi=pi`.  On the zero-onsite unit odd cycle
   this is not merely a leading-order accident: `H(pi)=-U H(0) U^*` for a
   diagonal vertex gauge `U`, while complex conjugation reverses `K` and
   fixes diagonal preparations and observations.  Hence the entire
   node-population channel cannot distinguish those two fluxes.  Nonconstant
   onsite potentials generally break this last `H -> -H` comparison.
6. **What is convention independent.**  Equations `(1)` and `(4)` are raw
   algebraic identities.  Laplace-memory elimination or a time-local slow
   normal form adds terms assembled from shorter words at the same asymptotic
   order.  Those additions are phase blind at the first girth order, so they
   cannot make Wilson data appear earlier.  An arbitrary phase-dependent
   near-identity change of slow coordinates can change even the apparent
   first order, because those coordinates are no longer literal node
   populations.  In the fixed population coordinates the invariant claims
   are the first possible order, the sine/cosine parity, and the raw
   coefficient `(4)`; no unique all-orders Markov generator is asserted.
7. **Nonuniform dephasing.**  The closed-walk/girth threshold survives when
   different coherences decay at different positive rates, but inverse
   decay weights are then interleaved between the commutators.  They weight
   different left/right step interlacings differently, so the universal
   binomial circulant `(4)` is special to uniform off-diagonal decay.

The practical lesson is sharper than “phases return on cycles”: the first
phase sidecar is a high cyclic finite difference `(I-S)^g`.  Longer coherent
cycles are therefore both later in `1/Gamma` and increasingly oscillatory in
the node coordinate.  That is the precise information erased when the
quantum walk collapses to its resistor network.
