# The underlying object is the support Dirichlet profile, not one reciprocal sum

*codex-2026-07-21. Companion to THM-2000/2005, correcting THM-1985/1990
through MISTAKE-209 and extending Downey--Ong--Sellers.*

The owner's phrase "the reciprocal of an integer sequence is a subset of the
harmonic numbers" changes the object in a precise way.  An indexed sequence
can repeat values; a subset cannot.  Once this distinction is respected, the
right hierarchy is

```text
indexed sequence
  -> positive support plus multiplicities
  -> support Dirichlet/counting profile
  -> block occupancies and singularity type
  -> the single scalar value D_A(1).
```

The last arrow is a severe quotient.  The triangular numbers and
`{1,2,3,6}` both have mass two.  The support `{2^j:j>=0}` also has mass two, and
independent Egyptian refinements generate continuum many further examples.
Those refinements increase `D_A(s)` below `s=1`, preserve it at one, and
decrease it above one.  The scalar is not a fingerprint: it is the
conservation hyperplane of a profile-valued flow.

## What the archived paper was really computing

Downey--Ong--Sellers use Maclaurin series, integration by parts, real partial
fractions, and trigonometric endpoint limits to evaluate even polygonal rows.
The shorter invariant is a pair of harmonic clocks:

```text
1/P_s(n)=2/(s-4)[1/(n-1+2/(s-2))-1/n].
```

Their finite trigonometric formula is Gauss's rational-digamma formula in
coordinates.  This reveals four facts that the paper's presentation hides:

1. odd polygons require no new argument;
2. the square/Basel row is where the two clocks coalesce, so digamma becomes
   trigamma;
3. polygonal support is one axis of the master array
   `(s-2)C(m+d-2,d)+C(m+d-2,d-1)`, whose other axis is the simplex ladder.
4. centered polygonal numbers form a second digamma clock, with the centered
   octagonal row confluent at `pi^2/8`.

The two-dimensional beta integral for that array is the first unified object
containing both the owner's triangular `2` and the paper's polygonal table.

## The exact boundary is an occupancy law

For support counting function `N_A(x)`, Abel summation says

```text
sum_(a<=X,a in A)1/a=N_A(X)/X+integral_1^X N_A(t)/t^2 dt.
```

On dyadic blocks, convergence is equivalent to summability of

```text
#(A intersect [2^k,2^(k+1)))/2^k.
```

This is the missing structural coordinate behind the vague word "growth":

| support type | relative block occupancy | verdict |
|---|---:|---|
| positive density | constant | diverges |
| `n log n` scale | `1/k` | diverges |
| `n(log n)^p` scale | `1/k^p` | converges iff `p>1` |
| quadratic / polygonal | `2^{-k/2}` | converges |
| uniformly bounded block occupancy (for example powers of two) | geometric `O(2^-k)` | converges rapidly |

The Abel--Dini lift makes the boundary recursive rather than merely
classificatory.  Starting with any divergent support, multiply its `k`th
denominator by a power of its accumulated harmonic mass and round upward.
The exponent one is critical; iterating it builds the full Bertrand tower of
integer supports.

## Sequence atlas after the correction

| family/thread | support-harmonic statement | status |
|---|---|---|
| simplex `C(n,k)` | `k/(k-1)` with exact finite remainder | proved |
| polygonal `P_s(n)` | digamma difference; square is `zeta(2)` | proved, all `s>=3` |
| centered polygonal `1+kn(n-1)/2` | digamma difference; `k=8` is `pi^2/8` | proved |
| master figurate `(s,d)` | beta double integral, mass in `(1,d/(d-1)]` | proved |
| max tournament `c3` | `75/4-24 log 2` | proved |
| reducibility ceilings `M_3(n-1)` | same mass; cumulative product of THM-2016 condensation ratios | proved |
| condensation defects `3/(1-tau_c(n))` | word `6,5,8,7,...`; support `{5,6,...}` and profile `zeta(s)-H_4^(s)` | proved |
| Forcade-order tournament arcs | `2E-2`, `E` Erdős--Borwein | proved |
| labeled tournament denominators | Ramanujan `psi_R(1/2)` product | proved |
| switching denominators | same support; indexed row has one collision | proved |
| Sylvester / `Phi_6` | mass `1`, exact next-term remainder | proved in HYP-3724 |
| factorial / Fibonacci / Catalan | subtract duplicated initial `1` | exact collision correction |
| Euler pentagonal theorem | signed reciprocal powers, not polygonal reciprocals | separated |
| Moser circle / Vandermonde truncations | positive reciprocal defect | numerical atlas; tail certification next |
| fibbinary / Moser--de Bruijn | Mahler equations, dimensions, exact block tails | proved |
| primitive residues mod `q` | Möbius Dirichlet profile; THM-819 measure | proved |
| A000568 | Burnside prefix through `n=20`, tail width `<3.11e-44` | arithmetic nature unknown |
| other census rows | start-index-sensitive numerical prefixes | tail certification open |
| maximum-H census | known prefixes/lower bounds only | tail interval needed |
| full Hamiltonian-path value set | all odds except `7,21` remains conjectural | cannot support a theorem |
| LRC witness/modulus sequences | plausible Bertrand-scale candidates | growth law open |

The exact defect identity for any positive truncation `Q_n=P_n-Delta_n` is

```text
1/Q_n-1/P_n=Delta_n/(P_n Q_n).
```

So the reciprocal shadow of a Vandermonde tail is itself a positive weighted
tail, not merely a changed decimal.  This is the clean bridge back to the
repo's truncation and tail machinery.

The condensation word supplies the sharpest guardrail on the support thesis.
The defect denominators `6,5,8,7,...` and their sorted version `5,6,7,8,...`
have exactly the same support Dirichlet profile.  But applying the ordered
prefix-product transform gives partition sums `75/4-24 log 2` and `2`, whose
difference is the positive parity-shuffle tax `67/4-24 log 2`.  Thus the
support profile is the complete additive reciprocal object, not a universal
sufficient statistic for order-sensitive transforms.  The missing coordinate
there is the hazard word (equivalently, its prefix-survival path).

## Tournament viewpoint and challenged vertices

On a finite family of distinct supports, order first by scalar mass and break
ties by the least element of the symmetric difference (then by name if
needed).  The resulting total order produces only a transitive tournament: it
has no cycles, one Hamiltonian path, and no structural content beyond a
ranking.  Scalar mass alone cannot orient its many ties.  The assumption that
sequence values should be vertices is as lossy as assuming runners must be
vertices in LRC work.

We compared terms, distinct values, multiplicity collisions, dyadic blocks,
counting events, residue-class clocks, automaton states, Egyptian moves,
condensation hazards, and proof obligations.  Bounded-ratio occupancy blocks
are a faithful vertex choice that preserves the convergence predicate.  The
Dirichlet profile is their sidecar; chronological block scale is the tie
Hamiltonian path.  Projecting to `D_A(1)` destroys exactly the information
needed to distinguish Bertrand tails, automatic dimensions, and
Egyptian-equivalent supports; forgetting the hazard order additionally
destroys prefix-product partition functions.

## What remains worth computing

1. Normalize the remaining census prefixes printed in THM-1985/1990 by offset
   and deduplication, then derive tail bounds analogous to the completed
   A000568 certificate.
2. Full Dirichlet profiles for Moser-circle, shallow-diagonal, and the
   remaining automatic/tournament-census supports.
3. Singularity types at `s=1` for LRC modulus/witness sequences; this is the
   honest test of whether any live sequence inhabits the Bertrand boundary.
4. Egyptian-refinement components with block/counting sidecars: characterize
   which structural invariants survive mass-preserving splits.
5. Extend the existing sorry-free `SupportHarmonicFigurate.lean` kernel beyond
   its master/ordinary-polygonal, maximum-`c3`, and block algebra to collision
   taxes, centered-polygonal algebra, and automatic block recurrences.
   Digamma and infinite-series endpoints can remain layered above that core.

The conceptual gain is recursive: once a scalar is recognized as an endpoint
of a profile, the next question is no longer "what decimal is it?" but "which
quotient preserved the theorem predicate, and which tail information did it
erase?"
