# Product-Gamma depth-three selector holotopy and jet compression

**Status:** research reflection.  The depth-at-most-three statements are
proved candidates in THM-3149 with exact companion evidence.  The depth-four
numbers below are a fresh finite-exact scout signal and are not yet a proved
dependency.

## Inheritance pass

- **Closest proved mechanism:** THM-3144 replaces pointwise selector existence
  by the cumulative section `C_11=intersection_(N=5)^11 P_N` and kills that
  section by four exact upset rows.
- **Canonical hostile:** the THM-3137 same-mean laws show that scalar pole means
  do not determine Hasse positivity; selector variance and higher type data
  remain visible.
- **Corrected near miss:** the THM-3144 four-row coefficients do not extend
  automatically to depth three.  They have one positive coordinate, the
  legal state `(1,1,2)`.
- **Least-used sidecar:** THM-3147 transports the joint profile by partition
  length and number of singleton parts.  It turns out to see every active
  THM-3149 facet exactly.

## Anchor / niche / wildcard

- **Anchor:** determine whether a stochastic pole-prefix construction can
  give a lawful positive partition current strong enough to contribute to
  NC2/GMC.
- **Niche:** classify the small dual face supporting the finite persistence
  obstruction rather than enumerate the whole Hasse cone again.
- **Wildcard:** regard selector laws across degree and prefix depth as a finite
  convex analogue of a section/homotopy problem, and ask for its obstruction
  class rather than one preferred selector.

## Exact result and its first failed extension

At support `(1,3),I2`, the physical pole multiset has `8`, `33`, and `93`
states in depths one, two, and three.  The THM-3144 coefficient vector is
negative on 133 of the resulting 134 states and positive only at

```text
(1,1,2): 744310486097342352.
```

Thus the old proof fails at exactly one coordinate.  Solving the active dual
system again does not introduce a new upset.  The same four rows admit a new
primitive positive coefficient vector, and its combination is strictly
negative on all 134 states.  The maximum is attained exactly at

```text
(1), (7), (1,1,2), (1,2,2),
```

while the minimum occurs at `(5,5,6)`.

This is stronger information than “the LP is infeasible.”  Prefix depth adds
a single escape, which hands the active face from `(1,1,2)` to `(1,2,2)`, but
the combinatorial support of the dual obstruction does not move.

At support `(1,2),I2`, all `5+13+24=42` states through depth three are killed
by two rows.  The coefficient pair is

```text
(14903617223,641337),
```

and the maximum equality pair is `(1),(3,4,5)`.  Hence mixing stopping depths
does not repair the earlier portability wall.

## The surprising THM-3147 compression

All six active upset rows factor through `(ell(lambda),m_1(lambda))`:

```text
degree 8 singleton complement:       ell<=7,
degree 10 singleton complement:      ell<=9,
degree 11 singleton complement:      ell<=10,
degree 11 three-shape complement:    ell<=8 or (ell,m_1)=(9,8),

degree 6 singleton complement:       ell<=5,
degree 9 three-shape complement:     ell<=6 or ell=7,m_1<=5.
```

The last two are exactly the two possible ways a singleton coordinate chooses
between the pair of partitions in the next-to-bottom length layer.  Thus the
THM-3147 endpoint kernel can transport the entire THM-3149 dual carrier.

This suggests a more useful next object than another full upset enumeration:
derive the pole-subtraction recurrence for the six upset masses directly from
the bivariate endpoint jet, then search for a positive dual section uniformly
in prefix depth.  A success would turn huge Farkas integers into a small
matrix cocycle.  A failure would identify the first extra partition statistic
needed after singleton count.

## Holotopy picture

Let `S_<=d` be the physical state bank through pole depth `d`, and let

```text
C_D^(d)={probability laws on S_<=d satisfying every upset at N<=D}.
```

For fixed `d`, the degree fibres `P_N^(d)` need not be nested, so their common
section can die while a later fibre remains nonempty, as in THM-3144.  For
fixed `D`, adjoining depth gives an inclusion of simplices by zero extension,
so `C_D^(d)` can only grow.  A dual certificate is a finite obstruction to a
section over the degree diagram; a new prefix state can kill that particular
dual class without producing a global section.

That is the right sense in which “holotopy” is useful here.  It is not a claim
of a topological invariant.  It is a discipline for tracking:

1. the changing fibres;
2. the common-section space;
3. the active dual carrier;
4. the new states on which the old carrier changes sign; and
5. the sidecar needed to continue or replace the carrier.

## Depth-four frontier signal

A separate exact scan of all 200 multiplicity-valid four-pole states tested
the THM-3149 four-row combination.  It is positive on exactly two states,

```text
(1,1,1,1),
(1,1,1,2),
```

with respective values

```text
8436435455549610720197534651041779150185883913024444286544,
6866903242237987020516836110255775595251299607671856528832.
```

It is negative on the other 198 states, with minimum

```text
-954397559119488445253787199314430240350430392204559768800384.
```

This makes THM-3149 sharp as stated.  It does **not** show that the full
depth-four common section is nonempty or empty.  The cheapest positive test
is the 136-state bank obtained by adjoining the two crossing states; failure
there would still not exclude a law using other depth-four states.  The
decisive test is the full 334-state LP followed by exactification of either a
law or a new sparse dual.

## Connection specification

| connection | source | target | map | preserved | destroyed | missing sidecar | cheapest test |
|---|---|---|---|---|---|---|---|
| endpoint jet to selector wall | THM-3147 `(ell,m_1)` profile | THM-3149 six upset rows | sum profile coefficients over the cells above | exact upset mass and pole-subtraction typing | all finer shape coordinates | positive dual recurrence | symbolically evolve the six-row cone one pole |
| depth inclusion | `S_<=3` | `S_<=4` | zero-extension of laws | old feasible laws and old row values | no control on new coordinates | depth-four dual/law | full 334-state `C_11` LP |
| support portability | `(1,3),I2` state bank | `(1,2),I2` state bank | reuse invariant/bank and physical submultisets | current construction and upset duality | pole alphabet and state counts | support-natural selector | compare the two compressed jet cones |
| finite selector to GMC | averaged virtual-prefix current | original response | currently absent | scalar zero-mass normalization | sequential provenance and original-response identity | lawful stopping decomposition | solve the prefix-mixture identity before any asymptotic claim |

## Live concept board

1. **Four-facet carrier:** stable from depth two to depth three after changing
   coefficients.
2. **Two shallow depth-four crossers:** the entire next escape locus for the
   current carrier.
3. **Length-singleton jet:** exact low-dimensional representation of every
   active row.
4. **Cross-support two-row wall:** evidence against a support-independent
   stochastic selector.
5. **Original-response realization:** still missing; a convex selector theorem
   alone does not prove NC2.
6. **Sequential coherence:** unordered state laws need not be stopping laws of
   a pole-removal Markov process.

The best next move is not to guess another selector.  It is to compute the
small jet-level dual cocycle and compare it with the full depth-four LP.  If
both choose the same carrier, one may have a depth-uniform obstruction.  If
they diverge, the first new active statistic tells us exactly what the
length-singleton compression forgot.
