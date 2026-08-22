# LRC norm phase: physical factorization and mixed-fibre boundary

**Status: FINITE-EXACT SCOUT + SELECTED TYPE AUDIT, promoted after independent
hostile reconstruction as THM-3267.**  This note records an exact
factorization/no-go boundary. It is not a physical owner current, a row
exclusion, or progress closing `LRC(14)`.

## Inheritance pass and bounded question

The closest proved algebraic mechanism is **THM-3255, twelve-balance
multiplicative Singer rank defect and phase-marker boundary**, whose chosen
model

```text
F_169=F_13[u]/(u^2-2),  alpha=1+2u,  Norm(alpha)=6
```

identifies `j mod 12` with the twelve field-norm values.  The canonical
hostile is **THM-3234, Singer owner compactification and pointed Heisenberg
carrier gate**: its exact map `j -> alpha^j` destroys interval geometry,
phase endpoints, `q`-values, and physical ancestry.  The corrected near miss
is **THM-2806, literal fixed-sheet central allocation scalar law and
endpoint-translation no-go**, where retaining a literal sheet still gives a
raw-flat allocation square.  The least-used relevant sidecar is **THM-2791,
full arm-orbit transfer and lower central chord**, which retains an actual
same-ancestry partial rail-sheet germ.

The bounded question was whether the twelve-valued norm phase factors
through, or is preserved by, any already constructed physical owner
placement or endpoint-current map near THM-3246/3252/3253.

## Targeted map audit

The closest chain separates sharply.

- **THM-3246, all-dilation second-owner seam stabilization and sign word**,
  has actual interval cells `j=0,...,167` and exact physical cell masses, but
  supplies no canonical endpoint current or owner-to-root map.
- **THM-3252, Singer-compactified owner Hodge word universal charged
  cyclicity**, places `q_j` at `alpha^(b+aj)` only by an external coefficient
  relabelling.  It explicitly supplies no owner-to-endpoint map.
- **THM-3253, positive owner-mass Newton cyclicity and maximal common
  Heisenberg module**, uses genuine positive cell masses, but their placement
  on the Singer plane and central slice is an abstract relocation, not the
  canonical endpoint current.
- **THM-2791** is the closest actual physical ancestry object: a literal
  same-ancestry partial rail-sheet translation germ.  It explicitly stops
  before endpoint allocation and names the missing map from the fixed rail
  sheet to an endpoint-origin current.
- **THM-3247, Heisenberg central-Fourier decomposition and canonical current
  cyclicity**, is the closest canonical coefficient current:

  ```text
  J_q(R)=P_(R+q) Q_R,  q != 0.
  ```

  Its connection contract retains an abstract `(v,w)` identification and
  `H_13` action, while leaving ancestry, high address, section compatibility,
  physical sign, and the physical intertwiner open.

Therefore there is no proved physical map in this closest chain on which to
build the requested commutative square.  The square that would be needed is

```text
same-ancestry rail sheet S  -- physical germ -->  S
          | Phi                                  | Phi
          v                                      v
allocated endpoint (q,R)    -- current action --> (q',R').
```

The vertical map `Phi` is precisely what is absent.  By contrast, the bare
abstract square `j -> j+1`, `alpha^j -> alpha^(j+1)` commutes, but THM-3234
already records the physical information it destroys.

## Exact factorization ladder

The new exact census finds a useful three-level boundary.

### 1. Full increment retains all twelve phases

Once a full nonzero vector `q=alpha^j` and the chosen quadratic-field
structure are retained,

```text
phase(q)=log_6 Norm(q)=j mod 12.
```

Thus the full `q`-field is sufficient in the abstract Singer model.  This is
not yet physical: the current canon has no same-ancestry owner-to-`q` map,
and the field multiplication/norm is additional structure on the endpoint
plane.

### 2. Projective direction retains exactly parity

A projective direction is one exponent class modulo `14`.  Its twelve
vectors have exponents `r+14k`, so their norm phases are

```text
r+2k mod 12,  0<=k<12.
```

Consequently every direction contains six phases, each twice: the six even
phases when `r` is even and the six odd phases when `r` is odd.  The largest
quotient of `C_12` that factors through projective direction is exactly
`C_2` phase parity; the internal `C_6` coordinate is lost.  This is the
sharp survivor because `gcd(12,14)=2`.

The companion exhausts all `8,064` primitive Singer gauges and all
`112,896=8,064*14` projective fibres.  Both the source phase `j mod 12` and
the placed target norm phase have the same six-phases-twice pattern in every
fibre.

### 3. Determinant fibre retains no phase

For the THM-3247 coordinate

```text
v=det(q,R),
```

each fixed `q != 0` has exactly thirteen origins `R` over every `v`.  Across
all increments, each of the thirteen determinant fibres therefore has
`2,184` states and all twelve norm phases exactly `182` times.  The
determinant coordinate is maximally phase-blind.

Even the pair `([q],v)` does not recover more than parity: each of its
`182=14*13` fibres has `156` states, comprising six phases exactly `26`
times each.  Thus adding the canonical determinant fibre to projective
direction leaves the same sharp `C_2` boundary.

## Affine/Heisenberg hostile

The full affine point `q` carries phase as a function, but phase is not an
invariant of the standard `H_13` point action.  Its translation
`(r,w)->(r+1,w)` gives explicit witnesses

```text
alpha^0=(1,0) -> (2,0)=alpha^70,  phase 0 -> 10;
alpha^1=(1,2) -> (2,2)=alpha^40,  phase 1 -> 4.
```

The second witness destroys even the parity that survives projectivization.
Among the `167` nonzero-to-nonzero transitions for this translation, only
`13` preserve phase; every nonzero phase difference occurs exactly `14`
times.  Aggregated over all `168` nonzero translations, phase difference
zero occurs `2,184` times and every nonzero difference occurs `2,352` times.

So an eventual physical bridge cannot treat norm phase as an
`H_13`-invariant label.  It must carry phase as a transforming sidecar, or
weaken to a quotient compatible with the actual physical action.  Full
affine covariance leaves no nonconstant invariant phase quotient.

## Seam-transversal structure

THM-3246's negative seam is

```text
N={0,1,2,3,4,5,162,163,164,165,166,167}.
```

It is simultaneously:

- one owner in each of the twelve norm phases; and
- twelve distinct projective directions in every primitive Singer gauge.

In the canonical gauge its phase-to-direction matching is

```text
0->0, 1->1, ..., 5->5,
6->8, 7->9, ..., 11->13,
```

with directions `6,7` missing.  This does not create a phase marker: it is a
complete phase transversal, not a selected phase, and its endpoint placement
is absent.  It is nevertheless the most concrete niche object exposed by
the scout: the physical seam supplies a twelve-edge partial matching between
phase and projective direction.  A future bridge should test this matching
before doing another rank census.

## Sharp type boundary and next decisive test

The first missing physical field is not another coefficient rank.  It is an
allocated full endpoint increment and origin

```text
(q,R)
```

on the same literal ancestry sheet retained by THM-2791.  To recover an
absolute twelve-valued phase, that record must additionally carry a
compatible `F_169` Singer structure (or equivalent norm gauge).  Retaining
only `[q]`, `v=det(q,R)`, or both cannot suffice.

The cheapest decisive next test is therefore small and typed: on the two
cylinders of the THM-2791 partial germ, construct one candidate allocated
`(q,R)` record while retaining the literal Boolean contributor labels; then
check the square above under the actual germ.  Controls must include the
six-phase projective fibres, the uniformly mixed determinant fibres, both
translation hostiles, and the twelve-edge seam matching.  Failure should be
recorded at the first field that changes or becomes undefined.

## Reproduction

Run

```text
python3 04-computation/lrc_norm_phase_physical_factorization_no_go_scout_20260803.py
python3 -O 04-computation/lrc_norm_phase_physical_factorization_no_go_scout_20260803.py
```

and compare with
`05-knowledge/results/lrc_norm_phase_physical_factorization_no_go_scout_20260803.out`.
The script uses exact finite-field arithmetic, exhausts all primitive gauges,
replays the promoted THM-3246 companion, pins the THM-3252/3253/3255 exact
companions, and gates the selected nearest connection contracts.  It contains
no `assert` statements, floating-point literals, randomness, or cached search.
