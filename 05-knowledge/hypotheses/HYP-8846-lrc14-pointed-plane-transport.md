---
id: HYP-8846
title: "LRC14 pointed transport on the rank-eleven relation planes"
status: >
  OPEN finite completion program after a proved large-direction theorem.
  THM-2052 reduces every hypothetical counterexample to a finite atlas of
  rational planes. THM-2053 now proves that every two-dimensional atlas cell
  has torus margin at least 1/13. Its exact safe gate is
  max_i|a z_i-b u_i|<=(a^2+b^2)/91; the round corollary is ||(a,b)||>=91L.
  Its adjacent-column refinement gives M(v)>=1/13-R/(2N), where N divides an
  exposed pair sum. Its exact residue deck D_N(m) eliminates whole N-fibers
  independently of longitudinal coefficients; THM-2055 then cuts the bad
  fibers into hull-owner tangent sectors. THM-2058 reduces each such
  intersection to one exact coprime interval in the longitudinal coordinate.
  The missing work is uniform atlas compression and discharge of the surviving
  rows, not an infinitary pointed transport theorem.
source: codex-2026-07-21-DC2-LRC14-termination
related:
  - THM-2050
  - THM-2052
  - THM-2053
  - THM-2054
  - THM-2055
  - THM-2056
  - THM-2057
  - THM-2058
  - THM-2059
  - HYP-2108
  - HYP-1842
  - HYP-1977
  - HYP-2647
  - HYP-2896
  - HYP-2986
  - HYP-3267
  - HYP-4336
  - HYP-4346
  - HYP-8841
  - HYP-8845
  - HYP-8850
  - HYP-8860
  - HYP-8871
  - HYP-8877
  - HYP-8878
  - HYP-8879
  - HYP-8885
  - HYP-8890
---

# HYP-8846 -- pointed transport, not an unpointed plane escape

The original proposal in this file had a real quantifier gap:

```text
known:  target v lies in a rational plane ker(W);
known:  rank rigidity supplies some escaping direction;
missing: "some direction is safe" does not imply the specified v is safe.
```

THM-2053 repairs the quantifier on all sufficiently long directions. If
`u,z` is a saturated integer basis of the plane, put

```text
L=max_i sqrt(u_i^2+z_i^2),
v(a,b)_i=a u_i+b z_i.
```

For every primitive `(a,b)`, the corresponding subcircle is
`1/(2sqrt(a^2+b^2))`-dense in the parameter two-torus. The phase-height
objective is `L`-Lipschitz. Keeping the actual normal displacement gives the
sharper error

```text
E(a,b)=max_i |a z_i-b u_i|/(2(a^2+b^2)).
```

Normalize each coefficient column `c_i` by the positive target value `v_i`.
Two adjacent normalized columns have a sum
parallel to no column. A perpendicular integral covector is therefore
full-support and has a repeated absolute speed, so settled LRC for at most
twelve distinct speeds gives two-torus margin at least `1/13`. Therefore

```text
max_i |a z_i-b u_i| <= (a^2+b^2)/91
                         ==> M(v(a,b))>=1/14.           (1)
```

In particular, `sqrt(a^2+b^2)>=91L` implies (1). This is pointed: it applies
to the specified `(a,b)`, not merely to a different direction in the same
plane.

The same adjacent pair gives a sharper intrinsic coordinate. If
`c_i+c_j=g w_0`, extend primitive `w_0` to a lattice basis and write
`c_k=a_k w_0+m_k eta`. Put

```text
N=phi(w_0),       R=max_k|m_k|.
```

The actual orbit samples the complete `N`-grid on the transverse circle, so

```text
M(v)>=1/13-R/(2N).                                    (2)
```

Hence `N>=91R` is safe, while every unresolved row has `N<91R`. Moreover
`N|(v_i+v_j)` and `v_k==m_kM (mod N)` for one unit `M`; this is a pointed
pair-sum-ruler cell, not merely a Euclidean disk.

The exact first filter is

```text
D_N(m)=max_(ell mod N) min_k |ell m_k|_N/N.             (3)
```

Every row in the cell satisfies `M(v)>=D_N(m)`, with no dependence on its
longitudinal coefficients or unit `M`. Bad `N` form a divisibility down-set.
The rational phase-height components of the template give a sharp conductor:
an interval `[alpha,beta]` meets the `N`-grid exactly when
`ceil(N alpha)<=floor(N beta)`, so HYP-3456's floor-count operation resolves
the conductor residue class by residue class. The generic `N<91R` region is
therefore only a crude outer box.

THM-2058 resolves the internal arithmetic of this filter. The safe packet has
the exact reduced-order split

```text
S_N(m)=disjoint_union_(d|N) (N/d)S_d^prim(m),
L_N=sum_(d|N)p_d,
```

and the entire labelled packet transports across the longitudinal fiber by
`M^(-1)`. For fixed bad `N` and hull owner, the positivity, owner-cone, and
determinant inequalities cut `M` to one explicit interval; primitive rows are
the coprime integers in that interval, minus explicit speed-collision points.
Non-hull representatives of a repeated normalized-column value must remain in
the deck sidecar even though THM-2055 may delete them as determinant owners.

## The revised terminal

The rank-eleven program now has two finite terminals rather than one:

```text
rank 12       -> maximal-minor box (THM-2052),
rank 11/plane -> parameter disk ||(a,b)||<91L,
                  sharpened to N<91R on an exposed pair-sum ruler (THM-2053).
```

The intended per-target routing remains useful inside those finite cells:

1. some peel has the HYP-2108 endpoint functional `P_w>=0`;
2. active owner data supplies a bounded relation outside `W`; or
3. an exact resolved phase/topological certificate discharges the row.

There is no longer a need to postulate an infinite wall-crossing transport to
close large directions. Wall transport is a finite-cell accelerator.

## Concrete completion algorithm

For each rank-eleven bounded triple code from THM-2052:

1. compute a saturated kernel basis by Smith/Hermite normal form;
2. lattice-reduce the basis; in each positive target chamber, order the
   normalized values and retain a divisor-covering sweep of labelled adjacent-
   pair representatives, including `(g,N,m)` component data;
3. compute the primitive phase packets and Ehrhart conductor; discard every
   whole `N`-fiber with `D_N(m)>=1/14` and transport all owner labels by
   `M^(-1)`;
4. compute the signed column hull and THM-2055 normal fan, and intersect each
   remaining fiber with its active-owner tangent disk;
5. use THM-2058 to enumerate coprime longitudinal coordinates in the resulting
   interval, then delete the explicit nonpositive and collision walls;
6. run the exact pair-sum phase-height test, beginning with the exposed ruler
   `N|(v_i+v_j)` and its template residues;
7. retain `P_w`, endpoint-owner, peel-tax, and relation-rank labels for any
   survivors, then quotient only by transformations preserving those fields.

Each tangent-disk boundary contains the origin, but the actual carrier is the
intersection of its owner sector with a labelled bad phase packet. The packet
must be computed before determinant hull deletion.

This is structural compression, not the first proof of finiteness: THM-763
already supplies the global ceiling `sum v_i<=91^12`. The potential gain is
that the plane gate discards whole projective tails before raw speed
enumeration and retains the geometry needed for resonance fans.

The finite atlas is astronomically large if generated as all bounded triple
rows. The theorem-facing next step is an atlas-compression lemma on rank-two
column configurations with determinant magnitudes and labelled normalized-
fiber representatives. A sign-only matroid is insufficient.

## What the unrelated repo work contributes

- HYP-4342 supplied the `(1,N)` Lipschitz/net mechanism; THM-2053 extends it to
  every primitive direction and the adjacent-column argument identifies an
  exact transverse `N`-grid on the specified row.
- THM-2058 imports HYP-3456's Beatty clock at count level and THM-685's grid
  discrepancy at primitive-order level. Its Mobius packets are phase strata,
  not THM-2041 Fourier projectors; unit Frobenius creates no safe seed.
- It also makes the repo's cusp diagnostic exact without modular forms. If the
  template height is strictly above `1/14`, every sufficiently large prime
  grid has a primitive safe phase. At equality, primitive support is finite
  and lies on explicit denominators `q` with `14|q` and `q|(x+y)<=2R`. Below
  threshold all packets vanish. This bulk/boundary/null trichotomy retains the
  signed pointwise predicate that MISTAKE-226 shows the level-14 newform lacks.
- The orbit-product argument from
  `THM-1605-tnc-proved-monodromy-transitivity.md` does transfer exactly to a
  primitive packet. If `H` is its unit stabilizer, there are
  `phi(N)/|H|` longitudinal packet images, each primitive phase has uniform
  incidence `|A|/|H|`, and their commutative product is a full-unit norm.
  Unlike TNC, there is no nonconstant `ct` factor to contradict that norm.
  Since `-1 in H`, the quotient also proves that signed orientation is a
  necessary sidecar rather than optional bookkeeping.
- HYP-4346 supplied the rank-two algebra but also exposed the wrong-quantifier
  trap. It is now optional acceleration, not the bridge.
- HYP-2896 is the scale-one model of the finite-disk fan. THM-2057 now closes
  the entire integral plane `span((1,...,13),e_12)`: the scaled walls are
  `12a|w` and `14a|w`; their cells carry a central-unit witness on the `12a`
  or `14a` clock, while simultaneous killing forces `84a|w` and the scaled
  affine binding phase `(35m+2)/(a(84m+5))`. In the basis
  `u=(1,...,13),z=e_12`, writing `w=12a+b`, THM-2053's determinant is
  `max(13|b|,|a-12b|)`, and THM-2056 Kelvin-inverts its gate to a rational
  polar parallelogram. The clock theorem, not the determinant gate, certifies
  every residual direction.
- The same THM-2057 sieve closes `span((1,...,13),e_13)`, equivalently every
  row `{a,2a,...,12a,w}`. Its missing clocks are `13,14`; simultaneous killing
  forces `182a|w`, and the explicit deep-well phase
  `14m/[a(182m+1)]` is strict. Thus both AP one-tail coordinate planes now
  have symbolic clock/binding closures.
- HYP-2986 supplies the faithful three-state terminal: open tope, boundary
  cocircuit, or forbidden wall packet.
- HYP-2647 supplies the addressed wall-transport matrix for moving between
  neighboring parameter cells. Its scalar signed total is not enough.
- HYP-1842 suggests private endpoint pivots as relation-rank suppliers, but a
  binary private pivot does not yet imply a rational relation outside `W`.
- HYP-1977 proves that even observer-pointed tournament classes can mix safe
  and unsafe cells; endpoint lengths/owners must remain in the quotient.
- HYP-3267's contact holonomy is a useful connection coordinate, but its
  empty/full collision prevents it from replacing the endpoint cell.
- HYP-8845/HYP-8850 parity duplicates a first survivor into its mirror. It
  cannot create that first survivor and therefore belongs after the pointed
  gate.
- THM-2054 is the proved analytic complement: bounded scalar resonances that
  all lift to vector-character relations give a whole-product Fejer
  decorrelation estimate without summing the relation lattice. Its
  model-specific full-torus identification debt remains explicit.
- HYP-8860's Paley-prime table usefully assigns roles to moduli `3,7,11`
  (resonance atom, period-14 apex, rank scale). It is a modulus-selection lens,
  not a carrier: Paley orientation discards the signed coefficients and
  endpoint owners needed to decide a tangent-disk point.
- THM-2055 replaces the raw tangent-disk union by the signed column polygon's
  normal fan. Only hull vertices own the determinant maximum; each owner cone
  has one disk and an owner-local radius. THM-2056 Kelvin-inverts those disks
  to one polar polygon and proves an acute-unimodular defect certificate for
  whole Farey cones. HYP-8871 now couples this exact address to scaled safe
  clocks, killed-clock divisibility, affine binding, and Euler sidecars.
  MISTAKE-225 records why Heegner form classes cannot replace the carrier.
- THM-2058 adds the exact phase-order packet and longitudinal-interval carrier
  inside each bad transverse deck. This is complementary to THM-2057: the
  former decomposes a fixed bad denominator, while the latter removes whole
  one-tail families by missing-clock divisibility before interval enumeration.
- THM-2059 proves the missing core/tail composition: their safe packets meet
  by an exact CRT reduction-histogram dot product for arbitrary `N`, with a
  zero-mode/finite-Fourier fluctuation split. THM-2057 is its rigid `N<=14`
  nonemptiness specialization; THM-2058 supplies primitive-order and owner-
  interval structure inside the resulting packet.
- HYP-8877/8878's GMC confinement, HYP-8885's cusp frame, and HYP-8890's
  saddle route are schedulers, not proof transfers. A unique antipodal
  maximizer can still have height `2/29<1/14`, and exact Beatty counting already
  replaces asymptotic saddle counting on one LRC packet. THM-2058's CRT,
  spanning-tree, mask-walk, fixed-jet, and unique-maximizer controls falsify
  any attempt to use resonance, uniqueness, or genus without a height/owner
  sidecar.
- Two further archive controls sharpen that warning. Tournament `c_3` is
  already spectral (`tr(A^3)=3c_3`), so it cannot recover a hidden cospectral
  sidecar. The tight AP row `{1,...,13}` has six ordinary two-wall transverse
  maxima, not a higher-multiplicity Jacobian-style cusp. Here “boundary” means
  zero phase measure plus exact denominator data, not singularity multiplicity.
- MISTAKE-227 gives the matching Fourier control: HYP-8879's integer-kernel
  integral is a legitimate strict-bulk observable after convergence control,
  but it vanishes on both unsafe rows and the tight AP boundary. Its finite
  Sidon/AP ratios neither reduce LRC to AP cores nor replace THM-2058's exact
  denominator/height packet.

## Tournament analysis

Vertices are proof carriers, not runners or raw relation rows:

```text
torus_geodesic_disk
rank12_maximal_minor
resolved_phase_height
endpoint_owner_tope_cocircuit
addressed_wall_transport
private_pivot_rank_supplier
contact_holonomy
unpointed_rank_rigidity
raw_relation_count
```

Pairwise observable: `(pointed target retention, finite terminal, exact LRC
predicate, owner/magnitude retention, cost)`. The gauge prefers a carrier only
when it preserves the specified row and supplies either a finite domain or an
exact exit. The tie path begins

```text
torus_geodesic_disk
> rank12_maximal_minor
> resolved_phase_height
> endpoint_owner_tope_cocircuit
> addressed_wall_transport
```

and places unpointed rank rigidity and raw relation count last. This quotient
preserves the target-row LRC predicate and destroys raw relation-matrix identity;
the destroyed identity must be recoverable through the saturated kernel basis.

## Honest remaining statement

LRC(14) is still open. THM-2053 proves that the two-dimensional branch is
finite in parameter space, not that the finite transverse cells are empty of
counterexamples. The next decisive target is:

> Bound and classify reduced saturated bases and transverse templates of
> THM-2052's bounded support-three rank-eleven codes strongly enough that every
> THM-2058 coprime interval inside an exact bad deck and THM-2055 tangent owner
> sector is either empty, covered by the existing finite window, or splits into
> finitely many HYP-2896-style resonance fans with explicit phase-height
> certificates.
