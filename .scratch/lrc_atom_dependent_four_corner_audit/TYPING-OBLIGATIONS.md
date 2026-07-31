# Hostile typing ledger for the atom-dependent four-corner witness

**Scope:** independent pre-code audit of
`.scratch/lrc_atom_dependent_four_corner/INCOMING-POSITIVE-WITNESS.md`
against the common-atom contract in
`THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction`.

**Initial verdict before the corrected companion:** the displayed
finite-field square is a valid numerical positive control, but it is not yet
a typed common-atom witness.  The final replay resolution is recorded below.
The shortest missing datum is a commuting map

```text
{one labelled physical ancestry copy xi} x {0,1}_S x {0,1}_T
          |                                      |
          | physical carrier allocation          | endpoint chart
          v                                      v
source/target bare-carried states  ----->  ((L+epsilon*s,R+eta*t),address)
```

whose four coefficient evaluations are made before aggregation and before
Fourier/Radon pushforward.  In particular, the chosen steps `s,t` must be
the images of the two carrier toggles under this map.  Appending any
independent determinant-one pair after obtaining four amplitudes does not
meet THM-2772.

The following obligations distinguish that map from four marginals or four
post-hoc endpoint evaluations.

## Final resolution after the stable corrected replay

Stable witness source:

```text
SHA-256 8f064ff8922a787211cb66ad7f7a9cdfe5a87ded84100d7dbd39e264e387482f
```

The ordinary and independently launched optimized transcripts are
byte-identical (86 lines, 10,862 bytes, SHA-256
`a970776ed95128b5745c1fd370af768778b409d931a57b68006e04271e00813f`).
The obligations below resolve as follows.

- **A1--A5 pass.**  The multiplicity unfolds as the literal product
  `966606*28534=27581135604`; complete source/target inverse-label sets
  coincide; `(59162,26,56658)` is a retained path; the pullback is
  `(x,a,b,eprime)->(x+SHIFT,a,b,eprime)`; and exact chamber slacks keep the
  whole half-open interval away from every seam.
- **B1--B4 pass as a raw constructor.**  With fixed outer selector `xi`, the
  exact carrier-twist arrays are

  ```text
  B=w, P=delta_(a=0)w, Q=delta_(b=0)w,
  H=delta_(a=0)delta_(b=0)w.
  ```

  Multiplication by the two carrier masks commutes.  **B5 passes as a
  hostile and rejects promotion:** supports are `(169,13,13,1)` and the sole
  fourfold point has vector `(w,w,w,w)` and mixed face zero.  **B6 passes:**
  the moving-sheet marginals are computed separately and explicitly typed
  as an extrinsic orbit bank.
- **C1 and C5 pass.**  The fixed base lift and determinant sign convention
  are explicit.  **C2--C4 fail as an allocation lift:** the four
  determinant-one endpoint corners exist only for the separately moving
  sheet.  The literal fixed-sheet toggles do not produce those four
  addresses.
- **D1--D5 fail in the desired nondegenerate form.**  On fixed `xi`, present
  is exactly `1/13` times bare at every endpoint, so both induced steps are
  zero and their determinant is zero.  On the moving-sheet bank, exhaustive
  search over all `169` shifts, `169` characters, and a global scalar finds
  no bare/present covariance on either side in either field.  **D6 passes as
  the decisive hostile:** the nonzero endpoint samples are not carrier
  flags.
- **E1--E2 pass only for the retained fixed-sheet marginal, not for the
  nonexistent determinant-one allocation lift.**  The full 2,197-cell
  endpoint-array gauge and all 169 fixed-mask rows on each side descend when
  mask origin is retained.  **E3 passes:** fixing the carrier while changing
  representative gives 436 mismatches on each side.  **E4 passes:** the
  companion keeps the physical moving-sheet orbit distinct from gauge
  change.
- **F1 and the base-sheet part of F2--F3 pass.**  The simultaneous
  twelve-factor signatures are all true, the path label is common, and the
  source-carry-12/target-carry-6 unit values agree on the retained ancestry
  sheet.  **The four-translated-corner part of F2 and F4 remains
  unconstructed** because C2--D5 fail.  **F5 passes:** the final scope
  separates positive interval weight, nonzero finite-field coefficients,
  and physical current.
- **G1, G1a, G3, and G4 pass.**  Both coefficient mixed faces, the separate
  determinant-label face, Segre/Pluecker identities, and the nonzero-versus-
  positive distinction are printed exactly.  **G2 passes only as two
  specializations of the same pinned cyclotomic constructor;** no ordered
  characteristic-zero positivity or LRC consequence is inferred.
- **H fails, as it must.**  Common cells occur on only clocks `1,2,3`, with
  census `(0,81,56,56,0,0,0)`.  There is no seven-clock/root-deck map,
  physical Cech correction, row exclusion, or LRC(14) conclusion.

The first missing object is earlier than the originally anticipated
endpoint map.  For Boolean carrier values `C_S,C_T`,

```text
(B,P,Q,H)=w(1,C_S,C_T,C_SC_T),
B-P-Q+H=w(1-C_S)(1-C_T).
```

Hence a pure membership indicator has zero mixed face on every common-
present atom.  A future construction first needs a non-idempotent same-atom
amplitude (phase, sign, character, cocycle, overlapping translate, or
equivalent local-system datum); only then can D1--D5 ask whether its two
toggles map to independent endpoint translations.

## A. One physical copy, not one unit of an aggregate weight

- **A1 -- explicit fibre.**  Give a finite labelled set `Xi_I` over the
  weighted interval
  `[142004992589460,142005019034340)` and prove
  `#Xi_I=27581135604`.  A step-profile value is presently only an integer
  multiplicity.
- **A2 -- section/copy.**  Name one `xi in Xi_I`, or give a deterministic
  section selecting it.  Replacing the multiplicity by the scalar `1`
  produces a normalized summand, not automatically a physical subatom.
- **A3 -- source/target equality.**  Exhibit the pullback map on labelled
  copies and prove that the retained source copy and pulled target copy are
  literally the same `xi`, not two unnamed copies from equal-cardinality
  fibres.
- **A4 -- copywise weight.**  Show that this copy contributes weight one on
  every one of the four allocation corners.  The equality of aggregate
  source and target weights is insufficient.
- **A5 -- boundary convention.**  Pin half-open/open endpoint conventions
  and verify that `xi` remains in the interior under every required chart
  translation; no seam copy may be silently duplicated or lost.

The label carried by `xi` must include, rather than merely reproduce after
summation:

```text
rail and physical cell (s,t,e);
clock, delayed/coefficient clock, semantic word and owner;
source and target predecessor carry/root labels;
factor-owner and simultaneous pass/failure signature;
marked triangle, Abel-boundary component, and endpoint ordering.
```

## B. A literal Boolean allocation square

- **B1 -- four states at zero harmonics.**  Retain
  `(epsilon_S,epsilon_T)` as two explicit Boolean coordinates even when
  `k=l=0`.  Zero Fourier mode is not carrier absence.
- **B2 -- lawful toggles.**  Define source and target toggle maps
  `tau_S,tau_T` on the same `xi`; prove they commute and preserve every
  label not declared to change.
- **B3 -- same two factors.**  Define exactly two source values
  `P_0(xi),P_1(xi)` and exactly two target values
  `Q_0(xi),Q_1(xi)`.  The corner values must be constructed atomwise as

  ```text
  (P_0 Q_0, P_1 Q_0, P_0 Q_1, P_1 Q_1)
  ```

  before any sum, DFT, determinant/Radon transform, or endpoint-origin
  pushforward.
- **B4 -- amplitude ordering.**  State that the four displayed endpoint
  amplitudes are `(P_0,P_1,Q_0,Q_1)` in that order, and give their raw
  definitions.  “Endpoint amplitudes” alone does not type their roles.
- **B5 -- no virtual one-corner route.**  Print support bits for all four
  corners at this same pre-Fourier address.  A transformed vector
  `(0,0,0,h)` passes both the nonzero Hadamard gate and the Pluecker
  identity, so those identities cannot certify co-support.
- **B6 -- hostile marginals.**  Recompute the four values after independently
  forgetting the ancestry label in each corner and show why that route is
  rejected by a failed identity-of-labels check, even though it may retain
  the same numerical rank-one square.

## C. The exact pullback address at all four corners

- **C1 -- base lift.**  Print the complete base address
  `(r,k,l,L,R)` rather than only `k=l=0` and `L=R=(0,0)`.  Verify the
  distinguished-`e_0` representative condition for `r`.
- **C2 -- four lifts.**  Print the complete address at each allocation
  corner.  With the claimed conventions they must project to

  ```text
  00: L=(0,0), R=(0,0),  q=(0,0), Delta=0
  10: L=(0,1), R=(0,0),  q=(0,1), Delta=0
  01: L=(0,0), R=(12,0), q=(1,0), Delta=0
  11: L=(0,1), R=(12,0), q=(1,1), Delta=1.
  ```

  Since `q` changes, either list the four representatives `r_epsilon_eta`
  or give and verify the section constructing them.  One fixed `(r,k,l)`
  cannot represent all four corners.
- **C3 -- typed endpoint factors.**  Evaluate endpoint factors at those
  four endpoint pairs inside the raw atom constructor.  Relabelling four
  already-computed marginal amplitudes by these points is not a lift.
- **C4 -- retained origin.**  Retain the base endpoint origin before the
  determinant/Radon pushforward.  The four translated `R_eta` values must
  arise from one base `R`, not four independently selected origins.
- **C5 -- sign and order.**  Pin the source/target convention
  `q=L-R`, `L_epsilon=L+epsilon*s`,
  `R_eta=R+eta*t`; then check
  `q_epsilon_eta=q+epsilon*s-eta*t` and
  `det((0,1),(12,0))=-12=1 mod 13`.

## D. The load-bearing allocation-to-translation map

- **D1 -- derive, do not append, `s`.**  Show that applying `tau_S` to
  `xi` changes its endpoint data by the translation `L -> L+(0,1)`.
  THM-2772's bit is an absent/present **carrier flag**, so this requires a
  bare bank `P_0` and carried bank `P_1`, together with the fibrewise
  identity

  ```text
  P_1(L)=chi_S(xi,L) P_0(L+s)
  ```

  where any scalar/character `chi_S` is explicitly tracked and `s` is
  derived from the primal carried factor by Fourier translation.
- **D2 -- derive, do not append, `t`.**  Show that applying `tau_T` to
  `xi` changes its endpoint data by the translation `R -> R+(12,0)`.
  Likewise require

  ```text
  Q_1(R)=chi_T(xi,R) Q_0(R+t)
  ```

  fibrewise, with `t` derived from the primal target-carrier character.
- **D3 -- unchanged side.**  The source toggle must leave the target factor
  and `R` fixed; the target toggle must leave the source factor and `L`
  fixed.  Otherwise the factorization in B3 is not the physical square.
- **D4 -- commutation.**  Both orders of toggling must produce the identical
  labelled `11` corner, not only the same `q` and determinant.
- **D5 -- normalization is intrinsic.**  Independence and determinant one
  must follow for these derived steps.  It is not enough that some freely
  chosen pair in `F_13^2` has determinant one.
- **D6 -- samples are not flags.**  Four fitted samples of a single endpoint
  transform,

  ```text
  P(L), P(L+s), Q(R), Q(R+t),
  ```

  establish a nonzero endpoint parallelogram but do not establish the
  absent/present allocation square.  The companion must construct the
  bare/carried banks independently, prove the identities in D1--D2
  globally or on the declared physical fibre, and only then evaluate the
  four corners.  Defining `P_1:=P_0(L+s)` or selecting `s,t` after seeing
  nonzero samples is circular.

This is the most likely minimal missing map in the incoming note.  The note
currently says “use” `s=(0,1),t=(12,0)`; it does not yet identify these
vectors as images of the physical allocation toggles.

## E. Representative gauge and quotient descent

- **E1 -- exhaustive covariance.**  Over all `13^3=2197` representative
  cells, transport together the copy `xi`, carrier address, endpoint
  origin, four endpoint factors, allocation bits, and semantic labels.
- **E2 -- cornerwise descent.**  Check representative covariance for each
  of the four banks and their corner addresses, not only for the aggregate
  source and target carrier arrays.
- **E3 -- fixed-carrier hostile.**  Include the known hostile in which the
  representative changes but endpoint factors are held fixed, and require
  failure.  Otherwise an aggregate invariant may conceal the missing
  endpoint transport.
- **E4 -- chart versus gauge.**  Distinguish an actual adjacent physical
  chart translation from a change of representative.  Gauge covariance
  alone does not manufacture the physical toggle map D1--D4.

## F. Semantic and factor ancestry

- **F1 -- simultaneous factor signature.**  Define
  `common_after_all_factors` by a simultaneous pass vector.  An ordered
  first-failure bucket is not an intrinsic factor label when several
  factors fail at once.
- **F2 -- common after all factors.**  Verify that the selected `xi` passes
  every inherited source factor and every pulled target factor at all four
  corners, with the complete pass vector printed.
- **F3 -- exact cospan.**  Explain how source carry `12`, root `1`, and
  target carry `6`, root `1` are two faces of the same labelled ancestry
  copy.  Equality of their delayed coefficient
  `103478815440` is numerical evidence, not an identity of owners.
- **F4 -- semantic preservation.**  Verify the same clock, word/owner,
  marked triangle, and Abel-boundary component at all four corners.
- **F5 -- no chronology claim.**  Keep separate:
  positive physical interval mass, nonzero exact semantic coefficient,
  nonzero finite-field amplitudes, and a physical current compatible with
  chronology.  None implies the next without a typed map.

## G. Arithmetic and field specialization

- **G1 -- exact recomputation.**  Independently recompute products,
  Hadamard coordinates, Segre, Pluecker, contrasts, and determinant in both
  certified fields.
- **G1a -- separate the two mixed faces.**  There are two nonvanishing
  statements of different types:

  ```text
  coefficient mixed face:
    D3=v00-v10-v01+v11=(P0-P1)(Q0-Q1)

  geometric determinant mixed face:
    delta_Delta=Delta00-Delta10-Delta01+Delta11
               =det(s,t)=1 in F13.
  ```

  `D3` is the third Hadamard coordinate of the **current coefficients**.
  `delta_Delta` is the Mobius boundary of THM-2620's **determinant-sector
  label**.  Neither equality identifies these two objects.  A physical
  promotion must state the map/function that pairs or transports the
  coefficient face through the determinant sector; otherwise it has proved
  two adjacent nonzero controls, not one determinant current.
- **G2 -- same integral/cyclotomic source.**  Pin the characteristic-zero
  object specialized into both fields, the evaluation/root choices, and
  the ordering of embeddings.  Nonvanishing in two residue fields proves
  a characteristic-zero element is nonzero only when both values are
  certified images of that same element.
- **G3 -- do not overread Pluecker.**  Record that
  `D_0 D_3=D_1 D_2` is automatic for any numerical rank-one square and also
  passes the post-Fourier virtual-one-corner hostile.  Its value is supposed
  to be **zero**; it certifies decomposability and is not either nonzero
  mixed face in G1a.  The actual nonvanishing checks are `D3 != 0` in each
  certified field and `delta_Delta=1` in `F13`.
- **G4 -- physical versus arithmetic positivity.**  A finite-field element
  is nonzero, never positive.  “Positive witness” may refer only to the
  physical multiplicity/interval unless an ordered characteristic-zero
  coefficient is exhibited.

The displayed arithmetic itself checks:

```text
field 352341050142921841:
  Segre=0, Pluecker=0, D3=170114988031260853
  (P0-P1)(Q0-Q1)=170114988031260853

field 956354278959359281:
  Segre=0, Pluecker=0, D3=757304766188814060
  (P0-P1)(Q0-Q1)=757304766188814060

det((0,1),(12,0))=1 mod 13.
```

## H. Exact promotion boundary

Even a pass through A--G proves at most:

```text
one physical ancestry atom
  -> one pre-Fourier four-allocation coefficient square
  -> one nonzero normalized determinant mixed face at one clock.
```

It does **not** by itself prove:

- seven coherently oriented clock faces or `sum_j c_j=-7a`;
- a carrier-equivariant root-deck functional `pi_j:mu_j -> c_j`;
- a common physical Abel/Cech correction;
- semantic arrival, all-91-unit current, row exclusion, or LRC(14).

For those stronger claims the companion must additionally give the
seven-clock/root-deck map and audit the full THM-2772 invoice.

## Replay requirements for the pending companion

The standalone file must be dependency-pinned, assertion-free, deterministic,
and replay identically in ordinary and optimized Python modes.  The audit
will:

1. pin its LF-normalized hash and every dependency hash;
2. compare ordinary, `python -O`, repeated, and stored transcripts;
3. independently reconstruct the selected physical fibre and all four raw
   evaluations without importing witness conclusions;
4. run the virtual-one-corner, four-independent-marginal, fixed-carrier,
   and seam/boundary hostiles;
5. exhaust all `2197` representative gauges for every corner; and
6. report the first failed typing gate rather than promoting a downstream
   numerical consequence.
