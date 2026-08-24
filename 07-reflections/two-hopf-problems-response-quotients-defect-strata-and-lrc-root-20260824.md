# Two Hopf problems: response quotients, defect strata, and a polyhedral LRC frontier

**Long-session synthesis -- 2026-08-24.**  This session jointly examined

- Brendle--Hung,
  [*A metric on `S2 x S2` with positive sectional
  curvature*](https://arxiv.org/abs/2608.19068), arXiv:2608.19068v1;
- the unsigned, undated 108-page manuscript at
  [alpo.ge/s6.pdf](https://alpo.ge/s6.pdf); and
- the two different questions commonly called a Hopf problem in those
  sources.

The source statuses remain separate:

```text
Brendle--Hung headline       PREPRINT CLAIM / UNDER AUDIT
S6 headline                  MANUSCRIPT CLAIM / UNDER AUDIT
LRC(14)                      OPEN
targeted BH identities       VERIFIED / FINITE-EXACT / RELATIVE REPAIR below
THM-3990/3991/3993/3995      PROVED in their stated scopes.              (1)
```

The durable source records are the
[Brendle--Hung referee ledger](../05-knowledge/reference/CORE-PAPERS-BRENDLE-HUNG-S2XS2-2026-08-24.md)
and the [S6 referee ledger](../05-knowledge/reference/CORE-PAPERS-HOPF-S6-2026-08-24.md).
Nothing here promotes either global construction.

## 1. Inheritance pass

### Closest proved mechanisms

1. [THM-3990](../01-canon/theorems/THM-3990-componentwise-harmonic-obstruction-and-repair-quotient.md)
   computes the smooth Laplacian repair quotient and its finite/integral
   analogue.  It makes component means necessary and sufficient for a strict
   third-order Poisson repair and retains Smith torsion over `Z`.
2. [THM-3955](../01-canon/theorems/THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion.md)
   and [THM-3957](../01-canon/theorems/THM-3957-triple-normal-crossing-cotangent-conductor-kernel-and-normalization-cokernel.md)
   prove that normalization can erase conductor-induced cotangent quotient
   and torsion data at nodes and triple crossings.
3. [THM-3991](../01-canon/theorems/THM-3991-periodic-unimodular-toric-cusp-factorial-euler-obstruction.md)
   proves `chi(W)=d n!` for the periodic unimodular cusp grammar.  In the
   irreducible one-Euler-fibre sphere case it isolates `n=2`; it neither builds
   the S6 manuscript's analytic quotient nor recognizes its total space.
4. [THM-2047](../01-canon/theorems/THM-2047-phase-height-toric-arrangement-for-lrc.md)
   is the lossless local LRC carrier: owner, sign, selected side, phase, and
   height are retained.
5. [THM-2050](../01-canon/theorems/THM-2050-period14-top-germs-do-not-determine-global-loneliness.md)
   is the canonical globalization hostile: AP13 and a loose lift have equal
   germs at every period-14 top but different global maxima.
6. [THM-3910](../01-canon/theorems/THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response.md)
   is the nearest current global LRC response mechanism.  Its integer sheet
   carrier needs a cubic mixed moment after pairwise data fail.  Concurrent
   [THM-3995](../01-canon/theorems/THM-3995-scale-two-parity-hole-support-and-integer-variance-tariff.md)
   adds parity holes and a sharper integer variance tariff for the sole
   scale-two survivor, but does not close it.

### Canonical hostile examples

```text
smooth repair hostile       disconnected residual components with means +1,-1;
integral repair hostile     triangle Laplacian, coker_Z has Z/3;
conductor hostile           reverse one S6 paired-branch orientation, gain Z/2;
fibre-quotient hostile      conductor torsion kills an ambient-to-fibre image;
LRC local hostile           AP13 and loose lift have equal top germs,
                            maxima 1/14 and 1/12;
symbolic-state hostile      notebook variables Vbc and ra are overwritten. (2)
```

### Corrected near misses

- MISTAKE-471 forbids replacing the original nonnormal regular locus by the
  full normalization.
- MISTAKE-480 forbids treating the kernel of raw S6 twist coordinates as
  analytic moduli before quotienting presentation and geometric gauge.
- The phrase “two Hopf preprints” formerly obscured that the `S6` artifact is
  unsigned and undated.  Its host suggests an attribution but does not supply
  explicit bibliographic authorship.
- A local response formula is not a global LRC certificate.  THM-2050 is an
  exact failure witness, not merely a caution about proof style.

### Least-used relevant sidecar

The least-used shared sidecar is the **labelled incidence/owner graph at the
residual stratum**.  In the `S6` cusp it records which normalization branches
meet along which conductor curves and with which orientation.  In LRC it
records which signed owner walls bound a selected top cell.  Brendle--Hung has
no conductor graph: connectedness of its smooth `Sigma` makes
`ker(Delta_Sigma)` the constants, which is why one mean suffices there.

## 2. Anchor / Niche / Wildcard portfolio

| lane | question | cheapest decisive operation | outcome |
|---|---|---|---|
| **Anchor -- source audit** | do the two frontier claims survive their displayed exceptional-stratum calculations? | source-level dependency trace, fresh exact hostile computations | two BH notebook omissions localized; `h_a` independently repaired, `V_bc` symmetry-repaired relative to `V_ac`; both headlines remain under audit |
| **Niche -- LRC response** | what is the honest polyhedral analogue of smooth normal minimization? | eliminate one labelled rising/falling wall pair and attack AP13 | [THM-3993](../01-canon/theorems/THM-3993-labelled-two-wall-polyhedral-elimination-and-local-top-response.md) **PROVED**; no global LRC consequence |
| **Wildcard -- scale/dimension** | does the S6 fan grammar extend to higher torus rank? | Euler count by toric orbit stratification | inherited THM-3991 gives `chi=d n!`; the irreducible sphere grammar isolates complex dimension three |

The niche produced the new theorem.  The wildcard produced an obstruction,
not a construction.  The anchor improved the referee state without confusing
a repaired cell with the paper's theorem.

## 3. Seven-concept board

The live board was deliberately small enough for every new computation to be
compared against every lane.

| concept | representation | invariant/consumer | decisive loss |
|---|---|---|---|
| `Z`, the BH zero-plane bundle | smooth submanifold of `Gr_2(TM)` | minimum sectional curvature | a pointwise plane before normal minimization |
| `Sigma`, the residual torus | connected closed defect locus | cubic mean in `coker Delta` | open generic pieces are not closed components |
| `W sup D sup T`, the S6 cusp | normalization + conductor curves + triple points | topology and differential extension | normalization erases branch incidence |
| S6 clutch row `ell` | integral presentation matrix | Smith group and `p=12ell_0-4ell_1-3ell_2` | raw coordinates include gauge |
| LRC labelled wedge | two owner/sign walls | displaced local top | detaching signed slopes from offsets loses weights; forgetting owner loses deletion semantics |
| LRC global exit | full phase-height carrier | `max_t min_v ||vt||` | top germs omit remote maxima |
| target map at an exceptional locus | fixed-bundle scalar / ambient-to-fibre quotient | closure versus kernel | conductor torsion can erase the fibre image |

The board updates were:

1. The BH `h_a` computation showed that exact reconstruction is feasible at
   the level of one tensor branch, but does not reduce the cubic witness.
2. The `V_bc` overwrite first looked like a missing symbolic identity.  The
   `a<->b` symmetry changed it into a naturality problem, while the exact point
   computation supplied a nontrivial arithmetic control.
3. THM-3990 turned “take the mean” into an actual cokernel statement.  The S6
   triangle hostile then forced the coefficient ring to remain part of the
   board.
4. AP13 forced replacement of a smooth Hessian by a labelled wedge.  The
   wedge retained a Schur-type elimination law but destroyed any automatic
   Poisson corrector.
5. THM-2050 then moved “global exit” from a warning to a mandatory state
   coordinate: every invariant of the recorded unperturbed germs can agree
   while the target differs.
6. The concurrently frozen S6 centralizer companion changed the continuous
   address from a raw-parameter guess into a marked local quotient: punctured
   shifts are integral, the order-four filling permits only even shifts, and
   `kappa=[c0-c2] in C/(2Z)` survives locally.  The cusp/overlap Cech class is
   still missing, so this is not a classification of completed threefolds.
7. THM-3995 changed the LRC lane by shrinking the support of any scale-two
   failure through four parity holes.  This is orthogonal to THM-3993's local
   top displacement: one controls the global sheet-count carrier, the other a
   labelled germ, and no lawful map between their response variables is yet
   known.

## 4. What the two Hopf problems actually are

Brendle--Hung addresses the **Hopf product problem**: whether `S2 x S2`
admits positive sectional curvature.  The paper claims yes by a perturbation
of a nonnegatively curved Cheeger--Mueter metric.

The `S6` manuscript addresses the classical **Hopf complex-structure
problem**: whether `S6` admits an integrable complex structure.  It claims a
compact complex threefold fibred by two-dimensional complex tori over `P1`
and a diffeomorphism of the total space with `S6`.

Neither target implies the other:

```text
S2 x S2 versus S6,
Riemannian sectional curvature versus integrable complex charts,
metric perturbation versus analytic gluing and topology recognition.       (3)
```

The Hopf sign conjecture is a third statement; the claimed positive metric
would be compatible with it because `chi(S2 x S2)=4`.  The Hopf
fibration/invariant-one/division-algebra story is a fourth historical cluster.
No target-preserving map from that cluster to either construction was found;
shared nomenclature is not a research bridge.

## 5. Exact Brendle--Hung audit advance

The v1 notebook's Lemma 5.4 cell assigns `ra=P1[L[hc]]`, so it checks the
`h_c` branch twice and omits `h_a`.  The independent
[SymPy audit](../04-computation/brendle_hung_lemma54_independent_audit_20260824.py)
rebuilds the moving-frame metric, connection, Riemann tensor, linearized
curvature, and `h_a`.  On both generic components of the parametrized torus it proves

```text
P1(L(h_a))=(0,0,0,0),                                    (4)
```

with four exact point controls and a hostile tensor giving
`(0,0,-1/3,2/9)`.  This is a **VERIFIED EXACT SYMBOLIC IDENTITY relative to
the reconstructed formulas**, not a finite census, and is restricted to this
one omitted `P1(L(h_a))=0` identity on `Sigma intersect M_generic`.  Deducing `z(h_a)=0`
additionally inherits the source's generic Hessian-invertibility assertion,
which this script does not audit.

The notebook also forms the intended `V_bc` and then overwrites it by `Vac`
before asserting zero.  At the exact generic quaternion frame `(1,2,3,4)`,
the [independent point audit](../04-computation/brendle_hung_vbc_exact_point_audit_20260824.py)
finds

```text
581/20250 - 1064/50625 - 259/33750 = 0.                 (5)
```

Because all three terms are nonzero, (5) is a strong point control, but it is
not an identity on `M_generic`.

The generic-locus repair comes from the background isometry

```text
F(p_1,p_2)=(R p_1,-R p_2),
R e_x=e_y,       R e_y=-e_x,       R e_z=e_z.            (6)
```

It sends

```text
h_a -> -h_b,        h_c -> h_c,        h_ac -> -h_bc,    (7)
```

so naturality gives `V_ac o F^(-1)=-V_bc`.  Conditional on the notebook's
correctly saved `V_ac=0` calculation, this repairs `V_bc=0` throughout
`M_generic`.  The highest-value remaining exact task is still an independent
reconstruction of the claimed cubic mean coefficient

```text
[lambda_c lambda_d^2] integral_Sigma V^(3)
  =pi^2/(18 sqrt(3)).                                    (8)
```

The paper's TeX notation, metric, index, cross-term, label, and cutoff defects
also require the repairs recorded in the source ledger.

## 6. Connection contract A -- defect modulo lawful repairs

### Source

Brendle--Hung's residual third coefficient `w=V^(3)|Sigma` and the conformal
correction law `w -> w-Delta_Sigma chi`.

### Target

The `S6` conductor/clutch complexes and, more distantly, an owner-labelled LRC
residual complex.

### Map

```text
defect space E, allowed correction A:C->E
       |-> obstruction class [w] in coker(A).             (9)
```

### Preserved predicate

Whether a lawful correction can remove the non-obstructed part of the defect.

### Destroyed information

The coefficient ring, geometry of the repair, interface conditions,
integrality, branch orientation, and whether the corrected object satisfies
the final target.

### Needed sidecar

Domain/codomain, coefficient ring, component/incidence labels, adjoint kernel,
and a realization theorem for the allowed corrections.

### Cheapest decisive tests

- disconnected smooth components with means `+1,-1`;
- the triangle Laplacian's real zero-average solve but integral `Z/3` class;
- the S6 reversed-branch `Z/2` hostile;
- a proposed LRC corrector applied to AP13's six isolated components.

The output of this contract is a grammar, not a functor between the two Hopf
problems.

## 7. Connection contract B -- fixed scalar closure versus changing sheaves

Brendle--Hung extends a scalar curvature inequality from the generic base
points to the diagonal and antidiagonal inside one fixed compact smooth
Grassmann bundle.  Continuity and compactness are the relevant closure
mechanisms.

The `S6` cusp is nonnormal.  At the special fibre, the differential sheaves,
normalization branches, conductor torsion, and gluing cokernel change.  The
ambient cotangent section at issue can remain generically nonzero while its
image in fibre differentials becomes conductor-supported torsion and dies in
the torsion-free quotient.  The lossy operation is this ambient-to-fibre
quotient/normalization map, not ordinary dense-open restriction of a locally
free ambient sheaf.

Typed closure test:

```text
source        ambient restricted cotangents on the special fibre
target        fibre differentials / torsion-free quotient / normalization
map           ambient restriction -> fibre quotient -> normalization
question      what are this map's kernel, cokernel, and supports?
sidecar       conductor incidence and specialization complex.              (10)
```

Density is decisive for Brendle--Hung because the target is a continuous
scalar on one fixed compact smooth bundle.  It does not repair a lossy
ambient-to-fibre quotient.  This distinction generated a new meta-pattern
rather than a false transfer to the S6 manuscript.

## 8. Connection contract C -- smooth and polyhedral elimination

Brendle--Hung eliminates a smooth normal plane coordinate.  If the normal
Hessian is invertible, the minimizer moves with the perturbation and produces
an effective Schur-complement response.

[THM-3993](../01-canon/theorems/THM-3993-labelled-two-wall-polyhedral-elimination-and-local-top-response.md)
gives the exact polyhedral analogue.  For

```text
ell_+=M+A(s)+p(s)h,       ell_-=M+B(s)-q(s)h,
p,q>0,                                                     (11)
```

with every other wall slack, the unique local top is

```text
h_*=(B-A)/(p+q),
M_loc=M+(qA+pB)/(p+q).                                    (12)
```

At the AP13 top `t=1/14`, `(p,q)=(1,13)`, so the weights are `13:1`.
The map preserves local extremum and the sign of the first effective
coefficient.  It loses smooth tubular geometry, a connected defect locus, a
lawful integer-row deformation, and all remote tops.

THM-2050 is decisive: AP13 and the loose lift have identical complete
unperturbed germs at all six period-14 top points, but global maxima `1/14`
and `1/12`.  Hence no invariant computed only from those germs can prove
LRC(14) without a global first-exit or gluing sidecar; the theorem does not
compare arbitrary deformation families.

The positive toy

```text
min(h+s^2 y^2+s^3,-h+s^2 y^2+s^3)=s^2y^2+s^3-|h|        (13)
```

has `inf_y max_h=s^3>0`.  Thus nonsmoothness is not the obstruction; lawful
global structure is.

## 9. Connection contract D -- symmetry is target-bearing data

Brendle--Hung's perturbation tensors depend on `e_x,e_y,e_z`.  This symmetry
breaking is necessary: a positively curved metric on `S2 x S2` cannot retain
a nontrivial Killing field.  Yet one discrete symmetry, (6), is exactly what
repairs the `V_bc` notebook omission.

This gives a two-sided lesson:

```text
averaging over all symmetry can destroy the desired object;
retaining a typed discrete symmetry can prove a missing identity.          (14)
```

For S6 the target-bearing symmetry data are branch orientation and clutch
gauge.  For LRC they are wall owner/sign and paired deletion.  Quotienting
before checking the consumer can either erase the theorem or invent a false
modulus.

## 10. Views generated from the board

| object | representation | invariant | operation | quotient/symmetry | scale |
|---|---|---|---|---|---|
| BH zero planes | compact smooth bundle | minimum curvature | metric jet + normal minimization | diagonal `SO(3)` broken to typed discrete symmetries | plane -> residual torus -> global bundle |
| BH residual | function on `Sigma` | component mean | Poisson correction | quotient by `im Delta` | cubic coefficient |
| S6 cusp | fan quotient + conductor complex | homology, differential extension | normalization/gluing/specialization | integral Smith quotient, branch gauge | node -> triple point -> fibre -> total space |
| LRC top | labelled wall envelope | local height response | phase elimination / paired deletion | owner/sign cannot be forgotten | germ -> all tops -> global first exit |

The table exposes the only plausible transfer axis: **operation followed by
its quotient**.  The objects, symmetries, and scales themselves do not match.

## 11. Frontier consequences

### Brendle--Hung

1. Independently recover coefficient (8) with a fresh immutable exact
   calculation and hostile perturbations.
2. Rebuild all ten quadratic summands with the corrected source indices and
   cutoff.
3. Produce an explicit parameter pair, spectral Poisson solution, admissible
   `s`, and interval-certified positive curvature lower bound.

### S6 manuscript

1. Rebuild the infinite-fan free/proper action, compactness, and oriented
   conductor quotient.
2. Audit the half-weight square root and `O(-1)` torsor.
3. Compare collapse and nearby-cycle topology through a labelled integral
   map, including attaching orientations and cup products.
4. Globalize or refute the conductor-induced ambient-to-fibre cotangent kernel
   class with relative duality and base change.
5. Compute the global overlap/cusp Cech class of the centralizer generator
   `b=2`; the FINITE-EXACT marked address `kappa=[c0-c2] in C/(2Z)` is only a
   local filling invariant.
6. Compute the `c_0` Kodaira--Spencer class and the unresolved `h^(1,2)`.

### LRC(14)

1. Seek a **lawful** deformation or corrector on THM-2047's selected complex;
   do not import a graph Laplacian without an operation realizing its edges.
2. Couple THM-3993's local response to THM-3910/3995's integer `t`-sheet,
   parity-hole, and variance carrier, retaining all competing tops and first
   global exit.
3. Attack every proposal first with AP13/loose-lift germ equality and isolated
   equality atoms.
4. A success must close a named one of THM-3910's 17 conditional pair types
   or the `t<U` branch; a new local statistic alone is not progress on the
   headline.

## 12. Stopping reasons for tempting transfers

1. **Direct theorem transfer:** different manifolds and predicates; stop.
2. **“Both have tori”:** residual torus, torus fibre, and phase circle have
   different roles; stop unless an operation intertwiner is constructed.
3. **Dense-generic proof for S6:** the ambient-to-torsion-free-fibre map has a
   conductor kernel; scalar density does not repair that quotient.
4. **Real ranks for integral gluing:** Smith torsion is target-bearing; stop.
5. **Local LRC jets:** THM-2050 gives an exact collision; retain first exit.
6. **Hopf invariant/fibration import:** no map to curvature positivity or
   complex integrability was found; stop at historical nomenclature.

## 13. Final truth boundary

### PROVED

- THM-3990's componentwise smooth/graph repair quotient;
- THM-3991's factorial Euler law in its periodic unimodular grammar;
- THM-3993's labelled two-wall local response and AP13 globalization hostile;
- THM-3995's scale-two parity-hole tariff in its conditional LRC slice;
- THM-3955/3957's local conductor differential mechanisms.

### VERIFIED -- exact symbolic identity

- the omitted Brendle--Hung `P1(L(h_a))=0` identity on both components of
  `Sigma intersect M_generic`, relative to the reconstructed formulas.

### FINITE-EXACT / targeted controls

- one nontrivial exact generic-point cancellation for `V_bc`;
- the S6 displayed-matrix, oriented-chain, centralizer, and marked filling
  audits in their stated conditional universes.

### RELATIVE REPAIR

- `V_bc=0` throughout `M_generic` follows by the typed isometry (6)--(7) from
  the notebook's correctly saved `V_ac=0` calculation.

### OPEN / UNDER AUDIT

- Brendle--Hung's cubic coefficient, full notebook replay, and headline
  positive metric;
- the S6 manuscript's analytic construction, topology realization,
  integrable complex structure, and CDP conflict;
- every global map from the repair-quotient grammar to LRC;
- LRC(14).

The session's main mathematical gain is a three-way separation:

```text
local elimination can be exact,
repair quotients can be universal in form,
global realization remains problem-specific.             (15)
```
