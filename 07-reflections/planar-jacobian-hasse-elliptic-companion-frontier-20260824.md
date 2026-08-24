# Planar Jacobian continuation: Hasse layers, elliptic reduction, and the companion ledger

**Research reflection, 2026-08-24.** This is a synthesis and task generator,
not an additional theorem. Canonical claims are only those cited below with
their statuses. The planar Jacobian conjecture remains **OPEN**.

## Session outcome

The inherited reduced `(2,3)` cusp-pole cell from
THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual had one
large unknown:

```text
G=C^2-A^3+(3a^2/4)A+a^3/4
 =gamma*u+(3a/(2gamma))*p+R(p,y),       R in (p^2,y).
```

Three distinct operations produced three proved advances.

1. **Change rows, not depth.** Passing from Laurent powers of `tau` to
   source-normal powers of `t` exposes sharp degree caps. The first two
   normal diagonals force `[p^2]R` and `[p^3]R`; `[y]R` and then `[py]R`
   are the surviving scalars. The zero-residual cell is empty. This is
   THM-3997, **PROVED + VERIFIED-EXACT + independently audited**.
2. **Freeze support, not coefficient degree.** The source grading closes the
   first cell wider than THM-3973's two-by-two theorem: three forced weights
   in `A`, two weights plus at most one arbitrary extra weight in `C`, with
   coefficient degrees in `u` completely unbounded. This is THM-3998,
   **PROVED + VERIFIED-EXACT + independently audited**.
3. **Divide before completing.** The quotient `Q=G/t` extends to the smooth
   completion. Its boundary order, endpoint scheme, and total divisor class
   are exact, while factor ownership is deliberately left unresolved. This
   is THM-3999, **PROVED + VERIFIED-EXACT + independently audited**.

The gains are compatible but logically independent. THM-3998 applies for
every residual `R`; THM-3997 says `R` cannot vanish; THM-3999 explains which
part of a nonzero `R` is visible at the completion boundary.

## Inheritance pass and corrected near misses

The closest proved mechanism was THM-3992's centered nodal residual. The
canonical hostile was THM-3980's disconnected formal algebraization: local
or formal solvability does not imply an element of `B_2`. The corrected near
miss was THM-3979: all finite formal Darboux jets lift, so asking for another
unstructured jet cannot be the obstruction. The least-used relevant sidecar
was the distinction between a normalization address, its global component
owner, and a completion-boundary endpoint.

That inheritance dictated the successful moves:

```text
unstructured formal t-jets
  --add Laurent-derived degree caps--> structured source diagonals;

nodal clutch roots
  --add complete-owner sidecar--> THM-3996 address balance;

regular equation Q with ord_D(Q)>0
  --remove fixed boundary component--> strict endpoint polynomial.
```

Two tempting implications failed and were repaired.

- `ord_D(Q)=2` does **not** imply that the strict companion meets `D`.
  The `R=0` control has `div(Q_0)=C_0+2D` with `C_0` disjoint from `D`.
- A two-edge graph critical group `Z/2` is **not** THM-3994's local `A_1`
  class group `Z/2`. The former can occur inside smooth `A2`; the latter is
  attached to one singular local ring. There is no canonical map between
  them.

## Live concept board

### 1. Hasse or `p`-Taylor repair layers

For a finite nonnegative Laurent tail

```text
Delta=sum_n d_n(s)tau^n,                 tau=p-s^2,
H_m=sum_(n>=m) binom(n,m)d_n(s)(-s^2)^(n-m),
```

THM-3997 proves

```text
Delta in k[p,y]  iff  deg_s H_m<=m for every m.
```

This is a lossless ordinalization of the repair data. The natural number
`m` is the layer; the actual monomials in that layer are secondary. In the
same spirit as replacing the `m`-th odd square by the index `m`, the useful
coordinate is membership in the selected layer, not the ambient numerical
size of its representatives. Unlike a heuristic ranking, this one is an
exact direct-sum decomposition.

### 2. Source-grading support

The grading `wt(x)=1, wt(t)=-2` retains both boundary colors:

```text
(B_2)_(-h)=x^-h u^ceil(h/2)(1+u)^ceil(h/3) k[u].
```

THM-3998 shows why this representation sees something the Laurent depth
does not. One homogeneous weight can contain arbitrarily many Laurent rows
as `deg_u` grows. Thus a finite number of weights is not a finite-jet or
finite-degree statement.

### 3. Elliptic generic fibres

Under the hostile specialization `R=0`, the source and target generic fibres
over `K=k(q)` become explicit elliptic curves. The Keller field inclusion
would make them isogenous. At `q=infinity`, the target has potential good
reduction and the source has negative `j`-valuation. This gives a second,
global proof that `R=0` is impossible.

The important methodological point is not merely “compute `j`.” Unequal
`j`-invariants are compatible with isogeny. The preserved predicate is
potential good reduction, and the proof must first establish the exact
source function field, the finite target inclusion, and the projective
genus-one morphism.

### 4. Companion divisor and endpoint polynomial

THM-3999 converts the mixed residual into

```text
R(p,y) |-> Q=G/t |-> (2, gamma-R(0,y), -2[D]).
```

The endpoint polynomial sees exactly the `p=0` shadow `R(0,y)`. It forgets
all terms divisible by `p`, all factor ownership, and all node addresses not
lying on the selected source line.

### 5. Node-address balance

THM-3996 is the correct graph sidecar. Vertices are normalization-cover
owners, edges are complete source addresses over the target node, and the
orientation is plus branch to minus branch. Finite-etale degree gives exact
in/out balance. A graph constructed only from the two visible roots on
`t=0` is incomplete and supports no cycle obstruction.

The incoming strengthening of THM-3996 adds a decisive global invoice. A
proper fibre has cardinality equal to the Keller field degree; degrees one
and two are excluded here, so it has at least three addresses. Consequently
the visible two-edge cycle can be a connected packet but cannot be the whole
proper fibre. Either another packet exists or the forced node is a Jelonek
value. This does not locate the extra address, but it makes the census task
nonvacuous.

### 6. Tournament inheritance

The useful tournament-era idea is not “force a tournament.” It is the older
discipline of identifying the intrinsic binary observable and retaining its
gauge. Here an address has an intrinsic ordered pair of branch owners only
after the target normalization labels plus and minus branches. Reversing
that label reverses every edge. Ties are not present, but missing addresses
are missing data, not ties.

This explains why the tournament ancestry mattered:

```text
vertices: normalization-cover owners;
observable: which owner contains each lifted branch;
orientation gauge: o_+ to o_-;
preserved target: degree balance and directed cycles;
lost data: incomplete node fibre and nonproper escaping addresses;
sidecar: THM-3996 properness/completeness audit.
```

The graph becomes legitimate exactly when those fields are declared. Before
then, sections, endpoints, residues, or proof obligations are better vertices
than source points.

### 7. Repair quotients across problems

The Hasse transform, source-grading coefficient ideals, and boundary quotient
`Q=G/t` are three versions of the same research move:

1. identify the apparent defect;
2. divide by the forced carrier;
3. characterize the exact image after division;
4. retain the missing coordinate as a sidecar.

For Hasse layers the carrier is `p-s^2`; for negative weights it is
`u^ceil(h/2)(1+u)^ceil(h/3)`; for the companion it is `t` followed by the
fixed boundary factor `x^2`. This is evidence for the existing repair-
quotient meta-pattern, not a new cross-problem theorem.

## Pairwise effects on the board

- **Hasse layers -> grading support:** `[y]R` identifies the first scalar a
  wider support cell must realize; THM-3998 shows the smallest three-weight
  `A` support still cannot realize any residual.
- **Hasse layers -> endpoint polynomial:** the first free scalar is the
  linear coefficient of `R(0,y)`, hence the first infinitesimal motion of a
  strict endpoint toward the boundary.
- **Elliptic fibres -> Hasse layers:** two unrelated invariants exclude the
  same `R=0` cell. This makes the local resonance calculation much less
  likely to be a gauge artifact.
- **Endpoint polynomial -> address graph:** boundary endpoints of companion
  closures are not target-node addresses. They live on different divisors
  and must not be merged into one vertex set.
- **Grading support -> companion factorization:** a factor of `Q` with narrow
  source support would induce constrained support in `A,C`, suggesting a
  factor-first hostile search.
- **Tournament discipline -> class ledger:** graph homology and divisor class
  groups may share presentations, but their generators and quotient maps must
  be typed before torsion is compared.

## Generated next tasks

### Anchor: continue the reduced `(2,3)` residual

1. **All-diagonal recurrence.** Derive the general source-normal elimination
   of the `(m+1)`-st jet and express the pure coefficient `[p^(m+1)]R` in
   terms of the earlier mixed coefficients. Cheapest test: compute rows four
   and five symbolically at generic `a,gamma,beta,delta`, then search for a
   triangular recurrence rather than another isolated formula.
2. **Pure-`p` residual lane.** Set `R=R(p)` but do not bound its degree. The
   source generic fibre remains a double cover with an explicit branch
   polynomial. Compute genus and reduction at `q=infinity`; compare with the
   target elliptic fibre. This is the closest extension of the successful
   elliptic proof.
3. **`R(0,y)=0` lane.** Here THM-3999 forces the strict companion to avoid
   `D`, while THM-3997 forces `R!=0`, hence `R in (p^2,py)` is a genuinely
   nonzero interior correction. Ask whether the class `-2[D]`, disjointness,
   and etale pullback force reducibility or a missing address.
4. **Endpoint multiplicity lane.** If `gamma-R(0,y)` has a repeated root,
   compute the first transverse coefficient of `Q/x^2`. Decide whether the
   repeated endpoint produces tangency, two components, or a smooth higher
   contact. THM-3994 warns that scalar multiplicity alone cannot decide.

### Niche: widen support by the smallest honest amount

5. **One extra `A` weight.** THM-3998 leaves this as the first unclosed
   grading cell. Enumerate collision weights against `C:{3,1}` before solving
   coefficients. The hostile probe should retain unbounded `u`-degree and use
   endpoint valuations to kill extreme rows.
6. **Two extra `C` weights.** Build the collision graph among bracket weights.
   Search first for cycles that can alter all four core rows; isolated
   collisions cannot repair the degree contradiction.
7. **Factor-aware support.** Assume `Q=Q_1Q_2` and distribute total class
   `-2[D]` and endpoint degrees among factors. Check whether either factor is
   forced into one of THM-3973/3998's closed support cells.

### Wildcard: reduction and monodromy

8. **Residual Newton polygon at infinity.** For sparse families such as
   `R=beta*y+alpha_2*p^2`, compute the generic source-fibre model and its
   stable-reduction polygon. A mismatch of potential reduction type is useful;
   a mere mismatch of `j` is not.
9. **Hasse associated graded.** Treat the degree condition `deg H_m<=m` as a
   filtered image problem. Compute the first obstruction module after fixing
   the forced layers. This may turn the infinite diagonal ladder into a finite
   generation question.
10. **Complete-address computation.** For bounded residual families, solve
    `A(x,t)=-a/2, C(x,t)=0` without restricting to `t=0`. Record every point,
    branch orientation, owner, and escape flag. This is the cheapest direct
    test of THM-3996's extra-address/Jelonek dichotomy.

## Stopping boundaries

The session does not establish any of the following:

- a contradiction for every nonzero coordinated residual;
- extension of the free scalars through all Laurent rows;
- irreducibility of the general companion `Q`;
- completeness or properness of the two visible node addresses;
- an obstruction from the two-edge cycle by itself;
- `JC(2)` or `DC(2)`.

The new frontier is nevertheless smaller and better typed. The formerly
featureless ideal `R in (p^2,y)` now has three independent coordinate systems:
Hasse layers, source-grading support, and the boundary endpoint polynomial.
The next breakthrough should make at least two of those systems constrain the
same nonzero residual, rather than pushing any one representation in isolation.

## Addendum: the first two invariant rows close the live seven-piece seam

[THM-4005](../01-canon/theorems/THM-4005-reduced-two-three-live-seam-invariant-support-atlas.md)
performs the promised comparison between Hasse layers and source grading. On
the continued seam `gamma=-a^3/2`, divide by the residual characters and use
the quotient coordinates

```text
A5=a^5,               b=[y](R/gamma),
d=[py](R/gamma).
```

The first two source-normal rows then force both outputs to carry at least
four retained weights. This closes `3x4/4x3` inside the oriented live `2:3`
seam, and a no-import `111`-gate audit shows that the exclusion survives the
specific THM-3992 linear shear. The stronger fixed-gauge `4x5` invoice does
not survive that transfer: an expanded normalized `C` may contain a complete
copy of `A`, which the inverse shear can cancel. That hostile is the reason
the theorem records the seven-piece exclusion but not a pre-gauge `4x5`
floor.

The concept board changes as follows.

- **Hasse layers x grading support:** `b` and `d` are no longer merely free
  coefficients; their vanishing pattern is exactly the support-stratum
  selector `5x6,5x6,4x5`.
- **Endpoint x support:** `b` is simultaneously the first endpoint derivative
  and the first coefficient whose nonvanishing widens both support invoices.
- **Clutch symmetry x next row:** on `b=d=0`, the known two-clutch packet is
  exchanged by `x->-x` through `t^2`; the missing
  `c21=[p^2y](R/gamma)` is the first source-odd sidecar that can break it.
- **Gauge x ownership:** a support claim can transfer through the recorded
  linear shear without making the finite Weierstrass factors global owners.
  The support and owner ledgers remain distinct.

The cheapest continuation is now sharply localized. Compute the next
source-normal/Hasse row and its residual triple

```text
(c40,c21,c02)=([p^4],[p^2y],[y^2])(R/gamma).
```

Here `c02` moves the endpoint after `b`, while `c21` breaks the odd symmetry
on the only fixed-gauge minimal-support stratum. A contradiction must use one
of those native effects or a later row; merely refining the already complete
`t^0,t^1,t^2` companion packet cannot close the cell.

## Addendum: the 6--11 weighted-face and companion campaign

**Research synthesis, 2026-08-24; JC(2) remains OPEN.** The requested
`6,7,8 -> 9,10,11` comparison became structural rather than numerological.
Four independently typed advances now replace the earlier task list.

1. [THM-4007](../01-canon/theorems/THM-4007-live-two-three-third-normal-row-five-weight-floor.md),
   **PROVED + VERIFIED-EXACT**, computes the third source-normal row. On
   `b=d=0`, it forces a nonzero `p^3` term and the affine lock

   ```text
   [p^4]Rtilde+(6/(7A5))[y^2]Rtilde=-11392/(105A5^4).
   ```

   The fixed-gauge support floor is now `5x5`; a determinant shortcut that
   appeared to give a contradiction was refuted by a direct fourth-jet
   hostile before promotion.
2. [THM-4008](../01-canon/theorems/THM-4008-pure-p-residual-totally-degenerate-generic-fibre-no-go.md),
   **PROVED + independently geometric-audited**, excludes the entire lane
   `R=f(p)`. Its genus-`m` source fibre degenerates to a rational curve with
   `m` nodes, while the target retains good `j=0` reduction; degree
   conservation forbids the map. This is an all-coefficient statement, not a
   truncation claim.
3. [THM-4011](../01-canon/theorems/THM-4011-companion-observer-kernel-etale-log-rh-and-endpoint-gate.md),
   **PROVED + independently audited**, identifies a large ambient kernel of
   the endpoint/class/clutch/finite-row observer:

   ```text
   G -> G(1+T(p^3-y^2)),          T in p*k[p,y].
   ```

   The prime insertions `1+p^M(p^3-y^2)` have growing genus but class zero
   and can start beyond any finite row packet. Actual etaleness recovers a
   log-Riemann--Hurwitz gate: odd `M` is impossible, while even `M` costs
   degree at least `M+3` and two distinct finite normalization-side exits.
   The boundary map, not the scalar endpoint polynomial, restores the target
   node label.
4. [THM-4012](../01-canon/theorems/THM-4012-weighted-leading-face-good-elliptic-factor-observer.md),
   **PROVED + independently geometric-audited**, supplies the weighted
   stable-reduction observer. Under its explicit face-stable hypothesis, a
   highest face must have a positive-genus Jacobian factor mapping to
   `E_0:Y^2=X^3+1`. Singleton faces fail. The exact forced max-weight-six
   model is proved face-stable without that hypothesis and is excluded by a
   two-prime nontorsion certificate. Thus on `b=d=0` the actual total weight
   satisfies the new unconditional floor `M>=7`; the stronger `M>=9` remains
   conditional on face-stability at weights seven and eight.

A parallel arithmetic control must not be conflated with item 4.
[THM-4016](../01-canon/theorems/THM-4016-sharp-five-by-five-elliptic-attachment-nontorsion.md)
proves nontorsion for the distinct sharp-support point
`(X^3,Y^2)=(43/84,127/84)`, where `p^4!=0`; its stable-specialization
application remains conditional. THM-4012's max-six point is
`(43/224,267/224)`, where `p^4=0`, and its boundary model is proved.

### Why six through eleven form one object

Put `U=SP`. A top monomial `p^i y^j` becomes `P^iU^j`, with weights
`wt(P)=2, wt(U)=3`. The cusp core has weight six:

```text
P^3+U^2.
```

The non-singleton supports from weights six through eleven are its smallest
semigroup translates:

| weight | normalized two-term face | genus | proved/candidate obstruction |
|---:|---|---:|---|
| 6 | `P^3+U^2=1` | 1 | `j=0` survives Hom, but the forced attachment is nontorsion |
| 7 | only singleton `P^2U=1` | 0 | Hom obstruction under face-stability |
| 8 | `P(P^3+U^2)=1` | 2 | Bolza; `Jac ~ E_8000^2`, no Hom to `E_0` |
| 9 | `U(P^3+U^2)=1` | 3 | primitive `Q(zeta_9)` CM candidate |
| 10 | `P^2(P^3+U^2)=1` | 2 | primitive `Q(zeta_5)` CM candidate |
| 11 | `PU(P^3+U^2)=1` | 5 | primitive `Q(zeta_11)` CM candidate |

Weight seven is the semigroup hole: weight one has no monomial translate of
the weight-six core, so only the isolated `P^2U` face occurs. Weight twelve
is the first support with three monomials,

```text
P^6,                    P^3U^2,                    U^4,
```

and is therefore the first face with a genuine coefficient modulus rather
than one coefficient ratio removable by diagonal scaling. This is the new
natural stopping boundary.

### The exact max-six torsion separator

The third-row lock forces, after residual normalization,

```text
epsilon=2752/(135A5^3),       kappa=-5696/(45A5^3).
```

An attachment on the normalized target `E_0` satisfies

```text
X^3=43/224,                    Y^2=267/224.
```

At primes above `11` and `17`, it reduces respectively to `(2,3)` of exact
order six and `(7,2)` of exact order nine. If its characteristic-zero order
were `N`, good reduction would force

```text
N=6*11^a=9*17^b,
```

which is impossible already 2-adically. The remaining geometric invoice is
not assumed: the weighted boundary chart has a persistent `A_5` point whose
normalization `r=z^3w` has unit quadratic discriminant, while the six finite
attachments are exact `A_35` models `UV=unit*rho^36`. Consequently the
special fibre has one elliptic component and otherwise only rational
components. Equality of the six attachment images forces torsion, producing
the contradiction.

### Provisional cyclotomic sieve at 9, 10, and 11

**VERIFIED-EXACT / PROVISIONAL, not yet a canon theorem.** The independent
THM-4012 certificate checks the Newton interiors, diagonal automorphisms, and
CM-type stabilizers. For a toric interior point `(i,j)`, the differential

```text
omega_ij=P^(i-1)U^(j-1)dP/F_U
```

has diagonal eigencharacter `lambda^i mu^j`. The three holomorphic spectra
are

```text
M=9:  {1,2,5} mod 9,
M=10: {1,2} mod 5,
M=11: {4,5,8,9,10} mod 11.
```

Together with their conjugates they exhaust the roots of
`Phi_9,Phi_5,Phi_11`. Each type has trivial multiplicative stabilizer, hence
is primitive; the standard primitive-CM classification then makes the
Jacobians simple of dimensions `3,2,5`. The imported theorem and its exact
scope are pinned in
[CORE-PAPERS](../05-knowledge/reference/CORE-PAPERS.md#florit--smith--milne--low-weight-jacobian-factors-and-primitive-cm-types).
Therefore none can have an elliptic quotient, in particular not `E_0`.

This would upgrade the conditional face-stable alternative from `M>=9` to
`M>=12`. It is not yet unconditional: the weighted boundary inventory for
weights `7,...,11` must be proved rather than inferred from the affine face.

### Revised Anchor / Niche / Wildcard board

- **Anchor — toroidal boundary inventory for weights 7--11.** For each
  singleton or cyclotomic face, normalize every point of `Z=0`, prove the
  closure flat, and list every component created by base change. The cheapest
  hostile is a coefficient-sum resonance that moves all finite attachments
  to the boundary. Success promotes the conditional `M>=12` floor.
- **Niche — actual companion factors.** THM-4011 shows every finite observer
  packet has an ambient class-zero kernel. The next useful invariant must see
  Darboux realizability: restrict a proposed factor to the target nodal
  normalization, compute its punctures, and apply log Riemann--Hurwitz before
  attempting factorization in `B_2`.
- **Wildcard — weight twelve as a moduli/Prym problem.** Classify
  `c_0P^6+c_1P^3U^2+c_2U^4=1` by its cross-ratio or discriminant. Determine
  exactly when its Jacobian contains `E_0`; then attach the node-orbit
  compatibility, since a Hom survivor is only necessary and never a lift.

The concept board now has a clear causal chain:

```text
Hasse row lock
 -> exact top-face coefficients
 -> proper weighted degeneration
 -> good elliptic-factor invoice
 -> attachment orbit
 -> arithmetic torsion separator.
```

The destroyed information at each arrow is also explicit: later residual
terms, boundary tails, the dual graph, attachment images, and the isogeny
kernel. The next session should restore exactly one of those sidecars rather
than compute another untyped row.
