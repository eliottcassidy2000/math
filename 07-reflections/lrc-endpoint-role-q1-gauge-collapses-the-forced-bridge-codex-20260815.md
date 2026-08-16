# Coarse target landing collapses the forced `H--q5` bridge

**Status: research synthesis around a FINITE-EXACT unnumbered sidecar; not a
truth source and not an LRC(14) result.  The audited statements incorporated
below now live in promoted THM-3479.**  Reproduce the exact finite claims with
`04-computation/lrc_endpoint_role_q1_gauge_quotient_obstruction_20260815.py`.

## 1. Inheritance and the decisive type split

The closest proved mechanism is
[THM-2334](../01-canon/theorems/THM-2334-relation-residue-current-and-character-twist-pushforward.md):
relation-residue pushforward is exactly a coordinate-twist Fourier transform.
The relevant hostile is MISTAKE-261: coarse target landing need not retain a
refined graph observable.  The corrected near miss is MISTAKE-409: an exact
graph coboundary can still have a nonzero weighted tree determinant, so zero
absolute graph `H^1` is not the present obstruction.  The least-used sidecar
is the q1-gauged role chart in
`lrc-relation-role-chart-is-b1-not-h1-and-the-d5-type-boundary-codex-20260815.md`.

Four objects must be kept separate.  Fix one of the two tuples `w=U_full` or
`w=U_clock` in current order

```text
(H,q1,q2,q3,q4,q5,c1,c2,c3).
```

1. The vector

   ```text
   a=(20110798,-41,-27,-27,-27,38,-27,-27,-27)
   ```

   is an exact **relation address** in
   `Lambda(w)=ker(w:Z^9 -> Z)`.  The scalar `C(a;X,m)` is the grouped
   exact-address orbit coefficient of THM-2334.  It is not computed by the
   present package.
2. With `K_w=ker(w mod 13)` and the six-row owner packet `L_w subset K_w`,
   the character-side endpoint current is

   ```text
   gamma_w = H_w : Ghat_w=L_w^perp/<w> -> K_cyc,       (1)
   ```

   where `K_cyc` denotes the relevant cyclotomic coefficient field (or a
   certified finite-field image used only to prove exact nonvanishing).
3. The address-side target response is the normalized inverse transform

   ```text
   A_w(q)=1/169 sum_(ell in Ghat_w)
                zeta_13^(-<ell,q>) gamma_w(ell),
   q in G_w=K_w/L_w.                                  (2)
   ```

   It is an unrestricted aggregate over every exact relation address landing
   at `q`.  It is neither `C(a;X,m)` nor the all-`91`-unit projector `B(q)`.
4. A **physical current** would be an actual Boolean/common-ancestry endpoint
   observable, or a lawful realization of a THM-2512 response contraction.
   Neither (1), (2), nor a graph gradient supplies that realization.

In particular, `gamma_w(e_i)` is generally ill-typed: a coordinate basis
vector need not lie in `L_w^perp/<w>`.  Coordinate labels live naturally on
the address side of (2), after a lawful gauge into `K_w`.

## 2. The cheapest lawful coordinate response

Since `w_q1` is a unit modulo thirteen, every retained coordinate `i` has a
canonical q1-gauged residue

```text
r_i=e_i-(w_i/w_q1)e_q1 in K_w.                       (3)
```

The cheapest type-correct eight-potential candidate is therefore

```text
P_i=A_w([r_i] in G_w),
i in {c1,c2,c3,H,q2,q3,q4,q5}.                       (4)
```

This q1 gauge must not be conflated with the rational deletion map in the
relation-role sidecar.  Equation (3) gauges an **address residue** using the
scalar word `w`; the sidecar deletes q1 from the **speed tuple** using the
fixed address `a`.  They have different domains.

Use the exact dual basis already present in the THM-3479 companion,

```text
v1=-e_q1+e_c2,
v2=-e_q2+e_c3.                                      (5)
```

Pairing (3) with (5) gives the same profile for both tuples:

| role label | class in `G_w ~= F_13^2` |
|---|---:|
| `c1` | `(0,0)` |
| `c2` | `(1,0)` |
| `c3` | `(0,1)` |
| `H`  | `(1,0)` |
| `q2` | `(1,12)` |
| `q3` | `(1,0)` |
| `q4` | `(1,0)` |
| `q5` | `(1,0)` |

The fixed relation address itself lands at `(1,0)`.  More importantly,

```text
[r_H]=[r_q5].                                       (6)
```

The mechanism is visible before any Fourier computation.  The guard row of
the owner packet is a nonzero scalar multiple of

```text
e_q5-e_H,
```

so the quotient by `L_w` deliberately identifies the two labels.  Equation
(6) therefore holds for every function on the coarse target quotient, not
only for the endpoint response (2).

## 3. Why every coarse role determinant is zero

The role contract sends `H` to the unique degree-seven hub and `q5` to the
unique repair leaf.  Their edge is the unique bridge joining the leaf to the
two-`K4` carrier.  From (4) and (6),

```text
P_H-P_q5=0.                                          (7)
```

Every spanning tree contains the unique bridge.  Hence its bridge weight is a
factor of the signed weighted reduced-Laplacian determinant, and (7) forces

```text
det L_red(dP)=0                                      (8)
```

for all 72 role charts, for both tuples, and for every possible endpoint
function `A_w`.  This conclusion does not depend on whether the provisional
U_full bank has `169/169` nonzero values.  Nonzero vertex values can coincide,
and their difference on a forced edge is still zero.

Thus the canonical coarse construction (4) does factor through eight labelled
slots, but it is spectrally dead.  Existing data supply no alternate
coordinate-resolved endpoint response that avoids the quotient identification.

## 4. The dimension-minimal repair is one guard character

Let `L_w^-` be the span of the five owner-packet rows obtained by deleting the
named guard row.  Then

```text
dim L_w^-=5,
dim(K_w/L_w^-)=3.                                    (9)
```

Among the six named one-row deletions, guard deletion is the only one for
which `r_H-r_q5` leaves the retained row span.  This is a statement about that
named deletion menu, not uniqueness among all hyperplanes of `L_w`.

The additional dual character can be chosen as simply as

```text
v3=e_H.                                             (10)
```

It annihilates all five retained rows and separates `r_H` from `r_q5`.  The
refined coordinate profile in the dual basis `(v1,v2,v3)` is

| role label | class in `K_w/L_w^- ~= F_13^3` |
|---|---:|
| `c1` | `(0,0,0)` |
| `c2` | `(1,0,0)` |
| `c3` | `(0,1,0)` |
| `H`  | `(1,0,1)` |
| `q2` | `(1,12,0)` |
| `q3` | `(1,0,0)` |
| `q4` | `(1,0,0)` |
| `q5` | `(1,0,0)` |

The exact companion assigns synthetic values `2^j` to the five distinct
classes as a positive structural control.  All 72 refined determinants are
nonzero; their minimum absolute value is `64,000`, with 18 distinct values and
sign census `(32,40)`.  Therefore the three-dimensional quotient removes the
**forced** zero.  It says nothing about the actual endpoint values.

The incoming
`lrc-nine-relation-slots-as-k4-edges-xor-axes-and-star-triangle-selector-codex-20260815.md`
gives this lost coordinate another interpretation.  Its four-vertex chart
places the six `-27` slots on `K4` edges and `{H,q1,q5}` on the three matching
axes.  The coarse owner packet kills the contrast between two of those axes,
`H-q5`; the character `e_H` restores exactly one matching-axis contrast.  It
does not restore the sidecar's three V4 edge characters or its Archimedean
star/triangle parity.  Indeed U_full and U_clock have identical coarse and
refined class **geometry** here, while that sidecar distinguishes them as
triangle and star.  Any endpoint separation must therefore enter through the
actual response values (and ultimately physical realization), not through the
mod-thirteen class labels alone.

## 5. The minimal implementable endpoint experiment

Define the refined character bank

```text
gamma_w^-(alpha,beta,tau)
 =H_w(alpha v1+beta v2+tau e_H),
(alpha,beta,tau) in F_13^3,                          (11)
```

and let `A_w^-` be its `13^3` inverse DFT.  THM-2334's general residue-current
transform makes (11) lawful.  The existing U_full interval engine already
evaluates `H_w(ell)` for a supplied coordinate twist; it must be extended from
`169` to `2,197` twists.

The smallest decisive computation is:

1. evaluate the `2,197` values in (11) under the existing certified embedding;
2. inverse-transform only the five distinct refined role classes in the table;
3. form their 13 edge differences;
4. test the bridge and both signed `K4` tree factors in all 72 charts.

The `tau=0` character slice must reproduce the old coarse transform, and
summing `A_w^-` over the new fibre must recover `A_w`.  The stored zero/v2
endpoint values, DFT reconstruction, a flat-response zero control, and
ordinary/optimized byte equality give cheap hostile and positive controls.

Even a nonzero refined determinant would remain an endpoint-aggregate result,
not yet a physical current.  A second sidecar would still have to prove that
the 13 edge differences are realized on one lawful Boolean/common-ancestry
row (and, for the original all-unit goal, survive the coupled mod-91
projector).  No grouped `C(a;X,m)` nonvanishing, all-unit `B(q)`, ancestry,
bispectrum, scalar-row exclusion, or LRC(14) conclusion follows here.

## 6. Connection contract

| field | exact content |
|---|---|
| source | THM-2334 coarse target response `A_w:G_w->K_cyc` |
| attempted map | q1-gauge each role coordinate by (3), evaluate (4), take graph differences |
| target | signed edge weights and reduced-Laplacian determinant on the two-`K4` carrier |
| preserved | target class, exact Fourier typing, role labels, graph automorphism orbit |
| destroyed | the guard-versus-q5 packet coordinate, exact address, all-unit mask, clock, ancestry, physical event meaning |
| hostile | the forced bridge has identical endpoint classes and zero weight |
| minimal restoration | refine `L_w` to `L_w^-` and retain the `e_H` character |
| next test | the `2,197`-twist experiment (11), followed by a separate physical-realization audit |

The useful lesson is a quotient rule: if the graph carrier has a forced edge,
check whether the source quotient identifies its endpoint labels before doing
any expensive spectral calculation.  Here that one check both kills the
coarse bridge and identifies the cheapest lawful refinement.

## 7. Exact companion

Run from the repository root:

```bash
python -B 04-computation/lrc_endpoint_role_q1_gauge_quotient_obstruction_20260815.py
python -B -O 04-computation/lrc_endpoint_role_q1_gauge_quotient_obstruction_20260815.py
```

Both LF-normalized transcripts equal
`05-knowledge/results/lrc_endpoint_role_q1_gauge_quotient_obstruction_20260815.out`
byte for byte.  The immutable LF hashes are

```text
script: 02d6cf3553edf3412da6d8eb99c7937f841d9441cc9599bc01319f72596b5887
output: 14e9483fd1ec42cc0cdffb156dc520a491bc95d0816a85828732c965858ec5d7
semantic: fe9fe2d4b1d98aba3b60e1a7f6823bcccabed04d7a319f3eb20c83d25c799422
```

## 8. 2026-08-16 update: the actual U_full guard refinement is nonzero

**Status: FINITE-EXACT unrestricted mod-thirteen endpoint-aggregate theorem
component, independently audited and incorporated in promoted THM-3479.**
It does not prove a physical current, a grouped exact-address coefficient, an
all-`91`-unit aggregate, a U_clock statement, a scalar-row exclusion, or
LRC(14).

The minimal experiment proposed in Section 5 has now been completed for
`U_full`.  Delete the named `H` guard row from the six-row algebraic packet,
retain the five other rows, and use the lawful dual basis below.  This does
not delete the Boolean guard: `PATTERN_E` and `PATTERN_QA` retain
`guard_safe`.

```text
(v1,v2,e_H).
```

The exact universe is the complete bank

```text
gamma^-(alpha,beta,tau)
 =H(alpha v1+beta v2+tau e_H),
(alpha,beta,tau) in F_13^3,                         (12)
```

of `2,197` twists.  The pinned THM-3479 periodic-subtraction endpoint engine
processed `71,070,080` exact interval components.  Its finite-field target is

```text
p=572252886246508880869,
h=279936,
ord_p(h)=81750412320929840124.                       (13)
```

The companion checks a complete Lucas certificate for the primality of `p`
and checks that `h` has the exact order in (13).  Therefore a nonzero image is
an exact nonvanishing certificate for the corresponding cyclotomic endpoint
quantity; the printed residues are not being identified with the original
cyclotomic numbers.

### 8.1 Recovery controls

The `tau=0` slice reproduces the earlier coarse `169`-twist bank exactly:

```text
coarse gamma digest:
  afcdb043eb1bf8095c313473a3d3bdcf4ce027f86b01b5f3cecbc7c87e6484b3,
coarse unnormalized target digest:
  726423df1b9e1c93b356966e5c3c386669e2f6b19da0bf8818606204eb2e9ee5.
                                                               (14)
```

The new full gamma digest is

```text
1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682. (15)
```

The inverse transform uses the normalized factors `13^-3` and `13^-2` in
the refined and coarse quotients respectively.  For each of the four distinct
coarse role bases

```text
(0,0), (0,1), (1,0), (1,12),                         (16)
```

the direct sum of all thirteen refined fibre values equals the corresponding
coarse response.  This is checked from the computed bank, rather than merely
invoked as a formal orthogonality identity.  A flat role response is the
hostile control: its bridge, both K4 tree factors, and their product vanish in
all `72` charts.

### 8.2 Actual role response and mechanism

At the five distinct refined role classes, the normalized inverse-transform
images are

```text
(0,0,0)  -> 405336876493642499425,
(0,1,0)  -> 518539850465495448196,
(1,0,0)  -> 503604956476841920373,
(1,0,1)  -> 320618948602619577408,
(1,12,0) ->  15703541686881447885.                   (17)
```

The entries in (17) are pairwise distinct in the certified image.  In
particular, the restored guard coordinate separates `H` from `q5`:

```text
P_H-P_q5=389266878372286537904 mod p !=0.             (18)
```

The labels `c2,q3,q4,q5` still share class `(1,0,0)`, so the refinement does
not restore every coordinate contrast.  Nevertheless those remaining edge
zeros do not kill either complete-`K4` tree polynomial.  Across the full
declared labelled chart orbit the exact zero census is

```text
(bridge, first K4, second K4, product)=(0,0,0,0)      (19)
```

where each entry counts zero values among `72` charts.  Each K4 factor takes
nine nonzero finite-field values, and the bridge value (18) is constant over
the orbit.  The complete chart-record digest is

```text
b7d8c2c9860e4f1aa542b1c85fdb7b65cf4985aba5a81a84ff3a324834d51c51. (20)
```

Thus the minimal repair has a precise mechanism: the extra `e_H` character
restores exactly the forced bridge contrast, while the complete K4 tree sums
tolerate the remaining within-wing class collision.  This proves nonzero
**endpoint-aggregate graph factors** for U_full under the refined quotient.
It does not realize those factors as one Boolean/common-ancestry physical
current.

### 8.3 Independent graph audit and scope ledger

A cheap independent audit imports the older frozen labelled role-chart engine,
feeds it only the five frozen residues (17), and reconstructs all `72` factor
rows.  It independently obtains (18)--(20) and the flat hostile.  The full
endpoint semantic digest is

```text
f06a258b85212ba981de4c18027147a8875c67fc35d1f55d9128c17d746994a9. (21)
```

The scope ledger is:

| object | status after this experiment |
|---|---|
| U_full refined `13^3` endpoint bank | FINITE-EXACT under one certified primitive embedding |
| `tau=0` coarse recovery | FINITE-EXACT, exact stored digests |
| refined-to-coarse fibre sums | FINITE-EXACT at all four role bases |
| bridge and two K4 factors | FINITE-EXACT nonzero in all `72` charts |
| flat-response hostile | FINITE-EXACT zero in all `72` charts |
| U_clock refined bank | OPEN / NOT COMPUTED |
| lawful Boolean/common-ancestry realization | OPEN |
| grouped `C(a;X,m)` or all-unit `B(q)` | OPEN |
| scalar-row exclusion or LRC(14) | OPEN |

U_clock was not attempted: its recorded exact-boundary interval counts are
materially larger than U_full's, and U_full already gives the decisive test
of whether the minimal quotient repair can work.  No inference about U_clock
is made.

Reproduce from the repository root with

```bash
PYTHONHASHSEED=0 python3 -B 04-computation/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.py
python3 -B 04-computation/lrc14_guard_deleted_refined_endpoint_role_graph_audit_20260816.py
```

The first command took about `2,078` seconds with four workers on the session
host.  The checked artifacts have LF hashes

```text
primary script: ee2105742abee578a9c41ff7ec954a07ada324fccc2c643429e7ac6e6e6f8fc2
primary output: 10a98351cc59615a5b6d2b8f555e0936d1a39566d9906127edc2b0fbc3918e73
graph script:   b75f5ef933c18c07b4fd2c4812fea468a9faaa6c5847a5c6a11190cc04676261
graph output:   d4f5e4b12854bd5802842df2784d580abdd75d9a4fc58dd77db3a185a10403f1
```
