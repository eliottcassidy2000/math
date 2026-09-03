---
id: THM-4384
title: "LRC14 small-defect short-relation master formula and sharp sector atlas"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4373 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED; LRC(14) OPEN. For every primitive distinct positive odd three-unit
  signed relation m b=s h p+t q with m>=h+1, m+h<=13, and 3 not dividing
  hm, labeled owners force defect zero and the exact comb has one determinant
  coordinate. Parity leaves exactly ten coefficient patterns, with all sharp
  maxima classified. The sectors overlap and cannot be summed. No universal
  nonresonant, seam-entry, ledger, all-tail, or LRC(14) consequence is asserted.
source: root + small_defect_relations + clean-room referee / next-sharp session, 2026-09-03
depends_on:
  - THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound
related:
  - THM-4382-lrc14-signed-one-four-one-comb-exact-measure-and-sharp-maximum
  - THM-4383-lrc14-signed-two-five-one-comb-exact-measure-and-sharp-maximum
primary_script: 04-computation/lrc14_small_defect_short_relation_master_atlas_thm4384.py
primary_output: 05-knowledge/results/lrc14_small_defect_short_relation_master_atlas_thm4384.out
primary_script_sha256: 743f49a9f1f21cfabcdeb288d89547000ebe4cb1aa1e969d620bdd9199fcad60
primary_output_sha256: 352923f43e52218906bde8930fd2103947f1b6e38b814bd89f707338f4d1837f
boundary_script: 04-computation/lrc14_small_defect_short_relation_boundary_hostile_thm4384.py
boundary_output: 05-knowledge/results/lrc14_small_defect_short_relation_boundary_hostile_thm4384.out
boundary_script_sha256: 287a966917deef07212f6605518955fe4d46fa852e290dcf7d13737d5c0a207e
boundary_output_sha256: 8ab0876f2e95afec113b548756e474a9e722c44a4b108af77542ee00cd9d0d50
independent_referee_script: 04-computation/lrc14_small_defect_short_relation_independent_referee_thm4384.py
independent_referee_output: 05-knowledge/results/lrc14_small_defect_short_relation_independent_referee_thm4384.out
independent_referee_script_sha256: b7da6c447e47e61e060817d203297ff604e5e9d6faeddbad7b12e576481be29a
independent_referee_output_sha256: f1bf7d03b8bf7ee7f1dc60232c8607e3dabc33430ef23f7d844846f8c5c919a0
hash_basis: raw LF bytes
audit: >
  PASS. The clean-room verifier reconstructs the defect-zero proof, cyclic
  quotient, exact quadrature, ten parity patterns, all 747 exceptional
  presentations, the complete sharp atlas, duplicate-presentation and
  cross-sector incidence, and the first defect-three boundary. It compares
  every exceptional presentation with a literal six-sheet physical-circle
  union and uses explicit checks that remain live under optimized Python.
  Normal, optimized, and hash-seeded replays match the frozen stream. The
  primary's optimized replay is determinism evidence only because it uses
  Python assert.
---

# THM-4384 -- Small-defect signed short relations: one determinant and a ten-sector atlas

**PROVED ELEMENTARY RELATIVE TO THM-4373 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. ALL TEN PARITY-COMPATIBLE SECTORS IN THE DECLARED SMALL-DEFECT RANGE
HAVE SHARP MAXIMA. THE SECTORS OVERLAP. THE UNIVERSAL NONRESONANT COMB BOUND,
SEAM ENTRY, ALL-TAIL TRANSFER, AND LRC(14) REMAIN OPEN.**

## Outcome

Let `h,m` be positive integers satisfying

```text
m>=h+1,          m+h<=13,          gcd(3,hm)=1.          (1)
```

Let `p,b,q` be primitive, distinct, positive odd three-units, and suppose
that for signs `s,t in {+1,-1}`,

```text
m b=s h p+t q>0.                                        (2)
```

For the physical scale-three failure comb of THM-4373, put

```text
A=3|q-p|/(14m),        B=3(p+q)/(14m),
E(x)=sum_(k>=1,3 does not divide k)(x-k)_+-x^2/3.       (3)
```

Then the exact master formula is

```text
mu(F_(p,b,q))
 =2 sum_(k>=1,3 does not divide k)
      m/(pq)[(B-k)_+-(A-k)_+]
 =6/(49m)+2m/(pq)[E(B)-E(A)].                           (4)
```

The formula is independent of `h,s,t` after the relation presentation has
selected `p,b,q`; those data are nevertheless essential in the proof.
Since `-1/3<=E<=0`,

```text
mu(F_(p,b,q))<=6/(49m)+2m/(3pq).                        (5)
```

Parity of odd speeds leaves exactly ten nonempty coefficient patterns in
the range (1):

```text
(1,2,1), (1,4,1), (1,8,1), (1,10,1),
(2,5,1), (2,7,1), (2,11,1),
(4,5,1), (4,7,1), (5,8,1).                             (6)
```

The exact sharp atlas is:

| pattern `(h,m,1)` | sharp maximum | primitive maximizing speed set(s) | analytic cutoff `pq>=N` | finite presentations / triples |
|---|---:|---|---:|---:|
| `(1,2,1)` | `6/77` | `{1,5,11}` | 80 | 26 / 13 |
| `(1,4,1)` | `12/301` | `{1,11,43}` | 289 | 62 / 31 |
| `(1,8,1)` | `6/287` | `{1,5,41}` | 953 | 116 / 58 |
| `(1,10,1)` | `6/343` | `{1,5,49}` | 1271 | 130 / 65 |
| `(2,5,1)` | `12/371` | `{1,11,53}` | 425 | 75 / 75 |
| `(2,7,1)` | `12/469` | `{1,19,67}` | 577 | 75 / 74 |
| `(2,11,1)` | `6/371` | `{1,5,53}` | 1455 | 132 / 131 |
| `(4,5,1)` | `82/2597` | `{5,7,53}`, `{7,41,53}` | 471 | 84 / 84 |
| `(4,7,1)` | `12/497` | `{5,13,71}` | 702 | 97 / 97 |
| `(5,8,1)` | `58/2765` | `{5,13,79}` | 941 | 117 / 116 |

The first two rows recover THM-4373 and THM-4382. The fifth recovers THM-4383.
The other seven are new sharp sector closures.

## Inheritance pass and concept board

The closest proved mechanism is THM-4373's signed-`(1,2,1)` determinant
quadrature. THM-4382 shows that the same mechanism survives a larger bounded
middle defect once labeled owners are retained. The canonical hostile is the
general rank-two nearest-integer lattice of an arbitrary triple. The
corrected near miss is to infer defect zero from interval eligibility:
THM-4383's `y=1/53` control has all three tails eligible but defect `-1`
and colliding owners. The least-used sidecar is the signed integer defect
before owner scalarization.

The live concepts are: the bounded integer defect; the owner sum in
`F_3`; the cyclic determinant quotient; the period-three quadrature; the
coefficient-pattern incidence hypergraph; and the first defect-`3`
boundary.

## 1. Labeled owners annihilate the bounded defect

Put `r=3/14` and `y=3x mod 1`. For an eligible speed `w`, write

```text
n_w=nint(wy),       e_w=wy-n_w,       |e_w|<r,
o_w=-w^(-1)n_w mod 3.                                    (7)
```

Physical triple failure is exactly eligibility of all three speeds and
`{o_p,o_b,o_q}=F_3`. Define

```text
delta=m n_b-s h n_p-t n_q.                               (8)
```

Equation (2) gives

```text
delta=s h e_p+t e_q-m e_b,
|delta|<(m+h+1)r<=3.                                    (9)
```

The last inequality is strict even when `m+h=13`; hence `|delta|<3`.

Reduce (2) modulo three. The three quantities
`A_0=s h p`, `C_0=t q`, and `m b=A_0+C_0` are all nonzero in
`F_3`. Two nonzero ternary residues have nonzero sum only when they are
equal, so

```text
A_0=C_0.                                                (10)
```

Using `n_w=-w o_w`, equations (8)--(10) give

```text
delta=A_0(o_p+o_b+o_q) mod 3.                            (11)
```

Three distinct ternary owners have sum zero. Thus `3|delta`, and (9)
forces the exact integer identity

```text
delta=0.                                                (12)
```

Conversely, endpoint eligibility for `p,q` plus (12) gives

```text
e_b=(s h e_p+t e_q)/m,
|e_b|<(h+1)r/m<=r,                                      (13)
```

so the middle speed is automatically eligible. Conditions (1) have two
separate roles: `m+h<=13` kills the integer defect, while
`m>=h+1` makes the middle interval redundant after the defect is killed.

## 2. The complete lift lattice is cyclic

Primitivity and (2) imply

```text
gcd(p,b)=gcd(p,b,q)=1.                                  (14)
```

Indeed, a common divisor of `p,b` also divides
`q=t(mb-shp)`. On the defect-zero plane define

```text
k=b n_p-p n_b.                                         (15)
```

After recovering

```text
n_q=t(m n_b-s h n_p),
```

Bezout for `p,b` shows that (15) is onto `Z`, and its kernel is exactly
the physical translation vector `Z(p,b,q)`. Therefore

```text
{delta=0 integer lifts}/Z(p,b,q) isomorphic to Z.        (16)
```

This remains true when `gcd(p,q)>1`; endpoint coprimality was never used.
The endpoint determinant is

```text
q n_p-p n_q=t m k.                                     (17)
```

Moreover `o_b-o_p=k/(pb) mod 3`. Under (12), equation (11) says the
owner sum is zero, so the owners are all distinct exactly when

```text
3 does not divide k.                                   (18)
```

Thus one natural-number address `|k|`, with its sign and ternary residue
retained, is complete for this consumer. Encoding it by an ordinal loses
nothing; dropping the sign or the owner residue would.

## 3. Exact component length and quadrature

For fixed `k`, (17) separates the endpoint interval centers by
`m|k|/(pq)`; their radii are `r/p,r/q`. The exact intersection length is

```text
f(|k|)=m/(pq)[(B-|k|)_+-(A-|k|)_+].                    (19)
```

Equation (16) supplies exactly one component on the `y` circle for each
`k`. Nearest-integer uniqueness makes different classes disjoint.
Removing `3|k` and pairing signs proves the series in (4).

For `u>=0`, direct summation over residues one and two modulo three gives

```text
E(u+3)=E(u),                    -1/3<=E(u)<=0.           (20)
```

Also

```text
B^2-A^2=9pq/(49m^2).
```

Substitution into (19) gives the second equality in (4), and (20) gives
(5).

## 4. Why exactly ten patterns remain

For odd `p,b,q`, equation (2) requires

```text
m=h+1 mod 2.                                           (21)
```

Enumerating the positive integer pairs satisfying (1) and (21) gives exactly
(6). This is a coefficient enumeration, not a bounded speed census.

For each row, once the displayed maximum `M` is found, the
smallest integer `N` with

```text
6/(49m)+2m/(3N)<M                                      (22)
```

is the cutoff in the table. Above it, (5) is strictly smaller than `M`.
Below it, exact enumeration ranges over every ordered odd three-unit
`(p,q)`, every sign pair, and

```text
b=(s h p+t q)/m>0,
```

then applies distinctness and `gcd(p,b,q)=1`. The table records the complete
finite universes. Formula (4) and a definition-level physical comb
implementation agree on every retained presentation.

The `h=1` rows have two ordered presentations because their endpoint
coefficients agree. Three other rows contain one multiply-presented speed
set. The theorem is presentation-wise, and the exact formula agrees on those
fibres.

## 5. The defect-width boundary is sharp

Odd-speed parity makes `m+h` odd, so the first case outside (1) has
`m+h=15`. At `(h,m)=(2,13)`,

```text
13*1=-2*5+23.
```

For `(p,b,q)=(5,1,23)` and `y=34/161`, the nearest integers, errors,
owners, and defect are

```text
(n_p,n_b,n_q)=(1,0,5),
(e_p,e_b,e_q)=(9/161,34/161,-1/7),
(o_p,o_b,o_q)=(1,0,2),
delta=-3.                                               (23)
```

All errors are strictly below `3/14` and all owners are distinct, but the
defect is not zero. Thus the implication (11) only gives `3|delta`, and the
range `m+h<=13` is sharp for this zero-defect mechanism. Exact controls at
`(4,11)`, `(5,10)`, and `(7,8)` likewise realize defect `+3`.

This is the correct next boundary: larger patterns require an extra defect
state `delta/3`, not another unmodified one-determinant proof.

## 6. Overlap and scope firewall

After quotienting the automatic endpoint swap in the `h=1` rows, the exact
coefficient-vector audit finds intrinsic multiple presentations precisely on

```text
(2,7,1):  {1,5,17},
(2,11,1): {5,7,41},
(5,8,1):  {11,13,23}.
```

The ten relation sectors also overlap. The complete cross-sector incidence
list is

```text
(1,4)-(5,8):   {1,7,11}
(1,8)-(2,5):   {1,5,13}
(1,8)-(4,7):   {1,5,13}, {1,11,19}
(1,10)-(2,5):  {1,7,17}
(1,10)-(4,7):  {1,13,23}
(2,5)-(4,7):   {1,5,13}, {5,7,11}
(2,7)-(2,11):  {1,7,25}
(2,7)-(5,8):   {1,5,17}, {7,11,19}
(2,11)-(5,8):  {5,7,31}, {7,19,29}.
```

For example,

```text
{1,5,13} lies in (1,8,1), (2,5,1), and (4,7,1);
{1,7,11} lies in (1,4,1) and (5,8,1).
```

Therefore the atlas is not a disjoint partition, and its sector measures
cannot be added. A relation presentation is a certificate for formula (4),
not an intrinsic unique type unless a separate coefficient-incidence audit
says so.

Every row controls only a physical three-tail scale-three failure comb.
Nothing here forces an arbitrary nonresonant triple to have one of these
relations, supplies a lower bound for a body-safe set, synchronizes several
tail triples on one phase, retains the body-component address after
scalarization, or proves seam entry. The universal nonresonant bound,
lower-shell composition, and LRC(14) remain **OPEN**.

## Reproduction

The primary performs the exact ten-row enumeration and compares every finite
presentation with the direct physical comb inherited from the audited
THM-4382 implementation. The boundary sidecar independently scans exact owner
cells for the first four `m+h=15` coefficient pairs and finds the nonzero
defects in (23). The clean-room referee rederives the master mechanism and
checks all 747 normalized exceptional presentations against its own literal
six-sheet physical-circle union; its 30 targeted owner/defect controls contain
4,952 exact cells.

```powershell
python -B 04-computation/lrc14_small_defect_short_relation_master_atlas_thm4384.py
python -B 04-computation/lrc14_small_defect_short_relation_boundary_hostile_thm4384.py
python -B 04-computation/lrc14_small_defect_short_relation_independent_referee_thm4384.py
python -B -O 04-computation/lrc14_small_defect_short_relation_independent_referee_thm4384.py
```

The primary's bare assertions are disabled by `-O`, so only its normal replay
is a live verification. The referee uses explicit runtime checks and its
normal, optimized, and hash-seeded streams agree. Raw-LF hashes are

```text
primary script:    743f49a9f1f21cfabcdeb288d89547000ebe4cb1aa1e969d620bdd9199fcad60
primary output:    352923f43e52218906bde8930fd2103947f1b6e38b814bd89f707338f4d1837f
boundary script:   287a966917deef07212f6605518955fe4d46fa852e290dcf7d13737d5c0a207e
boundary output:   8ab0876f2e95afec113b548756e474a9e722c44a4b108af77542ee00cd9d0d50
referee script:    b7da6c447e47e61e060817d203297ff604e5e9d6faeddbad7b12e576481be29a
referee output:    f1bf7d03b8bf7ee7f1dc60232c8607e3dabc33430ef23f7d844846f8c5c919a0
referee semantic:  283276a308e10e8af7656b64541576cce8263087c989c3c06b18cdf7ec0b880e
```

This theorem controls relation presentations, not their union measure. It
does not classify arbitrary nonresonant triples or prove a body-safe-set,
synchronization, seam-entry, ledger, all-tail, or LRC(14) result.

QED.
