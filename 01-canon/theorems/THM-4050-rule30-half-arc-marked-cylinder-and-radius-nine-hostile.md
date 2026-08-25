---
id: THM-4050
title: "Rule 30 half-arc rebase, marked-cylinder law, and radius-nine hostile"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED, with a FINITE-EXACT
  single-seed census. The physical terminal-current arc at depth k shortens
  universally to floor(k/2)+1 transitions, with a sharp even midpoint. Its
  nearest-zero stopping address is exactly a marked spatial cylinder
  1 followed by 2r-2 zeros. Under Haar-random initial data this gives the
  complete geometric stopping law, including the finite-line infinity
  class and exact conditional center probabilities. For the single seed,
  the census through k=262143 has maximum address nine and the exact witness
  k=79883 refutes a universal address-eight bound. Haar-to-temporal transfer,
  boundedness of the addresses, and all three Rule 30 prizes remain open.
source: root + rule30_marked_address + audit_row_atlas, 2026-08-24
audit: >
  PASS after adding the initially omitted finite-line infinity class. The
  referee independently rederived every index and orientation, the half-arc
  and cylinder proofs, the complete Haar partition, and the radius-nine
  witness; a separate standard-coordinate evolution reproduced the full
  census and center split. The companion itself compares packed and sparse
  complete Rule 30 rows through 512, exhausts characteristic dependency words
  through radius eight, and checks the single seed through 262143. Normal and
  optimized streams byte-match the frozen output. Universal claims come from
  the proof, not extrapolation from those ranges.
depends_on:
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
  - THM-3481-rule30-cyclic-arc-norm-rank-and-marked-innovation-spectrum
related:
  - THM-3468-rule30-radial-green-fold-innovation-discrepancy-and-fixed-seed-carrier-boundaries
  - THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge
  - THM-4047-rule30-left-front-affine-monodromy-clock
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
script: 04-computation/rule30_marked_half_arc_zero_radius_probe_20260824.py
output: 05-knowledge/results/rule30_marked_half_arc_zero_radius_probe_20260824.out
report: 07-reflections/rule30-marked-half-arc-and-characteristic-cylinder-codex-20260824.md
script_sha256: 429f36d4a6593b3cdae39896d7b332a548cb011d63c1320b272d303d80726d21
output_sha256: b617cd7b18439eee40f075d3ee939b7bd8dfdf4ecf0bcd0dc88bd5dd1b144aaf
report_sha256: 4124c46ebbf58fa12fcc8a6ead073d16d8095432de8e7606648509f05791700f
hash_basis: raw LF bytes
---

# THM-4050 -- the marked half-arc and its stopping cylinder

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED; FINITE-EXACT where
stated.** The light-cone boundary shortens THM-3471's calibrated terminal arc
by almost half. Stopping at its nearest physical zero turns that variable arc
into one explicit characteristic cylinder. The corresponding Haar law is
exact, but it lives in a different probability space from the deterministic
single-seed temporal orbit.

## 1. Terminal coordinates and the universal half-arc

Let `a_t(j)` be Rule 30 from one seed, over `F_2`:

```text
a_0(0)=1,  a_0(j)=0 for j!=0,
a_(t+1)(j)=a_t(j-1)+(a_t(j) or a_t(j+1)).             (1)
```

For `k>=2` and `-k<=h<=0`, retain THM-3471's physical terminal line and
transition current

```text
T_k(h)=a_(k+h)(h),
Q_k(h)=T_k(h+1)+T_k(h),
c_k=T_k(0)=a_k(0).                                   (2)
```

Put

```text
m_k=floor(k/2)+1.                                    (3)
```

At `h=-m_k`, the cell has time `k-m_k` and spatial distance `m_k` from the
seed. Since `2m_k>k`, it lies strictly outside the light cone. Therefore

```text
T_k(-m_k)=0,
c_k=xor_(h=-m_k)^(-1) Q_k(h).                        (4)
```

The second identity is telescoping. It replaces the old universal length-`k`
arc by exactly `floor(k/2)+1` current values.

The boundary is sharp. If `k=2r`, then `(t,j)=(r,-r)` is the extreme-left
ray, which is identically one because its unique live predecessor
neighborhood is `001`. Hence

```text
T_(2r)(-r)=1,
c_(2r)=1+xor_(h=-r)^(-1)Q_(2r)(h).                   (5)
```

Thus the uniform zero endpoint cannot be moved one step inward at even
depths. As an operator corollary of THM-3481, the complete cyclic
`m_k`-arc profile is a lossless recoding of the phase current exactly when
`m_k` is odd, equivalently when `k=0` or `1 mod 4`. This operator statement
does not reconstruct the one marked scalar from unpointed data.

## 2. Nearest-zero address and exact marked cylinder

The zero in `(4)` makes the stopping address

```text
z_k=min{a in {1,...,m_k}:T_k(-a)=0}                  (6)
```

well-defined for every single-seed depth `k>=2`. Telescoping from its actual
zero endpoint gives the adaptive identity

```text
c_k=xor_(h=-z_k)^(-1)Q_k(h).                         (7)
```

For every `1<=r<=m_k`, the following are equivalent:

```text
z_k>r;
T_k(-1)=T_k(-2)=...=T_k(-r)=1;
T_k(-1)=1 and Q_k(-2)=...=Q_k(-r)=0;
(a_(k-r)(-r),...,a_(k-r)(r-2))=(1,0,...,0).          (8)
```

The current list is empty when `r=1`. The first three clauses follow from
the definition and `(2)`. For the last, use the following elementary
right-characteristic lemma. For any binary configuration `x`,

```text
(F^j x)(j)=1 for 0<=j<r
iff (x(0),...,x(2r-2))=(1,0,...,0).                  (9)
```

For `r=1`, this just says `x(0)=1`. Inductively, after the forced prefix
`1 0^(2j-2)`, the next characteristic value equals one exactly when its two
new right parents vanish. The extreme-ray dependence makes those conditions
exactly `x(2j-1)=x(2j)=0`. This adds two forced zeros per step and proves
both directions of `(9)`. Translating its base point to `(k-r,-r)` proves
the last clause in `(8)`.

At the final transition there is also the exact two-coordinate readout

```text
c_k=Q_k(-1)+1_[z_k>1].                               (10)
```

The stopping address is therefore not a fitted statistic. It is the
physical location of a marked `1 0^(2r-2)` cylinder, and together with the
retained current arc it exactly recovers the center.

## 3. Complete Haar stopping law

Now start Rule 30 from a Bernoulli-`1/2` iid row instead of the single seed.
For fixed `k`, define the random finite-line address

```text
Z_k=min{a in {1,...,k}:T_k(-a)=0},
Z_k=infinity if no such a exists.                    (11)
```

Left permutivity makes Rule 30 preserve Bernoulli Haar measure: on every
finite output cylinder, free right-boundary data and successive solution for
the left input give the same number of preimages to every output word. By
`(8)--(9)`, for `1<=r<=k`, the event `Z_k>r` is one specified word among the
`2^(2r-1)` words in its dependency cylinder. Thus

```text
P(Z_k>r)=2^(-(2r-1)),                                (12)
P(Z_k=1)=1/2,
P(Z_k=r)=3*2^(-(2r-1))                 (2<=r<=k),
P(Z_k=infinity)=2^(-(2k-1)).                         (13)
```

Appending the center asks for one more right-characteristic one, so

```text
P(c_k=1,Z_k>r)=2^(-(2r+1)).                          (14)
```

Taking successive differences, including the terminal class, gives the
complete joint and conditional laws

```text
P(c_k=1,Z_k=1)=3/8,
P(c_k=1 | Z_k=1)=3/4,

P(c_k=1,Z_k=r)=3*2^(-(2r+1)),
P(c_k=1 | Z_k=r)=1/4                    (2<=r<=k),

P(c_k=1,Z_k=infinity)=2^(-(2k+1)),
P(c_k=1 | Z_k=infinity)=1/4.                        (15)
```

The infinity class is load-bearing: without it the displayed center-one
masses miss `2^(-(2k+1))`. With it they sum to `1/2`, as Haar preservation
requires. MISTAKE-492 records this finite-stopping repair.

## 4. The exact finite single-seed hostile

The companion performs the complete single-seed census

```text
2<=k<=2^18-1=262143.                                 (16)
```

It finds

```text
z_k     1       2      3     4     5    6   7  8  9
count 131120  98428  24449  6074  1570 383  90 20  8. (17)
```

No row is unresolved under the declared search cap `64`. The exact first
radius-nine witness is

```text
(k,z_k,c_k)=(79883,9,0),                             (18)
T_k(-9),...,T_k(0)=(0,1,1,1,1,1,1,1,1,0).
```

Its `z_k>8` event is the literal marked cylinder `1 0^14` on the earlier
row. Thus `(18)` **REFUTES** the proposed universal bound `z_k<=8`; this
refutation needs only the one verified exact witness. The assertion that
the maximum is nine is `FINITE-EXACT` only on `(16)`. No boundedness or
unboundedness theorem for the full address sequence follows.

## 5. Exact transfer obstruction and next target

The probability spaces in Sections 3 and 4 are different:

```text
Haar law:       random initial row, fixed marked cylinder;
prize target:   one deterministic seed, temporal average at a moving mark. (19)
```

The finite temporal histogram in `(17)` is close to, but exactly unequal to,
the Haar proportions. Neither Haar invariance, odd support, phase balance,
full Walsh support, nor the fixed-depth clock of THM-4047 transfers a random
spatial cylinder law to the distinguished orbit.

For fixed `r`, the cylinder theorem identifies the lawful discrepancy target

```text
D_r(N)=#{k:r<=k<=N and
             (a_(k-r)(-r),...,a_(k-r)(r-2))=1 0^(2r-2)}
       -(N-r+1)/2^(2r-1).                            (20)
```

Bounding `D_r(N)` would control the marginal stopping cylinder. Center
balance additionally requires the signed joint discrepancy with `Q_k(-1)`,
because `(10)` proves that the marginal law of `z_k` alone loses one bit.
The nested Haar cylinders have conditional factor `1/4`, but no deterministic
martingale or cylinder-discrepancy estimate is proved.

The preservation ledger is:

| map | preserves | destroys | needed sidecar / boundary |
|---|---|---|---|
| full terminal arc -> half arc | marked center | older transitions | seed light-cone zero |
| half arc -> nearest-zero arc | marked center | current before the stop | address `z_k` |
| terminal run -> spatial cylinder | exact event `z_k>r` | other row bits | physical characteristic orientation |
| Haar cylinder -> scalar mass | exact random-row probability | distinguished seed and time correlations | temporal discrepancy theorem |
| address -> center | terminal run | final transition bit | `Q_k(-1)` via `(10)` |

All three Rule 30 prizes remain **OPEN**. This theorem proves no temporal
equidistribution, eventual periodicity or nonperiodicity, center balance,
address boundedness, or computational lower bound.

## 6. Exact verification

Run

```text
python3 -B 04-computation/rule30_marked_half_arc_zero_radius_probe_20260824.py
python3 -B -O 04-computation/rule30_marked_half_arc_zero_radius_probe_20260824.py
```

Both modes reproduce the frozen output byte-for-byte. The companion compares
two complete Rule 30 implementations through time `512`, checks `131,327`
transition-current identities, `511` full and shortened arcs, `256` sharp
even midpoints, and `66,047` cylinder equivalences. It exhausts `43,690`
dependency words through radius eight and performs the census `(16)` with
zero unresolved rows. The independent referee separately rederived the
universal proof and reproduced the census in standard coordinates.

Universal statements are the proofs in Sections 1--3. The finite checks are
implementation audits or carry the explicit `FINITE-EXACT` universe `(16)`.

**QED.**
