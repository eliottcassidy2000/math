---
id: THM-4029
title: "LRC(14) AP-cover twelve-owner rational tail"
status: >
  PROVED + FINITE-EXACT BASE/SELECTORS + VERIFIED-EXACT + INDEPENDENTLY
  REFEREED. From m=13 onward, the seven-sector AP noncover set is, modulo
  finitely many walls, exactly twelve rational-owner components. Their
  max-min arrival radii give an exact 60-phase rational deficit law with
  limit m(1-a(m))=127/35. The sequence is strictly increasing, converges to
  one, and is neither eventually quasipolynomial nor eventually C-finite in
  characteristic zero. The formula fails at m=12. This does not prove AP
  extremality for sparse sets or LRC(14).
source: root + sturmian_frontier + sturmian_referee / sequence continuation, 2026-08-24
audit: >
  PASS. Direct x-walls and an independent Farey/Sturmian cover-time engine
  agree through m=30. The m=13 component-owner bijection and every individual
  radius are asserted. All 60 printed phase formulas are proved selector-stable
  by exact linear inequalities and checked at two phase values. The owner
  formula agrees with the Farey engine through m=200. Normal and optimized
  theorem-audit streams are byte-identical.
depends_on:
  - THM-536-lrc-seven-sector-sturmian-partial-sum-reframe
related:
  - THM-535
  - HYP-2608
component_script: 04-computation/lrc14_ap_cover_components_thm4029.py
sequence_script: 04-computation/lrc14_ap_cover_sequence_engine_thm4029.py
formula_script: 04-computation/lrc14_ap_cover_finite_owner_formula_thm4029.py
audit_script: 04-computation/lrc14_ap_cover_twelve_owner_rational_tail_thm4029.py
audit_output: 05-knowledge/results/lrc14_ap_cover_twelve_owner_rational_tail_thm4029.out
component_script_sha256: 883f98f81d24c19c57532b80d2d3cfa97cc6edc572b619ba4a8c9026eaf81f83
sequence_script_sha256: f34b0879262f8cce8e6e55a89ded44597e76e3f125a773d75934278fc01abcd4
formula_script_sha256: d23909f3c5aee6ecc89976de6e944167cd07f88a69bb94856bc512c0970c071e
audit_script_sha256: 8ae58a5e0e4fa0ad72a4aae3642678706ed65e9b5b10ef834f970eff4265adbb
audit_output_sha256: 7800161e43615b5199bc3571d8737896f49f913b6fc760007d07b9f57db88745
hash_basis: raw LF bytes
---

# THM-4029 -- the AP-cover twelve-owner rational tail

**PROVED + FINITE-EXACT BASE/SELECTORS + VERIFIED-EXACT + INDEPENDENTLY
REFEREED.** This theorem turns the AP-cover scalar sequence into a finite
labelled geometric object. The owner, side, arrival track, phase, and affine
denominator are load-bearing sidecars. LRC(14) remains open.

## 1. Object and inheritance

For `m>=1`, put

```text
a(m)=meas{x in [0,1):
  {floor(7 e x) mod 7:0<=e<m}=Z/7},
D(m)=1-a(m).                                           (1)
```

After `theta=7x`, the residues are the partial sums of a mechanical word,
as in THM-536. The closest proved mechanism is THM-536's pointwise inclusion
and AP threshold table; its corrected hostile is the finite-table sentinel
recorded in MISTAKE-489. The least-used sidecar is the rational slope that
owns each connected noncover component.

The cover sets increase with `m`. They increase strictly from the first
possible cover:

```text
a(1)=...=a(6)=0,
a(m)>a(m-1) for every m>=7.                            (2)
```

For `m=7`, every `theta in (1,7/6)` gives the residues `0,...,6`, while
six times cannot cover seven residues. For `m>=8`, every

```text
theta in I_m=(6/(m-1),6/(m-2))
```

first covers at time `m`: before `e=m-1` the nondecreasing walk stays below
six, and at `e=m-1` it reaches six. Hence

```text
a(m)-a(m-1)>=6/[7(m-1)(m-2)]>0.                       (3)
```

## 2. The twelve persistent owners

The infinite rotation fails to cover all seven sectors, equivalently misses
at least one sector, exactly at

```text
B={p/q mod 1:gcd(p,q)=1,1<=q<=6}
 ={0,1/2,1/3,2/3,1/4,3/4,
   1/5,2/5,3/5,4/5,1/6,5/6}.                          (4)
```

An irrational orbit is dense. A reduced rational orbit of denominator
`q>=7` has mesh at most `1/7` and covers every half-open sector; at `q=7`
it has one point on each sector grid position, and one extra iterate gives a
cover stable under perturbation on either side. A denominator below seven has
at most six orbit points and cannot cover. Thus `(4)` is exact and has twelve
members. Compactness away from this finite set also proves `a(m)->1`.

## 3. Exact max-min radii

Let `m>=13`, `n=m-1`, and fix `p/q in B`, with `0<=p<q`. For `0<=s<q`, define

```text
A_s=7ps/q,
b_s=floor(A_s) mod 7,
f_s=frac(A_s),
E_s(n)=max{e<=n:e=s mod q}=n-((n-s) mod q).            (5)
```

Put `R={b_s}` and `M=Z/7\R`. For `r in M`, let `d^+_(s,r)` and
`d^-_(s,r)` be the representatives in `{1,...,6}` of `r-b_s` and `b_s-r`
modulo seven, and put

```text
c^+_(s,r)=d^+_(s,r)-f_s,
c^-_(s,r)=d^-_(s,r)-1+f_s,

rho^+_(p/q)(n)=max_(r in M) min_(0<=s<q) c^+_(s,r)/E_s(n),
rho^-_(p/q)(n)=max_(r in M) min_(0<=s<q) c^-_(s,r)/E_s(n).   (6)
```

Modulo finitely many rational walls, the noncover set in theta-space is the
disjoint circular union of the twelve owner components with centre `7p/q`
and radii `(rho^-,rho^+)`. The component at zero is interpreted across the
circle seam. In particular the exact measure identity is

```text
D(m)=(1/7) sum_(p/q in B)(rho^-_(p/q)(m-1)+rho^+_(p/q)(m-1)). (7)
```

### Proof of propagation from the finite base

Write `theta=7p/q+delta`. Along the track `e=s mod q`,

```text
floor(e theta) mod 7=floor(A_s+e delta) mod 7.         (8)
```

The exact base-component guards at `m=13` give `q|delta|<1` and preserve a
representative of every base occupied sector. Successive track values then
cannot skip an integer. For each missing sector, the largest time `E_s`
reaches it precisely at the corresponding cost in `(6)`; the inner minimum
chooses the first track and the outer maximum the last missing sector.
Nestedness of the noncover sets propagates this description to every later
`m`.

At a positive-drift wall, arrival occurs at equality. At a negative-drift
wall the floor convention can require strict passage. Two owners (`1/6` and
`5/6`) also have zero left radius at the sharp base. Therefore the set
identity is deliberately stated modulo finitely many walls, while `(7)` is
an exact equality of measures.

## 4. Sharp base and hostile predecessor

**FINITE-EXACT:** at `m=13` there are exactly twelve positive-length circular
components. Each contains exactly one owner in `(4)`, every owner occurs
once, and both individual radii equal `(6)`. The only zero track-start guard
margins are the right endpoints

```text
theta=16/3 for owner 3/4,
theta=17/4 for owner 3/5.                              (9)
```

At each wall, the left cell is noncovering while the wall and right cell are
covering, so the modulo-walls convention is correct. At `m=14`, all track
guards have strict margin at least `1/13`.

The immediate predecessor is sharp and hostile:

```text
m=12: 14 components,
D_direct(12)=72937/194040,
D_12-owner(12)=68647/194040.                           (10)
```

Thus `(7)` must not be interpolated backward.

## 5. Exact 60-phase rational law

Fix `n mod 60`, where `60=lcm(1,...,6)`. Every horizon in `(5)` has the form

```text
E_s(n)=n-c_s,                 0<=c_s<=5,               (11)
```

with a fixed shift on that phase. Every candidate radius is therefore
`C/(n-c)`. Comparisons between two candidates become linear inequalities
after cross multiplication. The companion audit proves, at the first
admissible `n` in each of the 60 phases and by the sign of every linear
coefficient, that no selector ever crosses later in that phase. It prints
the actual aggregated formula for every phase and checks it again one period
later. Consequently

```text
D(m)=(1/7) sum_j C_(r,j)/(n-c_(r,j)),
n=m-1, r=n mod 60, m>=13,                              (12)

C_(r,j)>=0, 0<=c_(r,j)<=5,
sum_j C_(r,j)=127/5 in every phase.                    (13)
```

The leading contribution grouped by owner denominator is

| denominator `q` | owners | contribution to `lim mD(m)` |
|---:|---:|---:|
| 1 | 1 | `11/7` |
| 2 | 1 | `9/14` |
| 3 | 2 | `2/3` |
| 4 | 2 | `5/14` |
| 5 | 4 | `12/35` |
| 6 | 2 | `1/21` |

They sum to `127/35`. Since every denominator in `(12)` lies between `n-5`
and `n`, `(13)` gives

```text
127/[35(m-1)] <= D(m) <= 127/[35(m-6)]     (m>=13),   (14)
D(m)=127/(35m)+O(m^-2),
lim_(m->infinity) mD(m)=127/35.                        (15)
```

## 6. Stopping theorems for scalar sequence models

The following stronger scalar descriptions are impossible.

1. A bounded eventual quasipolynomial is constant on each sufficiently large
   residue class. It is eventually periodic, and a nondecreasing eventually
   periodic sequence is eventually constant, contradicting `(2)`.
2. If `a(m)` satisfied an eventual homogeneous constant-coefficient linear
   recurrence over a characteristic-zero field, then so would `D(m)`. Its
   convergence to zero would leave only exponentially decaying characteristic
   roots, contradicting the polynomial lower bound in `(14)`.
3. A scalar unary finite-state output is eventually periodic. A standard
   finite-dimensional weighted unary linear representation over a
   characteristic-zero field is C-finite. Hence neither can represent this
   scalar sequence eventually.

The valid finite-state description retains the phase together with the
unbounded affine denominator `n-c`; discarding that sidecar destroys the
theorem.

## 7. LRC(14) consequence and boundary

THM-536 gives, for every `E subset {0,...,N}`,

```text
meas(S7(E))<=a(N+1).                                  (16)
```

Equations `(12)--(15)` explain the corrected AP threshold row

```text
N*(8),...,N*(13)=7,8,10,13,26,infinity.               (17)
```

For `k=12`, the last exact success is `a(27)<=6/7` and the next value fails;
for `k=13`, every finite AP row has positive deficit. When `a(N+1)>cap_k`,
the scalar inclusion `(16)` alone cannot certify an individual sparse set.
Progress must retain cardinality, gaps, owner labels, or another sparsity
sidecar.

This theorem does not prove that the AP maximizes `meas(S7(E))` among sparse
sets, does not supply the colored finite-denominator placement lift, and does
not prove LRC(14).

## 8. Replay

From the repository root:

```text
python3 -B 04-computation/lrc14_ap_cover_twelve_owner_rational_tail_thm4029.py
python3 -B -O 04-computation/lrc14_ap_cover_twelve_owner_rational_tail_thm4029.py
```

Both runs execute the same explicit `require` gates and reproduce the frozen
raw-LF stream. **QED.**
