---
id: THM-4041
title: "LRC(14) d=2 affine defect-edge boundary and exact spoiled-phase law"
status: >
  PROVED RELATIVE TO THM-4024/4004 + VERIFIED-EXACT + HOSTILE-AUDITED.
  At the remaining d=2,c_2=9 equality boundary, the two exceptional body
  speeds are odd and full spoilage of the two labelled lifts is equivalent
  to one exact affine defect congruence and one strict window. The centre
  complex has one edge and no nontrivial circuit; its existence test is
  exactly THM-4004's reduced-sum condition. The complete spoiled-phase set
  has an exact component and measure law. Certificate-negative rows close,
  while positive rows still require the divided eleven-pack safe set. This
  does not prove LRC(14).
source: root + d2_affine_boundary / sequence continuation, 2026-08-24
audit: >
  PASS. Literal rational wall-cell subdivision, bounded affine-centre
  enumeration, defect-congruence enumeration, the reduced-sum criterion,
  and the exact measure/component formula agree on all 780 unordered
  distinct positive odd pairs through 79. The audit executes 6,957 explicit
  semantic gates and checks primitive/nonprimitive, endpoint, odd-dilation,
  H10-core, correctly typed eleven-pack, and strict full-row controls. Normal
  and optimized outputs are byte-identical.
depends_on:
  - THM-4024-lrc14-complete-divisor-incidence-envelope
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
related:
  - THM-4049-lrc14-d2-two-phase-residue-firewall
  - THM-4030-lrc14-d4-affine-defect-lattice-boundary
  - THM-4032-lrc14-d3-affine-defect-lattice-boundary
  - THM-4025-lrc14-owner-residue-odd-dilation-semigroup
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
script: 04-computation/lrc14_d2_affine_defect_edge_boundary_thm4041.py
output: 05-knowledge/results/lrc14_d2_affine_defect_edge_boundary_thm4041.out
script_sha256: 586a92c076b47f1d2be695771cc8b5e380327ba46cc7bad4184268ae489e6139
output_sha256: 8b55e6d9f35b14bb14fce2843e270449af5044539723bb995e7bb234cbd4be04
hash_basis: raw LF bytes
---

# THM-4041 -- the d=2 affine defect-edge boundary

**PROVED RELATIVE TO THM-4024/4004 + VERIFIED-EXACT +
HOSTILE-AUDITED.** The theorem completes the affine description of the three
small equality moduli from THM-4024. Its new phase-set law supplies the exact
sidecar needed for the next intersection test, but it does not prove
LRC(14).

## 1. Typed boundary and inheritance

Retain THM-4024's `d=2,c_2=9` boundary in the rank-eleven `11+2` branch:

```text
n=s(u_1,...,u_11) direct-sum t(p,q),
s in {1,2}, gcd(s,t)=gcd(p,q)=1.                       (1)
```

Here `2|t`, so `gcd(s,t)=1` forces `s=1`. Exactly nine body speeds are even,
and the two remaining body speeds are distinct positive odd integers

```text
alpha,beta.                                            (2)
```

The nine divided body speeds together with the two divided pair speeds form
an **eleven-speed** pack `H`. The cardinality is load-bearing: a literal
`{1,...,10}` can be a core of a control, but it is not the inherited divided
pack at `c_2=9`. Cited `LRCUpTo13` supplies a phase of `H` with clearance at
least `1/12`. For any pack phase `y in R/Z`, both lifts

```text
x_j=(y+j)/2,                    j=0,1,                 (3)
```

preserve every clearance of `H`.

Put

```text
D_w={x in R/Z: ||wx||<1/14}.                           (4)
```

Because each exception is odd, its two phases differ by `1/2`; hence it has
quotient order two and can spoil at most one label. Call `y` fully spoiled
when each label in `(3)` is dangerous for at least one of `alpha,beta`.

The closest proved mechanism is THM-4004's exact two-exception count. Its
corrected hostile is `(3,5)`, not the formerly described `(3,7)` minimum.
The least-used sidecar is the whole labelled spoiled-phase set, rather than
only its nonemptiness.

## 2. Exact affine-defect equivalence

Let

```text
g=gcd(alpha,beta),       alpha=ga, beta=gb,            (5)
```

where `a,b` are coprime odd integers. The following are equivalent.

1. Some phase `y` is fully spoiled.
2. There are integers `A,B` such that, for

   ```text
   N=2 beta A-2 alpha B+alpha beta,                    (6)
   ```

   one has

   ```text
   7|N|<alpha+beta.                                   (7)
   ```

3. There is an integer `N` satisfying `(7)` and the exact lattice sidecar

   ```text
   N == alpha beta (mod 2g).                           (8)
   ```

   Equivalently, `g|N` and `N/g` is odd.
4. The reduced exceptions satisfy

   ```text
   a+b>7.                                              (9)
   ```

Thus certificate absence closes the row, but the scalar existence test `(9)`
is exactly the condition already isolated in THM-4004.

Every fully spoiled phase has one exception on each label. Rotate the labels,
if necessary, so that `alpha` is dangerous at a real lift `z`; the reverse
assignment becomes this one after replacing `z` by `z+1/2`. Full spoilage is
therefore equivalent to

```text
z in D_alpha intersect (D_beta-1/2).                  (10)
```

Lift selected components of these open sets to real intervals. Their centres
and radii are

```text
c_alpha=A/alpha,          R_alpha=1/(14 alpha),
c_beta =B/beta-1/2,       R_beta =1/(14 beta).         (11)
```

Direct subtraction gives

```text
c_alpha-c_beta=N/(2 alpha beta),                       (12)
```

so the two intervals overlap exactly under the strict inequality `(7)`.
Conversely, `y={2z}` reconstructs the two labelled lifts. This proves
`1 iff 2`.

Formula `(6)` gives `(8)`. Conversely, suppose `(8)` holds. Write `N=gn`
with `n` odd. Then

```text
K=(N-alpha beta)/(2g)                                  (13)
```

is an integer, and coprimality of `a,b` lets Bezout solve

```text
bA-aB=K.                                               (14)
```

This reconstructs `(6)`, proving `2 iff 3`. The permitted defect lattice is
exactly `g` times the odd integers, so its least absolute value is `g`.
Consequently `(7)--(8)` have a solution exactly when `(9)` holds.

There is no attainable endpoint equality. After division by `g`, the value
`a+b` is even, whereas `7|N|/g` is odd. Thus the open danger convention leaves
an automatic arithmetic gap at every possible defect wall.

## 3. Why there is no circuit

The selected-centre complex has two vertices and one edge. Its cycle space
has rank zero. After quotienting the simultaneous-translation gauge

```text
(A,B)->(A+alpha,B+beta),                               (15)
```

the single edge coordinate is governed completely by `(8)`. There is no
nontrivial defect-circuit equation analogous to the support-three circuits in
THM-4030 and THM-4032. Retaining both edge orientations adds only
antisymmetry, not another obstruction.

This is a stopping theorem for one tempting route: no further scalar
Diophantine gate can be obtained by searching for a missing `d=2` circuit.
Further progress must retain how the one-edge spoiled set meets the divided
eleven-pack safe set.

## 4. Exact spoiled-phase law

Let `Sigma_(alpha,beta)` be the open set of pack phases `y` that are fully
spoiled. Swap the reduced names so that `a<b`, and for every positive odd
integer `r` with `7r<a+b` put

```text
L_r=min(1/(7b), (a+b-7r)/(14ab)).                      (16)
```

If

```text
q=#{r>0: r odd and 7r<a+b},                            (17)
```

then

```text
meas(Sigma_(alpha,beta))=4 sum_r L_r,
# positive-length components of Sigma_(alpha,beta)=2gq. (18)
```

Indeed, for the reduced pair, signed odd defects index the disjoint
fixed-orientation component intersections. At defect magnitude `r`, the two
centre radii and centre distance give intersection length

```text
min(2/(14b), 1/(14a)+1/(14b)-r/(2ab))=L_r.            (19)
```

The other label orientation is the half-turn translate. The quotient
`y=2z` identifies those half-turn partners and preserves normalized Haar
measure, giving `(18)` for `g=1`. Restoring the common odd dilation pulls the
reduced phase set back under `y->gy`; this preserves measure and multiplies
the number of components by `g`.

This separates three ordinal shadows of odd dilation. The existence bit and
measure are invariant, while component multiplicity scales by `g`. Reducing
to the primitive pair is therefore lossless for `(9)` but lossy for phase
placement.

## 5. The exact remaining pack-side certificate

Write

```text
G(H)={y: ||hy||>=1/14 for every h in H}.               (20)
```

If some `y in G(H)` is not in `Sigma_(alpha,beta)`, then one label in `(3)`
is safe for both exceptions, while both labels preserve every `H` clearance;
that lift is a lonely phase for the full row. Hence a hypothetical non-lonely
row at this boundary must satisfy the strictly stronger condition

```text
G(H) subset Sigma_(alpha,beta).                        (21)
```

The affine existence test forgets `G(H)` and cannot decide `(21)`. Formula
`(18)` now makes the next test exact: construct both rational interval sets
and decide their inclusion. In particular,

```text
meas(G(H))>meas(Sigma_(alpha,beta))                    (22)
```

is a sufficient closure certificate, though no uniform inequality `(22)` is
claimed here.

Downstream THM-4049 supplies a different exact certificate: if every member
of `H` avoids ten named classes modulo `56`, the four lifts over
`y in {1/14,5/56}` close the row for every odd exception pair, with no height
bound. In particular it closes the explicit pack `{1,...,10,12}` uniformly.
A typed THM-3818 row can hit THM-4049's forbidden class `11`, so no uniform
physical entry is possible and `(21)` remains open on the genuine
complementary image.

## 6. Sharp and typed controls

The reduced pairs `(1,5)` and `(1,3)` are selector-safe. The nonprimitive
pair `(3,9)` remains safe after reduction. The smallest possible maximum
speed of a fully spoiled distinct pair is five:

```text
(alpha,beta)=(3,5),        z=3/16,
bad labels=({1},{0}).                                   (23)
```

Its spoiled-phase measure is `2/105` and it has two `y` components. The odd
dilation `(9,15)` is hostile at `z=1/16`; it has the same measure and six
components, providing a direct control for `(18)`.

For a correctly typed pack containing the requested consecutive ten-core,
take

```text
H_10={1,...,10},             H=H_10 union {12},
y=1/11.                                                (24)
```

The eleven-pack `H` has clearance `1/11`. It is realized by

```text
s=1, t=2, (p,q)=(1,12),
body=(4,6,8,10,12,14,16,18,20,1,11).                  (25)
```

The body is primitive, exactly nine body owners are even, and all thirteen
full-row speeds are distinct. On the lifts of `(24)`, exceptions `(1,11)`
have masks

```text
({0},{1}),                                             (26)
```

so the selected pack phase is fully spoiled. Yet the same full row is
strictly lonely at

```text
x=229/560,             clearance=5/56>1/14.           (27)
```

Thus `(24)--(27)` are a hostile to arbitrary pack-phase selection, not an
LRC counterexample. They are not claimed to survive every inherited
high-height or decoder filter.

## 7. Audit and correction lineage

The exact companion compares literal rational wall cells, affine-centre
enumeration, defect-lattice enumeration, `(9)`, and the measure/component law
on all `780` unordered distinct positive odd pairs through `79`. It checks
`6,957` explicit semantic gates. No truth gate uses Python `assert`, so the
optimized replay retains the audit.

The proof labels both inverse fibres rather than calling their union
bijective (MISTAKE-393), keeps quotient order separate from literal gcd
(MISTAKE-390), and tests walls as well as positive cells (MISTAKE-464). The
strict control `(27)` makes its non-counterexample status independent of
closed-wall topology.

This theorem is relative to THM-4004's cited eleven-speed input and
THM-4024's `d=2,c_2=9` reduction. It does not prove `(21)` uniformly, close
the other arbitrary-body branches, or prove LRC(14).

## 8. Replay

```text
python3 -B 04-computation/lrc14_d2_affine_defect_edge_boundary_thm4041.py
python3 -B -O 04-computation/lrc14_d2_affine_defect_edge_boundary_thm4041.py
```

Both runs reproduce the frozen raw-LF output. **QED.**
