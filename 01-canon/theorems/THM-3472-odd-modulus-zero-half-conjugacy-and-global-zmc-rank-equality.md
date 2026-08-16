---
id: THM-3472
title: "All-modulus fixed-zero-to-half transport and global ZMC-rank equality"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For every q>=2, the full
  zero-mode-cochain rank rho_ZMC(q) equals the literal half-twist cover rank
  rho_H(q), including infinity.  The proof uses the active-gcd divisor
  interface, all-modulus one-way fixed-zero-to-half cover transport with even
  self-opposite deletion, and divisor dilation; it does not assert augmented-
  primitivity preservation or even-modulus sheet conjugacy.  Consequently
  the complete cap-seven rank-priority atlas has minimal period
  14362718970600 and exact natural/harmonic coefficients.  The >7 stratum is
  not one asserted exact rank.  No endpoint-current, bispectrum, decrement,
  or LRC(14) conclusion follows.
source: codex-2026-08-15-all-modulus-layer-transport-repair
audit: >
  first independent audit found MISTAKE-407: doubled owners are not
  augmented-primitive and Q=8 blocks conjugacy but not cover transport;
  repaired self-contained active-gcd divisor formula and all-modulus
  one-way transport; exact 10766900-cell transport, 6478224-cell divisor,
  3434000-cell dilation, 9216-state all-modulus CRT, 576-state odd-subatlas,
  minimal-period, dependency, semantic, security, and normal/optimized
  replay gates; independent clean-room re-audit of the repaired proof,
  all counts, eleven witnesses, immutable hashes, and 2643-byte transcript
  passed
depends_on:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3415-zero-mode-cochain-global-rank-five-support
  - THM-3453-global-literal-half-twist-cap-seven-support-classification
related:
  - THM-3464-u-spine-q123-rank-eight-break-and-divisor-layer-certificate
  - THM-3469-three-times-p-half-twist-eight-owner-cover-boundary
script: 04-computation/lrc_odd_zero_half_conjugacy_global_rank_thm3472.py
output: 05-knowledge/results/lrc_odd_zero_half_conjugacy_global_rank_thm3472.out
script_sha256: 3fce197e22e143df5944eeb8814d43f08ce9e0d71ef73b0309077f2ba84e6ce7
output_sha256: 0d02b2a4dc67797a078dfbaefdc140778f78e33f72ee1f206583bc4daca6ed0c
semantic_sha256: 116818fa2bbc5a0cada41b425f08c4b7afb9a3051e17d804463f002b1027d81a
hash_basis: LF-normalized bytes
---

# THM-3472 -- all-modulus fixed-zero-to-half transport and global rank equality

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The first independent audit rejected the original primitivity claim and its
even-modulus boundary.  MISTAKE-407 records the exact failure.  The repaired
proof, deterministic companion, and fresh independent immutable-package
audit all pass.

## 1. The two ranks

For `Q>=2`, define fixed-zero and literal half-twist masks

```text
Z_(Q,s)={j in Z/QZ:14 dist_Q(sj,0)<Q},                 (1)

B_(Q,r)={ell in Z/QZ:
  14 dist_(2Q)(r(2ell+1),0)<2Q}.                       (2)
```

Owners are transverse, sign-normalized, and distinct.  Let `rho_H(q)` be
THM-3453's minimum number of masks `(2)` covering every sheet.  Let
`rho_ZMC(q)` be the minimum number of selected THM-3398 modes at one common
centre with complete cochain zero.  Both ranks take values in the positive
integers together with infinity.

Every literal half cover is a zero-cochain cover at common centre `1/(2q)`,
so for every `q>=2`,

```text
rho_ZMC(q)<=rho_H(q).                                   (3)
```

The reverse inequality requires two maps: first the active-gcd divisor
reduction, then a fixed-zero-to-half cover transport.  Neither map needs to
preserve augmented primitivity.

## 2. Deriving the divisor interface explicitly

Take any active zero-cochain cover at ambient modulus `q`.  In THM-3405's
notation write

```text
U=dV,       gcd(V)=1,       g=gcd(q,d),
q=gQ,       d=g d_0,        gcd(Q,d_0)=1.              (4)
```

After the common sheet relabelling of that theorem, its centre is

```text
c_epsilon=epsilon g/(2qd),       epsilon in {0,1}.     (5)
```

For an owner `u=dv` and sheet `ell`, direct substitution gives

```text
u(c_epsilon+ell/q)
 = epsilon v/(2Q)+d_0 v ell/Q.                         (6)
```

The right side depends only on `ell mod Q`.  Multiplication by `d_0` is a
sheet permutation modulo `Q`.  Thus:

- at `epsilon=0`, the ambient family is the `g`-fold pullback of a fixed-zero
  cover on `Q`;
- at `epsilon=1`, it is the `g`-fold pullback of a literal half cover on `Q`.

Transversality descends exactly because

```text
q|u  iff  Q|d_0v  iff Q|v.                             (7)
```

Coincident or sign-equivalent quotient owners may be deleted.  This derives
the divisor-minimum interface used in THM-3415 directly from THM-3405 and
pins the required layer and sheet permutations; it does not silently import
an augmented prime-breaker gate.

## 3. Fixed-zero cover transport for every modulus

Canonicalize a fixed-zero owner to `1<=s<=Q/2`; replacing `s` by `Q-s` does
not change `(1)`.  For every `s<Q/2` and every modulus `Q`,

```text
B_(Q,2s)(ell)=Z_(Q,s)(2ell+1 mod Q).                  (8)
```

Indeed

```text
dist_(2Q)(2a,0)=2 dist_Q(a,0),                         (9)
```

so the strict inequalities agree exactly.  The doubled owner satisfies
`0<2s<Q` and is transverse.

If `Q` is odd, `phi_Q(ell)=2ell+1` is a sheet permutation.  Hence `(8)` is a
literal mask conjugacy and preserves incidence cardinalities.

If `Q` is even, `phi_Q` is two-to-one onto the odd residue coset.  The only
canonical owner not covered by `(8)` is the self-opposite `s=Q/2`.  For odd
`j`, however,

```text
dist_Q((Q/2)j,0)=Q/2,
14(Q/2)<Q  is false.                                  (10)
```

Thus its fixed-zero mask is empty on the image of `phi_Q`.  If a fixed-zero
family covers all `Q` sheets, the remaining owners already cover every image
sheet.  Applying `(8)` therefore gives a transverse literal half cover of all
target sheets, with no more owners than the source family.

This is one-way cover transport, not an isomorphism of typed presentations.
It does not preserve augmented primitivity.  The minimal hostile is

```text
Q=15,  S=(1,2,3,4,5,7),       gcd(15,S)=1,
2S=(2,4,6,8,10,14),           gcd(30,2S)=2.            (11)
```

Both displayed families are valid covers in their respective layers.  Thus
the original “preserves primitivity” sentence was false, while the cover
transport needed for ranks survives.

## 4. Equality of the global ranks

Let a zero-cochain cover at `q` use `k` owners.  Equations `(4)--(7)` produce
a cover on some `Q|q` in one Boolean layer.

- Retain it if it is already a literal half cover.
- If it is fixed-zero, apply the all-modulus transport `(8)--(10)`.

In either case there is a literal half cover on `Q` with at most `k` owners.
Write `q=lambda Q`.  Directly from `(2)`,

```text
B_(lambda Q,lambda r)(ell)=B_(Q,r)(ell mod Q).         (12)
```

Hence dilation gives a literal half cover on `q` with at most `k` owners.
Therefore `rho_H(q)<=rho_ZMC(q)`.  Together with `(3)`,

```text
rho_ZMC(q)=rho_H(q) for every q>=2,                    (13)
```

including infinity.

Equation `(13)` is equality of minimum grades, not presentation ancestry.
Distinct primitive and inherited realizations can coexist, as THM-3464 shows
at `q=123`.  The transport also may delete the even self-opposite owner and
may merge sign-equivalent labels.

## 5. Complete all-modulus cap-seven support

THM-3453 gives the literal half-twist atoms

```text
rank 4:  8,9;
rank 5:  10,12;
rank 6:  11,15,23,25;
rank 7:  13,14,29,38,51,68,148.                       (14)
```

Equation `(13)` promotes that classification to the full zero-cochain rank.
For every `q>=2`, apply rank priority:

```text
rho_ZMC(q)=4 iff 8|q or 9|q;

rho_ZMC(q)=5 iff no rank-4 atom divides q and
                  (10|q or 12|q);

rho_ZMC(q)=6 iff no lower atom divides q and
                  one of 11,15,23,25 divides q;

rho_ZMC(q)=7 iff no lower atom divides q and
                  one of 13,14,29,38,51,68,148 divides q;

rho_ZMC(q)>7 otherwise.                               (15)
```

The final line is a cap-seven stratum, not one asserted exact rank.

## 6. All-natural-number and odd harmonic atlases

The word in `(15)`, extended to `q=1` by the `>7` symbol, has minimal period

```text
P=2^3*3^2*5^2*7*11*13*17*19*23*29*37
 =14,362,718,970,600.                                  (16)
```

A `9,216`-state weighted CRT product gives the exact census:

| stratum | count in one period | natural density / harmonic coefficient |
|---:|---:|---:|
| rank 4 | `3,191,715,326,800` | `2/9` |
| rank 5 | `1,276,686,130,720` | `4/45` |
| rank 6 | `1,734,627,895,000` | `25/207` |
| rank 7 | `1,531,452,347,040` | `165741596/1554406815` |
| rank `>7` | `6,628,237,271,040` | `717341696/1554406815` |

For every row with density `delta`,

```text
sum_(q<=N,q in row) 1/q=delta log N+O(1).              (17)
```

Every prime in `(16)` is essential.  In the format
`(prime,q,value,q+P/prime,value)`, changing-shift witnesses are

```text
(2,4,>7,7181359485304,4),
(3,3,>7,4787572990203,4),
(5,5,>7,2872543794125,6),
(7,2,>7,2051816995802,7),
(11,2,>7,1305701724602,6),
(13,5,>7,1104824536205,7),
(17,6,>7,844865821806,7),
(19,38,7,755932577438,>7),
(23,1,>7,624466042201,6),
(29,29,7,495266171429,>7),
(37,148,7,388181593948,>7).                           (18)
```

The original odd-only atlas remains a useful conditioned slice.  Mark even
indices by zero and classify odd indices using atoms `9`; `11,15,23,25`;
and `13,29,51`.  Its minimal period is `729,664,650`, its `576`-state census
for `(0,4,5,6,7,>7)` is

```text
(364832325,40536925,0,64859080,31171360,228264960),   (19)
```

and its ambient coefficients for ranks `4,6,7,>7` are

```text
1/18, 4/45, 283376/6633315, 691712/2211105.           (20)
```

Thus the repair strictly extends rather than discards the audited odd
harmonic computation.

## 7. Information contract and hostile boundary

At `Q=8`, `phi_Q` reaches only `{1,3,5,7}` and `Z_(8,4)` transports to the
empty half mask.  This is the sharp failure of sheet conjugacy and of
augmented-primitivity preservation.  It is not a failure of cover transport
or rank equality.

```text
source:      a fixed-zero or half cover on the active divisor Q|q
target:      a literal half cover on the ambient q
maps:        active-gcd quotient, ell->2ell+1, owner doubling, dilation
preserved:   strict incidence on the affine image, cover upper bound, grade
destroyed:   even off-image sheets, augmented primitivity, layer ancestry
sidecars:    Q,d_0,epsilon, canonical owner sign, self-opposite deletion
hostiles:    Q=15 augmented gcd; Q=8 nonbijective sheet map
```

The theorem gives no endpoint current, relation-residue coefficient,
bispectrum, physical LRC row, decrement, or LRC(14) conclusion.

## 8. Exact companion

Run from the repository root:

```bash
python 04-computation/lrc_odd_zero_half_conjugacy_global_rank_thm3472.py
python -O 04-computation/lrc_odd_zero_half_conjugacy_global_rank_thm3472.py
```

The standard-library companion checks `40,200` canonical owner rows and
`10,766,900` transport cells for every `2<=Q<=401`; the `Q=15` augmented-gcd
hostile and `200` even self-opposite deletions; `118,072` active-gcd divisor
rows and `6,478,224` cells through `q=80`; `15,150` dilation rows and
`3,434,000` cells at scales `2,3,5`; the exact `9,216`-state all-modulus and
`576`-state odd CRT censuses; all eleven period witnesses; dependency,
semantic, and security gates.  It uses explicit exceptions under `-O` and
performs no file write, dynamic evaluation, subprocess, or network action.
