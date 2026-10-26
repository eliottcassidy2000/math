---
id: THM-4402
title: "LRC14 norm-eighteen vanishing live-carrier gap"
status: >
  PROVED ELEMENTARY + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the fixed
  minimal norm-18 coefficient sector (1,1,16), zero-defect raw carrier
  components can have arbitrarily small positive length. Therefore no
  positive component quantum depending only on this discrete shell data
  exists. This refutes one uniform-gap proof route, not LRC(14).
source: root + independent audit / LRC14 and percolation-transfer session, 2026-09-03
depends_on:
  - THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence
related:
  - THM-4393-lrc14-minimal-ternary-unit-norm-eighteen-shell
  - THM-4396-lrc14-finite-dual-exact-pair-hybrid-certificate
primary_script: 04-computation/lrc14_norm18_vanishing_carrier_gap_thm4402.py
primary_output: 05-knowledge/results/lrc14_norm18_vanishing_carrier_gap_thm4402.out
primary_script_sha256: d1dccd6e832e8ad762e37b5bf34c1cfc5fb3a33d17a17ee90a8c7b2985c622fd
primary_output_sha256: fb470ef1a49bdca8abe3bae49b1cd7b475c5d0a7e2e9cbd0316c2e699e26fe77
independent_audit_script: 04-computation/lrc14_norm18_vanishing_carrier_gap_thm4402_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_norm18_vanishing_carrier_gap_thm4402_independent_audit.out
independent_audit_script_sha256: f8d1187fa5b109ff0f98adc82b5a9fccaee39cb80896d2af4ce8528bc3a33525
independent_audit_output_sha256: e583efb53af2042363f15da5a61f5474d5966765aeb28222a3d5bf26ea43e2b5
hash_basis: raw LF bytes
audit: >
  PASS. The primary exhausts the complete coefficient l1-ball through norm
  17 uniformly in the speed parameter and symbolically checks the six-term
  carrier formula. The dependency-free referee reconstructs both the raw
  formula and literal interval intersection, checks 21,084 live predicates
  on 1,002 speed values including two large hostile controls, and brute-forces
  the full relation lattice through norm 18. Normal, optimized, and fixed
  hash-seed runs byte-match the frozen LF outputs; neither script contains an
  assert statement.
---

# THM-4402 -- LRC14 norm-eighteen vanishing live-carrier gap

**PROVED ELEMENTARY + VERIFIED-EXACT + INDEPENDENTLY AUDITED. THIS IS A
NO-GO THEOREM FOR A UNIFORM COMPONENT GAP, NOT A COUNTEREXAMPLE TO OR PROOF OF
`LRC(14)`.**

## 1. Statement

Use the LRC14 triple-comb radius `r=3/14` and raw-carrier conventions of
THM-4386.  For every integer

```text
m=6t+5,                 t>=2,
```

put

```text
w_m=(1,m,16m-1),
c=C=(1,-16,1).                                             (1)
```

Then:

1. `w_m` is a primitive sorted triple of distinct positive odd integers, all
   prime to three.
2. `c dot w_m=0`, and `c` is a primitive full-support ternary-unit relation.
   In fact the least `l1` norm of *any* nonzero integer relation on `w_m` is
   exactly `18`.
3. The raw carrier `C` is physically live, has defect zero and owner
   permutation `(0,1,2)`, and

```text
L_(w_m)(C)=3/[7(16m-1)]>0.                               (2)
```

Consequently

```text
inf_m L_(w_m)(C)=0.                                      (3)
```

The infimum is already zero after fixing relation-magnitude pattern
`(1,1,16)`, defect `0`, and owner permutation `(0,1,2)`.  Thus there is no
positive lower quantum for a nonempty raw component which depends only on
those discrete data.

## 2. Admissibility and exact relation minimum

The congruence in `(1)` gives

```text
m=1 mod 2,       16m-1=1 mod 2,
m=2 mod 3,       16m-1=1 mod 3.
```

Also `1<m<16m-1`, and the leading coordinate `1` makes the speed triple
primitive.  Directly,

```text
1-16m+(16m-1)=0,
```

so `c` is a primitive ternary-unit relation of norm `18`.

For the sharp lower bound, let `u=(a,b,d)` be any integer relation.  Its
equation is

```text
(a-d)+(b+16d)m=0.                                       (4)
```

If `b+16d` is nonzero and `||u||_1<=16`, then

```text
m=|a-d|/|b+16d| <= |a|+|d| <=16,
```

contrary to `m>=17`.  If `b+16d=0`, equation `(4)` gives `a=d`; a nonzero
solution is `d(1,-16,1)` and has norm `18|d|`.  Finally, because all three
speeds are odd,

```text
||u||_1 congruent a+b+d congruent u dot w_m congruent 0 (mod 2).
```

No relation has odd norm, so there is none of norm `17`.  This proves the
claimed minimum.  At `m=17` another norm-18 vector exists; no uniqueness is
claimed or needed.

## 3. One explicit physical carrier

Choose the nearest-integer lift

```text
n=(0,1,16).
```

Then

```text
w_m cross n=(1,-16,1)=C,       c dot n=0.               (5)
```

Modulo three, the speeds are `(1,2,1)` and the lift is `(0,1,1)`.  Hence the
owners

```text
o_i=-w_i^(-1)n_i mod 3
```

are exactly `(0,1,2)`.  Equivalently every coordinate of `C` is nonzero
modulo three, as required by THM-4386's exact owner gate.

Write `W=16m-1`.  The three strict nearest-integer intervals belonging to
`n` are

```text
I_1=(-3/14,3/14),
I_2=(11/(14m),17/(14m)),
I_3=(221/(14W),227/(14W)).                              (6)
```

For `m>=17`, `I_3` lies strictly inside both other intervals.  After clearing
positive denominators, the two comparisons with `I_2` have numerators

```text
45m+11,       45m-17,
```

and the upper comparison with `I_1` has numerator `48m-230`.  All are
positive.  Thus the physical component is exactly `I_3`, and its length is

```text
227/(14W)-221/(14W)=3/(7W),
```

which proves `(2)`.  Letting `t` tend to infinity proves `(3)`.

## 4. What the result closes, and what it does not

The Kozma--Nitzan-inspired Walsh/Bessel scout asked for a positive atom size
so that an average error could be made smaller than every live packet.  This
theorem proves that no such uniform atom size can depend only on the fixed
norm-18 shell, coefficient pattern, defect, and owner permutation.  In
particular, any fixed-dimensional Walsh/Bessel argument whose final step is
solely `error < universal live-component gap` cannot cover this shell.

The conclusion does **not** rule out:

- an owner-conditioned or speed-weighted inequality whose error decays with
  `m`;
- aggregation across all carriers or a positive lower bound for total comb
  measure;
- another finite-dual certificate with no componentwise gap step; or
- `LRC(14)` itself, which remains open.

## 5. Exact audit and reproduction

The primary certificate expands the six carrier-length candidates, proves
their minimum is `(2)` for every `m>=17` by coefficient positivity after a
shift, and exhausts all `7,174` nonzero coefficient vectors in the `l1<=17`
ball as affine equations in `m`.

The independent audit imports no repository code and uses exact `Fraction`
arithmetic.  It reconstructs the literal interval intersection, checks the
full relation minimum through norm 18 on 1,002 members of the family, and
includes large-speed controls.  Its finite sweep is a hostile check; the
uniform proof is the affine identity `(4)` and interval containments above.

Replay from the repository root:

```text
python3 -B 04-computation/lrc14_norm18_vanishing_carrier_gap_thm4402.py
python3 -B -O 04-computation/lrc14_norm18_vanishing_carrier_gap_thm4402.py
python3 -B 04-computation/lrc14_norm18_vanishing_carrier_gap_thm4402_independent_audit.py
python3 -B -O 04-computation/lrc14_norm18_vanishing_carrier_gap_thm4402_independent_audit.py
```

Normal and optimized runs byte-match the frozen LF outputs.  All theorem
checks remain live under optimization.
