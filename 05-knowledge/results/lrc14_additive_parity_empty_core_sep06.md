# Sharp all-height bound for additive mixed-parity triples

**Status: PROVED ANALYTICALLY + FINITE-EXACT; independent proof and finite
head audit PASS.** For primitive positive ternary-unit triples with
`a<b<c=a+b`, both the selected complete network projection and the physical
scale-three failure mass are at most `6/55`, with equality exactly at
`(a,b,c)=(1,10,11)`. This is a statement about the additive family. The
`6/55` bound for arbitrary nonadditive mixed-parity triples remains **OPEN**;
no body Haar floor, synchronization, or LRC conclusion follows.

## 1. Exact scope, inherited objects and hostile boundary

Assume

```text
1<=a<b<c,     c=a+b,     gcd(a,b,c)=1,     3∤abc.
```

Thus `gcd(a,b)=1`, and exactly one speed is even. Retain the parity-free
complete carrier formulas from
[THM-4434 — universal scale-three network bound, Section 1](../../01-canon/theorems/THM-4434-lrc14-universal-scale-three-network-projection-bound.md):

```text
Lambda(w)={C in Z^3:C.w=0, 3∤C_i,
                     14|C_i|<3(w_j+w_k) for every i},
q=3/(7c),
e_i(C)=min(q,[3(w_j+w_k)-14|C_i|]/[14w_jw_k]),
E_i=sum_(C in Lambda) e_i(C),
mu(F_w)=sum_(C in Lambda) min_i e_i(C).
```

We use its displayed interval/carrier identity, not its odd-only ceiling.
The same identity is verified directly by literal sheets in both finite
head implementations. Its preservation data come from
[THM-4032 — affine defect boundary](../../01-canon/theorems/THM-4032-lrc14-d3-affine-defect-lattice-boundary.md).

The motivating hostile is the canonical `(1,10,11)` selector in
[THM-4004 — three-detuned divisor combs, Section 3](../../01-canon/theorems/THM-4004-lrc14-three-detuned-divisor-comb-profile.md).
The [bounded parity probe](lrc14_parity_empty_core_sep06.md) attached exact
mass `6/55` to that old object and refuted extending the odd `6/77` ceiling.
The near miss `(2,11,20)` from that probe is nonadditive and has two even
speeds; this theorem does not cover it. The useful sidecar is the complete
ray address, with every third multiplier removed, developed in
[THM-4425 — rank-one carrier closure](../../01-canon/theorems/THM-4425-lrc14-all-height-rank-one-carrier-closure.md).
Here we prove its applicability directly without transferring oddness.

## 2. The complete carrier set is one explicit ray

Write `C=(x,y,z)`. The relation condition becomes

```text
a(x+z)+b(y+z)=0.
```

Since `gcd(a,b)=1`, there is an integer `h` with
`x+z=bh`, `y+z=-ah`, hence `x-y=ch`. But the strict roofs imply

```text
|x-y|<=|x|+|y|<3[(b+c)+(a+c)]/14=9c/14<c.
```

Therefore `h=0` and `C=k(1,1,-1)`. The third roof is the tightest; the owner
residue condition is exactly `3∤k`. Consequently the complete set is

```text
Lambda(w)={k(1,1,-1): k in Z, 0<|k|<3c/14, 3∤k}.     (1)
```

This is an equality of full sets, including the empty case; there are no
unexamined rays or affine defects.

## 3. Normalized projection and physical profiles

Put `r=3/14` and `t=a/c`, so `0<t<1/2`. For `0<=s<r`, define

```text
f_a(s)=min(2r,[r(2-t)-s]/(1-t)),
f_b(s)=min(2r,[r(1+t)-s]/t),
f_c(s)=min(2r,[r-s]/[t(1-t)]).
```

Set all three functions to zero for `s>=r`, retaining the strict address
cutoff. In particular `f_a` and `f_b` have a downward jump at `r`.
Each function is nonnegative and nonincreasing, with value `2r` at zero.
Formula `(1)` gives exactly

```text
E_i=(2/c) sum_(k>=1, 3∤k) f_i(k/c),
mu(F_w)=(2/c) sum_(k>=1, 3∤k) f_phys(k/c),
f_phys=min(f_a,f_b,f_c).                               (2)
```

The difference between the uncapped `b` and `a` expressions is
`(1-2t)(r-s)/[t(1-t)]>=0`. Capping preserves this inequality, so
`f_b>=f_a` and `min_i E_i=min(E_a,E_c)`.

The plateau of `f_a` ends at `s=rt`; that of `f_c` ends at
`s=r(1-2t+2t²)`. Direct integration gives

```text
A(t)=(4/3) integral_0^r f_a(s) ds=(9+3t)/98,
C(t)=(4/3) integral_0^r f_c(s) ds=12(1-t+t²)/98.       (3)
```

The first is increasing and the second decreasing on `(0,1/2)`. Their
unique crossing there is `t=1/4`; therefore

```text
min(A(t),C(t))<=39/392.                                (4)
```

The physical profile follows `f_a` up to `s=r(1-t)` and `f_c` thereafter.
More explicitly, its three pieces are `2r` on `[0,rt]`, the uncapped
`a` expression on `[rt,r(1-t)]`, and the uncapped `c` expression on
`[r(1-t),r]`. Integration yields

```text
integral_0^r f_phys(s) ds=3r²/2,
(4/3) integral_0^r f_phys(s) ds=2r²=9/98.              (5)
```

Thus the physical leading term is independent of the speed ratio. The
physical and selected projection profiles still differ: the finite control
`(4,7,11)` has physical mass `215/2156` and selected projection `223/2156`.

## 4. Strict residue quadrature and all-height tails

For `T>0`, put

```text
R_<(T)=#{k in Z:1<=k<T, 3∤k}.
```

With `M=ceil(T)-1`, the exact formula and elementary bound are

```text
R_<(T)=M-floor(M/3)<=(2M+2)/3<2T/3+2/3.              (6)
```

Let `f` be any of the four profiles above. For every `0<y<2r`, its strict
superlevel set under `u -> f(u/c)` is `[0,T_y)` for some positive `T_y`.
This remains true for the downward jump at the address cutoff. Layer-cake
integration of `(6)` gives

```text
sum_(k>=1,3∤k) f(k/c)
  =integral_0^(2r) R_<(T_y) dy
  <(2c/3) integral_0^r f(s) ds+4r/3.
```

Multiply by `2/c` and use `8r/3=4/7`. Equations `(2)`–`(5)` imply

```text
min_i E_i <39/392+4/(7c),
mu(F_w)   <9/98+4/(7c).                               (7)
```

The first right side is below `6/55` for `c>=60`: the margin at 60 is
`1/12936`. Equivalently its real cutoff is `12320/207<60`. The physical
right side is already below `6/55` for `c>=34`, with margin `41/91630` at
34 and real cutoff `3080/93<34`. Both conclusions are strict on their tails.

## 5. Complete finite bases and equality

Every remaining network triple is obtained uniquely by

```text
1<=c<=59,  1<=a<c/2,  b=c-a,
gcd(a,c)=1,  3∤a(c-a)c.
```

There are exactly 136 rows. The corresponding independent physical head
`c<=33` has 42 rows. The
[new producer](../../04-computation/lrc14_additive_parity_empty_core_sep06.py)
uses the normalized formulas `(2)` and compares every row against both
the frozen parity probe's literal interval engine and its complete raw
modular carrier enumeration. All projections, physical masses and complete
carrier sets agree. The [frozen output](lrc14_additive_parity_empty_core_sep06.out)
has exactly one equality for either target:

```text
w=(1,10,11), Lambda={±(1,1,-1),±(2,2,-2)},
E=(6/55,12/77,23/154), mu(F_w)=6/55.
```

This finite equality statement, together with the strict tails in `(7)`,
proves the asserted all-height bound and its unique equality boundary.

Reproduce the primary finite base from the repository root:

```sh
python3 -B 04-computation/lrc14_additive_parity_empty_core_sep06.py
python3 -B -O 04-computation/lrc14_additive_parity_empty_core_sep06.py
```

The primary check has 1230 new explicit gates and 4792 inherited gates;
normal and `python -O` outputs are byte-identical. Its complete semantic
digest is
`4e0380e61e715e44f8329f5a80cc2eb42bfaa231207b2576dcaacffb6819d86c`.

The [independent audit](../../04-computation/lrc14_additive_independent_audit_empty_core_sep06.py)
and its [output](lrc14_additive_independent_audit_empty_core_sep06.out) use
the original September 5 one-ray native interval source, pinned by commit
and source hash, without importing the new producer. They separately check
the full raw support, all normalized profiles, piecewise trapezoid
integrals, both complete finite heads, equality, the physical/projection
distinction, and exact tail margins. All 1503 new gates and 4793 inherited
gates pass; normal and optimized audit outputs are byte-identical.
Root independently audited the complete-ray confinement, normalized caps,
both continuum integrals, physical crossover, strict layer-cake error,
cutoffs and finite equality basis: **PASS**.

```sh
python3 -B 04-computation/lrc14_additive_independent_audit_empty_core_sep06.py
python3 -B -O 04-computation/lrc14_additive_independent_audit_empty_core_sep06.py
```

Raw SHA-256 hashes:

```text
primary source a1c37cc674a72c324226dc2945dfdf86e13119ab5f8cbf0d00683c4ad48a1d2a
primary output a9038e29bd75b8db4ee49e6f3cb71727b1d06495de6d82cab61e39767cb3a6ab
audit source   be2080e4e95e4cec479f8bdbd69672db28de185f8dd19eaf5e25f376f59e3372
audit output   74268bcbb0feb7092421a1d4937ea3dd9f35a023b7383d6d92a04c887c66b59b
```

## 6. Consumer map and next boundary

The source is the complete affine carrier identity; the target is a sharp
Haar ceiling for additive tails. The map records the signed ray multiplier,
then evaluates the monotone capped roof profiles. It preserves every owner
residue, strict endpoint, projection capacity and physical overlap. Taking
only a projection or continuum integral loses the physical intersection or
the finite residue discrepancy; `(2)` and `(6)` are their required sidecars.

The [body-specific compact/open transfer](lrc14_universal_haar_consumer_empty_core_certificate_sep06.md)
therefore applies to this family with body Haar floor `6/55`. It still
requires that floor for the actual divided body. Common dilation preserves
physical Haar mass, so the physical statement also extends to nonprimitive
additive ternary-unit triples after dividing by their common gcd; physical
equality means primitive reduction `(1,10,11)`.

The next decisive test is the nonadditive mixed-parity coefficient universe,
including norm five and the all-parity norm-four finite base. This theorem
closes the complete norm-three family and does not substitute for that test.
The incoming `THM-4437`,
`01-canon/theorems/THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits.md`,
is currently **RESERVED / UNPROVED EMPTY STUB** on `origin/main` at
`e47c08b98d`; it proposes reducing the all-parity problem to signed patterns
`(1,1,1)`, `(1,1,2)` and `(1,2,2)`. This proved additive result is a candidate
input closing that proposed reduction's norm-three residual. The reservation
is not a dependency or a proved reduction.
