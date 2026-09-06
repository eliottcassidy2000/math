---
id: THM-4441
title: "LRC14 signed (1,2,2) sharp ray closure"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4414 + FINITE-EXACT + INDEPENDENTLY
  AUDITED. Every primitive sorted distinct positive ternary-unit triple with
  a signed coefficient-magnitude (1,2,2) relation satisfies sharp
  min_i E_i<=46/665 and physical failure mass<=51/770. Both are below 6/77.
  Concurrent synthesis closes (1,1,1) and (1,1,2); entry, synchronization,
  and LRC(14) remain open.
source: root low-circuit continuation + independent referee, 2026-09-06
depends_on:
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
related:
  - THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits
primary_script: 04-computation/lrc14_signed_122_sharp_ray_closure_thm4441.py
primary_output: 05-knowledge/results/lrc14_signed_122_sharp_ray_closure_thm4441.out
primary_script_sha256: ebc145f97c57a040027be964ded10972911759a1ab6027ea7a32eb6d9c31403b
primary_output_sha256: b6d4ffbb8d2bab3bd07d64be11e748aa56dc90a6719936859b16629c007cf500
literal_script: 04-computation/lrc14_signed_122_literal_thm4441.cpp
literal_output: 05-knowledge/results/lrc14_signed_122_literal_h611_thm4441.out
literal_script_sha256: fdd593cae3e96c88015ebff163554cd67df66baed94ff87d22c1b8e987f573ad
literal_output_sha256: 877f4a68c2d391a2a6be5f28fbf556ba08388fac0540373a950ff64e273b495a
independent_script: 04-computation/lrc14_signed_122_cleanroom_referee_thm4441.py
independent_output: 05-knowledge/results/lrc14_signed_122_cleanroom_referee_thm4441.out
independent_script_sha256: 37513e8b35cac8a0f907ba0b3f3c01cd0d47e516fe98a50141d276342f302fd5
independent_output_sha256: 083b684f3261acd3c730a0dbefd00d5a67dd1def42fe7ad04b12e9ff07e342e2
audit: 05-knowledge/results/lrc14_signed_122_sharp_ray_closure_thm4441_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4441 -- LRC14 signed (1,2,2) sharp ray closure

**PROVED ELEMENTARY RELATIVE TO THM-4414 + FINITE-EXACT + INDEPENDENTLY
AUDITED.** The norm-five signed circuit is not an actual `6/77` hostile.
The concurrent [parity-free closure and exact threshold classification](../../05-knowledge/results/overnight4_20260906_lrc_parityfree_native.md#5-combined-nonadditive-ceiling-and-exact-old-threshold-classification)
also closes the two lower circuit families and combines them with this theorem.
Chart entry, synchronization, and `LRC(14)` remain **OPEN**.

## 1. Sharp statement

Let

```text
w=(a,b,c),  1<=a<b<c,  gcd(a,b,c)=1,  3∤abc,
```

and assume that `w` has a primitive signed integer relation whose coefficient
magnitudes are `(1,2,2)`.  Use the complete carrier set, raw projections, and
physical scale-three failure mass

```text
Lambda(w)={C in Z^3:C.w=0, 3∤C_i,
                     14|C_i|<3(w_j+w_k) for every i},
q=3/(7c),
e_i(C)=min(q,[3(w_j+w_k)-14|C_i|]/[14w_jw_k]),
E_i(w)=sum_(C in Lambda) e_i(C),
mu(F_w)=sum_(C in Lambda) min_i e_i(C).
```

Then the sharp all-height bounds are

```text
min_i E_i(w) <= 46/665,       equality iff w=(2,19,20),
mu(F_w)      <= 51/770,       equality iff w=(1,11,20).
```

Both constants are strictly below `6/77`.  Thus `(1,2,2)` is a residual of
the expanded coefficient-box argument, but not an actual projection hostile.

The four disjoint sorted sign cones have sharper individual constants:

| family | relation / primitive ray `u` | sharp `min E` | equality rows | sharp physical mass | equality row |
|---|---|---:|---|---:|---|
| F1 | `c=2(a+b)`, `(2,2,-1)` | `5/77` | `(1,10,22)` | `5/77` | `(1,10,22)` |
| F2 | `c=2(b-a)`, `(2,-2,1)` | `51/770` | `(1,11,20)` | `51/770` | `(1,11,20)` |
| F3 | `2c=a+2b`, `(1,2,-2)` | `46/665` | `(2,19,20)` | `173/2660` | `(2,19,20)` |
| F4 | `2c=2a+b`, `(2,1,-2)` | `3/49` | `(10,14,17)`, `(13,14,20)` | `731/12740` | `(13,14,20)` |

The exact extremal packets, all with positive address set `{1,2}`, are

```text
(1,10,22): E=(5/77,6/77,6/77),          mu=5/77;
(1,11,20): E=(51/770,3/35,3/35),        mu=51/770;
(2,19,20): E=(48/665,11/140,46/665),    mu=173/2660;
(10,14,17):E=(3/49,113/1190,3/49),      mu=461/8330;
(13,14,20):E=(3/49,149/1820,3/49),      mu=731/12740.
```

## 2. Exhaustive sign and arithmetic classification

After global sign reversal, positivity and `a<b<c` leave exactly four
possibilities.  If the coefficient of `c` has magnitude one, either `c` is
alone on one side, giving F1, or `c+2a=2b`, giving F2.  If its magnitude is
two, placing the unique coefficient one on `a` or `b` gives F3 or F4.
Every other sign placement contradicts the ordering.

The lossless primitive parameterizations are:

```text
F1: c=2(a+b), a<b, gcd(a,b)=1, a≡b (mod 3).
F2: c=2(b-a), b>2a, gcd(a,b)=1, a≠b (mod 3).
F3: a=2d, c=b+d, b>2d, gcd(b,d)=1, b≡d (mod 3).
F4: b=2d, c=a+d, d<a<2d, gcd(a,d)=1, a≡d (mod 3).
```

All displayed congruences are between nonzero classes modulo three.  They
are exactly the condition that the third speed is also a ternary unit.  No
all-odd triple occurs: the dot product with `(1,2,2)` modulo two forces at
least one even speed, while primitivity leaves exactly one or two even
speeds.  In F1/F2 this is controlled by whether `a,b` have equal parity; in
F3 by the parity of `(d,b)`; and in F4 by the parity of `(d,a)`.

The only potential overlap of distinct proper cones is F2/F3, which forces a
multiple of `(2,5,6)` and hence violates the ternary-unit condition.  The
other intersections violate positivity or strict ordering.  Thus the four
cases are disjoint in the typed universe.

## 3. The complete carrier set is exactly one ray

This is the key all-height step.  Put `L=ker_Z(w)`, let `u` be the displayed
primitive full-support norm-five relation, and set `r=3/14`.  The standard
carrier-zonotope lift writes every carrier as

```text
C=-w cross e,       |e_i|<r.
```

Primitivity of `w` makes `w cross (-):Z^3 -> L` surjective: if `q.w=1`,
then `w cross (u cross q)=u`. Choose `eta in Z^3` with
`w cross eta=u`, and define the integral covector
`phi(C)=eta dot C` on `L`. Its kernel in `L` is `Z u`. Then

```text
|phi(C)|=|u.e|<r ||u||_1=15/14.                       (1)
```

For ternary-unit `w`, any full-support mod-three kernel vector `X` satisfies

```text
X_1w_1=X_2w_2=X_3w_3  (mod 3),                        (2)
```

because three nonzero elements of `F_3` sum to zero only when they are all
equal. Both `C` and `u` obey `(2)`, so `C=+/-u (mod 3)`. Therefore
`phi(C)=0 (mod 3)`. Combined with `(1)`, this forces `phi(C)=0`, hence
`C in Z u`.

Conversely every live multiple under the roofs is a carrier.  Thus, exactly,

```text
Lambda(w)={k u: k in Z, 3∤k,
  0<|k|<B},
B=min_i 3(w_j+w_k)/(14|u_i|).                         (3)
```

This proves completeness, including the empty case; no other rational
direction or affine layer remains. No eligible address lies exactly on a
roof: `3` divides its right side, while `14k|u_i|` is nonzero modulo three
when `3` does not divide `k`.

## 4. Exact normalized projection and physical formulas

Normalize `z_i=w_i/c`, so `z_3=1`, and for `0<=x<B/c` define

```text
f_i(x)=min(2r,[r(z_j+z_k)-|u_i|x]/[z_jz_k]),
f_phys(x)=min_i f_i(x),                               (4)
```

setting all profiles to zero at and beyond the strict cutoff.  Equations
`(3)--(4)` give the exact discrete formulas

```text
E_i(w)=(2/c) sum_(k>=1,3∤k) f_i(k/c),
mu(F_w)=(2/c) sum_(k>=1,3∤k) f_phys(k/c).             (5)
```

The four normalized cones and support cutoffs are

| family | normalized speeds `z` | parameter range | `B/c` |
|---|---|---|---:|
| F1 | `(t,1/2-t,1)` | `0<t<1/4` | `r/2` |
| F2 | `(t,t+1/2,1)` | `0<t<1/2` | `r(1+t)/2` |
| F3 | `(t,1-t/2,1)` | `0<t<2/3` | `r(1/2+t/4)` |
| F4 | `(1-t/2,t,1)` | `2/3<t<1` | `r(1/2+t/4)` |

Direct cap-and-trapezoid integration gives the useful projection integrals:

```text
F1: I_1/r^2=7/8+t/4,            I_3/r^2=1-t+2t^2;
F2: I_1/r^2=(7+16t)/(8(1+2t));
F3: I_1/r^2=7/8+9t/16,          I_3/r^2=1-t/2+t^2/2;
F4: I_1/r^2=1-t/16,             I_3/r^2=1-t/2+t^2/2,
```

where `I_i=integral f_i`.  Consequently

```text
F1: min(I_1,I_3) <= (29/32)r^2, crossing at t=1/8;
F2: I_1 < (15/16)r^2;
F3: min(I_1,I_3) <= (121/128)r^2, crossing at t=1/8;
F4: min(I_1,I_3) <= (121/128)r^2, crossing at t=7/8. (6)
```

There is also a uniform and unexpectedly clean physical bulk identity:

```text
integral f_phys = (7/8)r^2                              (7)
```

in **all four sign cones**.  The envelope switches are:

```text
F1: f_1 -> f_3 at x=r(1-t)/2;
F2: f_1 -> f_2 at x=r/2;
F3: if t<=1/2, f_1 -> f_3 at x=r(1-t)/2;
    if t>=1/2, f_2 -> f_3 at x=rt/2;
F4: f_1 -> f_3 at x=rt/2.
```

Subtracting the relevant affine pieces verifies that the omitted coordinate
is nowhere lower.  Integrating either side of each switch yields `(7)`.
Since `(4/3)(7/8)r^2=3/56`, every norm-five sign cone has the same physical
continuous leading term; only the deleted-every-third-address discrepancy
distinguishes their finite sharp leaders.

## 5. Strict deleted-third quadrature and finite cutoff

For `T>0`, let

```text
R_<(T)=#{k:1<=k<T, 3∤k}.
```

With `M=ceil(T)-1`, three-term blocks give

```text
R_<(T)=M-floor(M/3) <=(2M+2)/3 <2T/3+2/3.             (8)
```

Layer cake applied to any profile in `(4)` therefore gives the strict rule

```text
(2/c) sum_(k>=1,3∤k) f(k/c)
  <(4/3) integral f + 4/(7c).                         (9)
```

Combining `(6)--(9)` yields

```text
family   network tail constant       tail beats sharp leader for
F1       87/1568 + 4/(7c)            c>=61
F2       45/784  + 4/(7c)            c>=65
F3       363/6272+ 4/(7c)            c>=51
F4       363/6272+ 4/(7c)            c>=171.
```

For every family the physical tail is

```text
mu(F_w)<3/56+4/(7c).                                  (10)
```

It beats the F1/F2/F3/F4 physical leaders from heights `51,46,50,151`,
respectively.  Hence the single exact head `c<=170` proves every sharp
claim in Section 1.

## 6. Exact checks

The exact head contains 1,951 typed rows, split

```text
F1/F2/F3/F4 = 280/559/744/368.
```

The Python producer reconstructs `(3)--(5)`, then independently enumerates
the complete raw integer-relation box under all three strict roofs.  Every
carrier set, projection, and physical mass agrees.  It also integrates the
capped affine profiles from independently generated breakpoints and checks
all displayed integral formulas, including the universal physical identity.
It performs 204,876 optimization-live checks; normal and `python -O` output
are byte-identical.

The C++ replay imports the literal common-sheet interval engine but is
independent of the ray formulas. Ordinary and optimized H170 outputs are
byte-identical and reproduce every leader above.  The optimized extension to
H611 contains 24,876 rows, split

```text
3553/7103/9483/4737,
```

with the same leaders and zero failures at `6/77`.

Reproduce from the session worktree root:

```powershell
python -B 04-computation/lrc14_signed_122_sharp_ray_closure_thm4441.py
python -B -O 04-computation/lrc14_signed_122_sharp_ray_closure_thm4441.py
g++ -std=c++20 -O0 04-computation/lrc14_signed_122_literal_thm4441.cpp -o pattern122_O0.exe
g++ -std=c++20 -O3 -DNDEBUG 04-computation/lrc14_signed_122_literal_thm4441.cpp -o pattern122_O3.exe
./pattern122_O0.exe 170
./pattern122_O3.exe 170
./pattern122_O3.exe 611
```

Canonical source hashes are recorded in THM-4441. The independent
[audit note](lrc14_signed_122_sharp_ray_closure_thm4441_independent_audit.md)
supplies a separate raw-carrier implementation and proof-line review.

## 7. Connection and next consequence

The source is a short full-support speed relation; the target is a complete
carrier ray.  The map preserves the mod-three owner class and the exact
strict roof.  Merely remembering the ray direction loses the natural address
and its deleted-third phase; `(3)` is the required sidecar.  Passing to the
continuum loses the finite extremizers; `(8)--(9)` plus the exact head repair
that loss.

The constant physical integral `(7)` is the unexpected structural survivor.
It is the norm-five analogue of the ratio-independent additive-family
physical bulk, but with leading term `3/56` rather than `9/98`.  This suggests
classifying other one-ray coefficient norms by the envelope area before
attacking residue discrepancy.  For the current reduction, however, the
decisive consequence is immediate: `(1,2,2)` may be removed entirely from
the list of possible `6/77` failures.  Only the `(1,1,1)` and `(1,1,2)`
families can remain actual low-circuit hostiles.
