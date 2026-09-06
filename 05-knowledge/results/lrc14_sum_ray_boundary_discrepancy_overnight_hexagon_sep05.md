# The additive ray: retain the projection jump and recover a two-sided physical law

**Status: PROVED exact discrepancy identities and all-height parity
obstruction; FINITE-EXACT smaller verification head; independent audit PASS.**
The sharp additive-family theorem `min_i E_i<=6/55`, with unique primitive
equality `(1,10,11)`, is already proved in the incoming
[additive-parity note](lrc14_additive_parity_empty_core_sep06.md), commit
`910dad3281880a9ec940d28a24fb784892b66c76`. This note credits that theorem
and supplies a stronger arithmetic remainder, not a second novelty claim
for its ceiling or equality witness.

The new two-sided consequence is: **every** primitive positive ternary-unit
sum triple `(a,b,c=a+b)` with `c>=49` has actual spoiled Haar measure greater
than `6/77`. This is an unbounded obstruction to extending the odd ceiling,
not an LRC counterexample or an arbitrary-body synchronization theorem.

## 1. Inheritance, failure boundary, and objects retained

Assume `1<=a<b<c=a+b`, `gcd(a,b)=1`, and `3` divides none of `a,b,c`.
The complete component network and physical spoiled set are the parity-free
objects in
[THM-4434 — universal scale-three network projection bound](../../01-canon/theorems/THM-4434-lrc14-universal-scale-three-network-projection-bound.md),
Section 1. We use its exact carrier identity, not its odd-only ceiling.

The closest proved mechanism is the period-three quadrature in
[THM-4373 — signed-121 resonant triple comb](../../01-canon/theorems/THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound.md),
Section 4. Its periodic function `E(t)` is called `P(t)` below. Its frozen
output already records `(1,10,11)` with actual mass `6/55`; the incoming
additive note recovers the earlier THM-4004 selector provenance. Neither
the witness nor this periodic primitive is new. The light-relation lock in
[THM-972 — relation lock and mediant triple count](../../01-canon/theorems/THM-972-relation-lock-mediant-triple-deathstar-S51.md)
is also an inherited operation. Its physical-radius theorem is not silently
transferred to the larger effective radius; the argument below proves the
needed version directly.

The corrected near miss is a pure `O(c^-2)` error for **every projection**.
Projected edge cost can remain positive when physical overlap disappears.
The least-used sidecar is the score at the strict raw-support boundary.
The live board is: complete ray; deleted residue class; affine roof;
physical overlap; boundary jump; actual body phase. Our map keeps complete
multipliers and the boundary score. A continuous integral loses the residue
primitive and boundary jump, and any tail statistic still loses the actual
body phase.

## 2. Complete ray and exact roofs, without assuming oddness

Write a carrier as `C=(x,y,z)`. The relation and primitivity give

```text
a(x+z)+b(y+z)=0,
x+z=bh, y+z=-ah, x-y=ch, h in Z.
```

The strict roofs imply `|x-y|<9c/14<c`, so `h=0`. Conversely the third roof
is the tightest on `C=k(1,1,-1)`. Thus the full dictionary is

```text
Omega={k(1,1,-1): 0<|k|<B=3c/14, 3 does not divide k}.       (1)
```

This is the incoming additive note's complete-ray proof. It also follows
from light-relation locking: effective errors of magnitude below `3/14`
force the integer `n_a+n_b-n_c` to have magnitude below `9/14`, hence to
vanish. No selected-ray restriction is imposed.

Put `q=3/(7c)` and, on the closed interval `[0,B]`, define

```text
f_i(t)=min(q, [3(w_j+w_k)/14-t]/(w_jw_k)),
h_a=3a/14, h_b=3b/14, h_c=3(a^2+b^2)/(14c).
```

Each `h_i` ends the corresponding plateau. The exact statistics are

```text
E_i=2 sum_(0<k<B,3∤k) f_i(k),
H=mu(F_w)=2 sum_(0<k<B,3∤k) min_i f_i(k).                    (2)
```

The closed endpoint values are *left limits*, not extra carriers:

```text
f_a(B)=f_b(B)=q/2,       f_c(B)=0.                           (3)
```

If `B` is integral, then `14|c` and `3|B`. Thus the strict cutoff never
removes a live integer at `B`. Inclusive counting at `B` below retains
exactly the required endpoint convention.

## 3. Exact projection error: the jump is a first-order term

For `t>=0`, set

```text
N(t)=floor(t)-floor(t/3),
D(t)=N(t)-2t/3,
P(t)=integral_0^t D(s) ds.
```

The inherited THM-4373 formula, with `r=t mod3`, is

```text
P(t)= -r^2/3              (0<=r<=1),
      r-1-r^2/3          (1<=r<=2),
      -(3-r)^2/3         (2<=r<3).
P(t+3)=P(t),   -1/3<=P(t)<=0,   -2/3<=D(t)<=2/3.            (4)
```

Define `R_i=(4/3) integral_0^B f_i(t) dt`. Integration by parts against
`N`, retaining the value at `B`, gives the **exact** law

```text
E_i-R_i = 2f_i(B)D(B)
          + 2[P(B)-P(h_i)]/(w_jw_k).                        (5)
```

Indeed `f_i'` is zero before `h_i` and `-1/(w_jw_k)` afterwards. Formula
(5) is also an elementary finite-sum identity. The jump vanishes for `i=c`,
but not for the other two projections.

Let `t=a/c`. The incoming continuum calculations are

```text
R_a=(9+3t)/98,
R_b=(12-3t)/98,
R_c=6(1-t+t^2)/49.
```

These are one-dimensional ray integrals, not the earlier two-dimensional
projection-area benchmarks. Since `f_a<=f_b`, `min_i E_i=min(E_a,E_c)`.
The increasing `R_a` and decreasing `R_c` cross at `t=1/4`, with common
value `39/392`. On `t<=1/4`, use `bc>=3c^2/4`; on `t>=1/4`, use
`ab>=3c^2/16`. Equations (3)--(5) imply

```text
t<=1/4: min E <= E_a <=39/392+2/(7c)+8/(9c^2),
t>=1/4: min E <= E_c <=39/392+32/(9c^2).                     (6)
```

Both right sides decrease with `c`; at `c=33` they are respectively
`418639/3841992` and `394783/3841992`, strictly below `6/55`.
This reduces the incoming sharp theorem's network head from `c<60` to
`c<33`: exactly 42 primitive ternary-unit sum triples. It is a smaller
independent certificate for the same sharp theorem.

The omitted-jump claim fails at its own extremizer. At `(1,10,11)`,
`E_a=6/55`, `R_a=51/539`, and `E_a-R_a>2/(3bc)=1/165`.
The quadratic integral correction alone cannot control `E_a`. The first
failed implication replaces a truncated positive projection roof by a
continuous compactly supported one. Formula (5) repairs that exact step.

## 4. Exact physical law and an unbounded parity obstruction

The physical minimum follows `f_a` up to `h_b` and `f_c` thereafter. Its
compact three-hinge representation is

```text
f_phys(t)=(B-t)_+/(ab)-(h_b-t)_+/(ac)-(h_a-t)_+/(bc).         (7)
```

Despite the signed terms this is nonnegative: it equals `q` before `h_a`,
the `a` roof from `h_a` to `h_b`, the `c` roof from `h_b` to `B`, and zero
afterwards. Apply the inherited hinge quadrature to (7). The identity
`c^3-a^3-b^3=3abc` gives

```text
H = 9/98 + 2P(B)/(ab)-2P(h_b)/(ac)-2P(h_a)/(bc),
|H-9/98| <= 2/(3ab) <= 2/[3(c-1)].                         (8)
```

Use `1/(ac)+1/(bc)=1/(ab)` and the range of `P` for the first bound, and
`ab-(c-1)=(a-1)(b-1)>=0` for the second. Thus `H` converges to `9/98`
**uniformly over all primitive sum triples**, including highly unbalanced
ratios, as `c` tends to infinity. The error is quadratic in height only
when the aspect ratio stays away from zero; (8) suppresses no qualification.

Since

```text
9/98-2/(3*48)-6/77 = 1/38808 >0,
```

every such triple with `c>=49` has `H>6/77`. The odd ceiling fails on an
entire unbounded family at the actual-overlap level, so changing the
projection selector cannot repair it. Common dilation preserves Haar mass;
the same statement holds when the **primitive reduced** height is at least
49, not merely when an unreduced height is large.

The inherited sharp ceiling and equality are independently reproduced:

```text
w=(1,10,11), E=(6/55,12/77,23/154), H=6/55.
```

The 42-row head has no other equality for either `min E` or `H`; the strict
tail (6) excludes larger heights. Neither `6/55` nor `6/77` is a proved
floor for an arbitrary actual ten-body safe set.

## 5. Reusable complete-ray sidecar

More generally, suppose an actual full carrier dictionary has the form
`{k d:0<|k|<B,3∤k}`, with `d` primitive and every coordinate a unit modulo
three. Put

```text
m_i=|d_i|/(w_jw_k),
f_i(t)=min(q,[3(w_j+w_k)/14-|d_i|t]/(w_jw_k)),
B=min_i 3(w_j+w_k)/(14|d_i|),
h_i=min(B, [3(w_j+w_k)/14-q w_jw_k]/|d_i|).
```

Here `q=3/(7 max(w))`. The displayed `h_i` is nonnegative since both
other speeds are at most `max(w)`. Integral `B` is again a deleted multiple
of three. The same argument proves

```text
E_i-(4/3) integral_0^B f_i
 =2f_i(B)D(B)+2m_i[P(B)-P(h_i)].                            (9)
```

The physical roof `f=min_i f_i` is concave on `[0,B]`, initially constant,
and vanishes at `B`. Its successive nonnegative slope magnitudes increase
from zero to at most `max_i m_i`. Apply (4) on each affine segment and
telescope to obtain

```text
|H-(4/3) integral_0^B f| <= (2/3) max_i m_i.                (10)
```

The coefficient of `P(B)` is the final slope; the negative coefficients at
slope-change points sum to that same slope, proving both error directions.
This recovers THM-4373's `4/(3pq)` remainder and (8)'s `2/(3ab)`.
The lemma requires a **complete** one-ray support and never permits
discarding other live directions or affine defects.

## Reproduction and audit

```bash
python3 -B 04-computation/lrc14_sum_ray_boundary_discrepancy_overnight_hexagon_sep05.py
python3 -B -O 04-computation/lrc14_sum_ray_boundary_discrepancy_overnight_hexagon_sep05.py
```

The [source](../../04-computation/lrc14_sum_ray_boundary_discrepancy_overnight_hexagon_sep05.py)
compares a closed residue-count formula, periodic identities, full raw
carrier sets, and native six-sheet interval projections for all 42 rows in
the complete `c<33` head. The earlier frozen one-ray script is imported
only for its independent raw and literal engines, not its odd-only theorem.
Controls retain the first mixed hostile, extremizer, false pure-quadratic
projection estimate, and actual-overlap obstructions above height 49.
The [output](lrc14_sum_ray_boundary_discrepancy_overnight_hexagon_sep05.out)
records `13,906` primary and `1,062` inherited raw/literal gates. Normal and
optimized runs agree byte-for-byte.

Independent `observer_collision` audit of Sections 2--4 and the complete
head: **PASS**, including cutoff, residue, projection/physical, inequality,
and equality checks. The separate independent audit of Section 5 also
**PASSES**: the integral cutoff, clipped plateau, and concave-roof
telescoping retain the stated constant without doubling it. No scarce ID
or shared navigation is edited.

Root independently read Sections 1--5 and audited the complete dictionary,
strict deleted cutoff, boundary integration, uniform height-49 lower bound,
and general concave-roof telescoping: **PASS**. The incoming sharp additive
theorem remains explicitly credited as a proved predecessor.

Frozen raw-LF SHA256:

```text
source 386a6c1c509dcbb47d5c9eee23894ca536d57f14bdafc8e554023a5d2dee8d98
output 04897b95d4a348f7c8bea7ab3c052d24a6b683f536b029cac73e4e8a88ef8dd8
```
