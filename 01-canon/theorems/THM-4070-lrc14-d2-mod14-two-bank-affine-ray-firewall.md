---
id: THM-4070
title: "LRC(14) d=2 complete small-denominator sieve and affine-ray firewall"
status: >
  PROVED RELATIVE TO THM-4041 FOR THE TYPED RESIDUAL CONSEQUENCE +
  VERIFIED-EXACT + INDEPENDENTLY AUDITED. For every q=2,...,14, absence of
  a q-divisible divided-pack speed closes every even-pack/two-odd-exception
  row by an explicit bank of at most four times. Hence a nonlonely survivor
  on THM-4041's typed d=2 boundary must divisibility-cover every q in that
  range. The q=14 bank strictly enlarges THM-4049's allowed residue set from
  46 to 52 classes, closes its typed physical hostile, and closes complete
  odd affine rays. Exact primitive defect windows, odd pullback, and the
  necessary-and-sufficient pack-safe containment criterion are also proved.
  The same banks handle arbitrarily many odd exceptions for q<=7 and up to
  three for even q=8,10,12,14. This is a scoped certificate theorem, not an
  untyped d=3 reduction and not LRC(14).
source: codex-frontier-synthesis-creative-20260825c / LRC d=2 affine-intercept lane
audit: >
  PASS. The primary exact-Fraction path reconstructs 364 primitive odd
  profiles below 60, checks 88,056 literal wall/cell memberships, checks 636
  odd pullbacks over 106 primitive bases, and audits all 559 natural odd-pair
  sieve profiles plus all 1,264 even-q unordered triple profiles. The
  independent integer path imports none of the primary code, loads every
  allowed residue simultaneously for every q, and redundantly checks all 406
  unordered odd pairs modulo 56 at q=14.
  It reproduces the strict THM-4049 gain, equality controls, physical hostile,
  and q=15 grid obstruction. Normal and optimized outputs are byte-identical;
  both scripts have zero assert nodes and zero float literals.
depends_on:
  - THM-4041-lrc14-d2-affine-defect-edge-boundary
  - THM-4049-lrc14-d2-two-phase-residue-firewall
  - THM-4052-lrc14-affine-component-width-escape-cones
  - THM-4066-lrc14-diagonal-intercept-pullback-and-exact-affine-ray-closure
related:
  - THM-4062-lrc-divisor-star-affine-intercept-obstruction
  - THM-615-folding-identity-and-AP-even-part-confinement
script: 04-computation/lrc14_d2_mod14_two_bank_affine_ray_firewall_thm4070.py
output: 05-knowledge/results/lrc14_d2_mod14_two_bank_affine_ray_firewall_thm4070.out
independent_audit_script: 04-computation/lrc14_d2_mod14_two_bank_affine_ray_firewall_thm4070_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_d2_mod14_two_bank_affine_ray_firewall_thm4070_independent_audit.out
script_sha256: 57f756a199af010f0e4581e58895e7724f1f567a02d8196b316d39758da06345
output_sha256: a00ef1bff9c9150d54b4f8a947d99b9c287be051fc80d2b00f16039565a233a8
independent_audit_script_sha256: 1489822e2c89006a10ec7f08730e8de721a3e756f0e5269a8c14e30387e29abc
independent_audit_output_sha256: e87f558ba1c3aed21ae2fa452ed6c40348cac4daf5aef89c249b8bc8e2940e50
hash_basis: raw LF bytes
---

# THM-4070 -- complete small-denominator sieve and affine-ray firewall

**PROVED RELATIVE TO THM-4041 FOR THE TYPED RESIDUAL CONSEQUENCE +
VERIFIED-EXACT + INDEPENDENTLY AUDITED.** A complete residue-only sieve for
denominators `2<=q<=14` now forces every still-nonlonely typed `d=2` row to
carry a divided pack owner for every one of those divisors. At `q=14`, the
new four-time bank strictly enlarges THM-4049's firewall and closes its
canonical physical hostile. The same argument closes infinite odd-dilation
affine rays. LRC(14) remains open.

Throughout,

```text
||z||=min_{m in Z}|z-m|                              (1)
```

is distance to the nearest integer. Safety is the **closed** inequality
`||z||>=1/14`; danger is the complementary **strict** inequality
`||z||<1/14`. This convention is load-bearing at the boundary times below.

## 1. Inheritance and the retained coordinate

The closest proved mechanism is THM-4041's exact `d=2` defect edge, together
with THM-4066's unit-dilation conjugacy. The canonical hostile is THM-4049's
typed physical row whose divided pack contains residue `11 mod 56`. The
corrected near miss is MISTAKE-490: THM-4049's residue theorem survived, but
its physical projection had to retain exception parity and acknowledge the
existing, unexecuted THM-3818 finite box. The least-used relevant sidecar is
THM-615's finite argmax bank for the AP even pack; it suggests retaining the
actual rational pack phases, not only the measure of their safe set.

For a finite integer pack `H`, put

```text
G(H)={y in R/Z: ||r y||>=1/14 for every r in H}.       (2)
```

For an odd exception `w` and a chosen representative `y`, retain the labelled
danger mask

```text
D_w(y)={j in C_2: ||w(y+j)/2||<1/14}.                 (3)
```

For an odd pair `E`, its fully spoiled pack phases are

```text
Sigma_2(E)={y in R/Z: union_{w in E}D_w(y)=C_2}.      (4)
```

Changing the representative of `y` rotates both labels. Coverage in `(4)`
is gauge invariant, while the individual masks in `(3)` are not.

The connection proved below is

```text
source:    THM-4041's labelled d=2 lifts over G(H)
target:    primitive open defect arcs and finite q-grid danger masks
map:       y=2x, then (alpha,beta)=h(a,b) and evaluation at y=1/q,k/q
preserves: strict spoilage, pack safety, and existence of an actual safe lift
loss:      the h-sheet under y->hy, off-grid safe-set geometry, physical origin
sidecar:   h odd, exception parity, and the divisibility vector (q|r)_{r in H}
test:      exact wall pullback plus maximal-residue odd-pair exhaustion.    (5)
```

## 2. Exact primitive windows and odd pullback

Let `1<=a<b` be coprime odd integers. For every signed odd integer `n` with

```text
7|n|<a+b,                                              (6)
```

choose integers `A_n,B_n` satisfying

```text
2b A_n-2a B_n+ab=n.                                   (7)
```

Such a pair exists because `(n-ab)/2` is integral and `gcd(a,b)=1`. Any two
choices differ by the simultaneous translation `(A,B)->(A+aL,B+bL)`, which
does not change the following circle arc. Define the real open interval

```text
J_n=( max(A_n/a-1/(14a), B_n/b-1/2-1/(14b)),
      min(A_n/a+1/(14a), B_n/b-1/2+1/(14b)) )         (8)
```

and its doubled circle image

```text
I_n={2z mod 1:z in J_n}.                              (9)
```

Then the exact spoiled set is the disjoint union of open arcs

```text
Sigma_2({a,b})=disjoint_union_{n odd, 7|n|<a+b} I_n, (10)
```

and

```text
|I_n|=min(2/(7b),(a+b-7|n|)/(7ab)).                  (11)
```

In particular, `Sigma_2({a,b})` is empty exactly when `a+b<=7`; otherwise
its largest component has length

```text
W(a,b)=min(2/(7b),(a+b-7)/(7ab)).                    (12)
```

To prove `(8)--(11)`, select the lift `z` on which `a` is dangerous. The
opposite lift is dangerous for `b` exactly when two open real intervals with
centres and radii

```text
A/a,       1/(14a),
B/b-1/2,   1/(14b)                                   (13)
```

overlap. Their signed centre difference is `n/(2ab)` by `(7)`. The
intersection length is therefore

```text
min(1/(7b),(a+b-7|n|)/(14ab)),                       (14)
```

and doubling gives `(11)`. Each fully spoiled phase has a unique opposite-
label assignment after rotating the two labels, and unique danger teeth;
hence the arcs are disjoint and exhaustive. All intervals are open. There is
also no attainable equality `7|n|=a+b`, because its two sides have opposite
parity.

Now let `h` be a positive odd integer. Since `h=1 mod 2`, label by label,

```text
D_{ha}(y)=D_a(hy),        D_{hb}(y)=D_b(hy).          (15)
```

Thus, for the degree-`h` circle map `[h](y)=hy`,

```text
Sigma_2({ha,hb})=[h]^(-1)Sigma_2({a,b}).              (16)
```

Every primitive arc `I_n` has exactly `h` preimages, each of length
`|I_n|/h`. Consequently the largest scaled component has the exact width

```text
W_h(a,b)=(1/h)min(2/(7b),(a+b-7)/(7ab))              (17)
```

when `a+b>7`. Oddness of `h` is exactly the unit sidecar in `C_2`; even
dilation changes the label type and `(15)` is false.

## 3. Exact pack-safe containment criterion

For any finite positive pack `H`, coprime odd `a<b`, and positive odd `h`,
put

```text
V=2H union {ha,hb}.                                   (18)
```

Then

```text
V is 1/14-lonely
  iff G(H) \ [h]^(-1)Sigma_2({a,b}) is nonempty,      (19)

V is not 1/14-lonely
  iff G(H) subset [h]^(-1)Sigma_2({a,b}).             (20)
```

Indeed, for `y in G(H)` both lifts `x_j=(y+j)/2` preserve every pack
inequality:

```text
||2r x_j||=||r(y+j)||=||ry||.                        (21)
```

If `y` is outside the spoiled set, one label is safe for both exceptions and
that lift is safe for the full row. Conversely, from any full-row safe time
`x`, take `y=2x mod 1`. Then `y in G(H)` and the label containing `x` is not
spoiled. This proves both directions, including all boundary points. Formula
`(20)`, not merely the scalar component width `(17)`, is the exact remaining
affine-intercept obstruction.

## 4. Complete small-denominator sieve

Let `H` be any finite set of integers, let `alpha,beta` be odd integers with
repetition allowed, and fix

```text
2<=q<=14,       q does not divide r for every r in H. (22)
```

Then `2H union {alpha,beta}` is `1/14`-lonely.

For `2<=q<=7`, the single time

```text
x=1/(2q)                                               (23)
```

works. For every pack owner, `||2rx||=||r/q||>=1/q`.
An odd exception has a nonzero odd residue modulo `2q`, so
`||alpha x||,||beta x||>=1/(2q)>=1/14`.

For `8<=q<=14`, use the following coprime multiplier `k` and four-time bank:

```text
T_q={1/(2q),(1+q)/(2q),k/(2q),(k+q)/(2q)}.           (24)

q     k     D_1          D_(1+q)     D_k          D_(k+q)
8     3     +/-1         +/-7         +/-5         +/-3
9     2     +/-1         {9}          {9}          +/-5
10    3     +/-1         +/-9         +/-7         +/-3
11    2     +/-1         {11}         {11}         +/-5
12    5     +/-1         +/-11        +/-5         +/-7
13    2     +/-1         {13}         {13}         +/-7
14    3     +/-1         +/-13        +/-9         +/-5.             (25)
```

The four danger classes in `(25)` are odd residues modulo `2q`, listed in
the order of the four numerators in `(24)`. At every candidate time, pack
safety follows from

```text
||2r(n/(2q))||=||rn/q||>=1/q>=1/14,                  (26)
```

because `n` is `1` or `k` modulo `q`, both units, and `q` does not divide
`r`.

It remains to prove the exception table. At denominator `2q`, danger means
that the circular residue magnitude `m` satisfies

```text
7m<q.                                                 (27)
```

For even `q`, all four numerators are odd units. The product with an odd
exception is odd, and `(27)` forces `m=1`; modular inversion gives the four
pairwise disjoint `+/-` classes in `(25)`. Each exception can therefore spoil
at most one of the four times. Two exceptions leave at least two winners.

For odd `q`, the numerators `1+q` and `k=2` are even. Their products are
even, so `(27)` forces `m=0`, whose unique odd solution is the class `{q}`.
The other two numerators are odd units and give the displayed inverse
classes. If the first two lifts are both spoiled, one exception lies in
`+/-1` and the other in `{q}`. At the last lift both are outside the displayed
`+/-5` or `+/-7` class, so that lift is safe. If the first bank is not fully
spoiled, it already contains a winner. This proves `(22)--(25)`.

The boundary convention is explicit. At `q=7`, an exception congruent to
`+/-1 mod 14` has distance exactly `1/14` at the singleton time and is safe.
At `q=14`, a pack owner with `rn=+/-1 mod 14` likewise has distance exactly
`1/14`. Replacing the closed safety inequality by a strict one would destroy
these endpoint conclusions.

There is a useful capacity corollary, still only for these explicit banks.
For `2<=q<=7`, the singleton time `(23)` is safe for **any finite number** of
odd exceptions, because no odd residue is dangerous there. For even
`q in {8,10,12,14}`, the four masks in `(25)` are pairwise disjoint, so each
odd exception can spoil at most one time; consequently the same four-time
bank handles any family of at most three odd exceptions. This is an analytic
fixed-bank statement. It is not promoted to a `d=3` LRC residual reduction,
whose inherited runner count and exception typing have not been checked here.

## 5. The typed residual must cover every small divisor

Retain exactly THM-4041's rank-eleven `11+2`, `d=2,c_2=9` boundary. Its two
exceptions are odd, and its eleven divided pack speeds form `H`. If the
corresponding thirteen-speed row is still not `1/14`-lonely, then

```text
for every q in {2,3,...,14} there exists r_q in H with q|r_q.        (28)
```

Equivalently, among the original even pack speeds there is, for every such
`q`, a speed `2r_q` divisible by `2q`. In particular, some divided pack speed
is divisible by `14`, equivalently some original even-pack speed is divisible
by `28`.

This is a proved necessary condition for the actual nonlonely residual, not
only a restatement of a certificate hypothesis. If `(28)` failed for any
one `q`, Section 4 would supply an explicit full-row safe time, contradicting
nonloneliness. Different divisors may of course be owned by the same pack
speed. No sufficiency is claimed for `(28)`.

The irredundant form of `(28)` keeps only the maximal divisibility antichain

```text
{8,9,10,11,12,13,14}.
```

Requiring, for each member of this antichain, some divisible owner in `H` is
equivalent to `(28)`: divisibility by `8` supplies `2,4`, by `9` supplies
`3`, by `10` supplies `5`, by `12` supplies `6`, and by `14` supplies `7`.
This does not require seven different owners; one pack speed may discharge
several or all of the divisibility obligations.

## 6. Infinite affine rays and strict gains

Let `H` be any eleven-element set of distinct positive integers. If `(22)`
holds for at least one `q in {2,...,14}`, then for every coprime odd shape
`a<b` and every positive odd dilation `h`, the row

```text
V_h=2H union {ha,hb}                                  (29)
```

has thirteen distinct positive speeds and is `1/14`-lonely. Thus `(29)` is
an exact infinite `d=2` affine-equality ray closure, with no exception-height
threshold. Section 3 simultaneously gives its exact pullback windows and
the necessary-and-sufficient containment test.

A strict-gain ray is

```text
H_0={1,2,3,4,5,6,7,8,9,11,13},       (a,b)=(1,15).   (30)
```

No member of `H_0` is divisible by `14`, while class `11` violates
THM-4049's old firewall. Hence every odd-dilation row

```text
2H_0 union {h,15h},       h positive and odd,         (31)
```

closes by the new `q=14` bank. At `h=1`, the four new-bank clearances are

```text
(1/28,1/28,1/14,1/14),                                (32)
```

so the last two times close exactly at equality. THM-4049's old four times
instead have clearances

```text
(1/28,1/28,1/56,1/56),                                (33)
```

and are all defeated. Moreover, the five dilations

```text
h in {1,3,5,7,9}                                      (34)
```

lie strictly inside THM-4052's residual width wedge: with `M=13` and
`g=h`,

```text
15h<12M,       15h^2<6M(16h-7h).                     (35)
```

Thus `(34)` is also a strict gain over the proved component-width criterion,
not merely a high-height tail already closed there.

The least-used THM-615 sidecar meets the new theorem at the base AP pack:
`H={1,...,11}` has no owner divisible by `q=12`, so Section 4 recovers its
`1/14` closure for arbitrary odd exceptions. It does not recover THM-615's
sharper exact value `1/12` or replace its argmax arithmetic.

## 7. Strict enlargement of THM-4049 and its physical hostile

At `q=14`, Section 4 uses

```text
T_14={1/28,15/28,3/28,17/28}.                         (36)
```

Modulo `56`, the new forbidden pack set is only

```text
R_new={0,14,28,42}.                                   (37)
```

THM-4049 used

```text
R_old={0,11,14,22,23,28,33,34,42,45}.                (38)
```

Therefore the old allowed set is a strict subset of the new one:

```text
(Z/56Z)\R_old  proper_subset  (Z/56Z)\R_new,
46 classes                         52 classes,
newly allowed={11,22,23,33,34,45}.                  (39)
```

This comparison is literal set containment, not a density heuristic.

It closes THM-4049's typed physical hostile. Its divided pack and exceptions
are

```text
H_*=(2,3,4,5,6,7,8,9,11,2^44,3*2^44),
(alpha,beta)=(1,15).                                  (40)
```

Every member of `H_*` is nonzero modulo `14`, despite the physical class
`11`. The old-bank full-row clearances are

```text
(1/28,1/28,1/56,1/56),                               (41)
```

whereas the new-bank clearances are

```text
(1/28,1/28,1/14,1/14).                               (42)
```

Thus `x=3/28` and `x=17/28` close this actual typed row at the closed
boundary. The hostile was never an LRC counterexample; it remains outside
THM-4049's old hypothesis but is removed from the unresolved residual.

Finally, why do odd `alpha,beta` reduce to the finite audits? At every
`q=14` candidate, exception safety depends only on the odd residue modulo
`28`. There are `14` such classes, hence `14*15/2=105` unordered pairs with
repetition. Equivalently one may retain the redundant odd residues modulo
`56`: there are `28`, hence

```text
28*29/2=406                                             (43)
```

unordered pairs with repetition. Every actual unordered odd pair maps to
one of these classes, and repetition is necessary because distinct positive
exceptions can have the same residue. Adding `28` to an exception changes
`w n/28` by the integer `n`, so the 406-class audit factors exactly to the
105-class mod-28 audit. Nothing is sampled.

## 8. Hostiles, sharp scope, and exact audit

The hypotheses are load-bearing for these fixed banks.

- If `H` contains `q`, then pack speed `2q` is integral at every time in
  `T_q`; for `q=14`, `(28,1,15)` defeats the whole bank.
- With `H={1,...,11}`, even exceptions `(28,56)` defeat all four `q=14`
  times. Oddness cannot be erased. These are fixed-bank hostiles, not LRC
  counterexamples.
- The endpoint `q=14` is sharp for this maximal-residue grid mechanism. At
  `q=15`, take every nonzero residue owner `H_15={1,...,14}`. For any
  `1<=k<15`, if `k` is a unit choose `r=k^(-1) mod 15`, giving distance
  `1/15<1/14`; if not, choose `r=15/gcd(k,15)`, giving distance zero.
  Hence no pack phase `k/15` is safe for this maximal residue pack. This does
  not rule out a different `q=15` certificate or assert a nonlonely row.

The primary companion uses exact `Fraction` arithmetic. It reconstructs all
`364` coprime odd pairs `a<b<60`, of which `362` have positive windows, and
checks `(10)` against every literal danger wall and adjacent cell at `88,056`
points. It independently builds direct wall components for `106` primitive
bases `a<b<32` and verifies `(16)--(17)` at each
`h in {1,3,5,7,9,11}`, for `636` scaled profiles. It then checks all `559`
natural unordered odd-pair profiles across `q=2,...,14`, all `1,264`
unordered odd triples for even `q in {8,10,12,14}`, the strict ray gain, both
old-bank defeats, and the physical equality closure. It executes `107,257`
active requirements.

The independent companion imports none of the primary program. It uses only
integer circular residues, loads all allowed pack residue classes at once
for each `q`, reproduces the complete danger table, and checks the same `559`
pair and `1,264` even-q triple profiles. At `q=14` it additionally loads all
`52` allowed pack classes modulo `56` simultaneously and exhausts all `406`
redundant odd-pair classes. Their safe-time histogram is

```text
2 safe times: 96 pairs;  3 safe times: 232;  4 safe times: 78.       (44)
```

It executes `4,619` active requirements. Both scripts contain no Python
`assert` node and no float literal. Normal and optimized runs reproduce the
frozen raw-LF outputs byte for byte.

## 9. Replay

```text
python3 -B 04-computation/lrc14_d2_mod14_two_bank_affine_ray_firewall_thm4070.py
python3 -B -O 04-computation/lrc14_d2_mod14_two_bank_affine_ray_firewall_thm4070.py
python3 -B 04-computation/lrc14_d2_mod14_two_bank_affine_ray_firewall_thm4070_independent_audit.py
python3 -B -O 04-computation/lrc14_d2_mod14_two_bank_affine_ray_firewall_thm4070_independent_audit.py
```

All four runs reproduce the two frozen outputs. **QED.**
