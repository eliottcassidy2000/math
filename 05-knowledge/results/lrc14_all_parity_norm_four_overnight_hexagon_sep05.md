# All-parity norm four and the completed parity-dependent local ceiling

**Status: PROVED ELEMENTARY + FINITE-EXACT; independently audited.**
The signed coefficient-magnitude `(1,1,2)` family has exact equality between
its best complete network projection and actual spoiled Haar mass. Removing
oddness changes the sharp constant to `11/140`, uniquely at `(2,11,20)`.
Every other primitive member is at most `6/77`, with equality only at
`(1,5,11)`.

**Concurrent credit.** Incoming
[THM-4444 — signed-112 sharp one-ray classification](../../01-canon/theorems/THM-4444-lrc14-signed-112-sharp-one-ray-classification.md)
is now **PROVED** at `fd73e66105`, with the same sharp statement and a
finer three-cone classification. This note is an independent co-discovery,
proof/referee sidecar, not a competing theorem-ID claim. Its distinct
mechanism is the exact endpoint-product discrepancy and forty-row head,
including the ten endpoint-gcd-two rows. The incoming
[native synthesis](overnight4_20260906_lrc_parityfree_native.md) independently
closes the universal local ceiling and the complete old-threshold locus;
the parity table below is a corollary, not a second global closure claim.

Combining this family with the incoming **PROVED** generic reduction,
additive theorem and norm-five theorem gives the all-parity ceiling `6/55`
and the three sharp primitive parity constants below. No arbitrary-body
Haar floor, entry map, synchronization, or LRC(14) closure is asserted.

## 1. Inheritance and the lost endpoint gcd

The closest proved mechanism is
[THM-4373 — signed-121 resonant triple comb](../../01-canon/theorems/THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound.md),
especially its exact interval profile and period-three quadrature. Its
original finite endpoint head is odd-typed. The incoming
[THM-4437 — all-parity low-circuit reduction](../../01-canon/theorems/THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits.md)
keeps this norm-four family as an explicit residual.

The canonical hostile `(2,11,20)` comes from the
[incoming parity probe](lrc14_parity_empty_core_sep06.md), which correctly
records actual mass `11/140>6/77`. The corrected near miss is to keep the
old endpoint gcd-one universe after removing parity. Here the coefficient-one
endpoints are `2,20`, whose gcd is **two**, while the full triple is
primitive. That omitted finite case changes the all-height answer.

The least-used sidecar is the complete three-coordinate integer carrier,
not a pair endpoint presentation alone. The concept board is: full-row
primitivity; endpoint gcd; redundant addressed tooth; physical/network
equality; period-three quadrature; parity-specific consumer. The map to
two endpoints preserves interval lengths only when its full primitive
carrier multiplicity and coefficient-two address survive.

## 2. Exhaustive endpoint representation for every parity

Let `w` be sorted, positive, distinct, primitive, and a ternary-unit triple
admitting a signed `(1,1,2)` relation. Let `p<q` be the speeds attached to
the magnitude-one coefficients and let `r` be the magnitude-two speed.
Positivity gives exactly the two possibilities

```text
r=(q+p)/2   or   r=(q-p)/2.                            (1)
```

Thus `p,q` have the same parity and `q` is the largest speed. Keep the exact
conditions `gcd(p,q,r)=1`, distinctness and `3∤pqr`.
Because `gcd(p,q)` divides `2r`, full primitivity implies

```text
gcd(p,q) in {1,2}.                                    (2)
```

The first case has odd endpoints. The second has even endpoints and odd
`r`; it is not removed. Formula (1) plus these full-triple filters exhausts
the family. No `gcd(p,q)=1` restriction is added.

Choose the signed relation vector `d=(1,-2,1)` for the mean branch or
`d=(1,2,-1)` for the difference branch, in the unsorted order `(p,r,q)`.
The parity-free raw-carrier identity from
[THM-4434 — universal scale-three network projection bound](../../01-canon/theorems/THM-4434-lrc14-universal-scale-three-network-projection-bound.md)
supplies `C=w cross e`, `|e_i|<3/14`. Since `d dot w=0`,
`d cross C=delta w`, where primitivity of the **full** triple makes
`delta` integral. But

```text
|delta|=|d dot e|<12/14<1,
```

so `delta=0`. The full live dictionary is therefore

```text
Omega={k d:0<|k|<B=3(p+q)/28, 3∤k}.                  (3)
```

The raw roof at coordinate `r` is tight; in the difference branch one other
roof ties it. Every multiplier in (3) is admitted, by the full primitive
three-speed carrier identity. In particular no incorrect pair-gcd-one
argument is used to count the even-endpoint case.

## 3. One projection equals physical overlap exactly

On every addressed failure component, the corresponding integer lifts
satisfy the same relation as the speeds, because the effective error bound
above forces the integer relation defect to vanish. Thus

```text
e_r=(e_p+e_q)/2       in the mean branch,
e_r=(e_q-e_p)/2       in the difference branch.         (4)
```

Whenever both endpoint errors lie in `(-3/14,3/14)`, (4) puts `e_r` in
that same interval. Hence the entire addressed pair intersection lies in
the coefficient-two tooth. That third tooth is redundant on this component,
so its contact capacity equals the actual overlap length. Consequently

```text
E_r=mu(F_w)=min_i E_i.                                (5)
```

This is an equality of the actual complete network and physical statistic,
not an estimate obtained by ignoring the third tail. The integer address
in (4) is essential, especially when both endpoint speeds are even.

For carrier `k d`, the exact endpoint interval length is

```text
lambda_k=min(3/(7q), 3/(14p)+3/(14q)-2|k|/(pq)).       (6)
```

Set `A=3(q-p)/28`, `B=3(q+p)/28` and extend by zero outside support. Then

```text
lambda(t)=2/(pq) [(B-t)_+-(A-t)_+],
mu(F_w)=2 sum_(k>=1,3∤k) lambda(k).                   (7)
```

These are the inherited THM-4373 formulas, now justified from the full
primitive dictionary for both endpoint gcd types. Integral `B` is a deleted
multiple of three; the hinge itself also vanishes there.

## 4. Exact quadrature and the expanded forty-row head

Retain THM-4373's period-three primitive `P`, with `-1/3<=P<=0`. Formula
(7) gives

```text
mu(F_w)=3/49+4/(pq)[P(B)-P(A)]
       <=3/49+4/(3pq).                                (8)
```

Therefore `pq>=81` implies strict inequality below `6/77`: the upper
bound at 81 is `925/11907=6/77-31/130977`.

The complete remaining endpoint universe is

```text
1<=p<q, pq<81, p=q mod2,
r=(q-p)/2 or (q+p)/2,
sort(p,r,q) distinct, primitive, and ternary-unit.     (9)
```

Since `p<sqrt(81)`, it suffices to run `p=1,...,8` and the explicitly
bounded `q` interval. There are exactly 40 distinct triples: 30 with
endpoint gcd one and 10 with endpoint gcd two. By number of even
coordinates, the counts are 13 all-odd, 17 one-even and 10 two-even.
The original odd head was only the first group.

Every complete raw dictionary and all three projections agree with literal
six-sheet interval calculations. The only rows in this full head with
physical mass at least `6/77` are

```text
(1,5,11):  E=(6/77,6/77,6/77),             H=6/77;
(2,11,20): E=(131/1540,11/140,3/35),       H=11/140.
```

All other head rows are strict below `6/77`. The strict tail (8) proves
the full all-height claims:

```text
min_i E_i=H<=11/140, equality iff w=(2,11,20);
min_i E_i>6/77 iff w=(2,11,20);
min_i E_i=6/77 iff w=(1,5,11).                         (10)
```

The last two equivalences are within the norm-four family. No all-coordinate
`11/140` bound is claimed; the extremizer itself has larger other projections.

## 5. Synthesis with proved incoming low-circuit work

The complete generic reduction is
[THM-4437](../../01-canon/theorems/THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits.md):
outside the three low circuits, every best projection is strictly below
`6/77`. The additive circuit has the independently audited
[sharp `6/55` theorem](lrc14_additive_parity_empty_core_sep06.md), uniquely
at `(1,10,11)`. The norm-five circuit is now
[THM-4441 — signed-122 sharp ray closure](../../01-canon/theorems/THM-4441-lrc14-signed-122-sharp-ray-closure.md),
**PROVED** at `058a8ded98`, with `min E<=46/665<6/77`; the
[independent norm-five proof](lrc14_norm_five_overnight_hexagon_sep05.md)
is an alternative audited dependency. No reserved theorem is used.

Together with (10), these cover every primitive distinct positive
ternary-unit triple. Both `min_i E_i` and actual mass have the following
sharp maxima, with unique primitive equality in each row:

| Number of even coordinates | Sharp maximum of each statistic | Equality triple |
| --- | ---: | --- |
| 0 | `6/77` | `(1,5,11)` |
| 1 | `6/55` | `(1,10,11)` |
| 2 | `11/140` | `(2,11,20)` |

Indeed a primitive additive triple has exactly one even coordinate. Outside
that family only norm four can exceed `6/77`, and (10) identifies its single
two-even exception. The odd norm-four equality is `(1,5,11)`; norm five
and all generic minima are strict. The equality packets above and the
incoming additive packet have physical mass equal to their selected
projections, proving sharpness for both targets.

In particular the universal all-parity sharp ceiling is `6/55`, uniquely
at `(1,10,11)`. Among **nonadditive** triples the only value above `6/77`
is `(2,11,20)`, and the only equality is `(1,5,11)`. These statements are
about minima and physical mass, not about every projection coordinate.

Common ternary-unit dilation preserves physical Haar mass and each complete
network projection: every contact interval has the corresponding number of
pullback copies, each with its length and capacity divided by the dilation.
Reduce by the common gcd before reading the parity table or equality classification.
A physical Haar consumer with a body floor at the corresponding threshold
is valid, but no such arbitrary-body floor is supplied. The recorded
ten-body floor counterexamples remain unchanged. This completes a local
triple statistic, not LRC(14).

## Reproduction and audit

```bash
python3 -B 04-computation/lrc14_all_parity_norm_four_overnight_hexagon_sep05.py
python3 -B -O 04-computation/lrc14_all_parity_norm_four_overnight_hexagon_sep05.py
```

The [source](../../04-computation/lrc14_all_parity_norm_four_overnight_hexagon_sep05.py)
retains the full-triple gcd and both endpoint branches before filtering.
It compares complete raw carriers, all literal interval projections,
physical overlap, coefficient-two projection equality and the inherited
quadrature on every row of (9). Large-product controls keep both endpoint
gcd types, and literal common dilations are tested separately. The
[output](lrc14_all_parity_norm_four_overnight_hexagon_sep05.out) records the
finite universe and equality loci. The independent `observer_collision`
audit passed the full proof and head replay, including the endpoint gcd-two
sidecar, integer defect, addressed-tooth redundancy, exact hinge normalization,
strict tail, and the parity/equality synthesis. The referee also read the
actual incoming THM-4441 YAML at `058a8ded98`, verifying its **PROVED** status
and the exact dependency used here. The root's separate full written-proof
review also passed. No scarce ID or navigation is edited here.

Normal and optimized executions have identical output: 471 primary checks
and 895 independent raw/literal checks. Frozen source SHA256:
`6ba0b055f01c2281e8d26b31ccdf3c96773bf2ac4c897e063eb73ffce29716f1`;
output SHA256:
`dc82c9a184a1fc97ef07f999f3f940f8ab0b526d000d848674cc8756aae03c8d`.
