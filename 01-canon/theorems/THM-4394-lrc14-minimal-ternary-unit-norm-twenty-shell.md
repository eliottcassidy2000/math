# Minimal primitive ternary-unit `l1=20` relation shell for the LRC14 triple comb

**Status: PROVED ELEMENTARY RELATIVE TO THM-4393 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.**  The seven minimal full-support ternary-unit
coefficient sectors of norm twenty have complete sharp maxima and unique
unordered maximizers.  The primary verifier reuses the audited carrier
primitives of THM-4393; a standalone clean-room verifier independently
reconstructs every theorem mechanism and proof window.

This is a theorem about one three-speed scale-three failure comb.  It gives no
seam entry, synchronization, arbitrary-triple bound, all-tail transfer, or
proof of LRC(14).  **LRC(14) remains OPEN.**

## Outcome

Let `w=(w1,w2,w3)` be an unordered primitive triple of distinct positive odd
integers, each prime to three.  Call `c dot w=0` **ternary-unit full-support**
when `c` is primitive and every `c_i` is nonzero modulo three.  Suppose the
least `l1` norm among such relations is twenty.  Sorting coefficient
magnitudes gives exactly seven patterns:

```text
(1,2,17), (1,5,14), (1,8,11), (2,5,13),
(2,7,11), (4,5,11), (5,7,8).                           (1)
```

Their sharp sector atlas is:

| pattern | sharp measure | unique unordered maximizer | masses at `delta=-3,0,+3` |
|---|---:|---|---|
| `(1,2,17)` | `36/1295` | `{1,11,185}` | `12/1295, 12/1295, 12/1295` |
| `(1,5,14)` | `24/1043` | `{1,11,149}` | `6/1043, 12/1043, 6/1043` |
| `(1,8,11)` | `412/16093` | `{11,19,121}` | `92/16093, 12/847, 92/16093` |
| `(2,5,13)` | `166/6853` | `{7,11,89}` | `50/6853, 6/623, 50/6853` |
| `(2,7,11)` | `608/24563` | `{11,29,121}` | `130/24563, 12/847, 130/24563` |
| `(4,5,11)` | `894/41503` | `{11,49,121}` | `153/41503, 12/847, 153/41503` |
| `(5,7,8)` | `216/8855` | `{13,23,55}` | `8/1771, 136/8855, 8/1771` |

Consequently the union of the seven declared minimal sectors has sharp
maximum

```text
36/1295, uniquely at {1,11,185}.                        (2)
```

The winner is governed by the relation

```text
17*11=185+2*1.
```

Its twelve raw carriers split evenly over the three defect sheets.

## 1. Complete coefficient shell and genuine minimality

For odd speeds, reducing `c dot w=0` modulo two shows that `||c||_1` is even.
Full support at the owner gate requires all three coefficient magnitudes to
be nonzero modulo three.  Enumerating

```text
a<=b<=c,  a+b+c=20,  gcd(a,b,c)=1,  3 does not divide abc
```

gives precisely (1).  The same defining enumeration gives twenty-six
patterns of even norm at most eighteen:

```text
112,114,125,118,127,145,1110,147,255,
1211,158,257,455,1114,1213,1411,1510,178,277,457,
116,1413,1710,2511,477,558.                             (3)
```

The verifier generates every signed coordinate placement, modulo simultaneous
negation, of all twenty-six patterns in (3).  It retains a norm-twenty
presentation only if its speed triple has none of those shorter relations.
Odd norm nineteen is impossible by parity.  Hence “minimal” here is proved
from the definition, not inferred from a shell label.

Relations with a zero coefficient or a coefficient divisible by three are
outside this invariant.  They may exist on a retained triple and must not be
silently treated as owner-forcing ternary-unit relations.

## 2. Three defect sheets and the affine carrier address

At a physical failure phase put

```text
n_i=nint(w_i y),  e_i=w_i y-n_i,  |e_i|<r,  r=3/14,
o_i=-w_i^(-1)n_i mod 3,  {o_1,o_2,o_3}=F_3.
```

For a norm-twenty relation define `delta=c dot n=-c dot e`.  Strictness gives

```text
|delta|<20r=30/7<5.                                     (4)
```

The three nonzero residues `c_i w_i` sum to zero, so they are equal.  Distinct
owners sum to zero; therefore `3|delta`.  Combining this with (4) yields

```text
delta in {-3,0,3}.                                      (5)
```

The chart-independent component address is `C=w cross n` in

```text
Lambda_w={C in Z^3:C dot w=0}.
```

For primitive `w`, the map `[n] -> w cross n` is the integral isomorphism
`Z^3/Zw -> Lambda_w` proved in THM-4386.  Choose `v in Z^3` with `c dot v=1`
and put `C_delta=w cross (delta v)`.  Every carrier in one defect fibre is
uniquely

```text
C=C_delta+k c,       k in Z.                            (6)
```

A different Bezout section only translates `k`.  Since `delta=0 mod3`, the
line (6) has one all-zero carrier class and two all-nonzero carrier classes
modulo three.  The latter two are exactly the distinct-owner gate.

For any carrier define

```text
L_w(C)=max(0,min(
  2r/w1,2r/w2,2r/w3,
  r/w1+r/w2-|C3|/(w1w2),
  r/w1+r/w3-|C2|/(w1w3),
  r/w2+r/w3-|C1|/(w2w3))).                              (7)
```

This is the exact common length of the three nearest-integer intervals.
Thus

```text
mu(F_w)=sum_{delta=-3,0,3}
        sum_{k in Z, all (C_delta+kc)_i nonzero mod3}
          L_w(C_delta+kc).                              (8)
```

Positivity in (7) supplies a finite strict interval of integers `k` for each
sheet.

## 3. Cube slices and finite global cutoffs

For real `t` set `f_delta(t)=L_w(C_delta+t c)`.  Its superlevel sets are
intervals, so it is nonnegative, compactly supported, and unimodal.  The
affine error-cube parametrization gives

```text
I_delta(c)=integral_R f_delta(t)dt
 =area([-r,r]^3 intersect {c dot e=-delta})/||c||.       (9)
```

For a magnitude pattern `(a,b,c)`, three-box convolution yields

```text
I_delta(a,b,c)=1/(2abc) sum_(S subset {1,2,3}) (-1)^|S|
 (|delta|+r(a+b+c)-2r sum_(i in S)a_i)_+^2.             (10)
```

For a nonnegative unimodal function sampled on one translate of `3Z`,
left/right rectangle comparison gives

```text
sum_j f(a+3j)<=integral(f)/3+sup(f).
```

There are two live residue classes on each of the three sheets, while
`sup(f_delta)<=3/(7W)` for `W=max(w)`.  Hence

```text
mu(F_w)<=B_c+18/(7W),
B_c=(2/3)(I_0(c)+2I_3(c)).                              (11)
```

The exact constants and proof windows are:

| pattern | `I_0` | `I_3` | `B_c` | equality threshold | proof height |
|---|---:|---:|---:|---:|---:|
| `(1,2,17)` | `9/833` | `9/833` | `18/833` | `22015/53` | `415` |
| `(1,5,14)` | `9/686` | `9/1372` | `6/343` | `21903/47` | `466` |
| `(1,8,11)` | `9/539` | `45/8624` | `39/2156` | `1158696/3385` | `342` |
| `(2,5,13)` | `9/637` | `18/3185` | `54/3185` | `4009005/11332` | `353` |
| `(2,7,11)` | `9/539` | `18/3773` | `6/343` | `1547469/4369` | `354` |
| `(4,5,11)` | `9/539` | `81/21560` | `87/5390` | `118580/249` | `476` |
| `(5,7,8)` | `279/13720` | `81/27440` | `6/343` | `185955/499` | `372` |

The threshold is `18/[7(M-B_c)]`, where `M` is that row's sharp value.
For every integer height above the displayed proof height, (11) is strictly
below `M`; the finite searches therefore prove global sector maxima.

## 4. Exact proof-window census

Counts are for sorted primitive distinct positive odd three-unit speed
triples.  `all` is incidence before the shorter filter; `minimal` is after;
`rays` counts signed primitive relation directions modulo overall sign.

| pattern | height | all | minimal | rays | positive | grouped carriers |
|---|---:|---:|---:|---:|---:|---:|
| `(1,2,17)` | 415 | 1,029 | 1,012 | 1,013 | 1,012 | 21,652 |
| `(1,5,14)` | 466 | 1,567 | 1,549 | 1,550 | 1,549 | 37,274 |
| `(1,8,11)` | 342 | 1,071 | 1,056 | 1,057 | 1,056 | 21,946 |
| `(2,5,13)` | 353 | 977 | 964 | 965 | 963 | 18,460 |
| `(2,7,11)` | 354 | 1,160 | 1,147 | 1,147 | 1,147 | 24,512 |
| `(4,5,11)` | 476 | 2,092 | 2,079 | 2,080 | 2,078 | 55,166 |
| `(5,7,8)` | 372 | 1,549 | 1,538 | 1,538 | 1,538 | 32,424 |

The empty row in the two affected proof windows is the same physical triple
`{1,19,41}`, carrying patterns `(2,5,13)` and `(4,5,11)`.  On every relation
ray, the affine dictionary from (8) agrees carrier by carrier with a separate
literal three-list intersection on the physical circle.

## 5. Incidence, relation rays, and physical triples

At common sidecar height `500`, the exact totals are

```text
sector-pattern incidences:                  14,926
signed relation rays modulo overall sign:   14,931
distinct underlying speed triples:          14,918
positive physical combs:                    14,917.     (12)
```

There are exactly eight cross-pattern overlaps:

| speed triple | patterns | physical measure |
|---|---|---:|
| `{1,19,41}` | `(2,5,13)`, `(4,5,11)` | `0` |
| `{5,13,49}` | `(1,2,17)`, `(1,8,11)` | `32/4459` |
| `{7,11,65}` | `(1,2,17)`, `(1,8,11)` | `6/455` |
| `{11,29,41}` | `(2,7,11)`, `(5,7,8)` | `138/8323` |
| `{11,29,79}` | `(1,2,17)`, `(1,5,14)` | `360/16037` |
| `{13,23,41}` | `(1,5,14)`, `(1,8,11)` | `38/6601` |
| `{13,31,37}` | `(2,7,11)`, `(5,7,8)` | `142/8029` |
| `{17,37,53}` | `(1,5,14)`, `(1,8,11)` | `90/6307` |

Five further speed rays have two presentations within one pattern:

| speed triple | pattern | physical measure |
|---|---|---:|
| `{7,11,97}` | `(1,2,17)` | `18/679` |
| `{13,23,47}` | `(4,5,11)` | `102/7567` |
| `{13,23,67}` | `(1,5,14)` | `282/20033` |
| `{17,23,53}` | `(2,5,13)` | `2592/145061` |
| `{25,29,43}` | `(1,8,11)` | `120/8729` |

Pairing every two labelled coefficient rays, taking their cross product, and
then imposing positivity, parity, three-unit, distinctness, primitivity, and
the complete shorter filter proves this thirteen-ray list at every height.
It is not an observation truncated at `500`.  Twelve of the thirteen combs
are nonempty; `{1,19,41}` is the unique empty one.

The raw carrier `C` remains the intrinsic component address.  Relation-sector
labels and `(delta,k)` charts count presentations, not new physical pieces,
so overlapping sector values must never be added.

## Reproduction and evidence

The primary verifier uses exact `Fraction` arithmetic, imports only the
audited THM-4393 carrier implementation, and performs 20 frozen theorem gates
plus 2,096,197 inherited optimization-live arithmetic checks.  The standalone
clean-room referee performs 3,846,013 further checks, including a second
coefficient-first/triple-first enumeration, all carrier dictionaries, and all
7,875 pairs of the 126 labelled norm-twenty coefficient rays.

```powershell
python -B 04-computation/lrc14_minimal_ternary_unit_norm20_shell_thm4394.py
python -O -B 04-computation/lrc14_minimal_ternary_unit_norm20_shell_thm4394.py
$env:PYTHONHASHSEED='20'
python -B 04-computation/lrc14_minimal_ternary_unit_norm20_shell_thm4394.py
python -B 04-computation/lrc14_minimal_ternary_unit_norm20_shell_independent_referee_thm4394.py
python -O -B 04-computation/lrc14_minimal_ternary_unit_norm20_shell_independent_referee_thm4394.py
```

Normal, optimized, and fixed-hash-seed streams agree with their frozen
canonical outputs:

```text
primary script       a9336358e1a9d0b220aa52e7e7186b095839577d4fe70166e25d01576984dd6e
primary output       1c999e4153d21e3903c681dc1296dccea4d13ffbafa23e9fcbaf5d0676637e90
independent script   e005c4ddbe95fe686f98a49d121c24456d3f6fa892b3cecd369c41cc42b98387
independent output   0dd6042305d4b8ce67d7329a6277fa31f26d77e80664c73635eb68e80e19edf5
```

## Honest next frontier

The norm-twenty shell is another structured family, not a universal
nonresonant theorem.  The more consequential next split is by coefficient
residue type: a primitive relation on a three-unit speed triple either has
all three coefficients nonzero modulo three, as above, or exactly one
coefficient zero modulo three.  In the latter case the owner permutation
forces a nonzero defect class and retains an orientation bit.  That route,
and a signed finite-data Fourier tail bound, may turn shell enumeration into
a stopping theorem.  Neither is asserted here.
