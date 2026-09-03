# Minimal primitive ternary-unit `l1=18` relation shell for the LRC14 triple comb

**Status: PROVED ELEMENTARY + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**  The
six minimal full-support ternary-unit coefficient sectors of norm eighteen
have complete sharp maxima and unique unordered maximizers.  The primary and
dependency-free clean-room verifiers compare raw-carrier dictionaries with a
literal physical-circle construction on every triple in each rigorous finite
proof window.

This is a theorem about one three-speed scale-three failure comb.  It gives no
seam entry, synchronization, arbitrary-triple bound, all-tail transfer, or
proof of LRC(14).  LRC(14) remains **OPEN**.

## Outcome

Let `w=(w1,w2,w3)` be an unordered primitive triple of distinct positive odd
integers, each prime to three.  Call an integer relation `c dot w=0`
**ternary-unit full-support** when `c` is primitive and every `c_i` is nonzero
modulo three.  This packet studies triples whose least `l1` norm among such
relations is exactly eighteen.

Sorting coefficient magnitudes gives exactly six possible patterns:

```text
116, 1413, 1710, 2511, 477, 558.                       (1)
```

The sharp sector atlas is:

| pattern | sharp measure | unique unordered maximizer | masses at `delta=-3,0,+3` |
|---|---:|---|---|
| `116` | `36/1225` | `{1,11,175}` | `12/1225, 12/1225, 12/1225` |
| `1413` | `12/497` | `{5,7,71}` | `3/497, 6/497, 3/497` |
| `1710` | `12/413` | `{1,7,59}` | `3/413, 6/413, 3/413` |
| `2511` | `318/14399` | `{11,17,121}` | `57/14399, 12/847, 57/14399` |
| `477` | `144/5831` | `{11,17,49}` | `22/5831, 100/5831, 22/5831` |
| `558` | `172/7175` | `{1,25,41}` | `22/7175, 128/7175, 22/7175` |

Consequently the union of all six minimal norm-eighteen sectors has sharp
maximum

```text
36/1225, uniquely at {1,11,175}.                       (2)
```

This identifies the earlier finite-census hostile `{1,11,175}` structurally:

```text
16*11=175+1.
```

It is not nonresonant; its first full-support ternary-unit relation has norm
eighteen.  Its twelve raw carriers split into four components on each of the
three defect sheets, each sheet having mass `12/1225`.

## 1. Why there are six patterns and why the shell is genuinely minimal

For odd speeds, reducing a signed relation modulo two shows that its `l1`
norm is even.  A primitive full-support carrier which can satisfy the owner
gate has three nonzero coefficient magnitudes modulo three.  Enumerating

```text
a<=b<=c,  a+b+c=18,  3 does not divide abc,  gcd(a,b,c)=1
```

gives exactly (1).  The same derivation gives exactly twenty patterns through
norm sixteen:

```text
112,114,118,1110,1114,
125,127,1211,1213,
145,147,1411,
158,1510,178,
255,257,277,455,457.                                   (3)
```

The finite search does not infer minimality from the norm-eighteen label.  It
tests every candidate against every signed coordinate placement of all twenty
patterns in (3), modulo simultaneous sign.  A separate relation-generation
route and a brute-force triple route agree for every pattern through height
79.  Since norm seventeen is impossible by parity, retaining an eighteen-row
candidate precisely when it has no relation from (3) proves that its minimal
ternary-unit full-support norm is eighteen.

Relations with a zero coefficient or a coefficient divisible by three are
deliberately outside this invariant.  They do not provide the owner-forcing
predicate used here.

## 2. Exact defect layers in raw-carrier coordinates

At a physical failure phase put

```text
n_i=nint(w_i y),  e_i=w_i y-n_i,  |e_i|<r,  r=3/14,
o_i=-w_i^(-1)n_i mod 3,  {o_1,o_2,o_3}=F_3.
```

For a primitive relation `c dot w=0` of norm eighteen define

```text
delta=c dot n=-c dot e.
```

Strict eligibility gives

```text
|delta|<18r=27/7<4.                                    (4)
```

Modulo three the three nonzero values `c_i w_i` sum to zero, so all three are
equal.  Distinct owners sum to zero, hence `3|delta`.  Combining this with
(4) yields exactly

```text
delta in {-3,0,3}.                                     (5)
```

The raw carrier

```text
C=w cross n
```

is the canonical component coordinate.  For completeness, choose an integer
Bezout vector `u` with `u dot w=1`.  Every `C dot w=0` has the explicit lift
`n=C cross u`, proving

```text
Z^3/Zw  isomorphic to  Lambda_w={C in Z^3:C dot w=0}.   (6)
```

Now choose `v in Z^3` with `c dot v=1`, and set

```text
C_delta=w cross (delta v).
```

Every carrier in the defect-`delta` fibre is uniquely

```text
C_delta+k c,             k in Z.                       (7)
```

Indeed, differences of two lifts with equal defect lie with `w` in the plane
`c^perp`, so their cross product is parallel to primitive `c`; conversely
(6) lifts every integer multiple of `c`.  A different Bezout choice merely
translates `k`, leaving the affine carrier line unchanged.

For any raw carrier put

```text
L_w(C)=max(0,min(
  2r/w1,2r/w2,2r/w3,
  r/w1+r/w2-|C3|/(w1w2),
  r/w1+r/w3-|C2|/(w1w3),
  r/w2+r/w3-|C1|/(w2w3))).                             (8)
```

This is exactly the length of the intersection of the three nearest-integer
intervals.  Therefore

```text
mu(F_w)=sum_{delta=-3,0,3}
        sum_{k in Z, all coordinates of C_delta+kc nonzero mod 3}
          L_w(C_delta+kc).                             (9)
```

Formula (9) requires no endpoint coprimality and silently loses no gcd
torsion.  It is the raw-carrier replacement for presentation-specific
endpoint determinants.

Because `delta=0 mod 3`, reduction of (7) modulo three lies on the line
spanned by `c`.  Exactly one value of `k mod 3` makes all three carrier
coordinates zero; the other two make all three nonzero.  Thus each defect
roof is sampled on exactly two residue classes modulo three.

## 3. Slice integrals and the rigorous all-height cutoff

For real `t` define the defect roof

```text
f_delta(t)=L_w(C_delta+t c).
```

Every superlevel set is an intersection of intervals in `t`, so the roof is
nonnegative, compactly supported, and unimodal.  Choose an integral `h` with
`w cross h=c`.  The map

```text
(y,t) -> e=w y-delta v-t h
```

has area Jacobian `||c||` and maps the roof region to

```text
[-r,r]^3 intersect {c dot e=-delta}.
```

Hence

```text
I_delta(c)=integral_R f_delta(t)dt
 =area([-r,r]^3 intersect {c dot e=-delta})/||c||.      (10)
```

After coordinate sign changes this depends only on the magnitude pattern
`(a,b,c)`.  Three-box convolution gives the exact formula

```text
I_delta(a,b,c)=1/(2abc) sum_(S subset {1,2,3}) (-1)^|S|
 (|delta|+r(a+b+c)-2r sum_(i in S)a_i)_+^2.            (11)
```

For a nonnegative unimodal function sampled on one class `a+3Z`, elementary
left/right rectangle comparison gives

```text
sum_j f(a+3j)<=integral(f)/3+sup(f).                    (12)
```

There are two allowed residue classes on each of three defect sheets, and
`sup f_delta<=2r/W=3/(7W)` for `W=max(w)`.  Therefore

```text
mu(F_w)<=B_c+18/(7W),
B_c=(2/3)(I_0(c)+2I_3(c)).                             (13)
```

The exact constants and thresholds obtained from (11)--(13) are:

| pattern | `I_0` | `I_3` | `B_c` | equality threshold `18/[7(M-B_c)]` | proof height |
|---|---:|---:|---:|---:|---:|
| `116` | `9/784` | `9/784` | `9/392` | `400` | `400` |
| `1413` | `9/637` | `27/5096` | `3/182` | `3692/11` | `335` |
| `1710` | `9/490` | `27/6860` | `6/343` | `2891/13` | `222` |
| `2511` | `9/539` | `9/2695` | `6/385` | `10285/26` | `395` |
| `477` | `54/2401` | `9/4802` | `6/343` | `357` | `357` |
| `558` | `27/1225` | `9/4900` | `3/175` | `18450/49` | `376` |

For every integer height strictly above the displayed proof height, (13) is
strictly smaller than that row's claimed maximum.  Thus the exact finite
search proves a global result in each sector rather than merely reporting a
discovery box.  The verifier also integrates the six-piece affine roof
directly at every maximizing presentation and recovers (11) independently.

## 4. Exact finite proof windows

All counts below are for sorted primitive distinct positive odd three-unit
speed triples.  `all` means triples having the displayed norm-eighteen
pattern before the shorter-relation filter; `minimal` means after that
filter; `rays` counts signed primitive relation directions modulo simultaneous
negation.  Physical components are counted once per speed triple in that row,
not once per presentation.

| pattern | height | all | minimal | rays | positive | grouped carrier components |
|---|---:|---:|---:|---:|---:|---:|
| `116` | 400 | 501 | 497 | 497 | 497 | 10,730 |
| `1413` | 335 | 876 | 866 | 867 | 864 | 15,456 |
| `1710` | 222 | 495 | 485 | 486 | 485 | 6,962 |
| `2511` | 395 | 1,439 | 1,429 | 1,430 | 1,428 | 32,904 |
| `477` | 357 | 789 | 785 | 785 | 784 | 15,756 |
| `558` | 376 | 855 | 851 | 851 | 851 | 17,668 |

The row-wise proof windows contain `4,913` minimal sector triples, `4,916`
relation rays, and `99,476` grouped raw-carrier/`y` components.  On every
relation ray the raw-carrier dictionary from (9) agrees component by component
with a separate three-list intersection of the literal physical intervals on
`[0,1]`.  The independent referee also constructs the six shifted `x`-sheet
combs on all six winners and all ten overlaps.  Each grouped carrier has
exactly three lifted `x` arcs, so the corresponding count is `298,428`;
grouping restores the identical carrier dictionary and mass.

The empty triples in the proof windows are exactly

```text
1413: {7,11,47}, {7,25,29};
2511: {7,11,43};
477:  {7,25,29}.
```

The repeated `{7,25,29}` is one physical triple with two sector labels.

## 5. Presentation incidence is not a physical-triple count

At the common finite sidecar height `400`, the six patterns have:

| pattern | all triples | minimal triples | minimal relation rays |
|---|---:|---:|---:|
| `116` | 501 | 497 | 497 |
| `1413` | 1,237 | 1,227 | 1,228 |
| `1710` | 1,612 | 1,602 | 1,603 |
| `2511` | 1,463 | 1,453 | 1,454 |
| `477` | 981 | 977 | 977 |
| `558` | 972 | 968 | 968 |

Thus there are three different exact counts:

```text
sector-pattern incidences:                  6,724
signed relation rays modulo overall sign:   6,727
distinct underlying speed triples:          6,717
positive physical combs:                    6,714.      (14)
```

The difference between the first and third counts is caused by exactly seven
cross-pattern overlaps:

| speed triple | patterns | physical measure |
|---|---|---:|
| `{1,19,35}` | `116,477` | `64/4655` |
| `{1,25,41}` | `116,558` | `172/7175` |
| `{1,37,53}` | `116,1710` | `300/13727` |
| `{5,13,41}` | `1413,1710` | `22/3731` |
| `{7,25,29}` | `1413,477` | `0` |
| `{11,29,43}` | `1413,1710` | `92/8729` |
| `{13,19,47}` | `1413,2511` | `14/893` |

The additional three relation rays come from within-pattern double
presentations:

```text
{1,17,55},  pattern 1413: (1,-13,4), (13,-4,1), mu=2/385;
{13,23,31}, pattern 1710: (1,-10,7), (10,-7,1), mu=22/4991;
{1,17,37},  pattern 2511: (2,-11,5), (11,-5,2), mu=8/4403.
```

An exact coefficient cross-product audit proves these ten speed rays are the
complete multiple-relation list, not merely the ones below height 400.

Nine of the ten multiply-related triples have nonempty combs.  This does not
contradict the zero-defect overlap obstruction: on every surviving component,
the two independent relations do not both have defect zero.  At least one
relation sees `delta=+/-3`.  The norm-eighteen shell therefore gives a sharp
hostile to extending the earlier pointwise zero-defect lemma by erasing the
defect sidecar.

The raw carrier `C` is chart-independent.  The pair `(delta,k)` in (7) is a
useful local natural-number coordinate only after naming the relation chart;
different presentations relabel the same physical carrier.  This is why raw
carrier counts, relation-ray counts, sector incidences, and physical triples
must remain separate.

## Reproduction and frozen evidence

The primary verifier imports only Python's standard library and uses exact
`Fraction` arithmetic.  All `923,085` checks are explicit and remain active
under optimization.  A dependency-free clean-room referee performs 189,085
further optimization-live checks, including a second cube-slice construction
by rational polygon clipping and the six-sheet normalization audit.

```powershell
python -B 04-computation/lrc14_minimal_ternary_unit_norm18_shell_thm4393.py
python -B -O 04-computation/lrc14_minimal_ternary_unit_norm18_shell_thm4393.py
$env:PYTHONHASHSEED='18'
python -B -O 04-computation/lrc14_minimal_ternary_unit_norm18_shell_thm4393.py
python -B 04-computation/lrc14_minimal_ternary_unit_norm18_shell_independent_referee_thm4393.py
python -B -O 04-computation/lrc14_minimal_ternary_unit_norm18_shell_independent_referee_thm4393.py
```

Normal, optimized, and fixed-hash-seed streams agree exactly.  Frozen raw-LF
hashes are:

```text
primary script       0d5fb1c0e65d5e8a048e01155da635e7cf5472c6d4e7bacbebd4c91108cb0d17
primary output       9176e308f9db189bc150b465681df246c22dbb77b6c8c24ce7e5d84589e0ee73
independent script   bd11c07be205126e90dd278a0b6b6dd82f01a1893324baf2b766785ca750f15c
independent output   481c2b5e9efb881395f027b19a5edb1e50d6867e8f74f38e0b3d3a2cf9b5fa60
```

Normal, optimized, and fixed-hash-seed replays match the two frozen canonical
outputs.  Discovery scripts are not part of the proof surface.

## Honest next frontier

This shell closes a natural structured family, not arbitrary triples.  The
best next uses of the result are:

1. carry (9)--(13) to the minimal norm-twenty shell, where the same three
   defect states persist because `20r=30/7<5`;
2. replace the coarse six-residue `18/(7W)` error by exact residue-roof
   discrepancy to shrink later proof windows;
3. express relation-chart changes on the ten multiply-presented rays as
   integral affine transformations of `(delta,k)`, with raw `C` as the common
   invariant; and
4. attack the genuinely relation-free rank-two carrier lattice rather than
   interpreting completion of one more coefficient shell as seam entry.

No claim beyond these six relation sectors is made.
