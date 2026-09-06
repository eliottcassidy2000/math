# Independent audit of the proposed THM-4437 all-parity reduction

## Verdict

**FAIL AS STATED; PASS AFTER A SHARP ONE-SYMBOL REPAIR.**  The proposed
conclusion that every non-low raw projection is *strictly* below `6/77` is
false.  There are exactly three primitive non-low projection equalities:

```text
w=(7,16,22):  E=(15/616, 6/77, 479/8624), equality in E_b;
w=(14,17,22): E=(3/119, 3/49, 6/77),       equality in E_c;
w=(4,19,22):  E=(39/1463, 23/308, 6/77),  equality in E_c.
```

None has a signed relation with coefficient magnitudes `(1,1,1)`,
`(1,1,2)`, or `(1,2,2)`.  Each instead has a primitive `(1,2,4)` carrier
ray, and its live carriers are the first two positive multiples and their
negatives.  Thus the finite intercept, not the asymptotic slope, creates the
boundary equality.

The following repaired all-height statement passes every analytic and exact
gate:

> If a primitive sorted distinct positive ternary-unit triple has no signed
> relation with magnitudes `(1,1,1)`, `(1,1,2)`, or `(1,2,2)`, then every
> raw projection is **at most** `6/77` and their minimum is **strictly below**
> `6/77`.  Individual equality occurs exactly at the three rows/coordinates
> displayed above.

This is enough for the advertised reduction of the minimum-projection target
to the three low circuit patterns.  It is not enough for the stronger
every-projection strictness written in the proposed statement.

## 1. Parity-free modulo-three dichotomy

Let every `w_i` be a unit in `F_3`.  A live residue `C` has all coordinates
nonzero and `C.w=0`.  The three products `C_i w_i` are each `+/-1`; three such
terms sum to zero in `F_3` exactly when all three agree.  Hence

```text
C_i = lambda / w_i, lambda in F_3^*.
```

For a primitive relation `v.w=0`, two zero coordinates modulo three would
force the third to vanish, which would contradict primitivity.  Thus `v` has
either zero or one zero residue.

- If all `v_i` are units, the same argument gives
  `v_i=mu/w_i`.  Therefore `v` and `C` are parallel modulo three,
  `v cross C=delta w` has `delta=0 mod 3`, and exactly two of the three
  longitudinal classes `C+t v` remain live.
- If exactly one `v_i` vanishes, `v` is not parallel to either live `C`.
  Then `delta` is nonzero modulo three.  The two moving coordinates vanish
  in two distinct `t` classes, leaving exactly one live longitudinal class.

This proves, without parity assumptions,

```text
F_v=(2/3) sum_k f(3k)                         (unit v),
F_v=(1/3)(sum_k f(k)-sum_k f(3k))             (one-zero v).
```

The independent finite-field audit exhausts all eight unit speed residues,
all 64 nonzero relation residues, and 128 live `(v,C)` observations.  It
finds 16 unit-relation and 48 one-zero-relation cases, with exactly the
claimed defect and longitudinal counts.

## 2. Slice, lattice discrepancy, and intercept

For a primitive relation `v`, set `S=||v||_1`, `M=max|v_i|`, and `r=3/14`.
Choose `v_i!=0` and map the error cube by

```text
e -> (v.e, (w cross e)_i/v_i).
```

The image is a centrally symmetric planar zonotope.  Its three generator-pair
determinants have absolute values `a,b,c`: for example, in the `i=1` chart
the minors are `-c,b,-a`, with the last identity using `v.w=0`.  Consequently

```text
I = integral f = 4r^2(a+b+c) <= 27c/49.
```

The closed section width is even and concave on the interior of its support.
Resetting a possibly positive endpoint section to zero matches the strict
carrier roofs and preserves monotonicity and the integral.  Rectangle
comparison therefore gives

```text
|sum_k f(hk)-I/h| <= f(0), h>0.
```

Using the modulo-three dichotomy in the correct upper/lower directions gives

```text
F_v <= 2I/9 + 2f(0)/3.
```

At defect zero, `w cross e=t v`, so in a largest coefficient coordinate

```text
f(0) <= 2r(w_j+w_k)/M <= 6c/(7M),
F_v/c <= 6/49+4/(7M).
```

For `M>=19`, the last expression is `142/931<15/98`.

For the intercept, a union of two live residue classes on one unit-relation
line has count strictly below `2f/3+4/3`; a single live residue class in the
one-zero case has count strictly below `f/3+1`.  Counting defects yields

```text
B_unit < 4S/21+4/3,
B_onezero < 2S/7+4/3,
N < F_v + B_v <= 15c/98+2S/7+4/3.
```

One minor all-odd sentence cannot be carried over literally: a one-zero
relation need not have even norm at least six.  In arbitrary parity the
one-zero defect set is empty for `S=1,2,3,4`; then `N=0` and the target is
automatic.  For `S>=5`, `+/-1` are available.  The clean-room arithmetic
script checks both intercept inequalities through `S=10000` and identifies
exactly this empty range.

Solving the last count bound against `N<2c/11` gives the same exact gate

```text
c >= (308/31)S + 4312/93.
```

## 3. Expanded coefficient box

With no parity restriction, the complete `M<=18` magnitude universe is

```text
0<=p1<=p2<=p3<=18,
at least two nonzero entries,
gcd(p1,p2,p3)=1,
at most one entry zero modulo three,
p!=(0,1,1).
```

The final pattern is impossible for distinct speeds.  There are exactly 750
patterns: 49 of support two and 701 of full support.  The three hypothesized
residuals are `(1,1,1)`, `(1,1,2)`, and `(1,2,2)`.

The clean-room compiler imports no repository code.  For every signed sector
it constructs the normalized speed polygon `{w in [0,1]^3:v.w=0}`.  For each
integer defect it intersects `v.e=delta` with all twelve exact cube edges and
evaluates the two extrema of `(w cross e)_i/v_i`.  Section width is convex in
`w`, so its sum is maximized on the speed-polygon vertices.  The audit covers
11,547 signed sectors, 32,235 speed vertices, and 427,180 explicit gates.

The low-pattern closed maxima are

```text
(1,1,1): 2/7,
(1,1,2): 2/7,
(1,2,2): 3/14.
```

Every other pattern has maximum at most `15/98`; equality is unique to
`(1,7,8)` and occurs at the equal-speed boundary `w=(1,1,1)`.  This proves
the required non-low slope bound.  Ordinary and optimized Python outputs are
byte-identical.

## 4. Minkowski cutoff

The projected relation lattice

```text
L={(x,y) in Z^2: ax+by=0 mod c}
```

has determinant `c`, because primitive `(a,b,c)` makes the map to `Z/c`
surjective.  The relation `l1` ball has exact area

```text
2c(ab+ac+bc)/((a+b)(a+c)(b+c)) L0^2 > (3/4)L0^2.
```

After setting `c=1`, strictness is the positive identity

```text
8(ab+a+b)-3(a+b)(a+1)(b+1)
=3a(1-a)(b+1)+3b(1-b)(a+1)+2a(1-b)+2b(1-a).
```

Planar Minkowski with `L0=4sqrt(c/3)` supplies a primitive relation with
`S<L0`.  The arbitrary-parity cutoff is:

- `S<=57`: the intercept threshold is at most
  `56980/93<613`;
- `S>=58`: `c>3S^2/16`, while at 58
  `3S^2/16=2523/4 > 57904/93`, and the difference is increasing.

Thus every non-low row with `c>=613` has `N<2c/11`, so **every** projection
is strictly below `6/77`.  Since `c` is a ternary unit, the exact residual
head is `c<=611`; `c=612` is inadmissible.

## 5. Definition-level height-611 head

The independent C++ engine does not include or import a project source.  It
enumerates the literal integer carriers from `C.w=0`, the three strict roof
inequalities, and the three nonzero residue conditions.  It solves one
Diophantine row only to avoid scanning empty lattice boxes.  For denominator
`D=14abc`, it accumulates the exact raw-projection numerators

```text
P_a += a min(p_a,6b),
P_b += b min(p_b,6a),
P_c += min(c p_c,6ab)
```

over all signed carriers, so `E_i=P_i/D`.  A separate quadratic-box
enumerator agrees on every one of the 1,205 eligible rows through height 31.
The explicit low-circuit equations agree with an independent signed-
permutation test on the same control universe.  Ordinary and optimized
height-199 transcripts are byte-identical.

The complete optimized height-611 result is

```text
eligible rows                                  9,720,930
non-low rows                                   9,653,396
low rows                                          67,534
positive carrier pairs                          182,289,317
non-low minimum >=6/77                                    0
non-low individual projection >6/77                       0
non-low individual projection =6/77                       3
```

The exact low-pattern cross-tab is

```text
pattern       rows     min>6/77   min=6/77
(1,1,1)      14,220       14,219            0
(1,1,2)      28,438            1            1
(1,2,2)      24,876            0            0
```

The lone safe `(1,1,1)` row is `(1,4,5)`, with
`E=(13/140,4/35,1/28)`.  The `(1,1,2)` failure is `(2,11,20)`, with

```text
E=(131/1540,11/140,3/35), min(E)-6/77=1/1540.
```

Its unique equality is the known `(1,5,11)`.  The totals happen to give
14,220 strict failures, but it would be wrong to identify that total with
the 14,220 additive-circuit rows: it is 14,219 additive failures plus the
single norm-four failure.

Among non-low rows, the largest minimum is `3/70` at `(2,11,40)`.  The three
individual equalities are exactly those in the verdict.  The complete-row
commutative digest is

```text
sum=11036112604704679290 xor=3479353455816149694.
```

## 6. Required hostile

`(2,5,7)` is a genuine excluded hostile.  Its only live carriers are
`+/-(1,1,-1)`, and direct exact summation gives

```text
E=(22/245,6/49,1/10),
```

so all three projections exceed `6/77`.  It is the first primitive sorted
ternary-unit minimum-projection failure.

## 7. Artifacts, hashes, and reproduction

The independent audit artifacts were promoted beside THM-4437. Hashes below
are for raw LF repository bytes after promotion.

```text
0caebf01c71e03e8f862e0e466b3a7c517359986e4a75c920664625afcb18782  lrc14_all_parity_analytic_chain_thm4437_independent.py
40238160e94ba9a35d71427e1472059ce1f32dcfcfd56761f1b098558e2e3329  lrc14_all_parity_analytic_chain_thm4437_independent.out
837dd8d7bdebdc9cc7535ea7845d21064498b925134950269575b6067f98548a  lrc14_all_parity_coefficient_box_thm4437_independent.py
734c6fddb23e3d8f732622a6974b2a639c4517e644743cf8f14674b3eb8b2aaa  lrc14_all_parity_coefficient_box_thm4437_independent.out
a931ecec254e5a2ae689d395a3c6b400608656aeebcebf0120440a291ba40c1b  lrc14_all_parity_raw_head_thm4437_independent.cpp
91ce462ddb039e59c239c1750b03792ad3f1b5f4d9645ff6255fcceb95d1ffde  lrc14_all_parity_raw_head_h199_thm4437_independent.out
1874b77966a59368a5760850ac8e253e8f6ea07b7a580bdb3be15523b15c566a  lrc14_all_parity_raw_head_h611_thm4437_independent.out
```

```powershell
python -B 04-computation/lrc14_all_parity_analytic_chain_thm4437_independent.py
python -B -O 04-computation/lrc14_all_parity_analytic_chain_thm4437_independent.py
python -B 04-computation/lrc14_all_parity_coefficient_box_thm4437_independent.py
python -B -O 04-computation/lrc14_all_parity_coefficient_box_thm4437_independent.py
g++ -std=c++17 -O0 04-computation/lrc14_all_parity_raw_head_thm4437_independent.cpp -o raw_O0
g++ -std=c++17 -O3 -DNDEBUG 04-computation/lrc14_all_parity_raw_head_thm4437_independent.cpp -o raw_O3
./raw_O0 199
./raw_O3 199
./raw_O3 611
```

## Scope

This audit proves the repaired raw-network reduction relative to the exact
THM-4414 carrier dictionary.  It does not prove low-circuit safety, physical
entry, synchronization, a ten-body Haar floor, or `LRC(14)`.
