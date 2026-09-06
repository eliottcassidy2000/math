# Independent audit of the universal odd ternary-unit local network bound

**Verdict: PASS at proof level, with explicit finite-evidence scope.** The
[global slope argument](lrc14_global_slope_empty_core_certificate_sep06.md)
and its referenced slice/count reductions have the correct quantifiers,
scales, residue densities, strictness, and cutoff. I independently recomputed
the **entire 308-pattern coefficient box** by a third width algorithm and
matched every exact pattern/defect/slope/intercept record. I audited the
native H601 source, universe, transcript, and frozen hashes, then replayed
a fresh complete H97 native head against a separate raw integer enumeration.
The full H601 interval computation was **not rerun** in this audit; its
checked-in FINITE-EXACT result remains the final finite premise.

The accepted assembled statement is precisely:

```text
primitive 1<=a<b<c, all speeds odd and nonzero modulo three:
min_i E_i(a,b,c)<=6/77,
equality iff (a,b,c)=(1,5,11).
```

Without a signed norm-four relation, all three projections are strictly
below `6/77`. The corresponding local physical triple-comb mass is at most
the minimum projection. Neither the local statement nor this audit proves
entry, synchronization, or LRC(14).

The current canonical reference is
[THM-4434 / lrc14-universal-scale-three-network-projection-bound](../../01-canon/theorems/THM-4434-lrc14-universal-scale-three-network-projection-bound.md),
verified directly at incoming revision `d2f64b809`. Its status is now
**PROVED ELEMENTARY RELATIVE TO THM-4414/4422 + FINITE-EXACT +
INDEPENDENTLY AUDITED**. It states the same odd ternary-unit bound and
explicitly leaves the quantitative body floor, arbitrary entry,
synchronization, and LRC(14) open. Its separate import-free coefficient
referee, introduced in `a60c1047d`, is additional independent evidence;
this audit uses its own third width algorithm below.

The separately proposed universal ten-body floor in the consumer note is
now **REFUTED** by recovered THM-530 Section A, as reconstructed in the
[current consumer audit](overnight3_20260906_lrc_consumers.md). This does
not affect the odd local theorem or its corollary with an explicit body
measure hypothesis. The global application needs a justified body subclass
or direct component noncontainment, in addition to its entry conditions.

For correction lineage, at the originally inspected revision
`2e1ef24c4aa3` the namespace was still **RESERVED / UNPROVED EMPTY STUB**,
despite the completed result-note chain. The incoming promotion resolves
that metadata boundary. This audit neither edited nor promoted the canon.

## 1. Inheritance and the exact consumer

The audited source chain is:

- [Full-support slice note](lrc14_empty_core_certificate_sep06.md), for the
  affine integer defect, complete scalar interval, residue word, and convex
  speed-polytope maximum.
- [Actual-zero-coordinate addendum](lrc14_pair_relation_empty_core_certificate_sep06.md),
  for support-two relations and their free error coordinate.
- [Complete coefficient box](lrc14_coefficient_box_empty_core_three_ray_sep06.md),
  for the finite normalized slope bound.
- [Global zonotope/short-relation reduction](lrc14_global_slope_empty_core_certificate_sep06.md).
- [Complete native finite head](lrc14_universal_literal_empty_core_sep06.md),
  whose source and saved H601 output are independently audited here.

The actual consumer is
[THM-4414 / six-separated-contact-capacity-collapse](../../01-canon/theorems/THM-4414-lrc14-six-separated-contact-capacity-collapse.md),
relative to THM-4409. It identifies the complete degree-zero sheet-network
capacities with

```text
E_i=sum_(C in Lambda(w)) min(q,p_i(C)/(14w_jw_k)),
q=3/(7c),     p_i(C)=3(w_j+w_k)-14|C_i|,
Lambda(w)={C in Z^3:C dot w=0, all C_i!=0 mod3, all p_i(C)>0}.
```

Thus `E_i<=qN`, where `N=|Lambda(w)|`, and
[THM-4422 / projection-deficit-and-beatty-row-reduction](../../01-canon/theorems/THM-4422-lrc14-projection-deficit-and-beatty-row-reduction.md)
supplies both the count gate `N<=2c/11` and the separate all-height signed
norm-four closure, including its unique equality. The physical mass uses
the minimum of all three roof terms **inside** the carrier sum. It is not
interchangeable with a chosen network projection.

The inherited hostiles are `(1,5,7)`, where two projections exceed the
target, and `(1,19,79)`, where network capacity exceeds physical mass. The
corrected near miss is dropping a live ray or a residue/affine layer. The
least-used sidecar is the entire slice integral together with its central
width; neither the average alone nor a small-direction census is enough.
The concept board was: raw address; integer defect; owner residue word;
zonotope section; lattice quadrature; shortest relation; native contact
scaling. Every stage below keeps the data consumed by the next one.

## 2. Affine address, owner densities, and normalization

For a primitive relation `v dot w=0`, choose integer lifts `n` with
`C=w cross n`, unique modulo `Zw`. The cross-product map is onto the
integral kernel for primitive `w`: if `z dot w=1`, then
`n=C cross z` satisfies `w cross n=C`. Its kernel is exactly `Zw`.

The three strict roof inequalities are equivalent to pairwise strict
intersection of the three allowed real phase intervals. For intervals on
a line, this implies a common interior point, hence an error vector in
`(-3/14,3/14)^3`. Therefore every carrier has an integer defect

```text
delta=v dot n=v dot e,       |delta|<(3/14)||v||_1,
v cross C=delta w.
```

At each integer defect, primitivity of `v` provides an integer lift, and
the integer carrier line is exactly `C_delta+Z v`: two integral points on
the same real line differ by an integer multiple of primitive `v`. No
additional lattice index is dropped.

Modulo three, there are exactly two primitive relation types. For a unit
relation, live carriers require `delta=0 mod3` and occupy two multiplier
classes out of three. For a relation with one zero residue, they require
`delta!=0 mod3` and occupy one class out of three. Two zero residues would
force the third to vanish because the speeds are units, contradicting
primitivity. This includes an actual zero integer coefficient. The new
verifier checks all **192** affine residue fibers over all ternary-unit
speed residue triples, with exactly these owner counts.

Let `T_delta` be the width of the complete projected error slice in units
of `v`. Live carriers occupy an **open** interval. One residue class has
count `<T_delta/3+1`; two residue classes have count
`<2T_delta/3+4/3`, since their alternating gaps are one and two. Summing
over the complete allowed integer defect list gives

```text
N<F_v(w)+B_v.
```

For fixed `v,delta`, the width is convex and homogeneous in `w`, so it
is enough to maximize on vertices of `W_v=[0,1]^3 intersect v^perp`.
Actual normalization gives `T_delta(v,w)=c T_delta(v,w/c)` while keeping
the same integer defect list; one does not need an integer lift at a
normalized polytope vertex. Equal or zero speeds at its boundary enlarge
the optimization domain and are legitimate for an upper bound.

## 3. Zonotope area and both discrepancy signs

For a nonzero pivot `v_i`, the map
`e -> (v dot e,(w cross e)_i/v_i)` takes the error cube to a planar
zonotope. With pivot one, its generator pair determinants are exactly
`-c,b,-a`. The last equality uses `v dot w=0`; omitting that relation
would give the wrong area and units. Hence

```text
I=integral f(delta)d delta=4(3/14)^2(a+b+c)
 =9(a+b+c)/49<=27c/49,
f(delta)=T_delta.
```

The section length is even and concave on its support, hence decreasing
on the positive half-line. At strictly excluded endpoints it is set to
zero. An actual zero coefficient can leave a vertical endpoint edge in
the closed zonotope, but replacing that endpoint width by zero changes
neither area nor monotonicity. The independent endpoint hostile
`v=(13,-1,0), w=(1,13,19), delta=3` has closed width `3/7`, explicitly
showing why endpoint continuity must not be assumed.

Rectangle comparison gives

```text
|sum_k f(hk)-I/h|<=f(0).
```

In the unit case, `F=(2/3)sum_k f(3k)`. In the one-zero case,
`F=(1/3)(sum_k f(k)-sum_k f(3k))`; its upper bound must use the upper
error on the first sum **and the lower error on the subtracted sum**.
Both cases therefore give

```text
F<=2I/9+2f(0)/3.
```

At defect zero, the cross product is a scalar multiple of `v`. Using a
coordinate of maximum coefficient `M` gives
`f(0)<=2(3/14)(w_j+w_k)/M<=6c/(7M)`, whence

```text
F/c<=6/49+4/(7M).
```

For `M>=19`, this is at most `142/931=15/98-1/1862`, strictly below
the desired slope. The remainder is exactly the coefficient box `M<=18`,
not a norm-twenty or speed-height extrapolation.

## 4. Complete coefficient box: independent third algorithm

The independent verifier generates the full ordered Cartesian cube
`{0,...,18}^3`, applies the primitive/even-norm/support/residue filters,
then sorts magnitudes to quotient coordinate permutations. It removes
only `(0,1,1)`, impossible for distinct positive speeds, and the inherited
norm-four exception `(1,1,2)`. It obtains exactly 308 patterns:
293 full-support and 15 support-two.

For the width computation, it **does not enumerate or clip error-slice
vertices**. Change each constrained error variable to
`x_i=sign(v_i)e_i` and then `z_i=x_i+3/14`. The slice equation becomes

```text
sum_i |v_i| z_i=delta+(3/14)||v||_1,
0<=z_i<=3/7.
```

A linear scalar objective is maximized by filling the variables in
descending objective-per-resource ratio, and minimized in ascending order.
The exchange argument for continuous knapsack proves this exactly: moving
resource from a smaller ratio to an available larger ratio improves the
objective until one bound is tight. A zero coefficient is a free variable
and contributes its exact independent interval width. This yields the
whole scalar width without either prior error-polygon implementation.

Speed-polytope vertices still arise from two active cube bounds. Only
nonzero elimination pivots are used. Three isolated coefficient signs
suffice: global negation preserves widths and the symmetric defect list,
and coordinate permutations preserve the entire speed/error cube maximum.

All exact slope and intercept records agree with the frozen compiler:

```text
maximum slope: 15/98
unique maximizing magnitude pattern: (1,7,8)
maximum |delta| used: 11
complete pattern/defect/slope/intercept semantic SHA-256:
3be81c2522a1df6a146e50634754620103f9d2d8840d17f34c9e9a4954e849f7
```

The norm-four hostile separately has slope `2/7`. Pattern equality here
is a closed normalized width maximum, not network equality at a speed
triple. The only wording correction identified in the global note is its
statement that eligible full-support norms can reach 54: 54 is the ambient
box bound, while the filtered universe actually reaches **52**, at
`(17,17,18)`. This does not alter the allowed defect maximum 11 or any bound.

## 5. Intercept, shortest relation, and cutoff 601

For `S=||v||_1`, the unit defect list has cardinality `<S/7+1`; its
intercept is `4|D|/3<2S/7+4/3`. The one-zero residue-deleted list likewise
has `|D|<2S/7+4/3`. Thus, outside norm four,

```text
N<15c/98+2S/7+4/3,
c>=(308/31)S+4312/93  ==>  N<2c/11.                  (2)
```

The relation lattice projected to `(x,y)` has determinant `c` because
`(x,y)->ax+by mod c` is onto. The projected norm ball has the six vertices
listed in the global note. Direct shoelace computation gives

```text
area(K_L0)=2c(ab+ac+bc)L0^2/[(a+b)(a+c)(b+c)]
          >(3/4)L0^2.
```

The stated positive remainder after scaling `c=1` is algebraically correct
and strictly positive for `0<a,b<1`. The new verifier checks its polynomial
identity exactly, as well as independent six-vertex areas on named rows.

Apply the fundamental-domain overlap argument to the **open interior**
at `L0=4sqrt(c/3)`. It has area greater than `4c`; its half-body has area
greater than a lattice fundamental domain. Two distinct lattice translates
overlap after reduction, producing a nonzero lattice difference in the
original open body. Hence a nonzero relation has
`S<4sqrt(c/3)`. Dividing by its content preserves the relation and decreases
the norm. Since all speeds are odd, the resulting primitive relation has
even norm and falls into the audited residue types.

The last cutoff uses that evenness exactly. If `S<=56`, the right side
of (2)'s cutoff is at most `56056/93<603`. If `S>=58`, then

```text
g(S)=3S^2/16-(308/31)S-4312/93,
g(58)=3023/372>0,
g'(S)=3S/8-308/31>0 for S>=58.
```

The short-relation bound gives `c>3S^2/16`, so (2) again applies.
Thus all `c>=603` are covered; `c=602` is inadmissible by parity, leaving
the entire finite head `c<=601`. There is no rounding gap at the cutoff.

If the selected relation is norm four, the inherited THM-4422 theorem
applies directly. Otherwise all three projections are strict by (2).
The finite head and the unique norm-four equality consequently imply
the global equality exactly at `(1,5,11)`.

## 6. Native finite-head integrity and the fresh replay

The H601 C++ source evaluates every eligible speed triple before any
classification. I independently counted **1,317,935** primitive odd
ternary-unit triples, including exactly **9,201** satisfying a signed
norm-four identity. The identities are exactly
`c=2a+b`, `c=a+2b`, or `2b=a+c`; they are classifications after evaluation,
not carrier or coefficient filters.

The native engine uses denominator `42abc`, radius `3abc/w_i`, and sorted
centers from residue `(-w_i s) mod3` with step `42abc/w_i`. These are the
literal shifted sheet intervals, not a raw congruence formula. All six
distinct sheet assignments are retained. The three-pointer scan tests
strict positive triple intersection and advances every interval tied at
the smallest right endpoint. Disjoint intervals within each sheet make
this a complete contact scan without duplicates or skipped overlaps.

At each contact the omitted-sheet capacity is the minimum of its whole
interval length and the intersection length of the other two. THM-4414's
six-separation result identifies this with the actual local network load.
Each raw carrier has three translated native contacts, each contributing
one third of its raw term; the native count is therefore `3N`, while the
summed projection and mass retain the raw normalization.

The cut-at-zero convention causes no hidden cap loss. If the label-zero
half interval of speed `d` meets a nonzero-label interval of speed `s`,
then `11/(42s)<1/(14d)`, so `s>11d/3`. Each other interval is consequently
shorter than the half interval. Thus, when the wrapped sheet is omitted,
the other pair's length is already below the half capacity; when it belongs
to the pair, the nonwrapping partner makes clipping irrelevant. The right
endpoint is symmetric.

All endpoints and sums fit signed 64-bit arithmetic in the H601 universe;
ratio comparisons use signed 128-bit cross-products. Explicit runtime
failure checks survive `-DNDEBUG`. The inspected source and saved H601
output have exactly the LF SHA-256 values frozen in the native report.
The saved result checks every minimum, every projection outside norm four,
and the sole equality, rather than relying on a scalar checksum.

For a fresh replay I compiled the unmodified native source at `H=97` and
requested its exact per-row TSV. A new raw implementation enumerates the
second and third integer coordinates and solves the first. With denominator
`D=42abc`, it independently sums raw integer terms

```text
min(18ab,3w_i p_i(C)).
```

All **5,409** rows agree exactly in all three projections, physical mass,
and `native_contacts=3N`. The complete expected universe is independently
generated, with no missing or repeated row. It has 241 norm-four rows and
sole equality `(1,5,11)`. This corroborates the H601 compiler without
pretending to execute its entire larger head again.

## 7. Dependency and continuation boundary

The all-height argument needs the complete coefficient box and the exact
H601 projection head. The earlier coefficient-dependent norm-twenty speed
heads are not extra premises of the universal reduction: their reusable
geometric and residue lemmas are the inherited input. The older raw H601
multi-ray census is corroboration; the native head covers one-ray and empty
rows as well.

The oddness assumption cannot be dropped. The separate
[additive physical audit](overnight3_20260906_lrc_additive_audit.md)
rederives the parity-free raw physical address and shows a cofinite additive
family with mass above `6/77`. It also retains the nonadditive mixed-parity
norm-four hostile `(2,11,20)`. Those facts do not damage the audited odd
theorem. They identify a necessary boundary before using it in an entry or
Haar argument. No additive-only description of all mixed-parity failures
is claimed.

The full source-to-consumer contract is: complete carriers map to integer
defect lines; exact scalar widths and owner densities bound their count;
the count bounds all network projections by the common cap. The lost data
are individual roof values, retained only for the norm-four exception and
the exact finite head. The shortest relation supplies a global norm scale;
the integral and peak bound its complete defect sum. This is a valid chain
of sufficient inequalities, not an equivalence or a proof of chart entry.

## Reproduction and frozen evidence

```text
g++ -O2 -DNDEBUG -std=c++17 04-computation/lrc14_universal_literal_empty_core_sep06.cpp -o C:/w/overnight3_20260906_lrc_universal_audit_native.exe
C:/w/overnight3_20260906_lrc_universal_audit_native.exe 97 C:/w/overnight3_20260906_lrc_universal_audit_h97.tsv
python 04-computation/overnight3_20260906_lrc_universal_audit.py --native-tsv C:/w/overnight3_20260906_lrc_universal_audit_h97.tsv
python -O 04-computation/overnight3_20260906_lrc_universal_audit.py --native-tsv C:/w/overnight3_20260906_lrc_universal_audit_h97.tsv
```

The new verifier imports no repository mathematics and passes **132,655**
explicit gates. Normal and optimized output byte-match. The optional TSV
is generated scratch evidence, not an additional checked-in dependency.
Its LF SHA-256 is
`832b38af1c947e3537b7c2c868b53ca497eeddeafb0eb113a30f8037a6261038`.

Companions:
[independent verifier](../../04-computation/overnight3_20260906_lrc_universal_audit.py),
[audit output](overnight3_20260906_lrc_universal_audit.out),
[fresh native output](overnight3_20260906_lrc_universal_audit_native_h97.out).

```text
SHA-256, LF bytes:
source 7caa2c00caec9d96623f37aa061423e63d54bd8bb99ec5338a732c2d76447f49
output 603aef0e457a842e6c2010a5aa1a16bf8fe53768aea1a5ee5bbd273764de3442
native H97 output 60968807376af17200bb8cd7c66a933b2af9e7a3230993a343844482e651ce70
```
