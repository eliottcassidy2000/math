---
id: THM-4425
title: "LRC14 all-height rank-one carrier closure"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4414/4422 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. Every eligible three-speed
  row whose complete live carrier set is contained in one rational line
  satisfies the degree-zero raw projection target at every height. Outside
  the three signed norm-four directions, the carrier count itself is below
  the automatic deficit threshold from height 35 onward; the remaining
  finite base is inherited from THM-4422. In addition, every row whose
  shortest primitive relation has l1 norm at most fourteen satisfies the
  projection target: the one-zero sector is automatic except for the strict
  A2 control (19,23,29), while the full-support sector reduces to rank one.
  Equality is unique at (1,5,11).
  This is not a rank-two carrier theorem, chart-entry or synchronization
  result, or a proof of LRC(14).
source: root + dense_geometry_dichotomy / LRC14 continuation session, 2026-09-05
depends_on:
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
  - THM-4422-lrc14-projection-deficit-and-beatty-row-reduction
related:
  - THM-4420-lrc14-near-doubling-ray-network-closure
  - THM-4423-lrc14-height-499-atlas-and-ninety-carrier-global-closure
script: 04-computation/lrc14_rank_one_carrier_closure_thm4425.py
output: 05-knowledge/results/lrc14_rank_one_carrier_closure_thm4425.out
script_sha256: f001d624822d86cdd9b62ccd7ca6159cadc0c5ca1bb9b7343cae094620c9c571
output_sha256: ba61466cc8f90f0f68924c85997da5e9e82aeb0c080ad8168ac22a58023e805d
independent_script: 04-computation/lrc14_rank_one_carrier_closure_thm4425_independent_referee.py
independent_output: 05-knowledge/results/lrc14_rank_one_carrier_closure_thm4425_independent_referee.out
independent_script_sha256: 0123a7b27c39c960a59185eafa16668f4b78f0d6ffcb763ab2406eaf3a9986d9
independent_output_sha256: a581163955173acffd67bc1655dd94c7a688f1b813a574fadc3207d6d9bb15c5
geometry_script: 04-computation/lrc14_carrier_zonotope_width_additive_dichotomy_referee.py
geometry_output: 05-knowledge/results/lrc14_carrier_zonotope_width_additive_dichotomy_referee.out
geometry_script_sha256: 7df1d419aeee70f77f738e51cf20989575ea30c7cb6ff1fac2b51631a6078dbf
geometry_output_sha256: c5c91e52cbc0e2c48f0df356b49450891498c4205c582b93d003315b55a6dba1
one_zero_script: 04-computation/lrc14_one_zero_shortest_width_le8_automatic_count_referee.py
one_zero_output: 05-knowledge/results/lrc14_one_zero_shortest_width_le8_automatic_count_referee.out
one_zero_script_sha256: babbe01472a6b3aae3630179ed461d345fcc2bca51a026890981e5c9a85b11fc
one_zero_output_sha256: d972435f2b8b6f29cce94b86b1e726b2dc3e35616a763541ad252d8831a6fe32
one_zero_wide_script: 04-computation/lrc14_one_zero_shortest_width_10_14_layer_count_referee.py
one_zero_wide_output: 05-knowledge/results/lrc14_one_zero_shortest_width_10_14_layer_count_referee.out
one_zero_wide_script_sha256: d2176e34ceb79d10af7f8e74936b8a696898437b6d1fd1ddd6fbd71174f79ece
one_zero_wide_output_sha256: 9973ed3d57944dc2310209ec8b027dcbb985e7f1b83f1268cc0f4252174aff81
hash_basis: raw LF bytes
audit: >
  PASS. A standalone definition-level lattice enumerator checks all 165
  primitive eligible triples below height 35, including 158 rank-at-most-one
  rows and the unique equality. It independently enumerates the small
  coefficient patterns, verifies every complete line-address equality, and
  checks the exact height threshold. Its 1,069
  explicit gates remain active under optimization; normal and optimized
  outputs byte-match the frozen LF artifact. A no-import referee repeats the
  complete height-31 base, extends line dictionaries through height 79, and
  audits the exact all-height arithmetic in 29,689 further gates. Its normal
  and optimized outputs also byte-match. A second no-import geometry referee
  checks the support, covolume, residue, width, layer and half-body identities,
  including their sharp controls, in 9,178 optimization-live gates. A third
  standalone exact referee exhausts the norm-six and norm-eight one-zero
  relation shapes, proves the long-height count bounds, and checks the finite
  base in 3,506 optimization-live gates; its normal and optimized outputs
  byte-match the frozen LF artifact. Its four-layer companion proves the
  exact norm-ten/twelve/fourteen chord bound, reconstructs the finite base,
  and isolates the unique nonautomatic A2 row in 28,434 further gates.
  Normal and optimized outputs again byte-match the frozen LF artifact.
---

# THM-4425 -- LRC14 all-height rank-one carrier closure

**PROVED ELEMENTARY RELATIVE TO THM-4414/4422 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.** Every one-direction live
carrier configuration satisfies the local raw projection target, as does
every row of shortest dual `l1` width at most fourteen. The proof does not
construct physical entry or synchronize the other eleven speeds; `LRC(14)`
remains **OPEN**.

## 1. Statement

Let

```text
w=(a,b,c),       1<=a<b<c,                              (1)
```

be primitive, odd, and nonzero modulo three. Retain THM-4422's complete
live carrier set

```text
Lambda(w)={C in Z^3:C dot w=0,
  every C_i nonzero mod 3,
  14|C_i|<3(w_j+w_k) for {i,j,k}={1,2,3}}.             (2)
```

Assume `Lambda(w)` is contained in a one-dimensional rational subspace of
`w^perp`. Then its three THM-4414 raw projection capacities satisfy

```text
min_i S_i(w)<=6/77.                                     (3)
```

Equality in `(3)` occurs only for `w=(1,5,11)`.

Independently of the rank-one hypothesis, the same conclusion holds whenever
the shortest primitive integer relation `n dot w=0` has `||n||_1<=14`.
The only non-count-automatic row in that added sector is the strict
`A_2` control `(19,23,29)`.

## 2. Exact natural address on a carrier line

The empty case is immediate. Otherwise choose the primitive integer
generator `u=(u_1,u_2,u_3)` of the carrier line. Every live carrier has the
form

```text
C=k u.                                                  (4)
```

One live value in `(4)` shows that every `u_i` and its address `k` are
nonzero modulo three. Put

```text
B=min_i 3(w_j+w_k)/(14|u_i|).                           (5)
```

Conversely, every multiple satisfying the strict roofs and residue gate is
a member of `(2)`. Hence the complete carrier line is exactly

```text
|k|<B,       3 does not divide k.                       (6)
```

```text
N=|Lambda(w)|=2 R_<(B),
R_<(B)=#{k in Z:1<=k<B, 3 does not divide k}.           (7)
```

This is the lossless part of the natural-number model: global sign supplies
the two orientations, while a positive address is one ordinary integer with
every third value deleted. Direction `u` remains a necessary sidecar; the
ordinal count alone does not determine a projection.

If `M=ceil(B)-1`, then `M<B` and three-term blocks give

```text
R_<(B)=M-floor(M/3)<=(2M+2)/3<(2B+2)/3.                (8)
```

## 3. Only norm four can evade the count bound

Since all speeds are odd and `u dot w=0`, reduction modulo two gives

```text
u_1+u_2+u_3=0 mod 2.                                   (9)
```

Suppose every `|u_i|<=2`. Primitivity, the ternary-unit condition and `(9)`
leave exactly the magnitude multiset

```text
{|u_1|,|u_2|,|u_3|}={1,1,2}.                           (10)
```

Indeed `(1,1,1)` and `(1,2,2)` have odd coordinate sum, while `(2,2,2)` is
not primitive. Equation `(10)` is precisely one of the three sorted signed
norm-four directions closed at every height by THM-4422.

Otherwise some `|u_i|>=4`. Because each sum of two speeds is strictly less
than `2c`, equations `(5)` and `(8)` give

```text
B<3c/28,
N<4B/3+4/3<c/7+4/3.                                   (11)
```

For `c>=35`,

```text
c/7+4/3<=2c/11                                         (12)
```

because `(12)` is equivalent to `c>=308/9`. Hence THM-4422's exact deficit
identity makes every non-norm-four rank-one row with `c>=35` automatic:

```text
N<=2c/11  implies  min_i S_i<=6/77.                    (13)
```

The eligible heights below 35 end at 31 and lie inside THM-4422's complete
height-79 exact universe. That finite base has no hostile. This proves `(3)`.

The inequalities in `(11)` are strict, so a non-norm-four row cannot attain
equality. THM-4422 proves that the norm-four families have unique equality
at `(1,5,11)`, completing the equality statement.

## 4. Exact support width and the short-relation closure

The rank-one hypothesis itself has a useful intrinsic supplier. Put
`r=3/14`, let `L=ker_Z(w)`, and consider the open support zonotope

```text
K=-w cross (-r,r)^3 subset w^perp.                     (14)
```

This is an equality, not merely containment. The fibre of
`e -> -w cross e` is a line parallel to `w`; the three coordinate conditions
on `e` are intervals on this line. Their pairwise-intersection conditions are
exactly the three carrier roofs, so the interval Helly property gives a common
preimage.

For the three generating segments `g_i=-w cross e_i`, the complementary
minor has absolute area `w_k||w||_2`. The planar zonotope formula therefore
gives

```text
area(K)=9||w||_2(a+b+c)/49,       covol(L)=||w||_2.    (15)
```

The live residue gate selects exactly two of the nine cosets of `3L`, namely
the classes `+/-w` modulo three. Thus its exact translation-mean, or
continuous bulk term, is

```text
N_bulk=2(a+b+c)/49.                                    (16)
```

This is an area/covolume term, not the exact count at the distinguished
translate; short directions can create linear discrepancy. Since distinct
sorted odd speeds obey `a+b+c<=3c-6`, every nonautomatic row has the explicit
linear excess

```text
N-N_bulk>(32c+132)/539.                                (17)
```

That discrepancy has an exact width coordinate. The quarter-turn map

```text
J_w:L^* -> L,       eta -> w cross eta                 (18)
```

is an integral isomorphism. If `n=J_w(eta)` is primitive, then the width of
`K` in the primitive `3L`-dual direction `eta/3` is exactly

```text
width_(eta/3)(K)=||n||_1/7.                            (19)
```

Indeed, for `C=-w cross e`, the layer coordinate is

```text
t=eta dot C=n dot e,       |t|<(3/14)||n||_1.          (20)
```

Now suppose `n` is full-support modulo three. Both `n` and every live `C`
reduce to one of the two full-support kernel classes, so they are parallel
modulo three. Since `eta dot n=0`, equation `(20)` has

```text
t=0 mod 3.                                              (21)
```

If `||n||_1<=14`, strictness in `(20)` gives `|t|<3`, hence `t=0`.
The primitive kernel of `eta` inside `L` is `Z n`, so every live carrier lies
on that line. The theorem therefore has the all-height corollary

```text
primitive full-support relation n with ||n||_1<=14
  implies min_i S_i<=6/77.                             (22)
```

### 4.1. Every one-zero width through fourteen is closed

Suppose instead that primitive `n` has exactly one zero coordinate modulo
three. The live classes then have

```text
t=+/-1 mod 3.                                           (22a)
```

Parity forces `||n||_1` to be even. Norm two repeats one of the distinct
speeds, and norm four has no feasible one-zero relation. At norm six, the
absolute shapes are exactly

```text
(0,1,5), (1,2,3).                                       (22b)
```

The ten sorted feasible signed placements each cut a live half-class into
one residue class modulo three along an open interval of length at most
`c/7`. Consequently its positive half-count obeys

```text
m<=ceil(c/21).                                          (22c)
```

This gives `N=2m<=2c/11` for every eligible `c>=11`. The sole smaller
eligible row admitting such a relation is `(1,5,7)`; it has a shorter
full-support norm-four relation and

```text
N=2,     (S_1,S_2,S_3)=(8/245,6/49,4/35).              (22d)
```

At norm eight the exhaustive absolute shapes are

```text
(0,1,7), (1,1,6), (1,3,4).                              (22e)
```

Every placement has a coefficient of magnitude at least four, so one roof
gives

```text
m<c/14+1<=c/11                 for c>=53.               (22f)
```

The complete exact base through `c=49` closes all remaining placements. Thus

```text
shortest primitive relation n one-zero mod 3,
||n||_1<=8   implies   N<=2c/11 and min_i S_i<=6/77.    (22g)
```

The word `shortest` removes the only count-hostile `(1,5,7)` because its
shortest relation is full-support norm four.

For norms ten, twelve, and fourteen, strictness in `(20)` permits exactly the
four possible layers

```text
t in {-2,-1,1,2}.                                       (22h)
```

Fix a Bezout vector `q` with `n dot q=1` and put `m=w cross q`. The complete
carrier dictionary is the disjoint union

```text
beta*m+k*n,       beta in {-2,-1,1,2},                  (22i)
```

subject to the original roofs, with one allowed residue class of `k` modulo
three on each affine line. The positive owner class uses `beta=1,2`. Let
`ell_1,ell_2` be the open chord lengths cut from those two lines.

Every feasible signed relation pattern of norm ten, twelve, or fourteen
satisfies the exact bound

```text
ell_1+ell_2<=3c/14.                                     (22j)
```

For maximum coefficient at least eight, two single-roof bounds prove
`(22j)`. The remaining 49 feasible small-coefficient signed patterns reduce
to one-dimensional rational speed cones; exact endpoint comparison gives
the same bound, sharply for the norm-twelve shape `(1,5,6)`. Since each chord
retains one residue class modulo three, the positive half-count satisfies

```text
m<(ell_1+ell_2)/3+2<=c/14+2<=c/11       for c>=103.    (22k)
```

The complete base through `c=101` has 5,937 eligible rows. It contains 373,
1,288, and 578 rows whose shortest one-zero relation has norm ten, twelve,
and fourteen respectively. The sole nonautomatic row is

```text
w=(19,23,29),       n=(3,-5,2),       N=6,
(S_1,S_2,S_3)=(156/4669,192/3857,3840/88711),          (22l)
```

whose minimum is strictly below `6/77`. This is the first four-layer `A_2`
carrier hexagon, not a count-automatic row.

Finally, for ternary-unit `w`, a primitive relation modulo three cannot have
two zero coordinates, while three zero coordinates would contradict
primitivity. Hence every shortest relation is either full-support or
one-zero. Combining `(22)` with `(22g)--(22l)` proves the intrinsic closure

```text
mu_1(w)<=14       implies       min_i S_i(w)<=6/77.     (22m)
```

Every remaining projection hostile therefore has shortest dual l1 width at
least sixteen.

## 5. A half-body additive boundary dichotomy

The same two live cosets supply a rigorous rank-versus-boundary statement.
Let

```text
A=K intersect (w+3L),       m=|A|=N/2,
P=A intersect (1/2)K,       h=|P|.                    (23)
```

For `x,y in P`, convexity gives `x+y in K`, while its residue is
`2w=-w mod 3`. Hence

```text
P+P subset -A,       |P+P|<=m.                        (24)
```

The torsion-free sumset bound gives `|P+P|>=2h-1`. If `P` has affine
dimension two, the elementary Freiman dimension lemma sharpens this to
`|P+P|>=3h-3`. Therefore

```text
h<=(m+1)/2;
and either P is collinear or h<=(m+3)/3.               (25)
```

Thus at least about half of a live class is outside the half-zonotope, and
about two thirds is outside when the inner core is genuinely two-dimensional.
This is not yet THM-4422's weighted deficit inequality: `(25)` forgets which
coordinate strip owns each outer carrier and how far it lies beyond that
projection's hinge.

The dimensional bound is sharp on a structured sparse control:

```text
w=(85,97,107),       N=12,
P={(-17,16,-1),(-2,-17,17),(19,1,-16)},
sum P=0,             A=P union (-2P),       P+P=-A.    (26)
```

This two-shell `A_2` circuit refutes the stronger guess that every half-body
core is collinear while identifying the equality objects that a labelled
boundary-depth refinement must preserve.

## 6. What has actually been removed

Every counterexample to the remaining THM-4414 projection inequality must
now contain two nonparallel live carrier directions and have shortest dual
`l1` width at least sixteen. This is stronger than removing a list of guessed
speed rays: it removes every rank-one carrier configuration, including
directions of arbitrarily large coefficient norm, and every narrow
multi-layer configuration through width fourteen.

The first finite two-dimensional dense control remains THM-4422's additive
hexagon

```text
w=(19,23,29),
Lambda(w)=+/-{u,v,u+v},                                (27)
```

with `u=(1,8,-7)` and `v=(10,-7,-1)`. It satisfies `(3)` but is not covered
by the present theorem. The sharp next object is therefore not another
one-ray coefficient atlas. It is the affine two-coset lattice polygon formed
by the positive natural-address class, with additive circuits such as `(27)`
retained rather than quotiented to a carrier count.

## 7. Reproduction

Run from the repository root:

```powershell
python -B 04-computation/lrc14_rank_one_carrier_closure_thm4425.py
python -B -O 04-computation/lrc14_rank_one_carrier_closure_thm4425.py
python -B 04-computation/lrc14_rank_one_carrier_closure_thm4425_independent_referee.py
python -B -O 04-computation/lrc14_rank_one_carrier_closure_thm4425_independent_referee.py
python -B 04-computation/lrc14_carrier_zonotope_width_additive_dichotomy_referee.py
python -B -O 04-computation/lrc14_carrier_zonotope_width_additive_dichotomy_referee.py
python -B 04-computation/lrc14_one_zero_shortest_width_le8_automatic_count_referee.py
python -B -O 04-computation/lrc14_one_zero_shortest_width_le8_automatic_count_referee.py
python -B 04-computation/lrc14_one_zero_shortest_width_10_14_layer_count_referee.py
python -B -O 04-computation/lrc14_one_zero_shortest_width_10_14_layer_count_referee.py
```

The primary rank-one verifier performs 1,069 checks, its clean-room referee
performs 29,689, and the geometry referee performs 9,178. Ordinary and
optimized outputs byte-match their respective frozen raw-LF artifacts. The
one-zero-width referee adds 3,506 exact optimization-live checks.
The four-layer referee adds 28,434 more.
