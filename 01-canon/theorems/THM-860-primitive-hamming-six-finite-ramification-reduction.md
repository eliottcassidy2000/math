---
id: THM-860
title: Primitive Hamming-six finite ramification reduction
status: PROVED STRUCTURAL + FINITE-EXACT — every primitive proper AP-centred H6 packet at or below 1/13 has c<=840 and belongs to a finite exact component recursion; no emptiness claim
source: codex-2026-07-15-S14/S15 primitive-H6 transport and independent audit
depends_on: [THM-765, THM-810, THM-815, THM-857, THM-858, THM-859]
related: [THM-823, THM-840, THM-856, THM-861, THM-862, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_six_primitive_ramification_verifier_codex_S15.py
  - 05-knowledge/results/lrc13_hamming_six_primitive_ramification_verifier_codex_S15.out
  - 04-computation/lrc13_hamming_six_joint_ramification_cap_codex_S10.py
  - 05-knowledge/results/lrc13_hamming_six_joint_ramification_cap_codex_S10.out
---

# THM-860 — primitive Hamming-six is finite-ramified

Put `F=F_13^*` and `[12]={1,...,12}`.  Let `R subset F` have six
elements, put `P=[12] minus R`, and consider the AP-centred packet

```text
A=cP union {w_i:i in R},
w_i=cr_i (mod 13),        w_i>0,      w_i!=cr_i.         (1)
```

Assume that `A` is primitive and

```text
M(A)<=1/13.                                             (2)
```

Thus every one of the six named coordinates is genuinely replaced, and the
theorem applies both to a primitive tight proper packet and to any hypothetical
primitive proper counterexample below the conjectured floor.  Since
every member of (1) would be divisible by thirteen if `13|c`, primitivity
implies `13` does not divide `c`.

For each replacement put

```text
g_i=gcd(c,w_i),       D_i=c/g_i,       u_i=w_i/g_i.     (3)
```

Then

```text
gcd(u_i,D_i)=1,             u_i=D_i r_i (mod 13).       (4)
```

The `D_i` are the effective orders and literal ramification degrees of
THM-859.

## Theorem

### A. Hereditary lcm and the six-colour relative cuts

Every five-order subword already carries the full common scale:

```text
lcm(D_j:j!=i)=c                              for every i. (5)
```

In particular `lcm_i D_i=c`, and every maximal prime power is carried by at
least two colours.

For every nonempty colour set `S`, put

```text
T={1,...,6} minus S,
L_T=lcm(D_j:j in T),                 lcm(empty)=1,
m_i=D_i/gcd(D_i,L_T),                i in S,
rho(m)=ceil(2m/13)/m.                                (6)
```

Then the complete relative-cut family holds:

```text
sum_(i in S) rho(m_i)>=1                 for every S!=empty. (7)
```

For every `S` with `|S|>=2` there is a label/unit-sensitive strengthening
which will matter in the remaining finite bank.  In the complement-missed
sheet fibre `F` used to prove (7), put `A_i=E_i(o) intersect F` for `i in S`.
Then

```text
sum_i |A_i|/|F| - max_(spanning trees tau on S)
  sum_(ij in tau)|A_i intersect A_j|/|F| >=1.           (7a)
```

This **fibre-Hunter cut** is exact for the actual affine sheet masks.  The
scalar cut (7) is its relaxation obtained by upper-bounding each singleton by
`rho(m_i)` and deleting the nonnegative tree-overlap term.

The singleton quantifier in (7) is not a formal six-colour copy of THM-858.
Its five-colour complement need not be known to miss a sheet.  It is instead
supplied by (5), which makes `m_i=1`.  For `|S|>=2`, the complement has at
most four colours and THM-858's missed-fibre argument does transport.

Six-colour attenuation gives the first pivot, and the all-six instance of (7)
sharpens it:

```text
sum_i 1/D_i>1/13,                         min_i D_i<=77,
sum_i rho(D_i)>=1,                        min_i D_i<=72. (8)
```

### B. Optimized valuation ranges and the common-scale cap

The full upper-set cuts in (7), not merely separate adjacent-gap maxima, give

```text
range_i v_2(D_i)<=5,
range_i v_3(D_i)<=3,
range_i v_5(D_i)<=1,
range_i v_7(D_i)<=1,
range_i v_p(D_i)=0             for p>=11, p!=13.        (9)
```

Choose a numerically least order `D_*`.  Since `c=lcm_i D_i`, (9) implies

```text
c/D_* divides 2^5*3^3*5*7=30,240.                     (10)
```

Together with the sharpened pivot in (8), this first gives the coarse product
bound

```text
                         c<=2,177,280.
```

Jointly imposing all `63` cuts is far stronger.  After dividing the six
orders by their common gcd, the complete exact arithmetic relaxation has
`8,449` normalized order words and `18,405` admissible common multipliers.
Its sharp order-only cap is `1,120`, attained by exactly four sorted words:

```text
(20,35,56,80,1120,1120),      (20,35,80,224,280,1120),
(20,35,80,280,1120,1120),     (20,35,112,160,280,1120).
```

None of these four words, the twelve arithmetic words at scale `1,008`, or the
twenty-three words at scale `882` admits even the unit-independent scalar owner
capacity required by common-sheet coverage.  The next arithmetic scale is
`840`.  Consequently every packet in the theorem satisfies the much smaller
exact bound

```text
                              c<=840.                   (11)
```

The last line of (9) is a **valuation-range** statement, not a prime
exclusion.  A prime at least eleven may occur to the same exponent in all six
orders and is then already contained in `D_*`.  For example the abstract
word `(11,11,11,11,11,11)` obeys (5) and every arithmetic cut (7): every
proper-complement relative order is one, while the all-colour cut is
`6rho(11)=12/11`.  No `{2,3,5,7}`-smoothness claim is made for the individual
orders.

### C. Every fixed ramification language has a finite exact metric tree

Fix `c`, `R`, the orders `D_i`, and units `e_i in (Z/D_i Z)^*`.  Let `u_i^0`
be the least positive CRT representative satisfying

```text
u_i^0=D_i r_i (mod 13),              u_i^0=e_i (mod D_i). (12)
```

Every positive metric lift in this language is exactly

```text
w_i=(c/D_i)u_i^0+13c k_i,                   k_i>=0,     (13)
```

with the unchanged zero-height speed omitted when the replacement is
required to be proper.

Numerically order the six replacements.  After `j` insertions, let `L_j` be
the longest component of the literal strict-safe set and let `r=6-j`.  Every
tight continuation satisfies THM-815's exact cap

```text
x_(j+1)<=floor(22r/[13(13-2r)L_j]).                    (14)
```

For `j<6` the prefix has `6+j<=11` speeds, so the settled lower-runner theorem
makes its strict-safe set nonempty at threshold `1/13`.  Hence `L_j>0`, and
(14) leaves finitely many members of each progression (13).  Literal residual
components, the remaining labelled progressions, and the last speed are
insertion-Markov, so induction gives a finite exact decision tree for every
fixed presentation.

At the root, dilation gives

```text
L_0(cP)=L_0(P)/c.                                      (15)
```

The exact census of all `C(12,6)=924` cores has the unique minimum

```text
L_0(P)=31/1430             at P={1,4,7,9,10,11}.        (16)
```

Thus the largest normalized real cap is

```text
B=132/[13L_0(P)]=14520/31.                             (17)
```

On a fixed ray, (13) and (17) give

```text
u_i^0/D_i+13k_i<=B,
k_i<B/13=14520/403=36.029... .                         (18)
```

There are therefore at most `37` legal first heights **per labelled ray**,
or at most `222` crude first children over all six labels.  At `c=2`, the
base is at least `1/2>12/31=B-468`, so the sharper per-ray count is `36`.

Equations (8), (11)--(14), and the finite label/unit alphabets prove that the
entire primitive AP-centred Hamming-six branch is finite-decidable.  At
`c=1`, all orders are one and THM-857 leaves only the nonprimitive equality
`2[12]`; hence any primitive row satisfying (2) has

```text
2<=c<=840.                                             (19)
```

This is a finite reduction, not a zero-row verdict.

### D. The first ramified sheet and its signed cycle

At `c=2`, every order is one or two.  The set of order-two colours is
nonempty: if all six orders were one, all replacements and all six core speeds
would be even, contradicting primitivity.  An order-one colour fills both sheets
at its own owner and none at a different owner.  An order-two colour fills
one sheet at its own owner; a distinct provider `r` fills its other sheet at
owner `o` exactly when

```text
r -> o       iff       o/r in {2,-2} in F_13.           (20)
```

Every order-two owner needs an incoming edge from another order-two colour.
Such a finite positive-indegree digraph contains a directed cycle.  A cycle
of length `m` would give

```text
(+/-2)^m=1 (mod 13).                                   (21)
```

For `2<=m<=5`, the possible unsigned residues are respectively
`+/-4,+/-8,+/-3,+/-6`, never one.  Since the order-two set is nonempty, its
cycle must therefore have length six. Thus all six colours have order two and
the cycle is Hamiltonian.  Since `2^6=-1 mod 13`, it closes exactly when the
six edge signs contain an odd number of negatives.  Any extra chord would
create a cycle of length at most five, so the signed cycle is unique.

There are `12*2^5` rooted odd-sign cycles and six roots per label set.  Hence
there are exactly

```text
12*32/6=64                                             (22)
```

common-sheet label sets.  The direct affine-mask census agrees, with six
multiplicative orbits of sizes

```text
12,12,12,4,12,12.                                      (23)
```

Every survivor has one SCC of size six, six directed edges, one directed
cycle, and six Hamiltonian paths.  The negative-edge histogram is

```text
{1:12,3:40,5:12}.                                      (24)
```

This relation is a sparse signed SCC, not a tournament.  The pair observable
is the signed ratio `o/r in {+2,-2}` and the switch exchanges the two signs.
If absent pairs are forcibly completed, increasing label supplies a tie
Hamiltonian path, but that completion destroys the sheet-cover predicate.
The signed cycle is the useful pair object; the metric verdict remains in the
component/progression incidence state.

THM-861 subsequently evaluates that complete incidence state over all 64
roots and every unbounded step-26 height. Its unique cover is the ordinary AP
`[12]`; hence the sporadic `c=2` face is empty.

### E. Common sheets do not imply common scale

The following `c=3` presentation is an exact guardrail:

```text
R={1,2,3,5,8,12},
(w_1,w_2,w_3,w_5,w_8,w_12)=(16,45,48,28,37,10),
A={10,12,16,18,21,27,28,30,33,37,45,48}.               (25)
```

The orders are three on `{1,5,8,12}` and one on `{2,3}`.  The quartet uses
the all-one opposite-pair unit word.  In provider order
`(1,5,8,12,2,3)`, the exact masks are

```text
owner 1 : 2,1,0,-,-,-
owner 2 : -,0,1,-,012,-
owner 3 : 1,-,-,2,-,012
owner 5 : 2,1,-,0,-,-
owner 8 : 0,-,1,2,-,-
owner 12: -,1,2,0,-,-.                                 (26)
```

Every owner covers all three sheets, although the packet has no common
factor three.  Exact breakpoint evaluation gives

```text
M(A)=5/29,                  witnesses t=1/58,57/58.     (27)
```

Thus it is loose, not an LRC counterexample, while decisively refuting the
shortcut “common sheets imply common scale.”

## Proof

### 1. Leave-one-out primitivity and (5)

Fix `i` and put `A_i=A minus {w_i}`.  The settled eleven-speed theorem gives
`M(A_i)>=1/12>1/13`.  Apply THM-765(B) at `L=1/13`: assumption (2) forces
`gcd(A_i)|w_i`.  That gcd then divides all of `A`, so primitivity makes

```text
gcd(A_i)=1.                                             (28)
```

Let `L_i=lcm(D_j:j!=i)` and `q=c/L_i`.  Since `D_j|L_i`,

```text
q divides c/D_j=gcd(c,w_j)          for every j!=i.     (29)
```

It also divides every core speed in `cP`; hence `q|gcd(A_i)=1`.  Therefore
`L_i=c`, proving (5).

### 2. Common-sheet fibres and (7)--(8)

The oriented core-safe germ argument of THM-810/823 applies to (2).  With
`u_i` reduced modulo `13D_i`, let `z_i(ell,o)` be the centred representative
of `u_i(o^(-1)+13ell)` modulo `13D_i`.  Colour `i` covers at owner `o` the
sheets

```text
E_i(o)={ell in Z/cZ:
 -D_i<z_i(ell,o)<=D_i}.                                  (30)
```

Their union is `Z/cZ` at every owner.

If `S={i}`, (5) gives `L_T=c` and `m_i=1`, so (7) is equality.  Suppose now
`|S|>=2`, hence `|T|<=4`.  Augment `T` to four colours.  If those four fail
scalar coverage at one of their own owners, `T` has a strict deficit there.
Otherwise THM-810 makes them either all order one or the order-three coset.
At either of the two external owners their total contribution is respectively
zero or at most `2/3`.  Thus in every case there is an owner at which `T`
misses at least one residue modulo `L_T`.

Restrict to the corresponding fibre in `Z/cZ`.  As in THM-858, colour `i`
samples a gcd coset of size `m_i` inside a half-open interval of length
`2m_i/13`; it covers at most `ceil(2m_i/13)` points.  The union bound on the
fibre proves (7).

For completeness, (7a) is pointwise.  At a sheet belonging to exactly `k>=1`
of the sets `A_i`, its contribution to the singleton sum is `k`.  A spanning
tree on `S` induces a forest on those `k` active vertices and hence contributes
at most `k-1` active edges.  The difference is at least one.  Sum over `F`,
then choose the maximum-overlap spanning tree.  This proves (7a) and makes the
connection to THM-856's Hunter/Kruskal functional exact; no metric interval
claim is being imported.

For each owner, THM-823 attenuation has

```text
C_(D_i)(r_i,o)-2/13<1/D_i.                              (31)
```

Six-colour scalar coverage and (31) give

```text
1/13<=sum_i(C_i-2/13)<sum_i1/D_i.                       (32)
```

If every order were at least `78`, the last sum would be at most `6/78=1/13`,
a contradiction, proving the first pivot in (8).

For the sharper pivot, take `S` to be all six colours in (7), so `T` is
empty and `m_i=D_i`.  For every `m>=73` with `13` not dividing `m`,

```text
rho(m)<1/6.                                             (32a)
```

Indeed `m=73,...,77` is a direct ceiling check, `m=78` is excluded by
`13` not dividing `c`, and for `m>=79`,

```text
rho(m)<2/13+1/m<2/13+1/78=1/6.
```

If all six orders were at least `73`, their all-colour cut would therefore
have sum strictly below one.  Hence `min_i D_i<=72`; note that
`rho(72)=1/6`, so this is the exact threshold supplied by that one cut.  This
completes (8).

### 3. Full upper-set optimization

For positive integers `q,k`,

```text
rho(qk)=ceil(k(2q/13))/(qk)<=ceil(2q/13)/q=rho(q).       (33)
```

Fix `p!=13`, sort `a_i=v_p(D_i)`, and consider a boundary after `a_k`.
Take `S` to be the strict upper set.  The complement lcm has valuation
`a_k`; hence each relative order in `S` is divisible by
`p^(a_i-a_k)`.  Equations (7) and (33) imply

```text
sum_(i>k) rho(p^(a_i-a_k))>=1.                          (34)
```

The maximum valuation occurs at least twice by (5).

For `p=2`, adjacent boundaries of upper sizes `5,4,3,2` have gaps at most
`3,3,1,1`, using

```text
rho(2)=1/2, rho(4)=rho(8)=1/4, rho(16)=3/16.            (35)
```

Write the four gaps as `x_5,x_4,x_3,x_2`.  If their sum were at least six,
then `x_5>0`.  When `x_5=1`, all remaining gaps are maximal and the size-four
cut is

```text
rho(8)+rho(16)+2rho(32)=3/4.                            (36)
```

When `x_5=2`, the five relative exponents dominate `(2,4,5,6,6)`; when
`x_5=3`, they dominate `(3,4,5,6,6)`.  Monotonicity from (33) bounds both
cut sums by

```text
1/4+3/16+5/32+2(5/32)=29/32.                           (37)
```

Thus the dyadic range is at most five.

For `p=3`, upper sizes `5,4,3` have gaps at most `2,1,1`.  Range four would
force normalized valuations `(0,2,3,4,4,4)`, whose first cut is

```text
rho(9)+rho(27)+3rho(81)=2/9+5/27+39/81=8/9.            (38)
```

For `p=5`, at most five upper colours and
`rho(5)=1/5`, `rho(25)=4/25` force a single gap of one.  For `p=7`, only
upper sizes five and four can occur, each with gap one.  Range two would
force `(0,1,2,2,2,2)`, but

```text
rho(7)+4rho(49)=2/7+32/49=46/49.                       (39)
```

Finally, if `p>=11`, `p!=13`, and `p|m`, then `rho(m)<1/5`.  For `m>=22`
this follows from `rho(m)<2/13+1/22<1/5`; the smaller possibilities
`11,17,19` are direct.  A nonconstant boundary has at most five upper
colours and therefore violates (7).  This proves (9).

For every prime,

```text
v_p(c/D_*)=max_i v_p(D_i)-v_p(D_*)
           <=range_i v_p(D_i).                          (40)
```

There is no thirteen-adic factor.  Multiplying (40) gives (10), and the
all-six refinement in (8) gives

```text
c<=72*30,240=2,177,280,
```

which is the preliminary product bound above.

The exact joint optimization is still finite and much smaller.  Put

```text
g=gcd(D_1,...,D_6),                    D_i=g d_i.        (40a)
```

Then `gcd(d_1,...,d_6)=1`.  Every prime at least eleven has valuation range
zero, so it cancels into `g`; the normalized orders are supported on
`{2,3,5,7}` within the ranges in (9).  For every proper colour set `S`, the
common factor cancels exactly:

```text
D_i/gcd(D_i,lcm_(j notin S)D_j)
 =d_i/gcd(d_i,lcm_(j notin S)d_j).                     (40b)
```

Only the all-six cut changes with `g`, while (8) leaves
`g min_i d_i<=72`.  Thus it is enough to enumerate the normalized valuation
profiles and then

```text
1<=g<=floor(72/min_i d_i),                  13 does not divide g. (40c)
```

For each prime, fix the dyadic profile by simultaneous colour relabelling and
align every distinct permutation of the later profiles.  If a partial-prime
word fails a cut it can be discarded permanently: adjoining another prime
multiplies each relative order by some integer, and (33) says its individual
`rho` contribution cannot increase.  The complete exact census is

| stage | raw distinct words | all-cut survivors |
|---|---:|---:|
| admissible `2` profiles | 38 | 38 |
| admissible `3` profiles | 10 | 10 |
| admissible `5` profiles | 2 | 2 |
| admissible `7` profiles | 3 | 3 |
| aligned through `2,3` | 3,509 | 1,221 |
| aligned through `2,3,5` | 6,977 | 2,913 |
| aligned through `2,3,5,7` | 47,915 | 8,449 |

The first alignment represents all `14,364` dyadic/ternary profile
alignments.  Applying (40c) gives `18,405` feasible `(d,g)` pairs.  Their
largest common scale is `1,120`, and the four words displayed before (11) are
all the maximizers.  Their all-six cut sums, in the same order, are

```text
281/280,       561/560,       1,       1.              (40d)
```

It remains legitimate to use label information, because common-sheet coverage
was already forced in the proof of (7).  For an order `D`, provider label `r`,
and owner `o`, define

```text
N_D(r,o)=#{z in Z:-D<z<=D and z=D r o^(-1) (mod 13)}.  (40e)
```

At common scale `c`, that provider's sheet mask has cardinality
`(c/D)N_D(r,o)`.  Hence the union bound forces, at every owner,

```text
sum_i (c/D_i)N_(D_i)(r_i,o)>=c.                        (40f)
```

Capacities depend only on label ratios.  Normalize the label on the first
order slot to one and exhaust the remaining `11P5=55,440` ordered distinct
labels.  For the four scale-`1,120` words, the largest possible minimum of the
left side of (40f) is respectively

```text
1046,       1047,       1045,       1049,              (40g)
```

all strictly below `1,120`.  Thus none has a common-sheet presentation.  The
twelve arithmetic words at scale `1,008` likewise have zero scalar-covering
assignments among `12*55,440=665,280` normalized label rows.  Their best
minimum owner capacity is only

```text
946/1008,                                               (40h)
```

with deficit `62`.  The next scale in the complete arithmetic census is
`882`.  Its twenty-three words likewise have zero scalar-covering assignments
among `23*55,440=1,275,120` normalized label rows, with best minimum owner
capacity

```text
817/882.                                                (40i)
```

The next arithmetic scale is `840`, proving (11).  Its `372` order rows have
not been scalar-owner scanned.

### 4. Exact finite verifier

The stored verifier independently performs four exact jobs:

1. enumerates every normalized sorted six-valuation word inside the coarse
   adjacent ranges and applies all upper cuts, reproducing the optimized
   ranges `5,3,1,1,0` and their cut-level extremal words;
2. reconstructs closed danger unions for all `924` six-cores, proving
   (16)--(18);
3. enumerates all primitive-compatible non-all-one `c=2` order/label rows by
   direct affine sheet masks and reproduces (22)--(24); and
4. replays every mask in (26) and computes (27) from exact breakpoint
   denominators.

The joint-cap verifier independently enumerates all `63` cuts, every normalized
prime-profile alignment and common multiplier in (40c), the four extremal
words, all `39*55,440=2,162,160` normalized label assignments at scales
`1,120`, `1,008`, and `882` in (40f), and the pairwise
Tournament Analysis loss audit.  Its two completed tournaments are transitive
on every extremal word.  The raw-order gauge has zero or one tie; after
conditioning a pair on the other four orders, all fifteen pairs tie.  The
fixed slot-order tie path gives score histogram `{0:1,...,5:1}`, zero directed
triangles, six singleton SCCs, and one Hamiltonian path, while flipping
fourteen or fifteen raw edges.  This is positive evidence that the decisive
object is the `63`-hyperedge cut system, not a completed order tournament.

Normal and optimized Python runs are byte-identical to the stored output.
The frozen integrity data are

```text
source SHA-256  d91e8e9be79b0339e80958ac3cedd1a6cc3efd7ccc25009673d0592978069a95
output SHA-256  15626665176d97acbed7e68b197d8c274c457b6965387c5e727705177e23245d
payload FNV-64  4ef6d591ed02b211
root FNV-64     321c8da774c59420
c=2 FNV-64      b2aad8bd497eb595.                       (41)
```

The additional joint certificate has

```text
joint source SHA-256  21d55646c21c4530e9fca1f00dec5c422aec4edc8993f9d3571e881f145e5b4f
joint output SHA-256  8d9968402a16e247c6ee6024b1dc573b69b4e3380faa2f48cbdc7393965d5ccb. (41a)
```

Reproduce with

```bash
python3 04-computation/lrc13_hamming_six_primitive_ramification_verifier_codex_S15.py |
  cmp - 05-knowledge/results/lrc13_hamming_six_primitive_ramification_verifier_codex_S15.out
python3 -O 04-computation/lrc13_hamming_six_primitive_ramification_verifier_codex_S15.py |
  cmp - 05-knowledge/results/lrc13_hamming_six_primitive_ramification_verifier_codex_S15.out
python3 04-computation/lrc13_hamming_six_joint_ramification_cap_codex_S10.py |
  cmp - 05-knowledge/results/lrc13_hamming_six_joint_ramification_cap_codex_S10.out
python3 -O 04-computation/lrc13_hamming_six_joint_ramification_cap_codex_S10.py |
  cmp - 05-knowledge/results/lrc13_hamming_six_joint_ramification_cap_codex_S10.out
```

This completes the proof. ∎

## Faithful carrier and scope guardrail

The useful object is the four-layer incidence system

```text
prime-power valuation upper sets
 -> common-sheet owner fibres
 -> labelled CRT progressions
 -> literal strict-safe components.                    (42)
```

The first two layers prove finiteness; the last two decide the metric fibre.
Bare runners, pairwise order comparisons, the root cap, or a completed
tournament destroy at least one required operation.  The `c=2` signed cycle
is informative precisely because it retains edge parity, but even there it
does not decide component coverage at arbitrary heights.

THM-860 makes every primitive proper AP-centred H6 language finite-decidable.
THM-861 subsequently proves that the only `c=2` cover is the ordinary AP
presentation `[12]`, and THM-862 classifies all `1,504` common-sheet contexts
at `c=3` without running their terminal metric recursion.  This theorem does
**not** close those `c=3` trees, the remaining `4<=c<=840` trees, the finite
ramified H5 metric bank, the seven-comb wall, or global `n=12`
sporadic-branch emptiness.
