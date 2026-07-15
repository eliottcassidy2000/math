---
id: HYP-6820
title: Uniformity audit for the LRC(14) q<=25 good-period claim and the n=12 sporadic branch
status: PARTIALLY RESOLVED — uniform q<=25 is DISPROVED; n=12 is uniformly finite and sheet-stratified, with the full AP-centred Hamming-one through Hamming-four stars uniformly loose, the scale-one Hamming-five chart reduced to two joint finite boxes with its complete height-at-most-two slice loose, scale-one Hamming six recursively finite, seven the first discrepancy-unbounded scale-one radius, and the arbitrary-scale five-deck interface reduced to a least-order pivot plus a bounded common-sheet survivor orbit; the two-sheet branch has finite bases, raw/parity support and anti-grid gates, an exact adaptive max-speed-cell/component selector, a pair-indexed uniform atomic-stalk factorization, an all-size two-radius factorization of the full symmetric return union at ratio `(13d,5d)`, and a sharp linear-satellite method limit; THM-829 gives related contragredient arithmetic-owner transport but is not yet joined to metric endpoint owners; branch emptiness remains OPEN in the Hamming-five boxes, unbounded common-sheet and survivor-strictness problems, all-scale descent, uniform violation or transport of the exact radius/sum-arc selectors, other odd ratios, the deep dyadic/collar residuals, and higher sheets
source: codex-2026-07-14-S3
renumber_note: reserved as HYP-6810 by codex-S3, which collided with opus-S298's earlier-pushed
  HYP-6810 claim (the assembly write-up); renumbered to HYP-6820 by opus-2026-07-14-S299 per the
  first-pusher protocol, as codex-S2's atlas itself requested; subsequently updated by codex-S3.
depends_on:
  - THM-758
  - THM-759
  - THM-762
  - THM-763
  - THM-764
  - THM-765
  - THM-766
  - THM-768
  - THM-769
  - THM-770
  - THM-772
  - THM-774
  - THM-775
  - THM-776
  - THM-780
  - THM-782
  - THM-789
  - THM-795
  - THM-797
  - THM-800
  - THM-803
  - THM-804
  - THM-806
  - THM-807
  - THM-810
  - THM-815
  - THM-816
  - THM-817
  - THM-820
  - THM-821
  - THM-822
  - THM-823
  - THM-824
  - HYP-6750
  - HYP-6775
related:
  - HYP-6780
  - HYP-6785
  - HYP-6800
  - HYP-6815
  - THM-771
  - THM-829
  - MISTAKE-143
---

# HYP-6820 — uniformity audit (formerly HYP-6810; renumbered after collision)

This audit began with two concrete proof obligations requested together:

1. determine whether every stated residual family has a rational lonely witness
   of denominator `q <= 25`, and separate exact banks from uniform claims;
2. prove, or reduce without quantifier loss, that the super-lonely-core (sporadic)
   branch of primitive tight 12-speed families is empty.

They now have different outcomes.  The first statement is false.  The second
has acquired several uniform reductions, but the requested emptiness theorem is
still open.

## A. The `q<=25` claim is disproved

THM-762 (independently restated in THM-764) proves the exact criterion visible
in the range `15<=q<=28`.  For a
covering family `S`, a `q`-witness exists exactly when:

1. no speed is divisible by `q`; and
2. the signed unit-pair deck
   `B_q(S)={[s] modulo +/-1 : gcd(s,q)=1}` misses a class of
   `(Z/qZ)^x/{+1,-1}`.

The proof is elementary: a nonunit multiplier reduces to a covered denominator
at most 14; for a unit multiplier and `q<=28`, the only strictly unsafe integer
residue distances are `0` and `1`.

Two primitive covering exact residuals refute uniform `q<=25`:

```
26*{1,...,12} union {339}                         first witness 2/27,
{81,91,131,151,157,196,258,274,313,328,330,339,348} first witness 3/26.
```

The second has no prime common to seven speeds, every leave-one-out gcd is one,
and exact `M=101/470`; coherence or near-tightness therefore cannot rescue the
claim.  Both pass the rational analogue of S312's band-residual predicate.
MISTAKE-143 corrects the universal wording.

The historical banks now have exact scope:

- S312: `120/120` sampled rows with `q<=25`, not an exhaustive class;
- S105: exact replay of its capped interval-core generator gives `91/8260`
  rows with no `q<=25` witness, but every row witnesses by `q<=38`.  The unique
  least-denominator-38 row is `{1,...,11,338,420}`, first witnessed by `3/38`.

THM-566 remains the global guardrail: no fixed raw denominator ladder can cover
all primitive covering families.

## B. The `n=12` branch is finite and tooth-stratified, not proved empty

Let `A={a_1<...<a_12}` be primitive and tight, `M(A)=1/13`, and let
`P=A\{a_12}`.  The sporadic branch is `M(P)>1/12`.

THM-763 gives the first quantifier-honest global finite bound:

```
sum A <= 78^11 = 650190514836423555072.
```

Thus uniform emptiness is a finite decision problem, conditional only on the
settled lower-dimensional LRC input, but a naive enumeration is infeasible.

THM-766 supplies a scale-normal reduction.  For a tight `n`-set, with
`m=a_1`, `b=a_(n-1)`, and `w=a_n`,

```
w >= n m,
b/m >= n^2/(n+2).
```

If `b<nm`, the first core-safe interval must fit one `w`-danger tooth, so for
some `1<=k<=n-1`,

```
((n+1)k-1)m <= w,
n w <= ((n+1)k+1)b.
```

At `n=12`, every candidate therefore has `a_11/a_1>=72/7`; the sub-12 span
lies in only eleven explicit projective cones.  At an exact core maximum
`t_0`, the stronger band

```
||a_12 t_0|| <= 1/13 - a_12(M(P)-1/13)/a_11
```

retains the residue alignment lost by THM-759's scalar ratio cap.

The faithful completion predicate is

```
G_P={t:min_(p in P)||pt||>1/13} subseteq
D_w={t:||wt||<=1/13}.
```

Every connected component of `G_P` must fit one `w`-tooth simultaneously.
THM-765 strengthens this component language with a gcd-deck theorem: every
leave-one-out core of a primitive tight twelve-set is itself primitive.  It
therefore removes all imprimitive-core sporadic candidates at every height,
not just in a box.
Total measure alone cannot decide this: for the nonextremal core
`P={1,...,10,12}`, exact `|G_P|=461/8190<2/13=|D_w|`.

An exact regression on the historical 77 nonextremal cores in `{1,...,13}`
uses the true THM-759 cap rather than an arbitrary killer window.  It leaves
790 primitive completions: 40 narrow-core and 750 wide-core.  THM-766's cone
test eliminates all 40 narrow candidates; exact pair-sum/difference/half-turn
evaluation finds zero tight completions among all 790, with bank minimum
`M=1/12`.  This is finite-exact for that box only.

The exact max-peel tooth atlas gives a further guardrail.  It exhausts all
primitive twelve-subsets `A` of `{1,...,20}` with `w=max(A)`,
`M(A\{w})>1/12`, and `M(A)<=1/10`, plus eleven tight AP deletion controls.
The component criterion agrees with an independent exact-maximum calculation
on all `2,464` rows.  Every one of the `2,453` escaping rows nevertheless has
cyclic nearest-tooth winding one, and `1,972` have pure endpoint owners.  The
first liar to “pure endpoints plus winding one implies cover” is
`{1,...,10,12,13}`; even requiring a transitive component-phase tournament
with one Hamiltonian path fails for `{1,...,9,11,12,15}`.  Their minimum
slacks are respectively `-11/26` and `-43/91`.  Thus winding is a checksum,
not the missing separator: the width term in
`sigma_J(w)=1/13-||wc_J||-w h_J` must survive every quotient.

## C. Binding-scale recursion: what the residue picture was missing

THM-768 eliminates one tempting deep configuration: if the largest speed is
the unique multiple of 13, an explicit perturbation of a shallow prime-grid
point makes all twelve runners strictly safer than `1/13`.  Thus a tight set's
largest speed cannot be its unique 13-multiple.

THM-769 gives the scale-normal statement at **every** rational global maximum.
Write a reduced maximizer as `p/(13s)`.  Its multiplied residues lie in the
packet `[s,12s]`, both endpoints occur, and endpoint owners are divisible by
`s`.  With
```text
E={v:s|v}=sU,                  F=A\E,
```

the familiar complete nonzero residue system is exactly the shallow `s=1`
branch.  If a 13-multiple blocks all shallow points, every maximizer is deep
and `F` must cover all `s` lifts of every point in
`G_U={tau:phi_U(tau)>1/13}`.  Putting

```text
D_w=s/gcd(w,s),
```

gives the necessary capacity inequality

```text
sum_(w in F) (floor(2D_w/13)+1)/D_w >= 1.
```

In particular `|F|>=2`.  Equality at two exceptions forces `s=2`, ten even
speeds and two odd speeds; tightness is then **equivalent** to persistent
opposite nearest-integer parity over all of `G_U`.  At three exceptions,
either a half-sheet tightener occurs or the equality edge is `s=3`, with nine
multiples of three and three nonmultiples persistently owning all three sheet
colours.  For `r=|F|<=6`, some exception satisfies
`D_w<=13r/(13-2r)`.  These are uniform reductions, not empirical patterns.

THM-772 makes the equality packets recursive rather than merely sparse.  In
the two-sheet branch

```text
A=2U union {x,y},       |U|=10,       x,y odd,
```

the quotient `U` is primitive, contains a multiple of every modulus
`2,...,12`, contains no multiple of `13`, and both exceptions are at most
`11 max(U)` (with a sharper one-exception and determinant tax).  The
three-sheet equality edge has the analogous primitive divisor transfer through
modulus `11`, an exact modulus-36 residual when `12` is missing, and the bound
`x,y,z<=10 max(U)`.  These conclusions use the simultaneous unit-fraction by
sheet-colour obligation hypergraph; no tournament on the moduli alone sees the
denominator-6/12 splice.

THM-774 identifies the exact metric object in the two-sheet branch.  With

```text
a=(x+y)/2,       b=(x-y)/2,
```

eligibility and opposite colour are together exactly

```text
||a tau||+||b tau||>=11/13.                               (C1)
```

Equivalently, if `B_(x,y)(tau)` is the larger of the two lift-clearance
minima, then

```text
B_(x,y)(tau)=(1-||a tau||-||b tau||)/2.                   (C2)
```

The folded diamond in (C1) has an exact half-grid formula and sharp universal
measure at most `8/117`, attained at reduced ratio `x:y=9:1`.  Hence tightness
forces

```text
G_U subset {tau:||a tau||+||b tau||>=11/13},
measure(G_U)<=8/117,
||a tau||,||b tau||>=9/26 on G_U.                         (C3)
```

The scalar cap does not close the branch: several explicit ten-speed loose
cores already have measure below `8/117`.  The theorem-facing datum is the
pointwise containment with component locations and widths.

The first complete folded slice is nevertheless now exact.  For every one of
the `binomial(19,10)=92,378` cores `U subset {1,...,19}`, the widest component
of `G_U` bounds any individually eligible odd runner by

```text
w<=floor(4/(13 ell_max(U))).
```

Testing all `767,700` intrinsically permitted `(U,w)` incidences finds no
single odd runner covering all of `G_U`, hence no two-runner colour cover.
This closes the `s=2`, `max(U)<=19` slice with odd exceptions unbounded.  It
does not bound `max(U)` globally.

THM-774 also makes the maximum deletion recursive.  If the largest speed is
odd, the sporadic inequality follows automatically from
`M(2U union {z})>=1/11`.  If it is `2R` with `R=max(U)`, put
`U^-=U\{R}` and `Q=||(x+y)tau/2||+||(x-y)tau/2||`.  The exact remaining
nine-core collar is

```text
E_(1/13)(U^-) subset D_R union {Q>=11/13},
E_(1/12)(U^-) intersect D_R intersect {Q<5/6} nonempty.
```

This isolates where the sporadic hypothesis genuinely adds information.

THM-775 proves the arithmetic descent that was previously only proposed. If a
one-element deletion of the two-sheet quotient is imprimitive, its gcd is two
and the deleted element is the unique odd member:

```text
U=2V union {u}.
```

The first four lifts carry an exact `2+1+1` ownership partition, and iteration
gives a finite dyadic chain of primitive divisor-complete quotients with binary
safe-child fibers, ending at a hereditarily primitive base. This does not
exclude the chain or its terminal base.

THM-776 supplies a complementary finite-exact landing region. For every
10-even/2-odd twelve-set with raw height at most 100, the odd-pair failure locus
atomizes into an incidence hypergraph whose transversal number is 12, while the
packet has only ten quotient speeds. Thus every such set has `M>1/13`, without
assuming primitivity or divisor completeness. Above height 100, no uniform
transversal lower bound has yet been proved.

THM-770 settles a very large but still bounded part of the shallow branch.
For the labelled packets

```text
W(k)={r+13k_r:1<=r<=12},       0<=k_r<=12,
```

an exact unique-owner CSP represents all `13^12=23,298,085,122,481` rows
without literal enumeration.  It has exactly thirteen zero-defect leaves,
`c*{1,...,12}` for `c=1,...,12,14`, and only the `c=1` leaf is primitive.
Consequently a primitive shallow tight set with `max A<=168` is the AP.
Its maximum deletion core is `{1,...,11}`, with `M=1/12`, so the bounded
shallow **sporadic** slice is empty.
Dilation gives the exact defect law `chi_13(cW)=c chi_13(W)`, so gcd descent is
lossless.  This result does not supply an unbounded height descent.

THM-795 supplies one scale-free piece of that descent.  For every `r` and
every `k>=1`,

```text
([12]\{r}) union {r+13k}
```

has a strict `1/13` witness.  More generally, if one coordinate `cr` of
`c[12]` is replaced by `w=cr mod 13`, with `13` not dividing `c`, tightness
forces `w=cr`.  The proof uses the `c` lifted missing-owner splice germs: all
their replacement phases would have to lie in one arc of length `2/13`, which
forces `c/gcd(c,w)=1`, followed by seven exact low atoms.

THM-800 closes the two-coordinate continuation.  The oriented core-safe germ
at a missing owner uses the half-open phase arc `(-1/13,1/13]`; exact deck
capacity then forces both replacement orders to one, hence common scale.  At
scale one every proper double lift has the sharp floor `M>=2/25`.  Thus the
entire residue-preserving Hamming-one/two star is uniformly loose.

THM-804 proves the next descent. If
three labelled coordinates of `c[12]` are replaced and the resulting packet
is tight, then all three effective orders

```text
D_i=c/gcd(c,w_i)
```

equal one.  The proof combines the same half-open germ with exact own/complement
capacities, a two-colour sublemma, and a finite case split on order-two and
order-at-least-three decks. Hence every arbitrary-scale Hamming-three packet
descends to a genuine scale-one triple lift. THM-806 closes that base. A
universal `1/156` owner collar and the exact half-open handoff bands force one
replacement into `[14,24]`; the settled ten- and eleven-speed bounds plus
connected two-comb geometry reduce the remaining ordered replacements to
`v<=262`, `w<=12v`; a larger legacy superset replay rejects all `5,713,539`
canonical rows. Thus the full residue-preserving Hamming-three AP star is
uniformly loose.

THM-810 identifies the first live recursion instead of leaving “radius at
least four” unstructured.  Four-colour oriented capacity forces either all
effective deck orders to one, or all four orders to three with labels in a
coset of `<5>={1,5,8,12}`.  Exact sheet overlap in the exceptional case is
possible precisely when the mod-three signs agree on complementary label
pairs.  After removing the common gcd, this becomes an `s=3` deep packet with
eight on-sheet speeds and four exceptions.  Hence Hamming four is the first
literal shallow/deep gluing interface.

THM-815 and THM-816 now close both sides of that interface uniformly.  In the
common-scale branch, the sharp interval-comb discrepancy recursively bounds
the four numerical lifts by

```text
x<=105,       v<=118,       w<=83,       z<=50,
```

after which a 35,640-row component-containment census and an independent
35,640-row endpoint-cell replay leave no tight row.  In the order-three branch,
the same discrepancy on a dynamic residual interval union reduces every coset
and parity pattern to 7,909 exact recursive states, with all 34 terminal
residuals nonempty; an independent 132,510-row endpoint implementation agrees.
The eight `q=39` equality-clock points from THM-810 lie on the boundary of the
strict core-safe set, so they are compatible with a surviving open component
elsewhere.  Thus the entire AP-centred, residue-preserving Hamming-four star is
loose at arbitrary height and scale.

The common-scale conclusion has a genuinely independent proof companion.  Its
collar cycle first reduces every all-large row to the band word `(2,2,2,5)`, a
least lift at most `245`, and spread at most four, rejecting `626,962` rows.
Otherwise one lift lies in `[14,24]` and the other three obey a recursive
doubling box; all `141,773` such rows also fail component containment.  The
combined 768,735-row C++ certificate agrees with the discrepancy proof while
using a different unbounded-to-finite reduction.

The recursion has a further exact consequence.  With `m<=6` proper scale-one
lifts, every prefix before the last has at most eleven speeds and hence a
positive strict-safe component by settled `LRC(<=13)`.  If the remaining
`r=m-j` combs cover a longest component of length `L_j`, their next speed is at
most

```text
floor(22r/[13(13-2r)L_j]).
```

Therefore radii five and six are finite exact decision trees.  Their first
component scans give caps `146` and `468`.  At radius seven the coefficient
`13-2r` changes sign, so this specific recursion supplies no initial cap; this
is a method boundary, not a proof that the chart is intrinsically infinite.

THM-820 gives the first uniform radius-five reduction.  For a proper
scale-one five-lift packet with ordered replacements `x<v<w<y<z`, every
hypothetical tight row satisfies

```text
x<=388.
```

Either `v<=2x` and recursive collar obligations give

```text
(x,v,w,y,z)<=(388,776,1552,3104,6208),
```

or `v>2x`, the top four labels form THM-815's `a{1,2,4,8}` cycle, and the
seven- and eight-speed reciprocal bounds combine to give

```text
x<=228,       v<=1986,       max_top<=7944.
```

These are THM-820's standalone boxes.  The same row is also subject to
THM-815's independent component-discrepancy cap `x<=146`.  Their intersection,
which is the sharp current scale-one Hamming-five frontier, is therefore

```text
(x,v,w,y,z)<=(146,292,584,1168,2336)
```

in the doubling branch, or

```text
x<=146,       v<=1986,       max_top<=7944
```

in the exceptional branch.

The least-coordinate sharpening is itself structural: when all five owners
are large, positive indegree gives a minimal cycle of length four or five.
The exact five-cycle band bank has five feasible types and forces a ratio at
least `3/2`; the four-cycle forces a factor-two edge.  This converts the
collar classification into the stronger reciprocal cap.

The theorem does not empty these boxes.  Its exact shallow census rejects all
`C(12,5)=792` height-one rows and finds the unique minimum

```text
2[11] union {11},               M=1/12.
```

Lifting the sixth odd label produces `2[12]`, the unique height-one
radius-six minimum at `1/13`.  Thus radius five approaches a literal
scale-changing AP face rather than an unrelated sporadic pattern.  THM-820
also supplies the method obstruction: the top-four SCC can persist at
unbounded scale, while two labelled five-cycles with identical tournament
telemetry have different integer centres and exact maxima.  The finite-box
closure must act on exact residual interval unions and remaining danger teeth,
not the collar graph alone.

THM-822 empties the complete height-at-most-two slice, all
`C(12,5)2^5=25,344` packets, with the same unique `1/12` minimum.  It also
audits the proposed three-face quotient exactly.  The labelled live relation
`H0` and its integer-centre decoration `H1` have identical fibres in this
ratio band, and their `111,006` ordered kernel rows include `3,810` fibres
mixing exact `M`; one shared relation has maxima `1/4` and `12/37`.  Literal
strict-safe component endpoint words `H2` make the bounded three-face join
injective, while the final residual alone still has fourteen collision pairs.
This is a static bounded reconstruction theorem.  Since every audited row is
loose, it says nothing nonvacuous about tightness purity, and it does not prove
that `H2` transports under arbitrary-height insertion or scale change.  The
next quotient target is an endpoint-owner/component-incidence codec between
`H1` and `H2`, minimized under the actual erosion action.

THM-823 proves that the remaining all-scale deck interface is not a finite
scalar-order problem.  Every five-colour scalar cover satisfies

```text
sum_i 1/D_i>3/13,                    hence min_i D_i<=21,
```

but `(1,2,3,5,10):(2,5,13q+1,2,13q+1)` scalar-covers for every `q>=3`.
At `q=3` all 1,024 unit words fail common-sheet coverage, with best minimum
`29/40`; scalar capacity and simultaneous sheets have separated completely.
In the exhaustive no-order-one bank `2<=D_i<=12`, only five of 2,190 scalar
presentations survive the sheet test.  They form one all-order-three orbit

```text
R=C union {b},              C=a{1,5,8,12},              b in 2C,
```

with THM-810 opposite-pair parity on `C` and a free parity on `b`.  All forty
least CRT packets are loose, with minimum `2/17`; arbitrary lifts retain six
mod-39 boundary witnesses and are only proved to have `M>=1/13`.  The
order-one branch is symbolically either all one or one plus the order-three
quartet, again with only the boundary clock at arbitrary height.  Remaining
work is unbounded common-sheet classification and strict metric erosion of
these survivor languages, not a stronger scalar cutoff.

## D. The two-sheet branch is now a folded dyadic cover

THM-772 first turns the `s=2` equality packet

```text
A=2U union {x,y},       |U|=10,       x,y odd
```

into a recursive arithmetic object.  The quotient `U` is primitive, contains
a multiple of every `m=2,...,12`, contains no 13-multiple, and, with
`B=max(U)`, satisfies

```text
x,y<=11B,       min(x,y)<=11B-36,
1/(xy)+2(M(U)-1/13)/B <= 2/(13x)+2/(13y).
```

THM-774 then removes the nearest-integer colour notation without losing the
predicate.  For

```text
a=(x+y)/2,       b=(x-y)/2,
```

the two odd exceptions cover the two lifts exactly at

```text
H_(x,y)={tau:||a tau||+||b tau||>=11/13}.
```

Thus tightness is precisely `G_U subset H_(x,y)`.  The folded diamond has the
sharp universal measure cap `8/117`, attained at reduced ratio `x:y=9:1`.
This is a useful necessary condition, not a noncoverage theorem: a small loose
set could still sit inside the diamond.

The same theorem closes the entire low-core slice `max(U)<=19` with `x,y`
unbounded: the widest strict component of `G_U` gives an intrinsic cap on any
individually eligible odd speed, and all 767,700 permitted incidences fail.
This is genuinely stronger than the measure filter: 52 of the 3,400
divisor-complete cores already have `measure(G_U)<=8/117`, with minimum
`41/858`.  The finite region is complementary to THM-776's larger
bounded-height box.

THM-775 identifies the recursive gcd pathology exactly.  If a deletion of
`U` is imprimitive, its gcd is two and

```text
U=2V union {u odd}.
```

Over `G_V` the four lifts carry a saturated `2+1+1` ownership partition: `2u`
owns one parity class and `x,y` own the two singleton sheets of the other
parity.  The quotient `V` is again primitive and divisor-complete.  Repeating
the argument produces a finite dyadic tower ending at a hereditarily primitive,
divisor-complete quotient.  At later levels the exact object is a canonical
assigned ownership tree; raw danger incidences can overlap and must not be
called a literal disjoint partition.  A depth-`r` tower has the explicit
normal form `U=2^r Q_r union {2^i h_i:0<=i<r}` with every `h_i` odd, so its
2-adic valuation profile contains exactly one speed at each level below `r`.

THM-776 supplies a substantial finite base.  Every twelve-set with exactly ten
even and two odd speeds and `max(A)<=100` has `M(A)>1/13`, without assuming
primitivity or divisor-completeness.  The exact proof reverses the tuple
quantifiers: for each of the 1,225 odd pairs, the bad atoms form a hitting-set
hypergraph whose transversal number is exactly twelve, while a packet core has
only ten quotient speeds.  The 6,876-point grid and full folded-diamond atlas
were independently regenerated.  This empties the entire two-sheet branch
through height 100, but it does not prove that the transversal deficit persists
when larger quotient teeth are admitted.

THM-782 adds a genuinely scale-free metric stalk.  Every ten-speed quotient
core has a `1/11`-deep time `t_0` and an anchored joint-phase return packet
`B` with

```text
|B|>=72^(-10),       ||u b||<1/72,
||u(t_0+b)||>1/13+1/10296       (u in U, b in B).
```

It also forces a component of `G_U` of length at least
`72^(-10)/(20 max(U))`.  If the packet were tight, both the translated return
packet and that component would have to lie in `H_(x,y)`.  This repairs the
earlier missing uniform measure-floor statement without closing the branch:
the constant is far below `8/117`, and the exact core-19 atlas has the explicit
profile `|G_U|=41/858<8/117=|H_(1,9)|`.  The sharpened target is therefore
structured phase-packet noncontainment, not scalar measure comparison.

THM-789 strengthens and delimits that stalk. Replacing one anchored translate
by the symmetric difference set `A-A` gives

```text
measure(G_U)>=2*72^(-10),
ell_max(G_U)>=72^(-10)/(5 max(U)).
```

Under tightness every safe time obeys the pointwise thickness tax

```text
||wt||+(w/max(U))(phi_U(t)-1/13)<=2/13,
```

and the deep set satisfies the exact erosion

```text
E_U subset H_(x,y) minus R_U,
R_U={d:max_u||ud||<2/143}.
```

This still does not close the branch. For
`U={1,2,3,5,7,8,9,10,11,12}`, `t_0=4/17`, `(x,y)=(13,9)`, the full return
set is exactly `R_U=(-1/858,1/858)` and its translate is contained in the
diamond. Every same-cell symmetric packet and literal extra-coordinate
refinement is trapped there too. Yet `14/19` is `2/19`-deep and lies outside
the diamond. The thickness tax also gives `w<=11 max(U)` for each odd
exception. Thus the new target is not stronger local refinement at an
arbitrary anchor; it is the global selection statement
`E_U not subset H_(x,y) minus R_U` among all deep components.

THM-797 turns odd divisors of the exceptions into such global selectors.  If
odd `q|x`, every deep unit class `p/q` is accepted by the other exception
exactly when the balanced residue `[yp]_q` is odd and lies in the shell
`1<=|[yp]_q|<=floor(2q/13)`.  Any deep class outside this shell proves the
full erosion failure.  At the mandatory prime `q=13`, write
`C=(Z/13)^*/{+/-1}` and let `S(U)` be the folded core support.  Containment
forces

```text
S(U)=C or S(U)=C\{+/-y}
```

when exactly `x` is 13-divisible, and forces `S(U)=C` when both exceptions
are.  Hence every support-at-most-four core, every misaligned five-class core,
and every non-full double-13 core is uniformly eliminated.

The gate is sharp but the global component method is stronger.  The aligned
row

```text
U=(1,2,3,4,7,9,10,11,12,16),       (x,y)=(13,5)
```

traps its only deep thirteenths `5/13,8/13`, yet `7/33` is `1/11`-deep and
has folded value `8/33<11/13`.  The associated full twelve-set

```text
(2,4,5,6,8,13,14,18,20,22,24,32)
```

is a genuine sporadic max-peel row—its maximum deletion has `M=1/8>1/12`—
but it is loose with `M=2/19>1/13`.  THM-803 corrects the status of this and
the companion signed-wall row: both are sharp for the prime-13 gate, but both
already escape at the mandatory quarter anti-grid `11/52`.

THM-803 supplies the strengthened exact carrier.  If a 13-divisible exception
is present, the half-divisor grid forces full **parity-twisted support**
`P(U)=C` in addition to THM-797's raw signed support.  Exactly the even
denominators `d=2,4,6` make the 13-divisible exception universally ineligible,
so the mandatory anti-shell ladder is `26,52,78` and no longer.  For the full
erosion predicate put

```text
K_U=E_U+closure(R_U).
```

Containment is equivalent to `K_U subset H_(x,y)` and is decided exactly at
the component endpoints of `K_U` together with the cusps `k/(2a),k/(2b)` of
the folded diamond.  THM-817 decomposes `closure(R_U)` exactly into `N_R`
signed cells, each inside one maximum-speed tooth with labelled left/right
endpoint owners.  Substituting those cells into the same proof gives the
adaptive scale-normal bound

```text
|Sigma(U;x,y)|<=2c_E N_R+2W-2g
               <=20 max(U)^2+22 max(U)-2g.
```

Thus the continuum selection problem is finite for each core, without a raw
height bound.  The factor `N_R` cannot be replaced by a constant using the
present arithmetic/scalar gates: THM-817's primitive divisor-complete signed-
complement family has `N_R=3+1440n=Theta(max(U))`.  It also has an explicit
central erosion escape, so it is a method limit rather than a tight candidate.
The current mandatory regression is

```text
U_*=(2,4,6,7,9,10,11,12,14,16),       (x,y)=(13,5).
```

It satisfies the exact signed complement, full parity-twisted support, every
odd exception-divisor test, and the `26/52/78` anti-shell ladder, while the
nonmaximal singleton component `7/22` escapes.  Consequently the uniform proof
must control the complete owner-labelled selector; neither global maximizers
nor any fixed short denominator ladder can replace it.

THM-807 and THM-817 now identify the selector's topology and arithmetic
exactly.  Thickening only by the mandatory central return component gives a
linear selector of size at most `42B-2`; it is equivalent to the full erosion
test when the closed return set is connected.  In general, every component of
`closure(R_U)` is one signed max-speed cell `(k+I_k)/B`, with its endpoint
owners read directly from the defining interval intersection.  If `N_R` cells
survive and `c_E` deep components are present, the exact adaptive bound is

```text
|Sigma(U;x,y)|<=2c_E N_R+2W-2g<=20B^2+22B-2g.
```

This is the right representation but not a new compactness theorem.  The
primitive divisor-complete exact signed-complement family
`B_n=506+360360n` has `N_R=3+1440n`.  Thus primitivity, divisor pins,
raw/parity support, and the known scalar taxes cannot force connected,
bounded, or sublinear return geometry.  Any uniform proof must act on the
signed cell/deep-component incidence and its margins.

THM-821 now identifies that incidence at atomic resolution for every positive
odd exception pair `x>y`. For every deep component `C` and signed return cell
`R_k`, the pair-indexed circular sum state `((x,y),C+R_k)` determines the
selector event, exact margin, and success verdict; global containment holds
exactly when every such atom succeeds.  Its 9,974-atom audit has 492 successes
and 9,482 failures and gives
explicit liar fibres showing that a return cell alone, a deep interval alone,
owners alone, widths, or selector-event shape can mix verdicts.  Exact input
intervals or the exact sum arc do not mix, as the factorization proves.  The
linear-satellite family has a uniformly failing central stalk at `(13,5)`, but
neither that failure nor the finite mixed-fibre table is an all-core/all-pair
row separator.

THM-824 proves that assembling **all** signed return cells restores enough
symmetry for an exact all-size compression at one fixed ratio.  Put

```text
H_(13,5)={t:||9t||+||4t||>=11/13}
        =B(5/13,2/169) union B(8/13,2/169).
```

For every nonempty compact `E` and compact `R=-R` with `0 in R`, its no-switch
lemma gives

```text
E+R subset H_(13,5)
iff max_(t in E) min(||t-5/13||,||t-8/13||)
   +max_(r in R)||r|| <=2/169.                              (C4)
```

For the common-dilate pair `(13d,5d)`, (C4) is intrinsically
`rho_(C,d)(E)+rho_d(R)<=2/169`.  It is equivalent to the phase inequality
`max_E||13dt||+13max_R||dr||<=2/13` only together with `dE subset H`; the
unguarded phase statement fails already for `E=R={0}`.  Thus the phase
inequality alone is the necessary THM-789 thickness tax.  Exact replay on 12,159 synthetic compact
packets and the 214 THM-821 cores has zero direct/radius mismatches; every
core fails the global budget.  This is a global-conjunction factorization,
not a factorization of an individual nonsymmetric `R_k`.  It makes present
evaluation linear in the exact component/cell lists, but does not prove that
owner ancestry may be discarded under descent or that other odd pairs have
the same two-cell geometry.

THM-829 proves a complementary all-size transport identity over primitive
arithmetic witnesses: if `bv=1` and `A in GL_2(Z)`, then the transported stalk
is `(Av,bA^{-1})`, with reflection handled by conjugating `A`.  This establishes
why an inverse residue without its denominator/Bezout row cannot transport.
It does not yet identify those witness owners with the endpoint owners of the
signed max-speed cells, so it is related transport infrastructure rather than
a new closure step in this audit.

The next all-pair object is therefore the target-component symmetric-switch
hypergraph together with its deck and metric sidecars, not a bare runner
tournament.

These results reveal the faithful carrier more precisely than “ten even plus
two odd.”  It is a folded bad-atom/core-tooth incidence hypergraph equipped
with a binary safe-child map, divisor unit columns, every owner-labelled deep
component, and its symmetric return packet, selector cusps, and escape margins. The runner
tournament records which odd runner owns a sheet but forgets the simultaneous
hitting number; the dyadic quotient chain records recursion but needs binary
sheet-fiber, incidence, component, and phase-cell sidecars.

## E. Remaining proof obligation

The uniform theorem now has two explicit residuals:

1. **Finish the finite scale-one radius-five action, enumerate radius six,
   and extend descent.**
   THM-795/800/804/806 close the complete AP-centred Hamming-one through
   Hamming-three stars.  THM-810 splits radius four into common scale and one
   order-three coset interface; THM-815 and THM-816 close those alternatives
   independently.  THM-820 has now classified the radius-five collar cycles
   and anchored recursion and reduced the full scale-one chart to the two
   explicit boxes above.  Empty those boxes by the commuting idempotent action
   `T_u(E)=E intersect Safe(u)` on literal residual interval unions, preferably
   through a minimized continuation automaton or a two-plus-three
   meet-in-the-middle cover join.  Audit every proposed codec through its
   operation-indexed kernel pair; literal labelled four-coordinate faces glue
   uniquely, but component counts, unordered endpoints, and tournament nodes
   do not.  Then prove the arbitrary-AP-scale five-deck descent, allowing for
   new sheet interfaces rather than assuming common scale. THM-815 Part C also
   makes scale-one radius six a finite exact recursion; enumerate it. At radius
   seven the discrepancy deficit `13-2m` changes sign, so a new potential must
   use overlap debt, owner diversity, or signed component/comb incidence rather
   than mean danger density alone. This is a method wall within the scale-one
   chart, not a claim that arbitrary higher-radius decks have been classified.
2. **Deep colour cover.**  At `s=2`, prove a scale-free transversal lower
   bound above ten for the folded bad-atom hypergraph beyond THM-774's
   `max(U)<=19` unbounded-odd slice and THM-776's full height-100 slice, or add
   quantitative bounds on every dyadic seam guard strong enough to place the
   reconstructed ten-core/full packet inside one of those finite bases, or
   prove the global erosion failure
   `E_U not subset H_(x,y) minus R_U`. THM-789 proves that fixed-anchor packet
   refinement alone cannot do this, while THM-797 proves it whenever an odd
   exception-divisor grid leaves the explicit acceptance shell. THM-803 adds
   raw/parity support and the complete `26/52/78` anti-shell ladder. THM-817
   makes the endpoint/cusp selector adaptive with size
   `2c_E N_R+2W-2g<=20B^2+22B-2g`, while proving that the current gates permit
   `N_R=Theta(B)`. THM-821 factors every positive-odd-pair selector into exact
   pair-indexed atomic sum arcs and proves that the coarser cell/component,
   width, owner, and event projections have mixed verdict fibres in the audited
   `(13,5)` bank.
   THM-824 proves that their full symmetric union at ratio `(13d,5d)` instead
   factors through the exact intrinsic radius budget
   `rho_(C,d)(E)+rho_d(R)<=2/169`.  Prove that every admissible core
   violates that budget, extend the no-switch geometry to other odd ratios,
   or prove a height-independent contradiction/decreasing invariant on the
   individual exact signed sum arcs with owner ancestry;
   the `U_*` row shows that all short gates and both global maximizers can pass
   while a nonmaximal component escapes, and THM-817's linear-satellite family
   rules out a connectedness shortcut.  At higher sheets, classify and rule
   out THM-769's remaining persistent colour covers; the special quartic
   order-three equality shell is already closed by THM-816.

Equivalently in the original top-peel language, prove that every primitive
nonextremal eleven-speed core below THM-763's finite height has at least one
safe component whose rational endpoint band is incompatible with every top
speed in its THM-759/766 cone.

For the shallow branch, let `kappa` be the number of components of the open
danger-tooth graph and `P_splice` the number of protected end/start
coincidences.  The settled twelve-speed LRC ensures that the open danger union
is not the whole circle, and exact endpoint sweep gives

```text
chi_13 = kappa-P_splice = number of open 1/13-safe components.
```

For full nonzero residues, the nominal complementary-pair capacity is
`P*=2 sum_(r=1)^6 gcd(w_r,w_(13-r))`, so

```text
chi_13=(kappa-P*)+(P*-P_splice).
```

This separates overlap-rank shortage from third-runner blocker debt.  The
height-one lift cube had 4,085 rank-short rows, 9 blocker-debt rows, and one
zero-defect row, the nonprimitive doubled AP; all 4,094 primitive rows had
`chi_13>=2`.  THM-770 supersedes that cube by an exact height-twelve theorem,
while THM-795/800/804/806 close the first three Hamming directions.  THM-810
isolates the radius-four ramification and THM-815/816 close both of its
branches.  THM-820 makes the scale-one radius-five warning finite and rejects
its full height-one face; THM-822 rejects the full height-at-most-two slice and
locates the static kernel boundary.  The two complete boxes and all-scale
descent remain the shallow residual; THM-815 also makes scale-one radius six
recursively finite and loses its discrepancy deficit at seven. Arbitrary
higher-radius decks and the remaining deep sheet packets are still
unclassified.

## F. Information-preservation / Tournament Analysis

The exact vertices are witness obligations `(q,a)`, safe components, endpoint
splices, or the `s` sheet fibres—not runners by default.  Modulus and runner
tournaments are telemetry: changing gauges flips many edges while the blocker
verdict stays fixed, and pairwise component compatibility does not imply
containment.  The component-phase atlas makes the same loss explicit: a
transitive tournament with singleton SCCs and one Hamiltonian path can still
have negative slack at every component because phase order forgets interval
width.  In THM-770 all 66 residue-obligation comparisons tie, while the
endpoint-owner hypergraph still distinguishes thirteen solutions.  The
deciding objects are therefore:

- the zero-owner/signed-pair blocker deck for small periods;
- the component-tooth incidence hypergraph with endpoint widths, divisor pins,
  and core-maximizer residues for tight completions;
- the off-sheet-runner by sheet incidence cover, with effective orders
  `s/gcd(w,s)` and persistent colour ownership over the quotient loose set;
- in the two-sheet equality packet, the folded half-frequency diamond together
  with the exact locations of all deep components, their pointwise escape
  margins, labelled phase-cell/anchor incidence, the symmetric return erosion
  of THM-782/789, THM-797's decorated divisor-class/exception-shell incidence,
  THM-803's parity-twisted support/anti-shells, and the complete owner-labelled
  endpoint/cusp selector; THM-807 separates its mandatory central return cell,
  while THM-817 represents every satellite as a signed maximum-speed cell and
  retains its left/right endpoint owners and incidence with each deep
  component. THM-821 makes each positive-odd-pair incidence atom the
  pair-indexed circular sum state `((x,y),C+R_k)`, whose endpoint/cusp minimum
  and signed margin determine its verdict; the fixed `(13,5)` mixed fibres
  prove that neither factor nor its coarse shadows can be ranked separately.
  THM-824 proves that after assembling the symmetric full return union, the
  same fixed-ratio conjunction factors through two extremal radii, while
  individual satellites and dynamic owner transport do not;
- for dyadic descent, the quotient chain and binary safe-child fibers with
  eligibility radii and divisor obligations;
- for bounded two-sheet truth, the inclusion-minimal bad-atom/quotient-speed
  incidence hypergraph, which can forget metric positions after atomization;
- in the shallow packet, missing-owner splice sheets joined to replacement
  teeth and then to strict-safe components; THM-795/800/804/806 close the first
  three radii, THM-810 identifies the radius-four scalar/parity ramification,
  and THM-815/816 close its two branches.  At radius five THM-820 proves that
  the owner-exit tournament supplies a finite reduction but loses integer
  centre and scale; THM-822 proves on the height-at-most-two bank that adding
  those centres still mixes exact `M`, whereas literal endpoint faces are
  statically injective.  THM-823 shows that five-colour scalar attenuation has
  infinite cones, while bounded common sheets retain only the forward flags
  of a three-coset cycle with a parity fibre and unresolved metric strictness.
  The exact action state is `(E_S,V)`, the literal residual
  interval union together with its remaining labelled operation bank; its
  continuation equivalence is equality of the terminal emptiness verdict
  after every legal future subset of `V`. Endpoint owners and widths decorate
  this operation-closed carrier, and its discrepancy potential stays coercive
  through six remaining combs.

These objects preserve the LRC predicate.  Their tournament quotients destroy
joint blocker ownership, multiplier identity, scale, ramification, and
simultaneous alignment.  THM-774's two scalar gauges on the finite exceptional
proof obligations flip 114 edges while remaining transitive; neither ranking
can decide the containment in (C3).
