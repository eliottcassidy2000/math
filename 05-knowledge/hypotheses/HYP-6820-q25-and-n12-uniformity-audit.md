---
id: HYP-6820
title: Uniformity audit for the LRC(14) q<=25 good-period claim and the n=12 sporadic branch
status: PARTIALLY RESOLVED — uniform q<=25 is DISPROVED; n=12 is uniformly finite and sheet-stratified, with the full AP-centred Hamming-one through Hamming-four stars, the proper scale-one Hamming-five chart, every effective-order-at-most-21 common-sheet H5 survivor language reduced to the already closed all-one/all-three/mixed lifts, the entire proper scale-one Hamming-six chart plus all common dilations closed except for the doubled-AP equality family `2c[12]`, and the first ramified face `c=2` closed except for its ordinary-AP equality presentation `[12]`; the remaining common-sheet H5 presentation bank is finite with `min D_i<=21<max D_i<=10,584`, is {2,3,7}-smooth, and has no private prime power; effective order is the exact H6 ramification degree, every H6 branch containing order one is finite-decidable, and THM-860 makes the entire primitive proper AP-centred H6 ramified branch finite-decidable with `2<=c<=2,177,280`; the operation boundary is exact—literal endpoints are Markov for monotone insertion, while H0/H1 are not and endpoints alone fail deletion; the two-sheet branch has exact sum-arc, target-switch, and thin owner-shell carriers, with the first shell-five congruence class excluded and both the universal single-column endpoint-grid template and every fixed two-column unit endpoint family refuted in the relaxed model for the other three; global emptiness remains OPEN in the finite smooth-ramified H5 bank and its metric languages, the primitive H6 metric context trees at `c>=3`, radius-seven endpoint/third-moment and correlated AP-window residuals, uniform radius/sum-arc exclusion, lift-dependent `s=5` deep-shell certificates, dyadic/collar residuals, and higher sheets
source: codex-2026-07-14-S3
progress_note: >-
  THM-857 closes all 924 proper scale-one H6 deletion roots by a
  580,919,164-node exact component recursion; the unique covering terminal is
  2[12].  This subsumes the former 903 `f<=2` primitive-core row residual and
  the odd-label mixed-parity branch at scale one. THM-859 conjugates the tree
  across common dilation, identifies D as the deck-mask ramification degree,
  and makes every H6 order-one gate finite-decidable. THM-860 supplies the
  hereditary five-order lcm identity, all six-colour relative cuts, the
  optimized scale bound c<=2,177,280, and a finite exact recursion for every
  primitive proper AP-centred ramification language. THM-861 evaluates its
  first ramified face: a 41,882,982-node exact tree leaves only the ordinary
  AP [12], so every sporadic c=2 packet is loose. THM-862 classifies the next
  `c=3` common-sheet stalk into 212 presentations and 1,504 arithmetic
  contexts, whose unit fibres are affine matching codes. Its 68 sheet orbits
  are not metric quotients; the exact first layer has 146,912 edges and the
  metric recursion remains open. THM-836 §§6B--6C also rule out
  a U-independent single-numerator endpoint-grid proof and every fixed pair of
  unit endpoint columns in the relaxed shell-admissible model, without closing
  the three structured shell-five classes.
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
  - THM-831
  - THM-836
  - THM-837
  - THM-840
  - THM-844
  - THM-845
  - THM-847
  - THM-856
  - THM-857
  - THM-858
  - THM-859
  - THM-860
  - THM-861
  - THM-862
  - HYP-6750
  - HYP-6775
related:
  - HYP-6780
  - HYP-6785
  - HYP-6800
  - HYP-6815
  - THM-771
  - THM-829
  - THM-831
  - THM-836
  - THM-837
  - THM-840
  - THM-844
  - THM-845
  - THM-847
  - THM-856
  - THM-857
  - THM-858
  - THM-859
  - THM-860
  - THM-861
  - THM-862
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

THM-820 alone does not empty these boxes.  Its exact shallow census rejects all
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

THM-845 supplies exactly that missing recursion.  At every prefix it chooses a
longest strict-safe component and applies THM-815's discrepancy bound to the
least remaining speed, while retaining THM-820's exhaustive doubling versus
exceptional-cycle dichotomy.  The exact state counts at depths one through
five are

```text
40,590, 612,221, 111,675, 7,255, 9.
```

The exceptional branch has `415,178,1,0` states from depth two through five
and therefore no completion.  Each of the nine doubling terminals has a
nonempty strict-safe interval; independent exact endpoint reconstruction gives
`M` strictly above `1/13` in every row.  Hence every proper scale-one
residue-preserving five-coordinate lift is loose, uniformly in all five lift
heights.  The former H5 boxes are closed, not merely searched to a height.

THM-822 empties the complete height-at-most-two slice, all
`C(12,5)2^5=25,344` packets, with the same unique `1/12` minimum.  It also
audits the proposed three-face quotient exactly.  The labelled live relation
`H0` and its integer-centre decoration `H1` have identical fibres in this
ratio band, and their `111,006` ordered kernel rows include `3,810` fibres
mixing exact `M`; one shared relation has maxima `1/4` and `12/37`.  Literal
strict-safe component endpoint words `H2` make the bounded three-face join
injective, while the final residual alone still has fourteen collision pairs.
This is a static bounded reconstruction theorem.  THM-845 subsequently shows
that the literal residual endpoint word is sufficient for arbitrary-height
scale-one recursion, while the coarser `H0=H1` quotient is not.  What remains
unproved is transport of that state across AP-scale ramification: the next
quotient target is an endpoint-owner/component-incidence codec between `H1`
and `H2` compatible with the oriented deck action.

THM-815 C.1 now closes the nonprimitive part of the adjacent scale-one H6
face.  If

```text
W=([12] minus R) union {r+13h_r:r in R},       |R|=6,
```

is nonprimitive, its six retained labels force `gcd(W)=2`,
`R={1,3,5,7,9,11}`, and all six heights odd.  Writing
`h_(2i-1)=2k_i+1` and dividing by two gives exactly

```text
{1,...,6} union {6+i+13k_i:1<=i<=6}.                    (B6)
```

Zero positive `k_i` is `[12]`; one through five positive coordinates are the
already closed H1--H5 charts.  The all-positive top-half H6 chamber is closed
by a frozen 136,288-prefix longest-component recursion with depth counts
`1,54,3612,130515,2104,2,0`, no covering prefix, and no depth-six state.
Both deepest leaves agree exactly with an independent closed-danger-union
reconstruction.  Therefore `2[12]` is the only possible tight nonprimitive
scale-one H6 row.  This argument by itself is not primitive H6 closure and does not transport
through arbitrary AP scale.  On the doubled-AP equality side, the core
`{1,...,6}` has twelve strict-safe components and the six danger combs
`D_7,...,D_12` cover them with 34 incidences; four components have zero
overlap debt and a unique owner.  This weighted component--comb incidence,
not pairwise comb overlap, is the predicate-preserving Kakeya-needle carrier.

THM-815 C.2 makes the first uniform cut inside primitive H6.  For a missing
row `R`, let `f(R)` count the full antipodal pairs `{r,13-r}` contained in
`R`.  The 923 rows whose retained core is primitive have distribution

```text
f=0,1,2,3:                 63,480,360,20.                (B6a)
```

Exactly `2f` nonzero points of the thirteenth grid lie in the retained-core
strict-safe set. At an oriented AP cusp the owner has an exact germ cap
`B_r(P)` and provider edge `r->s` has multiplier `c([as])/2`. If every owner
evades its cap, a selected provider digraph contains a directed cycle whose
multiplier product must be at most one. The twenty `f=3` rows have
minimum-cycle-product census `1/16:6, 1:2, 3/2:12`; the twelve expanding rows
and two equality rows therefore force a cap. In six of those fourteen rows no
proper lift satisfies the cap, so they close immediately. In the other eight,
the cap forces `u_5=18` and/or `u_6=19`; ten fixed-coordinate
longest-component trees close those slices in `3,699` exact states, with
aggregate depth counts `10,228,2176,1236,49,0` and no covering prefix. Hence
fourteen rows close by the AP-germ reduction.

The six product-`1/16` rows show exactly where that quotient stops: their
contracting cycles satisfy the one-step handoff constraints.  Restoring the
literal residual endpoint union and numerically ordering all six remaining
combs gives six exact longest-component trees with aggregate depth counts

```text
6,494,46813,2531670,74415,202,0,                       (B6b)
```

or `2,653,600` states.  Every prefix is noncovering, none reaches depth six,
and independent closed-danger complements reproduce all six roots and all 202
deepest dead leaves, for 208 crosschecks.  Thus all twenty `f=3` rows are
loose and the primitive-core row bank falls from 923 to the 903 rows with
`f<=2`.  The separate odd-label row with mixed height parity can be a primitive
packet although its retained core has gcd two, so the open scale-one H6
label-pattern count is 904.  The doubled AP `2[12]` is the separate
nonprimitive equality.  The cusp count is a routing grade, not a quotient: it
forgets pair signs, germ owners, exact endpoints, divisor obligations, heights,
and component overlap.  The six-root closure therefore supports, rather than
replaces, the incoming tropical defect-graph view: edge/cycle data routes the
search, while the residual endpoint state is the operation-complete carrier.

THM-857 now closes that entire scale-one residual by retaining exactly the
information the cusp quotient discarded.  It numerically orders the six
proper lifts and recursively carries the literal strict-safe component union,
remaining labelled residue progressions, last speed, and exact shortcut
witness.  THM-815's longest-component bound makes every next-speed range
finite without a height cutoff.  Over all 924 deletion roots the frozen tree
has

```text
nodes by depth = 924, 83,881, 8,906,315, 559,202,706,
                 12,671,505, 53,812, 21;
logical nodes = 580,919,164; candidate edges = 580,918,240.
```

Exactly twenty terminals retain a strict-safe interval and exactly one covers:
the missing odd labels lifted to `14,16,18,20,22,24`, namely `2[12]`.  Hence
all 923 primitive retained-core roots and every mixed-parity completion of the
exceptional odd-label root are loose.  A separate closed-danger-union replay
reconstructs each expanded prefix from scratch and hashes the labelled root
trees.  The remaining H6 obligation is metric evaluation of ramified
common-sheet languages, not the former 903-row plus exceptional-root scale-one
bank. THM-859 proves that common dilation is exactly conjugate and every branch
containing `D=1` is finite-decidable; THM-860 makes every primitive
AP-centred `D>1` language finite-decidable as well.

An independent C.3 replay crosschecks a nontrivial subforest using the earlier
residual-component implementation: all 52 `f=2` roots with exact first-speed
cap at most 156 close in `9,888,159` states, with zero covering prefix and
1,323 closed-danger reconstructions.  It is a subcertificate of THM-857, not
a new frontier count.

At seven remaining combs THM-815's single-comb discrepancy coefficient changes
sign.  THM-856 proves the correct second-order functional is the
Hunter--Kounias maximum spanning tree of restricted pair overlaps.  At ideal
densities its lower-bound coefficient is `(165-22m)/169`, positive for seven
combs and negative for eight.  The S14 live-pull referee corrects its original
pair formula.  For `x=ga,y=gb`, `(a,b)=1`,

```text
mu(D_x intersect D_y)
 =4/169+[Q(a+b)-Q(|a-b|)]/(169ab),
Q(c)=r(13-r),  r=c mod 13,                              (B7)
```

and if `E` has `c_E` components then restriction adds at most `2c_E/g` to
the projective-product discrepancy.  The sawtooth in (B7) has both signs:
`mu(D_6 intersect D_7)=2/91<4/169`.  Hence neither the original one-sided
overlap floor nor a deficit bounded only by raw `min(x,y)` is valid.  The
radius-seven state must retain the reduced pair ratio, common scale, mod-13
sawtooth, and endpoint discrepancy on each tree edge.  THM-856's exact pilots
show this functional certifies generic packets, not that all packets are
closed.

The S15 node/edge defect decomposition makes that state exact.  With
`e=mu(E)`, put

```text
s_i=mu(E intersect D_(x_i))-2e/13,
c_ij=e[Q(a_ij+b_ij)-Q(b_ij-a_ij)]/(169a_ijb_ij)+eta_ij,
```

where `(a_ij,b_ij)` is the reduced pair and `eta_ij` is its restricted
endpoint discrepancy.  Then for `m` remaining combs the Hunter lower bound is

```text
(165-22m)e/169-sum_i s_i+MST(c).                       (B8)
```

At seven combs the baseline is `11e/169`.  Kruskal rewrites `MST(c)` as the
sum of edge-credit thresholds times jumps in the graphic rank
`m-(number of threshold-graph components)`.  Thus the live finite quotient is
a vertex-period and edge-threshold incidence language.  An edge histogram or
average pair density is not sufficient.

The S15 refinement identifies the algebra behind those colours.  With
`p=2/13`, `F_E=1_E-e`, and `F_i=1_(D_(x_i))-p`,

```text
s_i=<F_E,F_i>,                 h_ij=<F_i,F_j>,
theta_(Eij)=integral F_E F_i F_j,
c_ij=e h_ij+p(s_i+s_j)+theta_(Eij).                       (B8a)
```

Thus the completed restricted matrix

```text
G^E_ij=c_ij-p(s_i+s_j),                         i!=j,
G^E_ii=22e/169+9s_i/13                                  (B8b)
```

is PSD for every packet.  The edge discrepancy is not free noise: after its
node contribution is removed, it is the centered prefix--comb--comb third
moment.  All principal minors of (B8b) constrain admissible node and edge
colours.

The global projective covariance has two further exact forms.  For arbitrary
gcd patterns, `h_ij=<F_i,F_j>` has a PSD completion with constant diagonal
`22/169`.  On `Z/13`, the sawtooth numerator is

```text
H_(r,s)=Q(r+s)-Q(r-s)=C_Q(R-I),                            (B8c)
```

a rank-six PSD kernel whose kernel is the seven-dimensional reflection-even
sector.  A common-scale packet with pairwise-coprime quotients therefore has
an explicit six-dimensional sine Gram realization.  The sharp scalar floor
is `h_ij>=-11/1014`, at reduced ratio `1:12`; six such edge floors exactly
cancel, but cannot cross, the ideal `11/169` margin when treated independently.
The complete graph forces strict ratio consistency.  At threshold
`-3/676`, the bad edges have exactly ten reduced ratio types, and their
directed multiplicative relation is triangle-free.  The good complement has
at most two components; using the second global edge level on the cross cut
gives

```text
MST(h)>=-237/7436,
11/169+MST(h)>=19/572.                                  (B8c1)
```

This is a strict projective margin for every seven-speed packet, independent
of its gcd pattern.

The tree functional is recursively closed by more than its present value.
Tropical deletion--contraction is exact, while insertion of a vertex `v` obeys

```text
MST(G+v)=max_(pi partition of V)
 [sum_(B in pi)MST(G[B])+sum_(B in pi)max_(i in B)c_(iv)]. (B8d)
```

Hence the induced-subtree/partition-defect profile is sufficient for the next
comb insertion.  Ordering graph edges by credit gives a transitive tournament;
its graphic-rank derivative is the Kruskal acceptance word.  The tournament
alone loses incidence.  After higher levels are contracted, a level edge is a
loop (in no maximum tree), a bridge (in every maximum tree), or an optional
nonloop (in some but not all), giving an exact three-colour edge classifier.

Finally, the strict tree margin gives an all-packet endpoint criterion.  Let
`Gamma_h(x)` be the least `sum_(ij in T)1/gcd(x_i,x_j)` over projective-
maximizing trees.  Then

```text
L_H>=19e/572-2c_E[sum_i 1/x_i+Gamma_h(x)].               (B8e)
```

If `x_i=g a_i` for arbitrary distinct positive `a_i`, the bracket is at most
`(H_7+6)/g=1203/(140g)`, so

```text
L_H>=19e/572-(1203/70)c_E/g.                             (B8f)
```

Every fixed common-dilate ray is therefore loose beyond an explicit scale,
with no pairwise-coprime hypothesis.  What remains is the small-gcd endpoint
field.  In the exact six-packet replay, all four Hunter successes have positive-
credit spanning trees.  Both consecutive failures instead have all 21
endpoint and total credits negative despite 14 and 18 positive projective
terms.  Their conditional-overlap tournaments are still transitive, locating
the failure in the joined third-moment edge field.

THM-840 now separates those operations exactly.  The THM-822 quantitative
liar pair has the same `H0=H1`, but inserting the same labelled speed `(6,19)`
creates different handoffs; neither coarse observation is insertion-Markov.
Literal exact endpoint words do update at every height under a monotone comb
addition by `E->E intersect Safe(u)`.  They are not deletion-Markov and do not
determine exact `M`: two thirteen-speed rows with empty strict residuals split
after deleting their common speed `13`.  Thus a pure-addition recursion may
forget owners after forming exact endpoints, while deletion, replacement, and
transport must retain the labelled tooth bank and ancestry.

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
quartet.  The all-one case normalizes to the scale-one chart and is now closed
by THM-845 (zero-height faces reduce to already closed lower Hamming radii).
THM-847 now closes every proper mixed one-plus-order-three lift; its order-one
speed is exactly `3(b+13h)`, `h>=1`, while `h=0` is THM-816's Hamming-four
face.  Thus every survivor language with all five effective orders at most
twelve is closed.  At that stage the apparent remaining work was unbounded
common-sheet classification, not a stronger scalar cutoff; THM-858 below turns
that apparent unboundedness into a finite ramification bank.

THM-837 performs that metric erosion on one of the 96 all-order-three
directed-flag/parity contexts, namely `C={1,5,8,12}`, `b=10`, bits
`(1,1,1)`.  The operation-complete active-run recursion has 75,371 states,
zero covering prefixes, and 57 nonempty terminals, proving this context loose
at every height.  Its raw and conditioned comb tournaments are both
transitive on all terminals while one through ten edges flip.  The proof is
in the residual endpoints plus remaining CRT progressions.

THM-844 strengthens the recursive inequality before extending the census.  A
longest residual component has length at least `L/K`, so its THM-815 cap is no
larger than THM-837's global `K/L` cap.  Applying this sharper cap separately
at every state closes all 96 contexts in 28,876 exact states: there is no
covering prefix and every branch dies by depth four, before any depth-five
terminal.  The THM-837 context itself falls from 75,371 states to 213.  The
binary tournaments remain transitive in both gauges across all contexts even
though 492 edges flip.  The exact carrier is the evolving incidence of literal
residual components with remaining labelled CRT comb obligations, decorated by
active endpoints and the last chosen speed.  This closes the bounded
all-order-three common-sheet language.

THM-847 applies the same carrier to the mixed language.  The exact recursion
visits `31,715` states over `96` contexts, finds no covering prefix, and dies
before depth five.  A standalone implementation reproduces the certificate,
and an independent literal sheet census finds exactly those contexts among all
`63,360` marked mixed presentations.  Both tournament gauges remain transitive
while 395 edges flip.  Together THM-844/845/847 empty THM-823's complete
effective-order-at-most-twelve common-sheet bank; before THM-858, common sheets
outside that bank and arbitrary deck descent remained.

THM-858 reaches the least-order pivot without asserting a false all-order
cutoff.  For every nonempty colour set `S`, put

```text
L=lcm(D_j:j notin S),        m_i=D_i/gcd(D_i,L).
```

The complement misses some sheet at some owner, and colour `i` meets the
resulting `L`-fibre in at most `ceil(2m_i/13)/m_i` of its points.  Hence every
common-sheet cover obeys

```text
sum_(i in S) ceil(2m_i/13)/m_i>=1.                     (B8)
```

Singleton cuts forbid private maximal prime powers.  In the exact new bank
through effective order 21 they reject 1,870 of 2,055 scalar rows; top-prime
cuts reject the remaining 185, and a separate literal census finds zero
covers among 51,360 unit words.  THM-823's uniform order-one classification
supplies the branch omitted by that source.  Therefore all common-sheet
languages through order 21 are the already closed all-one, all-three, and
mixed languages.  Any new common-sheet row must have

```text
min D_i<=21<max D_i<=10,584,
every D_i {2,3,7}-smooth,
D_i divides lcm(D_j:j!=i),
and all cuts (B8).                                      (B9)
```

The upper endpoint is structural.  Apply (B8) at each adjacent `p`-adic
valuation level.  For `p=2`, the only possible upper-set sizes are `2,3,4`,
with level gaps at most `1,1,3`.  The apparent sum five cannot occur: a
size-four initial gap three forces all four relative orders to equal eight,
while a gap two forces them into `{4,8}`; either equality case leaves total
dyadic spread at most three.  If the initial gap is one, the later size-three
and size-two gaps add at most one each.  For `p=3`, sizes `3,4` have gaps at
most one; for `p=7`, only size four is possible and its gap is one.  Thus

```text
range(v_2)<=3,             range(v_3)<=2,
range(v_7)<=1.                                            (B10)
```

Comparing every order primewise with THM-823's order `D_*<=21` gives
`D_i<=D_*2^3 3^2 7<=10,584`.  The common-sheet presentation bank is therefore
uniformly finite, though the arbitrary metric lifts in its surviving
languages are not decided by this arithmetic bound.

Raw and complement-conditioned capacity tournaments are transitive on all
185 residual rows despite 574 edge flips.  The rejecting operation is a
subset/complement cut, so the faithful carrier is instead the maximal-prime-
power hypergraph decorated by complement-lcm fibres and affine owner-sheet
intervals.  This is a genuine finite reduction of the H5 common-sheet object, not a
global H5 or `n=12` closure.

THM-859 identifies how this H5 base enters H6.  Common dilation conjugates the
entire THM-857 insertion tree, while an order-`D` insertion creates exactly
`D` distinct deck masks; the scale quotient is a continuation congruence
exactly at `D=1`.  Deleting an order-one H6 colour leaves an H5 common-sheet
presentation.  Hence two or more order-one colours give either the already
closed all-one dilation or exactly `336` two-order-one/four-order-three
contexts.  Exactly one order-one colour leaves the finite THM-858 H5 bank; in
its through-21 face this gives `96*7=672` marked contexts.  Every fixed context
has a finite longest-component recursion because its future speeds lie in
labelled progressions modulo `13 lcm(D_i)`.  The reduction is uniform and
finite, but no zero-cover verdict is yet claimed for those ramified trees.

THM-860 removes the remaining unbounded-scale quantifier for every primitive
proper AP-centred H6 packet.  THM-765 first makes every leave-one-out core primitive,
which is equivalent in effective-order language to

```text
lcm(D_j:j!=i)=c                         for every i.
```

For a nonempty colour set `S`, the complement has at most four colours unless
`S` is a singleton.  The former case inherits THM-858's missed-fibre proof;
the singleton case is equality because the displayed lcm identity makes its
relative order one.  Applying all resulting upper-set cuts, rather than only
separate adjacent-gap maxima, gives

```text
range(v_2,v_3,v_5,v_7)<=(5,3,1,1),
range(v_p)=0 for p>=11, p!=13,
c/min(D_i) divides 2^5*3^3*5*7=30,240.
```

Six-colour attenuation supplies `min D_i<=77`; the all-six relative cut
sharpens this to `min D_i<=72`, hence `c<=2,177,280`.  Range zero at a large prime means a common factor, not prime
absence; individual H6 orders are not asserted smooth.  Every fixed
order/unit presentation has six labelled progressions modulo `13c`, and the
THM-815 component cap makes its exact insertion tree finite.  The unique
worst root has longest component `31/1430`, so a labelled progression has at
most 37 first heights.  The first ramified scale `c=2` consists of exactly 64
signed-doubling Hamiltonian-cycle contexts. THM-861 evaluates their complete
unbounded step-26 metric trees in `41,882,982` logical nodes. The unique cover
replaces labels `7,...,12` by `1,3,...,11`, producing exactly `[12]`; every
sporadic `c=2` packet is loose. Thus “arbitrary-scale H6 transport” has become
“evaluate the finite primitive H6 metric context bank at `c>=3`,” beginning
with the THM-862 scale-three bank.

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

THM-831 completes that target-side extension sharply.  For every primitive
coprime opposite-parity half-frequency pair `(alpha,beta)`, it gives all folded
target components as exact Bezout-offset balls.  The nonempty symmetric-return
predicate has a no-switch radius factorization exactly for the sixteen rows
`4<=alpha<=9`; targets with `alpha<=3` are empty, and every `alpha>=10` has an
explicit strictly negative weighted switch triple.  The abstract theorem also
requires `4h<1` to exclude same-label winding; this is automatic for every
folded target but sharp for arbitrary balls.  Common gcd creates raw deck
switches, but multiplication by the gcd preserves the quotient-scaled radius
criterion.  Its 518-pair formula audit and 6,160 compact-packet replay are
exact.  Hence the remaining all-pair problem is no longer to guess which
targets scalarize: it is to force arithmetic violation on the legal sixteen
ratio types, or to transport the weighted ternary switch/owner carrier on all
larger types.  Individual nonsymmetric satellites remain outside the scalar
factorization.

THM-836 supplies the first all-size arithmetic exclusion inside the guarded
`(13d,5d)` row.  At the forced signed-wall centre, the two directional owner
exit tables are permutations of
`(9,20,31,42,53,64,75,86,97,108)`.  If `B<13d`, both directions must use the
coefficient-nine owners; their real interval and mod-13 separation empty
`2B-13d=1` and `3`.  Hence every configuration satisfying THM-836's ten-
distinct-speed, exact signed-complement, and guarded-containment hypotheses
obeys

```text
13d<=2B-5.                                                (C6)
```

At the next shell `s=2B-13d=5`, local feasibility is exact:

```text
d mod 52 in {11,15,37,41},
missing signed class={3,10},       owners={B-3,B-2}.      (C7)
```

The class `d=11 mod 52` is nevertheless uniformly impossible.  The mandatory
raw-exception grid `q=5d` and the unit numerator
`p=(45d-1)/26` put every allowed signed-complement lift at depth at least
`1/11`, while the folded `(13d,5d)` target value is `2/5<11/13`.
Thus even the guarded containment fails.  A separate exact 1,605,632-row
census closes the single value `d=15`: after divisor completeness and parity
support, 3,004 rows pass the `q=75` grid, four pass `q=195`, and none pass
both.  That finite statement does not extrapolate to `d=15 mod 52`.
Consequently the uniformly open shell-five classes are exactly
`15,37,41 mod 52`; their next proof must add a divisor/component/satellite
incidence beyond the local two-owner packing.

THM-836 §6B rules out the most direct proposed continuation uniformly.  Let
`P_d` contain every possible free raw lift up to `B`, together with the forced
owners `B-3,B-2,B`.  An exact 35-speed skeleton and a no-jump progression
argument localize the common deep set for `P_d` to four thin intervals around
`4/13` and `9/26`.  Neither `q=5d` nor `q=13d` has a unit numerator in those
intervals for any of the three open classes.  Thus no single endpoint-grid
column chosen from `d` alone can work for every allowed `U`.  This is a
method obstruction, not a shell closure: the next certificate must depend on
the actual lift set, combine several numerator columns, or leave the endpoint
grids.

THM-836 §6C also closes the first fixed multi-column repair in the relaxed
shell model.  Every unit `q=13d` column has free danger support at least three;
every nonforced unit `q=5d` column has support at least two, except
`d=41 mod 52`, `p=+/-(45d+1)/26`, whose support is exactly `{4}`.  A labelled
matching therefore constructs one relaxed shell-admissible lift set killing
any fixed pair of unit endpoint columns; the two singleton signs share the
same concrete killer `B-18`.  The 432-row affine proof has 429 generic rows,
three rows forming the exceptional family, and minimum strict margin
`1/6545`.  This does not survive automatically after imposing THM-772 divisor
completeness or THM-803 parity support.  The next multi-column proof must use
those structured restrictions, depend on `U`, use at least three columns, or
change denominator.

These results reveal the faithful carrier more precisely than “ten even plus
two odd.”  It is a folded bad-atom/core-tooth incidence hypergraph equipped
with a binary safe-child map, divisor unit columns, every owner-labelled deep
component, and its symmetric return packet, selector cusps, and escape margins. The runner
tournament records which odd runner owns a sheet but forgets the simultaneous
hitting number; the dyadic quotient chain records recursion but needs binary
sheet-fiber, incidence, component, and phase-cell sidecars.

## E. Remaining proof obligation

The uniform theorem now has two explicit residuals:

1. **Transport the closed scale-one shallow action across AP scale and common
   sheets.**
   THM-795/800/804/806 close the complete AP-centred Hamming-one through
   Hamming-three stars.  THM-810 splits radius four into common scale and one
   order-three coset interface; THM-815 and THM-816 close those alternatives
   independently.  THM-820 has now classified the radius-five collar cycles
   and anchored recursion and reduced the full scale-one chart to the two
   explicit branches above.  THM-845 empties them by the commuting idempotent
   action `T_u(E)=E intersect Safe(u)` on literal residual interval unions,
   with a row-wise component-length cap. THM-840 permits exact endpoints as the
   Markov state for this monotone insertion action, but requires the labelled
   tooth bank as soon as deletion or replacement is allowed. THM-837 first
   closes one all-order-three context with this state, and THM-844 closes all
   96 by the stronger longest-component recursion. THM-847 closes all 96 mixed
   contexts, and THM-845 supplies all order one. THM-858 extends this closure
   through effective order 21 and reduces every remaining common-sheet row to
   the finite `{2,3,7}`-smooth strip (B9). Classify and erode that shared-
   prime-power/complement-fibre language. Prove the remaining arbitrary-AP-scale
   descent without assuming common scale. THM-857 has now completed THM-815
   Part C's finite scale-one radius-six tree over all 924 roots: every proper
   row is loose except the genuine equality `2[12]`.  The former 903 `f<=2`
   primitive-root residual and exceptional mixed-parity branch are no longer
   open at scale one.  Transport the literal component/progression action and
   its exact shortcut witnesses across AP scale and common sheets. At radius
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
   violates that budget. Under its exact signed-complement and guarded-
   containment hypotheses, THM-836 reduces this fixed-ratio task to the shell
   `2B-13d>=5`; at equality it leaves exactly four classes modulo 52.  The
   explicit `q=5d` witness uniformly removes `d=11 mod 52`, and the single
   value `d=15` is finite-exact empty, but the classes `15,37,41 mod 52`
   remain open uniformly.  THM-836 §6B proves that no U-independent single
   unit numerator on either endpoint grid can close them.  Their next step
   must use U-dependent or multi-column incidence, a non-endpoint denominator,
   or another component/divisor/satellite obligation. THM-831 has
   classified the other fifteen viable folds in its coprime opposite-parity
   reduced family: repeat the
   arithmetic exclusion on the other fifteen primitive rows in (C5), and for
   `alpha>=10` retain the weighted ordered centre-switch three-hypergraph
   rather than a nonexistent radius quotient.  Alternatively prove a
   height-independent contradiction/decreasing invariant on the individual
   exact signed sum arcs with owner ancestry;
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
locates the static kernel boundary. THM-845 closes the two complete boxes at
arbitrary lift height. THM-840 identifies the operation-congruence boundary,
and THM-837 closes one of the 96 arbitrary-height order-three contexts by the
resulting active-endpoint recursion. THM-844 then closes all 96 with the
state-wise longest-component cap, while THM-847 closes all 96 mixed contexts.
THM-858 closes every new common-sheet row through order 21.  The shallow H5
residual is the finite smooth-ramified strip (B9) plus other all-scale deck
descent; THM-857 completes THM-815's scale-one radius-six recursion over every
root and leaves only `2[12]`. THM-859 conjugates this tree across every common
dilation and proves that the scale quotient fails precisely at effective order
`D>1`; branches containing `D=1` form a finite family of finite exact trees.
THM-860 then bounds every primitive proper AP-centred ramified scale by `2,177,280`
and transports each order/unit language to a finite exact component tree.
THM-861 closes its `c=2` face: the only cover is the ordinary AP `[12]`.
No proper primitive or mixed-parity unramified H6 chamber remains; the finite
ramified metric context bank at `c>=3` is open, and the discrepancy deficit is lost at seven. Arbitrary
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
  individual satellites and dynamic owner transport do not; THM-831
  proves that this primitive no-switch quotient exists exactly on sixteen
  reduced half-frequency types and replaces component tournaments outside
  them by the weighted ordered centre-switch three-hypergraph. THM-836's
  fixed-ratio deduction instead needs two directional owner minima and their
  mod-13 packing; its binary owner tournament loses the simultaneous
  obligation that empties the first two shells.  At shell five the faithful
  carrier is the unit-numerator by signed-lift by divisor/parity incidence:
  it uniformly removes `d=11 mod 52`, while the transitive four-class
  tournament does not distinguish the three open classes.  A column chosen
  independently of the actual free lifts is now ruled out on both endpoint
  grids, so the carrier must retain lift-dependent or multi-column data;
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
  statically injective. THM-840 proves the latter are Markov exactly for
  monotone insertion, while `H0/H1` fail even there and endpoint words fail
  under deletion without the labelled tooth bank. THM-823 shows that
  five-colour scalar attenuation has
  infinite cones, while bounded common sheets retain only the forward flags
  of a three-coset cycle with a parity fibre, together with the symbolic
  order-one alternatives. THM-837 closes one of the 96 all-order-three metric
  contexts, THM-844 strengthens the cap and closes all 96, THM-845 closes all
  order one, and THM-847 closes all 96 mixed contexts. THM-858 adds the
  complement-lcm fibre cut, closes every new order-through-21 row, and shows
  that the live language is a finite `{2,3,7}`-smooth bank with
  `max D_i<=10,584` and no private maximal prime power. On this edge the useful alternate vertices are prime-power
  carrier hyperedges; they must retain complement fibres and affine interval
  phase to preserve common-sheet coverage.
  The exact action state is `(E_S,V)`, the literal residual
  interval union together with its remaining labelled operation bank; its
  continuation equivalence is equality of the terminal emptiness verdict
  after every legal future subset of `V`. Endpoint owners and widths decorate
  this operation-closed carrier, and its discrepancy potential stays coercive
  through six remaining combs.  In nonprimitive H6, divisibility contracts
  the packet to the core `[6]` plus combs `7,...,12`; the resulting weighted
  component--comb incidence has four zero-overlap-debt unique-owner pins.
  THM-857 then retains the full literal action state over all 924 roots and
  closes every proper scale-one H6 row except the doubled AP equality. THM-859
  makes common dilation an exact action conjugacy and identifies the
  `Z/DZ` phase action as the missing sidecar at ramification. THM-860 adds the
  prime-power upper-set layer, bounds primitive proper AP-centred H6 scale by
  `2,177,280`, and proves that the resulting labelled progression/component
  bank is finite. THM-861 then evaluates the complete `c=2` fibre and finds
  only the ordinary AP `[12]`; its sparse signed cycle routes sheet parity but
  cannot replace the literal component/ray state that decides coverage. The
  earlier oriented AP cusps and weighted handoff cycles remain explanatory
  local vertices, but their intermediate 909-row frontier and eventual
  903-row plus exceptional-root residual were quotient artefacts rather than
  surviving scale-one branches.

These objects preserve the LRC predicate.  Their tournament quotients destroy
joint blocker ownership, multiplier identity, scale, ramification, and
simultaneous alignment.  THM-774's two scalar gauges on the finite exceptional
proof obligations flip 114 edges while remaining transitive; neither ranking
can decide the containment in (C3).
