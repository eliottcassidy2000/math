---
id: HYP-6820
title: Uniformity audit for the LRC(14) q<=25 good-period claim and the n=12 sporadic branch
status: PARTIALLY RESOLVED — uniform q<=25 is DISPROVED; n=12 is uniformly finite and sheet-stratified, shallow-exact through lift height twelve, and two-sheet-exact for max(U)<=19 with unbounded odd exceptions and for max(A)<=100, but branch emptiness remains OPEN
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
  - HYP-6750
  - HYP-6775
related:
  - HYP-6780
  - HYP-6785
  - HYP-6800
  - HYP-6815
  - THM-771
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

These results reveal the faithful carrier more precisely than “ten even plus
two odd.”  It is a folded bad-atom/core-tooth incidence hypergraph equipped
with a binary safe-child map and divisor unit columns.  The runner tournament
records which odd runner owns a sheet but forgets the simultaneous hitting
number; the sheet tree records recursion but needs the incidence sidecar.

## E. Remaining proof obligation

The uniform theorem now has two explicit residuals:

1. **Shallow descent.**  Prove that a primitive full-nonzero-residue packet
   with `chi_13=0` descends into THM-770's height-twelve box, or extend its
   owner-CSP by a scale-free coherence argument.
2. **Deep colour cover.**  At `s=2`, prove a scale-free transversal lower
   bound above ten for the folded bad-atom hypergraph beyond THM-774's
   `max(U)<=19` unbounded-odd slice and THM-776's full height-100 slice, or add
   quantitative bounds on every dyadic seam guard strong enough to place the
   reconstructed ten-core/full packet inside one of those finite bases.  At
   higher sheets, classify and rule out THM-769's persistent colour covers,
   beginning with the explicit `s=3` equality shells of THM-772.

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
but the same warning remains: bounded exactness is not the global coherence
lemma.

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
  `s/gcd(w,s)` and persistent colour ownership over the quotient loose set.
- in the two-sheet equality packet, the folded half-frequency diamond together
  with the exact locations of the quotient loose components.

These objects preserve the LRC predicate.  Their tournament quotients destroy
joint blocker ownership, multiplier identity, scale, ramification, and
simultaneous alignment.  THM-774's two scalar gauges on the finite exceptional
proof obligations flip 114 edges while remaining transitive; neither ranking
can decide the containment in (C3).
