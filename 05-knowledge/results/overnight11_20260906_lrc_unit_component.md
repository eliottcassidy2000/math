# Actual rank-eleven entry closes when the primitive eleven-component contains one

**Status: PROVED relative to inherited THM-3818 and cited lower-dimensional
LRC + FINITE-EXACT + INDEPENDENTLY AUDITED.** This is an actual-entry
subclass, conditional on the full decoder equality hypotheses below. It does
not infer those hypotheses from a pair atlas, pairwise cross heights, or the
rank of a selected relation family. LRC(14) remains open.

## 1. Inheritance, duplicate search, and the statement

The closest proved mechanism is the symmetric cyclic gluing theorem in
**THM-3818**, Section 6.5, (15aa), together with the forbidden crossing rows
and internal pair-height bound in Section 6.4, (15o)--(15q):
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`.
The lower-dimensional safe-arc supplier is the same cited LRCUpTo13 input
already used there. The arithmetic compiler below joins these two mechanisms.

Targeted searches covered current canon and result notes for unit/one/minimum
in an eleven-component, Euclidean division, quotient/remainder, bounded
crossing relations, and the constants Q(K+1). No exact prior statement of this
entry closure was found. This is a scoped search result, not an external
priority claim. The closest additional records read were:

* **THM-3878**, `01-canon/theorems/THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse.md`:
  core-scale absorption and remaining scale-one/two seams, without this unit
  compiler. The present proof does not need its seam census.
* **THM-4117** and **THM-4129**, respectively
  `01-canon/theorems/THM-4117-physical-eleven-plus-two-primitive-support-obstruction.md`
  and `01-canon/theorems/THM-4129-universal-two-speed-completion-of-the-eleven-speed-lrc14-body.md`:
  a genuine unit-core rank-eleven row and its fixed-body universal completion.
  Those results already close that particular body, not arbitrary unit shapes.
* **THM-4052**, `01-canon/theorems/THM-4052-lrc14-affine-component-width-escape-cones.md`:
  different affine width cones and a primitive nonunit component with huge
  maximum but small internal pair heights, used as the hostile below.
* The least-used sidecar `lrc14_unit_component_address_escape_probe_20260902`
  concerns selected tails that are units modulo a pack clock; its meaning of
  unit is different from the literal normalized coordinate one here.

The corrected near miss is the ninth-round h=M+1 family: large pairwise cross
heights and an atlas-admissible pair did not establish entry, because the
height-one crossing row 3h-42c-3=0 excluded its displayed W=V_dec split.
MISTAKE-486/490 and the current THM-3818 scope were checked; all-scale decoder
typing and the actual finite physical box are retained.

Let Q=91^6=567869252041. Consider a primitive thirteen-speed row in the actual
THM-3818 finite box, with positive distinct coordinates and sum at most Q^2.
Assume its decoder graph has exactly the two components of sizes eleven and
two, and its full support-at-most-three, height-at-most-Q relation span W
equals the rank-eleven decoder span V_dec. Write its primitive shapes and
scales as

    n = t V disjoint-union g(p,q),
    t,g positive integers, gcd(t,g)=1,
    |V|=11, gcd(V)=1, 1 in V,
    1<=p<q, gcd(p,q)=1,

where (p,q) belongs to THM-3818's 5,855-entry inert-prime atlas. In particular
p+q<=356, p<=177, and q<=355. The literal unit in V means that the minimum
of the physical eleven-component equals its gcd t. It does not require the
whole physical row to contain the speed one.

**Theorem.** Every such row n has a common phase of clearance at least 1/14.
Consequently an unresolved counterexample in this actual 11+2 equality branch
must have minimum at least two in its primitive eleven-component shape.

The portfolio is the LRC anchor alongside the contemporaneous jet/Boolean and
no-three-line work. The concept board has five objects: bounded crossing rows,
literal normalized unit, exact scale quotient, complete cyclic grids, and
protected safe arcs. The connection preserves the actual component split and
integral coefficient budget before passing to the scale inequality. The unit
and the coprime scale permutation are essential sidecars, not disposable
normalizations.

## 2. Bounded division and the all-parameter proof

First put K=max V. There are eleven distinct positive entries, so K>=11 and
the unit and maximum are different columns. THM-3818 (15q) gives
K<=Q, because the reduced internal pair (1,K) has height K. This is inherited
from the actual finite-box equality branch, not asserted from abstract rank
or from gcd(V)=1 alone.

The elementary bounded-division compiler is slightly stronger than ordinary
Euclidean division. For integers 1<=x<=Q(K+1), 2<=K<=Q, put

    a=min(floor(x/K),Q),       r=x-aK.                 (1)

Then 0<=a,r<=Q and x=aK+r. If the quotient is at most Q, its remainder is
less than K; otherwise a=Q and x<=QK+Q bounds the final remainder by Q.
The inclusive endpoint x=Q(K+1) uses a=r=Q. Thus it must not be omitted
by insisting that the last remainder remain strictly below K. This box is
sharp for the two fixed coordinates 1,K: |aK+r|<=Q(K+1) for all signed
coefficients of magnitude at most Q. No sharpness claim for the full
eleven-coordinate relation set is made.

For any coprime pair p<q, a phase with both clearances at least 1/3 exists.
An elementary choice is z=j/(p+q), where

    p j = floor((p+q)/2) modulo (p+q).

Both pair clearances are floor((p+q)/2)/(p+q)>=1/3. The pair therefore has a
closed 1/14-safe arc of length at least

    2(1/3-1/14)/q = 11/(21q).                        (2)

THM-3818's cyclic gluing, with the two components interchanged as necessary,
now closes whenever

    11t>=21q.                                       (3)

Otherwise 11t<21q, so the integer core scale satisfies

    t<=floor((21q-1)/11)<=677<Q.                      (4)

Set delta=gcd(t,p). Since gcd(t,g)=1, this also equals gcd(t,gp). The
integer c=t/delta is at most Q. If gp/delta<=Q(K+1), apply (1) to
x=gp/delta. It produces the literal relation

    c(gp) - a(tK) - r t = 0.                         (5)

Its support is at most three, its height is at most Q, and it meets both
prescribed components. Zero a or r merely drops the support. Its partial
sum on the pair component is the nonzero number cgp, so it lies outside
V_dec, whereas its support and height put it in W. This contradicts
W=V_dec. Thus an actual equality row in the small-t branch must obey

    gp/delta>Q(K+1),
    g>delta Q(K+1)/p>42K,                            (6)

where the final inequality follows from delta>=1, p<=177, and
Q>42*177. This is an exact all-parameter comparison, with no height census.

Cited LRC for eleven nonzero speeds supplies a phase y for V with every
clearance at least 1/12. Its closed arc of radius 1/(84K) is 1/14-safe,
so one component of G(V) has length at least

    2(1/12-1/14)/K = 1/(42K).                        (7)

Take the g physical lifts x_j=(z+j)/g of a safe pair phase. In the primitive
core clock they are

    t x_j = t z/g + t j/g modulo one.

Coprimality of t and g makes these a complete shifted g-grid. By (6), its
spacing is strictly smaller than (7), so one point lies in the core-safe
arc. That physical x_j is simultaneously safe for tV and g(p,q), proving
the theorem. This argument is exactly the symmetric (15aa) gluing, with the
previously missing quantitative scale separation supplied by (5).

The large-t branch only uses the weak closed-arc statement, including an
equality in (3). The theorem is therefore stated at the weak target 1/14;
it does not make an unnecessary all-case strictness claim. In the small-t
branch the strict spacing in (6) can select the interior of (7).

## 3. Physical unit, d=3 consumer, and the old virtual families

If the physical eleven-component itself contains one, then t=1 and delta=1.
The crossing compiler immediately gives gp>Q(K+1), and the second branch
alone closes. This specialization holds for every physical 11+2 shape
satisfying the actual entry assumptions, without any ternary parity condition.

In the d=3 chart write the eleven physical core as 3A union T, with |A|=8,
|T|=3, and one in that physical core. Put M=max A, E=max T and
K=max(3M,E). If the decoder pair is 3h(p,q), then g=3h and

    h>Q(K+1)/(3p)  implies  14h>29K.

Hence both 14h>=87M and 14h>=29E hold. The row therefore also enters
incoming **THM-4448**'s uniform attachment cone,
`01-canon/theorems/THM-4448-lrc14-general-shore-attachment-and-decoder-pair-cones.md`.
The direct cyclic proof above is more general and does not require this chart
or THM-4448's finite pair-component maximum.

For the ninth-round ratio-four virtual family, the proposed eleven physical
component is 3A union T with A the eight base frequencies. It contains one,
has K=42c, and the proposed pair is (3h,12h). Its prescribed split therefore
has the following exact scope:

* If K>Q, it fails the inherited internal pair-height condition of actual
  finite-box equality entry.
* If K<=Q and 3h<=Q(K+1), (1) gives a bounded crossing row using physical
  coordinates 3h,42c,1; the prescribed W=V_dec split is excluded.
* Any remaining actual equality entry has 3h>Q(K+1) and is already in the
  attachment regime just proved.

This covers the arbitrary-coprime-h virtual family as an entry audit, not
as an assertion that every member actually has that decoder graph. The
h=M+1 subfamily has the stronger height-one obstruction 3h-42c-3=0 already
recorded in the ninth audit. The old ratio-two pair (1,2) additionally fails
the inert-prime atlas itself. The constructive virtual-wall results remain
valid and retain their exact addresses and sharp margins, but they do not
leave an untreated live decoder subclass when the core contains one.

## 4. Nonvacuity and hostile normalization

The theorem does not say that a unit component makes actual entry impossible.
The canonical THM-4049/4117 row

    V=(1,4,6,8,10,12,14,15,16,18,22),
    t=1, g=2^45, (p,q)=(1,3)

is a genuine positive control. The companion independently reconstructs its
decoder graph with component sizes 11+2 and checks the finite box. Since
g>2Q max V, no bounded crossing row exists: two core terms cannot balance a
nonzero pair multiple of g, and two pair terms sum to a multiple of g that
one nonzero core term cannot balance. The connected internal decoder graphs
then give W=V_dec. Its safety was already proved for this specific core by
THM-4117/4129; here it prevents vacuous interpretation of the general theorem.

Conversely gcd(V)=1 cannot replace the literal normalized unit. THM-4052's
actual primitive eleven-component supplies the hostile. For

    R=(37,43,61,67,73,79,97,103,127), P=15 product R,
    V={2P/r:r in R} union {P/3,P/5},

one has gcd(V)=1, min V>1, and

    max V=237907127334685115>Q,
    max_(u,v in V) max(u,v)/gcd(u,v)=127.

The internal pair-height bound therefore coexists with a huge primitive
maximum when the unit column is absent. Dividing by the minimum would make
some coordinates nonintegral and destroy the bounded integer relation budget.
The first failed implication in that tempting extension is K<=Q, before
any gluing or owner argument is reached. Arbitrary primitive eleven-component
shapes remain outside this result.

## 5. Exact controls and reproduction

The standalone standard-library script uses **21,260 live gates**:
5,187 bounded-division controls including the inclusive endpoint; 3,865
cleared crossing controls at nontrivial component scales; all 5,855 inherited
pair atlas entries for the p<=177,q<=355,t<=677 inequalities; the canonical
actual-entry row; the primitive nonunit hostile; and 45 exact physical grid
witnesses across three core shapes and both proof branches. It never scans
to the large physical heights and does not infer actual entry from its generic
grid-witness bank. Both a complete residue permutation and literal physical
clearances are checked.

Run:

    python3 -B 04-computation/overnight11_20260906_lrc_unit_component.py
    python3 -B -O 04-computation/overnight11_20260906_lrc_unit_component.py

Normal and optimized outputs are frozen separately and must agree exactly.
The independent audit is coordinated with the orthogonal jet/Boolean lane;
its proof and literal control engine are separate from this producer.

Frozen source SHA256:
`d2fe7245e8d304557b036b72a03a76a46c06e871f60057e5a137892d1ce27bea`.
Frozen output and optimized-output SHA256:
`33ee43062b8f3035d2fa150972d9e817bfa72ec210f71b92db5f4fb282ea2608`.
The two output streams are byte-identical LF.

Independent audit: `overnight11_20260906_lrc_unit_component_audit.md` and its
separate standard-library source/output. All 131,958 referee gates pass with
byte-identical normal/optimized output; no mathematical repair was requested.
Referee source SHA256:
`37d2d41531770f047cb3979c214ea0a57b5a8ef8574d1c2b3f4570734b2ba51d`.
Referee output SHA256:
`22e866526bc150cff6aac2ce41b5610946f2ea799fdaf98f600c64855544dac7`.

**Filing:** root integrated these audited artifacts in the eleventh checkpoint;
reproduction commands above are relative to the repository root. Outside-worktree
locations preserve author provenance, not the present reproduction location.
