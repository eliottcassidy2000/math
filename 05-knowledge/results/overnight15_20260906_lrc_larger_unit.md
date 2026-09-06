# A unit in the larger primitive component closes five actual decoder split types

**Status: PROVED RELATIVE TO THE STATED PROVED AND CITED SUPPLIERS +
INDEPENDENTLY AUDITED; FINITE-EXACT controls.** The
[independent referee](overnight15_20260906_lrc_component_closure_audit.md)
accepts the full proof including the exact native refinement, with 29,288
always-active gates and a separate phase-selection construction. The root
also independently read and accepted the proof. This closes a subclass of
the actual two-component equality branch. It does not close LRC(14), and it
does not assume every primitive component contains a coordinate equal to one.

Let Q=91^6=567,869,252,041. For an actual primitive thirteen-speed row in the
physical box, assume its decoder has two components of sizes a+b=13, a<=b,
and W_(Q,3)=V_dec. Write the row as

    tV union gU,   gcd(V)=gcd(U)=gcd(t,g)=1,
    |V|=a, |U|=b, 1 in U, v=min V, L=max U.

**Theorem.** The row has a common time of clearance at least 1/14 whenever

    v <= Qa/[7(b+1)].                                      (1)

This condition is automatic for a=1,2,3,4,5. Thus, when the larger primitive
component contains one, the splits 1+12, 2+11, 3+10, 4+9, and 5+8 all close.
For 6+7 the sufficient condition is

    min V <= floor(3Q/28) = 60,843,134,147.                  (2)

The remaining high-minimum 6+7 subclass and larger shapes without a unit are
not eliminated here. In particular, the connected-shape bound alone does not
prove (2) for every primitive six-component.

## 1. Inheritance, hostile, and the exact bridge

The closest proved mechanism is **THM-3818**, Sections 6.1, 6.4, and 6.5,
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`:
the connected primitive bound (15e), the actual crossing prohibition and
internal pair-height bound (15o)--(15q), and the coprime phase-grid lemma
(15y)--(15aa). Its finite-box supplier is **THM-2052**,
`01-canon/theorems/THM-2052-finite-height-forces-high-rank-bounded-relation-code.md`.
These arguments were read, rather than inferred from their headlines.

The audited [eleventh unit-component theorem](overnight11_20260906_lrc_unit_component.md)
closed the unit-containing eleven-shape in the 11+2 branch. The
[twelfth signed-box theorem](overnight12_20260906_lrc_gcd_semigroup.md)
retained the primitive pair and its physical gcd. Here its simplest complete
interval, the pair (1,L), suffices; no general coefficient-minimization
compiler is needed for the proof.

The new second supplier is the completed, independently audited
`05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.md`.
It proves that any primitive strict counterexample has every seven-label
gcd at most 90, and every eight-label gcd at most 30. Its stronger table for
subset sizes 12,11,10,9,8,7 is 1,2,4,9,30,90. We use the completed proof note,
not a reserved theorem namespace. These are restrictions on hypothetical
failure, not necessary conditions on every safe equality entry.

The lower-dimensional phase supplier is **CITED**: Touch Sungkawichai and
Tanupat Trakulthongchai, *Eleven, twelve, and thirteen lonely runners*,
[arXiv:2604.23906v2, Theorem 1.3, printed page 2](https://arxiv.org/pdf/2604.23906v2).
For k<=12 nonzero speeds it supplies a common clearance 1/(k+1). Its
computer-assisted proof is an external supplier, not a computation repeated
by our finite controls. For the singleton component the phase 1/2 suffices.

The corrected near miss is replacing actual entry by a convenient partition
or a rank label. The [thirteenth exact entry decoder](overnight13_20260906_lrc_entry_decoder.md)
shows why both orientations of mixed triples and the physical box matter.
Its out-of-box equality example has internal pair height greater than Q.
The less-used sidecar is the exact good-arc length in THM-3818 (15aa), applied
here to the larger component. The live board is bounded three-label rows,
physical gcd scales, unit-pair intervals, connected-tree cofactors, and
coprime translated phase grids.

The source-to-target map is explicit: a forbidden small-support relation
forces a lower bound on the smaller component's scale t; that scale is the
number of lifts of one smaller-component safe phase, and coprimality turns
those lifts into a complete t-grid in the larger component's clock. The
preserved predicate is every smaller-component clearance. The initial phase
and its selected lift are not retained, so the larger component's safe-arc
length is the required sidecar. No large-support Bezout relation is silently
converted into support three, and no marginal measure bound is used as a
coherent phase certificate.

## 2. Native entry and the positive-box crossing

All thirteen physical speeds are distinct positive integers, their collective
gcd is one, and their sum is at most Q^2. The components are those of the
actual THM-3818 inert-cube decoder atlas; the primitive edge coefficients
have sum at most 356 and maximum at most 355. W_(Q,3) is the rational span of
all integral zero relations of support at most three and coefficient height
at most Q. V_dec is the rational span of the actual decoder edge rows.

Connected incidence identifies V_dec with the direct sum of the two internal
weighted kernels. Consequently a mixed zero relation of support at most
three lies outside V_dec: its partial sum on the component containing only
one occupied label cannot vanish. Entry W_(Q,3)=V_dec excludes every such
relation. The finite-box argument in THM-3818 (15q) gives primitive height
at most Q for each internal pair. In particular, the two distinct U labels
1 and L give

    b <= L <= Q.                                           (3)

Here b>=7, so the unit and maximum really are separate labels. The primitive
scales are coprime because the entire physical row is primitive.

Put

    delta=gcd(g,v)=gcd(g,tv), c=g/delta, x=tv/delta.

Suppose c<=Q. If x<=Q(L+1), define

    s=min(floor(x/L),Q), r=x-sL.

Then 0<=r,s<=Q. Indeed, before the cap is reached r<L<=Q; after the cap is
reached r=x-QL lies in [0,Q]. Thus

    c(tv)-r g-s(gL)=0                                    (4)

is a literal mixed relation on at most three physical labels, of height at
most Q. Its distinguished coefficient c is positive, so it is nonzero and
is excluded by actual equality entry. We conclude

    c<=Q  ==>  tv/delta > Q(L+1),
    hence t > delta Q(L+1)/v.                            (5)

This argument needs no unproved implication from a collective gcd to a pair
gcd. It also retains the outside coefficient: x alone is not enough if c>Q.

## 3. The phase-grid consequence and its exact native form

The k-speed supplier gives phases eta and zeta such that all V-clearances
at eta are at least 1/(a+1) and all U-clearances at zeta are at least 1/(b+1).
The U minimum-clearance function is L-Lipschitz on the circle. It is therefore
at least 1/14 on the closed arc centered at zeta of radius

    [1/(b+1)-1/14]/L = a/[14(b+1)L].

This arc has full length

    ell = a/[7(b+1)L].                                    (6)

The t physical phases x_j=(eta+j)/t, 0<=j<t, preserve every V-clearance.
Their U-clock images are g eta/t+g j/t modulo one. Since gcd(g,t)=1, these
are exactly a complete translated t-grid. Every closed arc of length 1/t
meets the grid. An arc of length strictly greater than 1/t contains a grid
point in its interior, where the Lipschitz estimate gives strict clearance
greater than 1/14. The smaller component has strict clearance as well since
1/(a+1)>1/14.

Combining (5) and (6) gives a useful stronger **native sufficient criterion**:

    g/gcd(g,v) <= Q,
    7(b+1)Lv <= a gcd(g,v) Q(L+1)                         (7)

imply a common time of strict clearance greater than 1/14. Equality in the
second inequality is allowed: the inequality forced in (5) is strict.
The simpler condition (1) implies its second line for every delta>=1.

To prove the stated theorem, suppose the row had no common weakly safe
time. Since U has b>=7 labels, the incoming joint-gcd theorem forces g<=90;
indeed g divides the gcd of any seven of its physical labels. Thus c<=90<Q
and (1) gives (7), a contradiction. Equivalently, if g>90, the incoming
theorem already excludes failure; if g<=90, the explicit argument above
constructs a strictly safe time. The overall theorem asserts weak safety,
because the g>90 supplier is only used with that conclusion.

## 4. Automatic split types and the remaining boundary

The connected-tree cofactor bound (15e) of THM-3818 gives

    min V <= max V <= 355^(a-1).

For a=1 the empty tree minor is one and the primitive singleton is V=(1).
For a>=2 each maximal minor is a product of a-1 positive edge coefficients;
dividing their common gcd produces the primitive shape and can only lower
the bound. The exact comparisons are:

| a+b | Connected maximum bound | Integer minimum cutoff from (1) |
|---|---:|---:|
| 1+12 | 1 | 6,240,321,451 |
| 2+11 | 355 | 13,520,696,477 |
| 3+10 | 126,025 | 22,124,776,053 |
| 4+9 | 44,738,875 | 32,449,671,545 |
| 5+8 | 15,882,300,625 | 45,068,988,257 |
| 6+7 | 5,638,216,721,875 | 60,843,134,147 |

The first five comparisons establish the automatic range. The last is a
failure of this sufficient comparison, not evidence of an unsafe row. It
leaves the exact native refinement (7) available in the balanced case. No
claim is made that a unit in the smaller component can replace the required
unit in U; that is a different placement and needs a different argument.

## 5. Nonvacuity, hostile controls, and reproducibility

The standalone verifier uses exact integers and rational phases. Its seven
actual equality entries have g=1 and t=2QL+1. Each has two connected actual
atlas components, a primitive distinct physical row inside Q^2, and no mixed
bounded row: in any mixed zero relation of support at most three, the smaller-component
partial sum would be a nonzero multiple of t, whereas the larger-component partial
sum has absolute value at most 2QL<t. It also independently checks every
mixed triple support in both orientations using exact signed-box membership.
Their counts for a=1,...,6 are 66,121,165,198,220,231. Each entry has an
explicit rational strictly safe physical phase printed in the transcript.

All six split types occur. The smaller shapes are (1), (2,3), (2,3,6),
(2,3,4,6), (2,3,4,6,9), and (2,3,4,6,9,12). Larger shapes are initial
segments of (1,4,6,8,10,12,14,15,16,18,22,24). These are actual graph
controls, not guessed partitions. They avoid obtaining nonvacuity solely
from an already excluded large-g branch.

The additional 5+8 stress uses the connected primitive unitless shape

    V=(1013861907,1036929995,1060507875,1084612635,1096868145).

It is the star with center P=179*181*183*185 and leaves
(356-q)P/q for q=179,181,183,185. Every displayed edge has reduced sum356,
which is in the actual atlas. With U=(1,4,6,8,10,12,14,15), g=1,
t=30Q+1, the full physical sum is

    90,168,220,083,627,413,785,737 < Q^2.

Thus the automatic criterion applies to a genuine equality entry whose
smaller primitive minimum exceeds one billion. Its literal core phases are
eta=1/2, zeta=1/9; the selected strict physical phase is printed in the output.

The finite universe consists of all 210 pairs 2<=L<=B<=21 and all 31,395
nonnegative points up to B(L+1); 8,325 coprime scale/coordinate controls; 3,773
exact arc-hit controls; and the seven full actual entries. Two explicit
hostiles retain the proof boundaries: B=2,L=6 has x=3 missing from its signed
box despite |x|<=B(L+1), and t=4,g=2 produces only {0,1/2}, missing the arc
[1/10,2/5] although its length exceeds 1/4. The first shows why an internal
height condition is needed; the second shows why a grid cannot be inferred
after dropping the coprime-scale condition.

[Source](../../04-computation/overnight15_20260906_lrc_larger_unit.py) and
[normal output](overnight15_20260906_lrc_larger_unit.out):

```
python -B 04-computation/overnight15_20260906_lrc_larger_unit.py
python -B -O 04-computation/overnight15_20260906_lrc_larger_unit.py
```

Both pass 54,930 always-active gates. Normal and optimized transcripts are
byte-identical LF. No finite bank is being used to infer the all-parameter
phase or entry statements; those are proved in Sections 2--4.

Frozen source SHA256:
`23f42086c09e367bf9ddc87958587901b338a8a27d712e34e6eaea64ba8ce935`.
Frozen normal/optimized output SHA256:
`3b8c091e3e387ce0cacc24e22690bd78ecc9c3d313fffd4fcc0f4f1127106b36`.
The independent source SHA256 is
`72dc83b007422f079d5663686fa7312013669acf89dfabb348b7f7c326ee2bf2`,
and its normal/optimized output SHA256 is
`409ec49be081fa208382d37d2e065e36bc93620180c5966501e21fb5509f7979`.

**Filing:** root read the complete proof and independent audit. All source and
output bytes are retained; reproduction commands are relative to repository root.
