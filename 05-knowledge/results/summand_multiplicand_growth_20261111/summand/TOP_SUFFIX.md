# Eventual complete saturation of strict synchronous seed growth

Date: 2026-11-11. **PROVED:** every finite seed set of at least two distinct
positive integers eventually fills its entire attainable closure below its
moving maximum. The transition has a finite constructive bound; no claim
is made that this bound is close to the minimum transition time.

## Exact statement and inherited input

Let S be a finite set of at least two distinct positive integers. Set

    S_0=S,
    S_(t+1)=S_t union {a+b:a,b∈S_t, a<b},
    C(S)=union_(t≥0) S_t,
    M_t=max S_t,   d=gcd(S).

Then there is a finite T such that for every t≥T,

    S_t=C(S) intersection [1,M_t],                     (1)
    M_(t+1)=2M_t−d.                                   (2)

The previously proved input is [CLOSURE.md](CLOSURE.md): the exact mixed-seed
ideal, cofinite normalized closure, a leaf-count arrival bound, and the
consecutive-block chase giving eventual top three spaced by d. The current
note closes the remaining gap between a complete lower interior and a
complete top suffix. It does not assume those two growing intervals meet.

Normalize by dividing all seed values by d for Sections 1–3, so gcd(S)=1.
The original scaling is restored in Section 4.

## 1. Top three force top four after finitely many steps

Suppose generation T3 has top three `B3−2,B3−1,B3`, where B3≥3. From
then onward these top-three values persist in moving coordinates and
`M_(t+1)=2M_t−1` by CLOSURE.md §5.

For an available vertex x below M_t, write its deficit as `q=M_t−x>0`.
The legal pair x,M_t gives a vertex in the next generation with deficit

    (2M_t−1)−(x+M_t)=q−1.                             (3)

Thus any available deficit q≥3 can be decreased by one at each step until
it reaches 3. The corresponding vertex, together with the persistent top
three, then gives top four.

There is always such a deficit after at most one step. Retain the old
vertex x=B3−2 in generation T3+1; its new deficit is
`(2B3−1)−(B3−2)=B3+1≥4`. Applying (3) another B3−2 times reaches deficit 3.
Therefore one may take

    T4≤T3+B3−1.                                      (4)

An existing top fourth vertex gives a much earlier T4. Only existence and
an explicit finite bound are needed. Write B=M_T4. At this stage
`[B−3,B]⊆S_T4`, B≥4, and for every j≥0,

    M_(T4+j)=(B−1)2^j+1.                             (5)

The transient can be long even with only three seeds. For
`S={A,A+1,A+2}`, A≥2, the first top-four generation is exactly A−1. At
generation 1 the smallest deficit above 2 is A+1. Before it reaches 3,
sums using only deficits 0,1,2 create no deficit above 2; a sum involving
any other deficit has smallest possible new deficit q−1, attained by
pairing deficit q with 0. Retained old vertices have new deficit at least
M_t−1>q−1. Thus the smallest deficit above 2 is exactly A+2−t for
`1≤t≤A−1`, proving the claim. The maximum when top four first appear is
`(A+1)2^(A−1)+1`. This is an exact transient, not a lower bound inferred
from a truncated graph.

## 2. Four consecutive top values fill a fixed-threshold suffix

Put

    L=B−2,  U=B−1,  I_h=[hL−1,hU+1] for h≥1.

The initial interval is `I_1=[B−3,B]`. We claim that for every j≥0,

    I_h⊆S_(T4+j) for every 1≤h≤2^j.                  (6)

The claim retains a complete family of integer intervals, rather than
only their union's density.

**Interval addition lemma.** If p,q≥1, the strict sums of values in I_p
and I_q contain I_(p+q).

If p=q, the strict sums of an integer interval `[a,b]` of length at least
3 form exactly `[2a+1,2b−1]`. Here this is
`[2pL−1,2pU+1]=I_(2p)`.

If p≠q, the two intervals have different lower endpoints and different
upper endpoints because L,U>0. Every integer strictly between the extreme
ordinary sums has at least two ordered representations, one with a value
from each interval: the feasible first coordinates form an integer
interval, whose size is at least two away from the two extreme sums.
At most one of those representations has equal coordinates. The two
extreme sums themselves have unequal parents because the respective
endpoints differ. Thus every ordinary sum is also a legal strict sum,
and their full interval is

    [(p+q)L−2,(p+q)U+2],

which contains I_(p+q). This proof applies even when the two input
intervals overlap. QED.

To prove (6), retain I_1 and, for 2≤h≤2^j, split h into its two balanced
positive parts p=floor(h/2), q=ceil(h/2). Both are at most 2^(j−1).
The interval addition lemma and induction complete the construction.

Adjacent intervals I_h and I_(h+1) overlap or are consecutive exactly when

    (h+1)L−1≤hU+2,
    h≥L−3=B−5.

Set

    h0=max(1,B−5),   K=h0(B−2)−1.                    (7)

For every j with 2^j≥h0, the union of I_h over h0≤h≤2^j is one integer
interval. Equations (5)–(6) therefore give

    [K,M_(T4+j)]⊆S_(T4+j).                           (8)

The lower endpoint K is fixed in time. This is stronger than a suffix
whose length is merely a positive fraction of the maximum. In particular
the permanent closure contains every integer at least K.

## 3. The finite lower remainder completes saturation

Only finitely many attainable values lie below K. Every one eventually
enters by the definition of C(S), so (8) already proves (1). A constructive
bound follows from the mixed-seed expression theorem in CLOSURE.md §1.

Let m0=min S. Every attainable nonseed n has a mixed expression with at
most floor(n/m0) seed occurrences, and hence a strict sequential genealogy
of depth at most floor(n/m0)−1. Seeds have depth zero. Consequently every
attainable value below K is present by

    T_low=max(0,floor((K−1)/m0)−1).

Thus the explicitly finite choice

    T=max{T4+ceil(log2 h0), T_low}                    (9)

ensures (1) for every t≥T. The inclusion S_t⊆C(S) is automatic. Values
below K are permanently present when attainable; all values from K to
the current maximum are present by (8). No claim about minimum genealogy
or minimum T has entered the argument.

This also supplies a finite alternative to the leaf-count bound: use the
exact Apéry membership compiler to list the attainable values below K,
construct one genealogy for each, and take the largest of their depths
and `T4+ceil(log2 h0)`.

## 4. Exact eventual cardinality and target arrival time

Return to the original seeds with gcd d. Let h be the finite number of
positive integers absent from the normalized closure C(S/d), computable
exactly from the Apéry state in CLOSURE.md §3. The specific T in (9)
satisfies `M_T/d≥K`, because `T≥T4+ceil(log2 h0)` in normalized coordinates.
Equation (8) puts every permanent normalized hole below K. Thus this
chosen maximum lies past all permanent holes; finite saturation alone
would not establish that. Therefore, writing
N_t=|S_t|, for every t≥T,

    N_t=M_t/d−h,
    N_(t+1)=2N_t+h−1,                               (10)
    N_(T+j)=(M_T/d−1)2^j+1−h.                       (11)

The distinction between normalized integers and all integers is essential:
when d>1, nonmultiples of d never enter and are not counted as the finite
hole number h.

For every attainable n>M_T, its exact first arrival generation is

    depth_S(n)=T+ceil(log2((n/d−1)/(M_T/d−1))).       (12)

Indeed before that generation the maximum is below n, and at that first
generation (1) guarantees n is present. Thus the earlier `log2 n+O_S(1)`
depth estimate eventually sharpens to an exact rounded-logarithm law.

## Scope, controls, and transfer contract

The closest inherited mechanism is the top-three shrinking-window chase.
The missing coordinate was a whole indexed interval family I_h: retaining
only a bottom fraction and a top fraction does not show that their gaps
close. The useful operation is balanced addition of intervals, including
an explicit repair for the forbidden equal-parent diagonal.

Source→target: a stable top-four interval→the family I_h→a suffix with
fixed lower threshold. Preserved: actual legal parent witnesses and exact
maximum. Lost: other vertices, representation counts, and minimum onset.
Sidecar: the interval index h and common stage bound h≤2^j.

Source→target: fixed-threshold suffix plus the finite attainable lower
remainder→eventual full saturation. Preserved: exact membership at every
later generation. The permanent holes are retained by the Apéry state;
they are not filled by the conclusion.

The single-seed case is excluded: no strict sum exists. The immediate
top-three claim is also insufficient: e.g. the state {100,101,102} does
not contain its fourth top value. Step (3) proves that it will eventually
appear; the interval-family argument begins only after that event.

Nor does the first occurrence of finite saturation imply permanence.
For seeds {2,3}, the states at generations 0,1,2 are already equal to
`C(S)∩[1,M_t]`, but generation 3 has maximum 15 and misses attainable 14.
Even adding stable top three to that premature criterion is insufficient:
seeds {5,6,7} are saturated at generations 0 and 1 and have top three
forever, yet generation 2 has maximum 25 and misses
`21=5+(5+(5+6))`, which arrives at generation 3. These examples are why
the proof retains the fixed suffix threshold and the finite lower remainder.

`python3 -B top_suffix_check.py` and its `-O` version give separate
finite controls: 832 interval additions; all 28 two-seed subsets of[1,8]
and ten named multiseed states through two generations beyond their safe
saturation certificate; and exact fixed-window top-four transients for
A=2,3,7,25,100. The last control never materializes the exponentially
large full state. The receipt does not substitute for the all-seed proof.
