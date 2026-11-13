# Guarded descent, transported endpoints, and precise failure boundaries

**PROVED within the explicit templates and operation languages below;
FINITE-EXACT certificate classifications. Chan–Pak remains OPEN.**
Packet labels D1–D3,E1–E2,P1 are local, not reserved global theorem IDs.

## D1. A genuine decreasing recursion on nine infinite seed charts

Use the inherited definitions from
[R0, seed recursion](../seed_recursion_births_20261111/PROOFS.md):
`delta=H-W(O)` and `U=sum_(j in O)(s_j+sum_i Q_ij)`.
The complete paired endpoint cover is the source of W and U. A unit seed
has delta=1; if U=c*2^j with c odd, its physical tail supplies
`H=1+c*2^s` only for s>=j+1.

There are78 free-port presentations in the70 inherited birth charts.
For each stored source, replace the two specified coupled representatives
by TT_m and the free representative by TT_z, m,z>=1. The selected vertices
remain singletons. Existing fixed twins stay present; m therefore counts
the replaced representative, not the total size of its maximal twin block.
The exact signed tensor gives delta(m,m,z)=1 in all78 cases.

Nine presentations have

    U(m,m,z)=a*4^m+2^b B(m)*2^z,
    a=d*2^alpha, d odd,

where B(m) is positive, odd and strictly increasing at every m>=1.
The receipt records its complete integer exponential coefficients and a
nonnegative Stirling expansion. In that expansion F_r(m)=r!S(m,r), F_1=1,
and F_r(m+1)=r(F_r(m)+F_(r-1)(m)); the strictly positive F_2 coefficient
proves strict increase, even when the exponential expansion has negatives.

Set

    z=2m+alpha-b+k,  k>=1.
    j=2m+alpha,       c=d+2^k B(m).                    (D1)

All these z are positive. Since d is odd and k>=1, v2(U)=j exactly.
This is a physical parametrization, not a formal normalization at clock0.
The state retains the family, actual source rows, selected set, named
coupled/free blocks, m,k, and j. It determines every endpoint incidence;
the D,U tensors are its finite seed readouts, not a claim to encode every
possible future graph operation.

Two inverse moves are now valid:

* If k>1, delete one free twin. Then c'=(c+d)/2<c, j'=j.
* If m>1, delete one twin from each coupled block and two free twins.
  Then c'=d+2^k B(m-1)<c, j'=j-2.

Each move preserves the unit condition on the whole positive parameter
domain. Each strictly lowers both the coefficient c and seed order.
No unstable contraction is made in isolation: the coupled deletions are
synchronized. The first move reaches k=1; the second reaches m=1.
Thus every recognized state terminates at one of nine explicit base
presentations. Their coefficients are5553,3153,14433,16653,14435,10859,
8193,9453,8195, in receipt order. This is a constructor and decoder for
these families, not a recognition theorem for all odd numbers.

Recognition within each chart is exact: compute k=v2(c-d), require k>=1,
then search the strictly increasing sequence B(m) for (c-d)/2^k. Exponential
bracketing followed by binary search supplies m or certifies nonmembership.
Thus the family guard is executable; it is not an unspecified oracle.
The k-decrement boundary is exact: at k=1 the putative predecessor d+B(m)
is even. For k>1 it remains odd. Retaining this parity boundary is essential.

**Certificate proof.** Acyclicity of the three variable ports bounds runs
by the number of fixed vertices plus2. The producer enumerates full
endpoint tables and convolves two-path covers. The independent audit
enumerates one-cut words for D and, for U, uses a different decomposition:
Hamilton paths of the literal two-vertex tail are divided at fixed vertices,
with an acyclic port chain in each gap. Its entire coefficient polynomial
equals D+2U. Equality is of complete tensors, not interpolation. Thirty-six
literal endpoint computations include unequal coupled lengths as hostile
controls. The remaining69 presentations are retained as30 more general
binary-clock profiles and39 profiles with a third free-clock mode.

## D2. Four local changes of core with decreasing inverse targets

Two parent presentations, `even:0` and `even:5`, each have two eligible
order13 children. They are genuine one-vertex extensions of their stored
order12 source. The coupled representatives and selected set are retained;
the free length increases by one to align valuations.

For `even:0` put

    B(m)=-25+100*2^m+200*4^m+300*6^m.

Its source rows are
`[32,5,33,7,15,986,415,287,543,3295,2559,511]`, O={1,3},
coupled representatives6,10, and free representative4, all zero-based.
The parent has U=48*4^m+2B(m)*2^z, so alpha=4,b=1,d=3.
Child `odd:0` adds a vertex which beats only old vertex0. Its U is

    96*4^m+(26/5)*B(m)*2^z.

For fixed m,k the child's free length is one larger, its j is one larger,
and its odd coefficient is `c_child=3+(13/5)*2^k B(m)`. Therefore

    13*c_parent-5*c_child=24,
    c_parent=(5*c_child+24)/13 < c_child.               (D2a)

Child `odd:1` beats precisely old vertices0,2. Its B is3 times the
parent's B and its other offsets are the same as `odd:0`. Thus

    c_child=3*c_parent-6,
    c_parent=(c_child+6)/3 < c_child.                   (D2b)

For `even:5`, replace the +200*4^m in B by -100*4^m and use its free
representative5. Its children `odd:15,odd:16` satisfy exactly the same
two laws. Their whole tensors, actual added arcs and literal parent
restriction are stored. The independent compiler verifies both sides.

The ratio13/5 itself is inherited: the earlier
[selected-interval theorem, section2](../five_stress_connection_atlas_20261112/SEEDS.md)
already gives5=2+3 and13=2^2+3^2 for an allowed first gap insertion.
Here the new result is its four explicit synchronized-chart lifts with
positive normalized anchor3, transported valuation j, and all-length
endpoint certificates. We do not claim a new arithmetic origin of13/5.

These operations escape the twelve-command experiment: they add a vertex
which beats only one or two old vertices, rather than one which dominates
every forgotten vertex. The actual exceptional old vertices must remain
in the construction state. Ordinary seventeen-query equality does not
license this operation on an arbitrary representative.

**Five terminal charts and confluence.** Combining D1 with D2 leaves exactly
five terminal family types: even:0,even:5,odd:2,odd:3,odd:17, with respective
base coefficients5553,3153,14435,10859,8195. Every core move leaves m,k
unchanged, so it commutes with either coordinate decrement; the two coordinate
decrements commute with each other. These are all overlapping pairs of rules.
Termination and these diamonds give a unique terminal chart for every state
in this grammar. This restricted confluence does not undo the earlier J2
counterexample for arbitrary stable contractions.

If states retain their nine source labels, their seed-order generating
function is `(3z^17+2z^18+4z^19)/((1-z)(1-z^4))`. Its first coefficients
are3,5,9,9,12,14,18,18,21,23,27,27 at orders17–28. These count formal
construction states, not unmarked tournaments or distinct integers.
An independent capped integer Hamilton DP verifies the literal order20
parent and order21 child at external exponent8, obtaining1,421,569 and
3,694,849. Saturation at the requested target plus1 prevents overflow from
masquerading as equality.

## D3. The tail threshold is an essential guard

Keep an external ray exponent s and require s>=j+1. The D1/D2 inverse
moves lower j or leave it fixed, so the predecessor still has a physical
tail. Its count1+c'*2^s is strictly smaller than1+c*2^s.
For the m move the tail grows by two vertices as the seed loses four;
for a D2 move it grows by one as the seed loses two. Thus full graph order
also decreases. The known endpoint data is transported with the template;
it is not asserted numerically invariant under a change of core.

The guard is not removable. At even:0,m=1,k=1 the coefficient is5553
and j=6. This seed supplies the ray only from s=7 onward. The recursion
does not pay the six smaller exponents. A scalar satisfying the congruence
in D2a need not belong to the stated child family; that membership sidecar
is part of the rule, not a consequence of divisibility.
For example c_child=29 makes (5c_child+24)/13=13 integral and odd, yet
neither coefficient belongs to the asserted parent/child charts. This
does not claim that13 or29 lack other seeds; it rejects that unsupported
use of this particular recursive rule.

Only O((log X)^2) states in these nine families have c<=X: B(m) has a
positive dominant base6 and c grows geometrically in k. Four cross-core
links do not change that bound. This proves that the displayed grammar
alone cannot cover all odd coefficients. No old repair is declared closed.

## E1. Modes of the complete endpoint observer

The inherited [universal two-port observer](../adaptive_state_induction_20261103/COMMON_OBSERVER.md)
has17 static queries and exact rational updates for all12 commands.
For a fixed actual source whose unnamed TT clock has length a, write its
full query row as X(a)=sum_b b^a X_b. A fixed command word w gives
X(a)A_w=sum_b b^a X_b A_w. Therefore an absent mode cannot be created by
any word in that language. The statement concerns the complete row.
Here the same fixed word is used uniformly in a. An unbounded word chosen
as a function of a is outside this mode-preservation conclusion.

Exactly39 of the78 stored source/free-clock presentations have full
observer support within{1,2};39 reach base4. Naming the first and last
vertices of a split free module TT_(a+2), with its middle TT_a unnamed,
gives maximum base3 in all78. These are statements at the stored coupled
lengths. They do not classify all synchronized m.
Consequently every fixed-word continuation of the split-clock source bank
misses every whole target profile requiring base4 or5. A mode exclusion is
not an exclusion of an individual target integer at a chosen clock length.

An initial-H-only version of this assertion was mistakenly stated in
progress and immediately corrected. The explicit hostile is even:8,
clock10, roots0,4: the initial H has only1,2, but its one-command query is
`-133+2036*2^a+9090*3^a`. The complete observer is the necessary coordinate.

## E2. Even the complete ordinary observer forgets the selected seed gate

The presentations even:0 and even:2 use the same order12 row array.
The first selects O={1,3}, clock4; the second selects O={3,4}, clock1.
Both use roots6,10. Their complete17-query clock profiles are identical.
Hence every word in the twelve-command language gives identical H counts
at every positive clock length a. Yet their selected slopes at coupled
length1 are respectively

    U_0(a)=192+5550*2^a,
    U_2(a)=-9+192*2^a+5550*3^a.                        (E2)

At a=1 both have H=2663, but U values11292 and17025. Their two-vertex
unit tails therefore have counts22585 and34051. The second U is odd for
every a>=1; the first has different valuation and lower-head obligations.
The failure is precisely that the selected set is not a mark of the
ordinary two-port query language. Retain O and its paired cover data when
transporting a seed. This does not contradict the observer's completeness
for its declared twelve operations.

Among all66 pairs of equal observers,16 are already isomorphic preserving
the two ordered roots and clock. The other50 are isomorphic only after
forgetting some of those marks. None has a different unmarked source core.
These collisions must not be advertised as new unmarked core reductions.

## P1. A finite spectral barrier for root-changing periodic repairs

Restrict to the six one-exit commands with incoming mask1 or2. They add
x with v->x and x->every other old vertex, and may retain or replace the
two names. They preserve strongness. The inherited rooted endpoint cone
[GRAMMAR G1–G2](../source_selection_closure_20261103/GRAMMAR.md)
gives q>=s_v, h>=s_v+q and F>=h+q-s_v. Therefore

    h_new=h+F-s_v >= h+(h-s_v) >= 3h/2.

Consider any actual strongly connected source, a fixed prefix and suffix
of these insertion operations, and a fixed period w repeated ell times.
If its output is A+B*2^(r ell)+C*3^(r ell), with A,B,C nonzero, its period matrix
must have eigenvalues1,2^r,3^r. Also (3/2)^|w|<=3^r. For strides r=1,2,3
the length bounds are2,5,8. Multiplying the exact17-dimensional matrices
by2 clears denominators. Reduction modulo101 excludes the required
eigenvalues for all but the following words:

| Stride | Periods examined | Survivors |
|---:|---:|---|
|1|42|5,8|
|2|9,330|all four words of length2 on{5,8}|
|3|2,015,538|all eight words of length3 on{5,8}|

Codes5,8 keep both roots fixed. Thus **no period that replaces a root**
can give those three nonzero target modes on any of these three strides.
Surviving static-root words still require physical-source and endpoint
tests; they are not asserted impossible or successful. Strides>=4,
aperiodic choices, other operations, or a different target mode set remain
outside this certificate. In particular this is not an all-language
extension of the earlier O/Z/J theorem.
