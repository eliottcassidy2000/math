# Split-head colors: the exact compatibility requirement

**PROVED semantic relaxation.** The old common head/tail shift condition is
sufficient but not necessary for the weak B invariant. It may safely be dropped
if head colors and tail shifts are each compared across the two output words.
This does not alter or weaken the validity of any previously frozen certificate.

## Typed input and ownership

Let input paths B_n and B_(n+1) have the same position parity for every common
label, with label1 at position1. Let q be odd, R=q², h=(R−1)/2. A source label
x at one-based input position j is lifted in offset c to

    L_c(x)=R*x+(-1)^j*c,     -h<=c<=h.

The two input representations give identical lifted values for every x<=n,
because their position parities agree. The exceptional new input label n+1
has position parity n+1: compare the counts of odd-position labels in the
length-n and length-(n+1) paths after canceling all common labels.
Thus selecting member i(c)=1 iff R+(-1)^(n+1)*c<=d gives exactly the top partial
band needed for output interval [1,R*n+d].

The band of x is [R*x−h,R*x+h]. These bands are pairwise disjoint; offset c
is also unique inside a band once the input parity of x is fixed. Low singleton
labels1..h fill the initial gap. Consequently every common output label has a
unique source type and source label. There is no identification between a head
in one offset and a tail in another.

## Necessary and sufficient profile comparison for these words

A prefix cut separates input label1 from the remaining source path. Its head is
R−c, in band1=[h+1,3h+1]. Every tail label is in a band x>=2, starting at3h+2.
Let H_(b,c) be the head's **actual output position parity** in output word b.
Let S_(b,c) be the shift such that each tail label x has output position parity

    input_position_parity(x) + S_(b,c) mod2.

Orientation and atom lengths determine S. The compatibility obligations are:

1. every low singleton has the same output parity in the two words;
2. H_(0,c)=H_(1,c) for every split head;
3. S_(0,c)=S_(1,c) for each tail or uncut whole block.

These conditions imply common-output-label position parity agreement by the
unique typed decomposition above. Conversely, wherever the corresponding
common source set is nonempty, agreement of the actual label parities implies
these profile equalities. The input sizes here are at least42, so tails have
many common labels and this nonempty qualification is automatic.

**There is no requirement H_(b,c)=1+S_(b,c) within a word.** That additional
identity says the separated head happens to retain the same shift as its tail.
When it fails, the head's own color still agrees between outputs; every other
label's color still agrees between outputs. The next recursive lift requires
only this per-label agreement and therefore remains valid.

This is a semantic condition on B pairs, not a claim that an arbitrary
underconstrained port model gives a Hamiltonian path. Exact atom ownership,
length, endpoint, connector, prefix and connected-traversal checks all remain
mandatory. `B/audit_independent.py` verifies the independent-head certificate family
without the within-word identity, while still requiring complete equality of
the head/singleton/tail profiles between outputs and independently expanding
both complete paths on a real compatible base. The producer's earlier generic
checker used the explicit `independent_head_colors: true` certificate flag.

The new Q25 family uses H=(1,3,6,19), protected offset8, and independent
shared head colors. Its prefix recurrence is

    P_next=(1,3,6,19) followed by (25*p_j+(-1)^j*8).

Its all-state existence is a separate certificate obligation, not a consequence
of this semantic relaxation lemma alone.

That separate obligation is discharged by the complete B catalogs and
`B/INDEPENDENT_AUDIT.md`; the resulting extension theorem is stated in
`TERNARY_SCHEDULE_THEOREM.md`.
