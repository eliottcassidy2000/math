# Independent transitive-profile audit

Status: PASS, including the finalized producer proof and its incoming-priority
qualification. The elementary propagation formulas below are independently
derived. This audit imports no producer code.

The seed subsets are classified by forbidding cyclic triangles, using a
least-vertex removal recurrence, rather than searching for a universal
source. Every root out-neighborhood is checked. Every seven-vertex
neighborhood Hamilton count is counted by all 7! literal orders. All five-
vertex motif classes are independently expanded over permutations, giving
a disjoint cover of all 2^10 labelled tournaments. The full global and
neighborhood five-motif histograms agree with the producer. Global H and
ordered two-path-cover counts are recomputed by forward subset DP.

Both seeds have the complete transitive-subset polynomial

    F = 1+15z+105z²+315z³+375z⁴+150z⁵+15z⁶,

and every out-neighborhood has polynomial

    G = 1+7z+21z²+25z³+10z⁴+z⁵.

Their neighborhood Hamilton counts are 101 versus 89. Every neighborhood
contains respectively zero versus one regular five-vertex tournament,
checked directly by five internal scores equal to two, without using
canonicalization. Thus the rooted source-plus-regular-five motif counts
are 0 versus 15, also verified by literal global six-subset enumeration. All four-vertex motif counts agree; the failure is not a
failure of every local motif statistic.

## All-size identities checked algebraically

Let L_r(T)=C3[T,{a},P_r], with T→a→P_r→T and P_r transitive. Put
P=(1+z)^r. A transitive set cannot meet all three blocks. Inclusion-exclusion
of the three two-block choices yields

    F_(L_r(T)) = (P+z)F_T + z(P−1).

For an old root v, its out-neighborhood is T[N⁺(v)] followed by the sink a,
so its polynomial is (1+z)G_v. The root a sees P_r, giving P. A vertex of
P_r with j later vertices sees their transitive block followed by T, giving
(1+z)^j F_T, for j=0,…,r−1. Therefore equality of F and the full multiset
of out-neighborhood polynomials is preserved by every r≥1.

These formulas have been checked directly on every subset of both lifted
seeds at orders 17,18,19,20, independently of the proof and producer's
smaller controls.

The same audit records a stronger, precisely typed companion: if C_v is
the generating polynomial for transitive subsets containing v, then old
roots have (P+z)C_v; root a has z(F_T+P−1); and every P_r root has
z(1+z)^(r−1)(F_T+z). Both seed roots have

    C_v=z+14z²+63z³+100z⁴+50z⁵+6z⁶.

These contained-root identities were also checked on all literal lifts.
They do not claim preservation of arbitrary associations between different
unlabelled root statistics: such a joint claim must retain that association
or use the homogeneity of these particular seeds.

Source padding has F→(1+z)F, preserves old out-neighborhood polynomials,
and adds the new root polynomial F. The universal source supplies no new
cyclic triangle. Thus the equality statements survive source padding.
The strong cyclic lift preserves exact E by the independently derived
incidence update in the propagation packet. None of this says a finite
transitive profile determines the full tournament or its Hamilton count.

## Files and reproduction

    python -B check.py --producer PATH
    python -O -B check.py --producer PATH

PATH contains `check_transitive_profile.py` and
`transitive_profile_receipt.json`, and `TRANSITIVE_PROFILE_OBSTRUCTION.md`.
Without PATH the workspace sibling
`session38-root` is used. Default checks compare the frozen receipt;
`--write` intentionally regenerates it. No external solver or unrecorded
input is used. `receipt.json` pins producer hashes and contains all seed
values, literal lift polynomials, root-profile hashes and check counts.
