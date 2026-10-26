# Antihydra valuations: exact coding, integer ancestry, and barrier debt

2026-10-26, Session64. **PROVED:** the elementary statements below for the
explicit arithmetic maps. **CITED:** the identification of the parameter
map with Antihydra's native machine. **FINITE-EXACT:** the accompanying
rational/integer checker and a 10,000-return macro prefix. **OPEN:** whether
the particular Antihydra blank-tape orbit ever hits its barrier; positive
Collatz convergence. Neither question is answered by a full-shift coding,
a generic measure statement, a decreasing real coordinate, or this census.

## 0. Inheritance and scope

The primary [arXiv:2610.22603v2, section 2.2](https://arxiv.org/html/2610.22603v2#S2.SS2)
gives the parameter transitions and blank configuration used below. The
repository's [source audit](../normseven_boundary_20261025/PAPERS_AND_TRANSPORT.md)
already flags the displayed negative native-clock value at the halt-family
boundary. We do not use that clock or assert that the machine formalization
has been rebuilt. The macros' arithmetic targets are our explicit input.
THM-4770/4771 were still RESERVED source stubs, not proved dependencies.

Closest proved mechanisms, not new discoveries:

- [Session56 valuation trie/LCA, sections 2–3 and 6](../session56_rmq_lca_connections_20261025/valuation-lca/REPORT.md):
  least-significant digit prefix depth is the valuation of a labelled
  difference; a rational marked path can measure finite repeated-word fuel.
- [Session56 chronological debt monoid, sections 2–5](../session56_rmq_lca_connections_20261025/debt-range/RESULT.md):
  preserve both total change and the minimum of every prefix; exact source
  cylinders require the final oddness bit and their own clock/guard data.
- [Collatz bugs, section 2](../collatz_bugs_20260925.md): a forward-meeting
  LCA is available on the known root basin after deleting its cycle, not
  automatically on all positive integers.

Canonical hostiles here are the singular numerator, the negative Collatz
fixed point, the nonintegral twisted fixed point, and distinct integer
inverse-branch patterns. The live board is: coding domain / valuation-prefix
metric / actual common future / chronological resource / integer realization /
measure versus one marked seed. The new contribution is their exact comparison
for these two maps, especially integer injectivity and the Antihydra barrier.

This is **not** the tournament counting [reverse valuation potential](../port_residual_integral_20261025/counting/REVERSE_POTENTIAL.md).
That potential is a minimum-cost gauge on a finite acyclic transfer graph.
The barrier recurrence below is an inherited ordered-debt specialization,
and the Haar harmonic function is a third, probabilistic object. Their shared
method is to retain the precise state and local valuation before reducing
coordinates; there is no map transporting tournament counts to these orbits.

## 1. Exact odd section and barrier (PROVED arithmetic; native bridge CITED)

The parameter rules, with nonnegative integer a,b, are

    (a,2b)       -> (a+2,3b+2),
    (a,2b+1)     -> (a-1,3b+3) if a>0,
    (0,2b+1)    -> HALT.

Set x=b+4, where b now denotes the actual second coordinate. Even x gives
x'=3x/2 and a'=a+2; odd x gives x'=(3x-1)/2 and a'=a-1, provided a>0.
The blank parameter state (0,4) takes three even macro steps to (a,x)=(6,27).
For odd x let k=v_2(3x-1). One odd step and k-1 even steps give

    H(x)=3^(k-1)(3x-1)/2^k,
    a'=a+2k-3.                                      (1)

The odd guard a>0 is sufficient for the whole block: after its first step
the counter is a-1>=0, and every subsequent step increases it by two.
The returned odd state is halting precisely when a'=0. Thus the unguarded
integer map H always continues, even when the actual parameter process has
already halted. That extension must not be confused with the machine orbit.

For valuations k_0,k_1,... and A_n=sum_(i<n) k_i,

    a_n=a_0+2A_n-3n.                                (2)

The blank process survives at its n-th odd section iff all preceding guards
and a_n>0 hold. To **execute** n blocks one needs the guards only at j<n;
to be safe also at the returned section one needs them at j<=n. Every barrier
statement below explicitly uses the latter, endpoint-inclusive convention.
The block clock is n odd returns and A_n parameter macros, not A_n native
Turing-machine transitions. Native clocks remain outside this packet.

## 2. Exact domains of full-shift coding (PROVED)

Let O=1+2 Z_2. Define ordinary accelerated Collatz on O away from -1/3 by

    C(x)=(3x+1)/2^k,  k=v_2(3x+1).

Define H by (1) away from +1/3. Both valuations are positive integers where
they are finite. Their inverse branches, for k>=1 and y in O, are

    c_k(y)=(2^k y-1)/3,
    h_k(y)=1/3+(2/3)^k y.                            (3)

Each branch maps O bijectively onto the class with exact first valuation k,
and contracts 2-adic distances by 2^-k. Let c_w,h_w mean compositions in
chronological valuation order, so the first letter is the outermost branch.
Remove the complete finite preimage trees of the singular points:

    E_C={c_w(-1/3): w any finite positive word},
    E_H={h_w(+1/3): w any finite positive word},
    Omega_C=O\E_C,   Omega_H=O\E_H.                  (4)

The empty word is included. These are countable exceptional sets, consisting
of trajectories whose next numerator eventually vanishes. A valuation of zero
is infinity, not a positive finite letter. The accelerated maps are not
endomorphisms of every odd 2-adic integer.

**Coding theorem.** Each itinerary map Omega_C or Omega_H -> N_{>=1}^N is
a homeomorphic conjugacy with the one-sided left shift. For an infinite word
k and its partial sums A_0=0,A_1,..., the inverse codes are

    X_C(k)=-sum_(j>=0) 2^A_j/3^(j+1),
    X_H(k)=(1/3) sum_(j>=0) (2/3)^A_j.              (5)

**Proof.** Iterating (3) gives the corresponding finite sum plus a tail
multiplied respectively by 2^A_n/3^n or (2/3)^A_n. Since A_n>=n, this tail
tends to zero 2-adically. The first summand is odd and all later summands are
even, so the limit belongs to O. Applying the inverse identities to each
infinite tail gives the prescribed finite valuation at every stage; hence the
limit avoids (4). Conversely every nonsingular orbit has these finite affine
identities, whose tails vanish 2-adically because their endpoints remain odd.
Thus the limit is unique, establishes surjectivity, and intertwines the maps.
The cylinder descriptions below prove continuity in both directions. QED.

A word w of length n and valuation sum A determines exactly one residue
class modulo 2^(A+1), not just modulo 2^A: its image is c_w(O) or h_w(O).
Its finite exceptional tail must still be removed for the infinite domain.
All odd rational integers (positive or negative) lie in both domains: integer iterates cannot hit
+/-1/3. Negative integers and noninteger odd-denominator rationals are not
extraneous; they are essential members of the coding spaces.

The usual compact binary-parity shift for an unaccelerated shortcut map is a
different coding assertion. Here the alphabet is countable and the punctured
odd domains are noncompact. No all-Z_2 accelerated full-shift claim is made.

## 3. A genuine isometry, but two different meanings of LCA (PROVED)

Suppose two unequal infinite valuation words share n letters with sum A,
and their next letters are p!=q. The difference of their next inverse-branch
values has valuation min(p,q): its two nonconstant terms have distinct
valuations and its common constant cancels. The shared prefix multiplies the
difference by a 2-adic number of valuation A. Therefore in **both** maps

    v_2(X(k)-X(l))=A+min(p,q).                      (6)

Coincident infinite words give the same point and difference valuation infinity.
Consequently X_H composed with X_C^(-1) is a 2-adic isometry, not merely a
topological conjugacy. It preserves every labelled residue-prefix collision
and therefore the Session56 binary-trie LCA depth. Counting shared valuation
letters alone is not the binary depth: their accumulated sizes and the first
differing letters in (6) are required. Finite precision supplies only a lower
bound if the first separating bit has not yet been observed.

A common **forward future**, in contrast, means C^i(x)=C^j(y) or H^i(x)=H^j(y)
for some i,j; in codes this is equality of two shifted infinite tails. The
conjugacy preserves that relation on its full domains, but common prefixes do
not characterize it. Collatz 3 and 13 have different first valuations, 1 and 3,
yet both map to 5. The points 1 and -5/3 have the same first valuation 2,
but map respectively to the distinct fixed points 1 and -1, so they never meet.
A rooted dynamical LCA additionally requires a specified common rooted basin,
a cycle convention and a chosen clock. None follows from a prefix trie.

## 4. Integer geometry is not preserved (PROVED all odd integers)

For constant valuation words,

    X_C(k^infinity)=1/(2^k-3),
    X_H(k^infinity)=3^(k-1)/(3^k-2^k).              (7)

Thus Collatz 1 (word 2^infinity) maps to twisted 3/5, and twisted 1
(word 1^infinity) corresponds to Collatz -1. The conjugacy does not preserve
positive integers, integer height or the native parameter condition x>=4.

**Integer inverse theorem.** If y=H(x) for an odd integer x, then

    r=v_3(y)=k-1,
    m=y/3^r,
    k=r+1,
    x=(2^(r+1)m+1)/3.                              (8)

There is such an odd integer predecessor iff m=(-1)^r modulo 3; when it
exists it is unique. Indeed 3x-1 is prime to 3, so (1) gives r=k-1. Solving
for x proves necessity, and the congruence makes the displayed x integral
and odd. Then 3x-1=2^(r+1)m has precisely the required 2-adic valuation,
proving sufficiency. The proof works for negative integers too.

H(1)=1. For odd integers x>=3,

    H(x)>= (3x-1)/2 > x;

for odd integers x<=-1, multiplication of the negative (3x-1)/2^k by
3^(k-1) gives H(x)<= (3x-1)/2 < x. Thus every positive integer other than 1
escapes to +infinity, every negative integer escapes to -infinity, and 1 is
the only integer cycle. The inverse uniqueness makes each nontrivial integer
component an unbranched ray: inverse steps strictly approach 1 or 0 and hence
terminate after finitely many integer steps. No two distinct branches merge.
This sharply differs from Collatz's integer sibling merges. In a twisted
integer component, two points meet forward only when one lies on the other's
ray; a nontrivial branching LCA cannot be transported from integer Collatz.

Every periodic H word gives a rational number in (1/3,1]. One proof sums
its geometric-series fixed point; another uses the real version of (5):
A_j>=j bounds the sum by 3, with a strict positive tail after 1/3. A finite
prefix of the maps h_k preserves that interval, so the same holds for every
eventually periodic word. Hence the H valuation word of every integer other
than 1 is **not eventually periodic**. The converse rational->eventually
periodic assertion is neither used nor claimed. Also, for a nonperiodic word,
its real sum and its 2-adic sum live in different completions; one may not
infer that the 2-adic integer encoded by (5) lies in that real interval.

## 5. Exact backward barrier and inherited range summaries (PROVED all words)

For a finite word w define delta_j=2k_j-3 and partial change D_j=2A_j-3j.
The endpoint-inclusive least integer initial resource is

    R(w)=1-min_(0<=j<=n) D_j
        =1+max_(0<=j<=n)(3j-2A_j).                 (9)

Equivalently,

    R(empty)=1,
    R(k followed by w)=max(1,R(w)+3-2k).          (10)

These follow directly from a+D_j>=1 at every retained section. They are
necessary and sufficient for this fixed word, not just sufficient drift
bounds. Keeping the empty prefix makes R>=1 and retains the initial guard.

This is exactly Session56's additive debt monoid with strict integer
threshold: summarize w by (Delta,m), where Delta=D_n and m=min_j D_j. Then

    (d,m)*(e,n)=(d+e,min(m,d+n)),  R=1-m,          (11)
    R(uv)=max(R(u),R(v)-Delta(u)).                (12)

Ordered disjoint range composition is valid; omitting intermediate minima
or replacing chronological composition by an overlapping RMQ rule is not.
The word (1,2) has changes (-1,+1), total zero and R=2; (2,1) has the same
total and R=1. Endpoint drift alone forgets an actual halt.

For the repeated word w^infinity, let Delta=Delta(w), m=m(w). The resource
needed to keep every section positive is R(w) if Delta>=0, and no finite
resource suffices if Delta<0. This is an all-word classification of the
**formal periodic barrier**, not an integer realization theorem. In particular
1^infinity has negative resource drift but twisted point 1 is not a native
x=b+4 state. All finite bad prefixes are nevertheless integer-realizable:

    x_0=1+2^(a+1),  word 1 repeated a times,
    H^a(x_0)=1+2*3^a.

Starting counter a>=1 reaches zero at that returned odd section. For native
x>=5 take a>=1; x_0>=5. Thus positive integer realization alone is no
universal barrier-survival principle. It says nothing adverse or favorable
about the particular marked source (6,27).

## 6. Exact affine height ledger and direction (PROVED)

For H, let L_n=(2/3)^A_n x_n. Finite inverse composition gives

    L_n=x_0-(1/3) sum_(j<n)(2/3)^A_j,
    L_(n+1)=L_n-(1/3)(2/3)^A_n.                   (13)

For x_0>1, the exact real inequalities are x_0-1<L_n<=x_0. Consequently
x_n>(x_0-1)(3/2)^A_n. This decreasing real coordinate is bounded away from
zero and is not well-founded. It proves growth of x_n, not termination or
resource survival. All comparisons can be performed as rational/integer
comparisons without logarithmic roundoff.

Ordinary Collatz instead has

    (2^A_n/3^n)x_n=x_0+sum_(j<n)2^A_j/3^(j+1).    (14)

If B_n=sum_(j<n)3^(n-1-j)2^A_j, exact n-step descent is

    (2^A_n-3^n)x_0 > B_n.                          (15)

The coefficient inequality 2^A_n>3^n alone is not sufficient at small
height: x_0=1 and k=2 has coefficient 3/4 but is a fixed point. Neither
these identities nor the coding law prove positive Collatz convergence.

## 7. Haar harmonic potential and exact ruin probability (PROVED measure scope)

Use normalized Haar measure on O. A length-n word has cylinder mass 2^-A_n;
the exceptional sets (4) have measure zero. Therefore valuation letters are
independent with P(k)=2^-k under either itinerary coding. This is a measure
statement on 2-adic starting points, not independence along a specified
ordinary integer orbit.

Let the resource begin at an integer a>=1 and evolve by increments 2k-3,
stopping on first arrival at zero. These increments are bounded below by -1,
so it cannot jump over zero from above. Put

    q=(sqrt(5)-1)/2,   0<q<1,   q^2+q=1.

The kernel is exactly

    E[q^(2k-3)] = 1/[q(2-q^2)] = 1.              (16)

Thus q^a is harmonic for the live counter. **The exact eventual ruin
probability is q^a.** Here is a bounded stopping proof, including the boundary
at infinity. Fix integer M>a and stop on reaching 0 or >=M. This stopping
time is almost surely finite: from any state 1,...,M-1, M consecutive k=1
letters force exit downward, and their probability is 2^-M independently of
the past. The stopped q^a process is bounded between 0 and 1. The finite-time
martingale identity and bounded convergence yield

    P(hit 0 before >=M) <= q^a
       <= P(hit 0 before >=M)+q^M.                (17)

On upward exit the counter is >=M even with overshoot. As M increases, the
lower events exhaust eventual ruin: any finite ruin path has finite maximum.
Taking M to infinity proves the claim, without assuming a law of large
numbers or excluding escaping paths. For resource 6 it is q^6=9-4sqrt(5).

The mean formal increment is 1. The exact multiplicative coefficient drifts
under the same Haar law have log means log(3/4) for C and 2log(3/2) for H;
these signs describe coefficient statistics only. Equation (15) still retains
the affine carry, and the native Antihydra barrier still retains its prefixes.
The particular source x=27 is a Haar-null point: neither its survival nor its
ruin follows from (16)–(17). No probability is assigned to the truth of its
halting question.

## 8. Closed exact controls and reproduction

Run in this directory:

    python3 -B check_valuation.py --output report.json > normal.json
    python3 -O -B check_valuation.py > optimized.json
    cmp normal.json optimized.json

All load-bearing checks raise explicit exceptions, so optimized Python does
not erase them. The standard-library implementation uses integers and exact
Fractions only. It checks every word of length 0..5 over valuations 1..4:
1,365 words, 38,232 branch-letter round trips, 2,730 cylinder checks, 2,730
finite singular prefixes, 2,728 periodic-code closures and 12,285 barrier
threshold comparisons. There are 510 cross-map exact isometry comparisons,
all word cuts for debt concatenation, and 2,000 signed integer inverse checks.
The brute inverse enumeration uses letters 1..12 and targets |y|<=999;
(8) independently bounds the only possible letter in that target universe.
Negative/noninteger periodic points, exact shared-future witnesses, same-prefix
disjoint cycles, resource-zero boundaries, and bad positive integer prefixes
are explicit hostile controls.

The finite-horizon Haar ruin computation aggregates all infinitely many
letters that cannot reach zero by the remaining horizon into a safe sink,
rather than silently dropping or renormalizing a finite alphabet. Its exact
masses at horizons 6,12,24,48 are recorded in `report.json`; each is checked
below 9-4sqrt(5) by squaring a rational inequality with its sign retained.

From the blank-derived odd section (6,27), 10,000 odd returns were independently
replayed both through (1) and through all 19,939 individual parameter macros.
The final resource is 9,884; its minimum over odd sections is 6 at the starting
section (`report.json` field `min_a` uses this odd-section convention). Over all
parameter macros from that first odd section, the minimum is 5, attained after
the first odd macro. The final x has 11,669 bits. The valuation CSV and final integer hex encoding have
SHA-256 receipts in `report.json`. The result is **FINITE_PREFIX_ONLY**, not
nonhalting. The source and proofs carry the declared unbounded statements;
the finite checks are not extrapolated beyond their universes.
