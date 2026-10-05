# Repair a supplied common-future family by cutting its child path

2026-10-04. **PROVED:** child-cut slope and guard laws, a sufficient repayment
test, and the same-pair continuation obstruction. **FINITE-EXACT:** the stated
supplied-certificate examples and bounded hostile probe. **OPEN:** universal
availability of a useful cut and universal grounded Collatz coverage.
No literature-priority claim.

## 1. Inheritance and the changed operation

The closest proved mechanism is the least-period family lift in
[inverse shields and lifted joins](collatz_join_shields_and_lifts_20261004.md).
Its canonical hostile is an actual smaller child that produces an expanding
lift: source 233, child 231. The corrected near miss is that moving a join
later is not the same operation as changing the child. The least-used
sidecar is the source guard's ternary depth, which changes when a prefix is
removed from the child route.

The [binary/ternary fusion](collatz_binary_ternary_guard_fusion_20261004.md)
already bounds the finitely eligible sibling addresses of a supplied
integer. It proves a natural density strictly below one for that single-stage
bank. Its dense ternary address space does not imply source coverage.
The [grounded observation union](adaptive_observation_union_20261004.md)
already keeps actual unfinished edges and available suffixes. This note
extracts a further parametric consequence from such supplied suffixes.

The live board is: the original source; a certified child; exact words;
family slope; mixed guard period; child suffix; first-hit root. The anchor
is source-directed proof reuse, the niche is removal of excess ternary guard
depth, and the wildcard is whether every expanding pair admits a repairing
cut. We do not increase an inverse-depth cap or treat a larger convergence
census as the result.

## 2. Exact child-cut theorem

Use \(U(n)=\operatorname{oddpart}(3n+1)\). Suppose checked actual words give

\[
 n\xrightarrow{w}J\xleftarrow{v}m,\qquad 0<m<n,
\]

where \(w,v\) have lengths \(r,s\) and total valuations \(A,D\).
Empty words are allowed. Retain a supplied first-hit root certificate for
\(m\). The inherited least simultaneous positive periods are

\[
 P=2^{A+1}3^{\max(s-r,0)},\qquad
 Q=2^{D+1}3^{\max(r-s,0)},\qquad
 \lambda=Q/P.
\]

They give the same exact words from \(n+Pt\) and \(m+Qt\) for every
integer \(t\ge0\). The child remains smaller than the original source at
every height when \(\lambda\le1\). For \(\lambda>1\), a successful seed
does not imply that uniform source-size decrease.

Choose a cut after \(p\) edges of the checked child prefix, \(0\le p\le s\).
Let its value be \(x\), valuation sum \(C\), and affine coefficient
\(\mu=3^p/2^C\). If \(x<n\), the retained child suffix supplies a new
diagram

\[
 n\xrightarrow{w}J\xleftarrow{v[p:]}x.
\]

**PROVED.** Its lift slope and least source period are

\[
 \lambda'=\lambda\mu,\qquad
 P'=2^{A+1}3^{\max(s-p-r,0)},\qquad
 \frac P{P'}=3^{\min(p,\max(s-r,0))}.                 \tag{1}
\]

The new child period is
\[
 Q'=2^{D-C+1}3^{\max(r-s+p,0)}.
\]

These formulas follow by substituting the remaining child length and cost
in the least-period theorem. In particular, \(P'\) divides \(P\). The
source progression of the new diagram contains the old progression:
the old parameter \(t\) becomes \((P/P')t\). Removing a child prefix thus
relaxes only a power of three in the source guard; the original source
word and its binary exactness are unchanged.

If \(p>0\) and \(x<m\), then

\[
 x=\mu m+\frac{B_{\mathrm{prefix}}}{2^C}>\mu m
\]

with positive carry. Therefore \(\mu<x/m<1\), and the cut strictly
improves both the child label and the lift slope. It need not improve the
slope far enough. The exact all-height test remains \(\lambda\mu\le1\).

A useful sufficient repayment test for a formerly expanding lift is

\[
 \lambda>1,\qquad x\le m/\lambda
 \quad\Longrightarrow\quad \lambda'<1.              \tag{2}
\]

It is sufficient, not necessary. It uses a checked value on the supplied
child path, so it needs no new Collatz query. The child certificate at
\(x\) is literally the suffix of the supplied certificate at \(m\).

The operation retains every eligible cut, each with its own child identity,
word, period and coefficient. A smaller label alone does not justify dropping
other alternatives. Formula (1) compares source applicability for this
particular fixed source word; it is not a dominance theorem over every
possible proof representation.

## 3. An expanding lift becomes a larger decreasing family

The inherited actual route

    231,347,521,391,587,881,661,31,47,71,107,161,121,91,137,103,155,233

has length 17 and valuation sum 27. Using the empty source word at 233
gives

\[
 \lambda=\frac{2^{27}}{3^{17}}>1,\qquad
 P=2\cdot3^{17}=258280326,\quad Q=2^{28}.
\]

The seed child is smaller by two, but the very next positive parameter
already has child larger than source.

Cut the child route at 31 after seven edges of valuation sum 14.
The remaining word has length 10 and sum 13. Formula (1) gives

\[
 \lambda'=\frac{8192}{59049}<1,\qquad
 n_t=233+118098t,\quad m_t=31+16384t,\quad t\ge0.     \tag{3}
\]

The child follows
\((1,1,1,1,2,2,1,2,1,1)\) to \(n_t\).
The source period is smaller by exactly \(3^7=2187\).
Consequently (3) supplies an all-height decreasing dependency on a strictly
larger source progression. It consumes a supplied certificate for \(m_t\);
the certificate for seed 31 does not certify every free parameter.

This repair changes the child. Keeping child 231 and moving the common
endpoint down the shared route never changes the expanding slope.

**Partial-repair hostile.** Source 75 and child 73 first share the root.
Cutting the child's route at 71 after five edges decreases the child label,
but the new slope is

\[
 \frac{18014398509481984}{16677181699666569}>1.
\]

It is an improvement, not an all-height proof that the child stays smaller.
The exact coefficient
test prevents that invalid promotion.

## 4. Why moving the same pair's endpoint cannot fix its lift

For a supplied positive odd first-hit root certificate define

\[
 \alpha(n)=3^{R(n)}/2^{A(n)},
\]

where \(R\) is its odd rank and \(A\) its total halving count.
At every common point of two such routes, cancellation of their shared
remaining suffix gives

\[
 \lambda(n,m)=\frac{\alpha(n)}{\alpha(m)}.           \tag{4}
\]

Equivalently, appending the same checked word to both sides of any join
cancels its length and cost from the period ratio. For fixed certified
sources \(n,m\), every common endpoint on their first-hit routes has the
same slope. An expanding pair therefore has an all-length obstruction
to being repaired merely by moving its endpoint later.

Unequal formal root padding can change the ratio: adding one root
valuation-two edge on the source side multiplies it by \(3/4\).
That requires the source already to have reached the root, and lies
outside strict first-hit words. It is not an extension algorithm for an
unresolved source. The program rejects such padding.

Nor is \(\alpha(n)\) an oracle for a supplied unresolved \(n\). Its definition
here requires a supplied completed route. The local form of (1) needs only
the checked source prefix and the supplied child route; after their join
has been checked, it can construct the source certificate. Equation (4)
then explains the fixed-pair limitation and helps organize learned families.

## 5. The stubborn supplied states and the role of 257

The inherited all-smaller-source audit proves that 7, 27 and 703 first
meet a smaller completed route at odd steps 4, 37 and 51. Their respective
endpoints are 5, 23 and 157. Our control instances retain children 3, 15
and 123, reaching those endpoints in 1, 1 and 4 steps. These are already
checked finite routes; the present construction does not discover an
earlier escape or change the inherited lower bounds. Their exact diagrams
do lift to decreasing arithmetic families.

The value 257 has a specific role in the grounded graph:

    171 -> 257 -> 193 -> 145 -> 109 -> 41,
    27 -> 41.

The smaller ancestor 171 does not become a root certificate merely because
171 is below 257. A proof attempt using 257 to ground 171 and 171 to ground
257 would be circular. Once the actual shared suffix at 41 is grounded,
the checked chain is reusable. For example,

    515 --(1,4)--> 145 <--(2,2)-- 257

is an exact half-child join with slope \(1/2\). Its child premise must
still be available. No modulus, primality or special threshold property
of the integer 257 is being asserted.

The inherited inverse shield already shows why arbitrary additional
inverse depth is not a universal remedy at the first reset-two endpoint.
The source threshold and original integer remain attached to every
candidate. The new cut operation improves some discovered diagrams;
it does not remove that shield or guarantee a discovered diagram.

## 6. Finite controls, interface and stopping reason

Program:
[supplied_state_join_extensions_20261004.py](../../04-computation/experiments/supplied_state_join_extensions_20261004.py).
Saved output:
[supplied_state_join_extensions_20261004.out](supplied_state_join_extensions_20261004.out).

    python 04-computation/experiments/supplied_state_join_extensions_20261004.py
    python -O 04-computation/experiments/supplied_state_join_extensions_20261004.py

The join object accepts an exact source prefix, the actual supplied child
certificate and the cut position on that certificate. Validation checks
the original source, endpoint equality and strict source-size decrease. The suffix
enumerator keeps every cut below the original source. The attachment
function constructs the source AST by retaining the child suffix and
prepending the checked source word; it never calls the source encoder.
Source expansions pass the supplied certificate's own conservative codec
bit bound. A one-edge certificate beyond 10000 bits is an explicit control
against accidentally imposing the codec's default export cap. This is
finite exact arithmetic, not a claim of constant-time access to huge integers.

Controls include the full 233/231 repair, its 2187-fold guard relaxation,
positive lift parameters through one million, all 30 common strict
endpoints of that pair, the incomplete 75/73 repair, the four supplied
state examples, and five rejected endpoint/type/root-padding cases.

A separate cheap hostile probe supplies completed routes for all 501 odd
integers through 1001, then examines every ordered pair \(1\le m<n\le1001\).
For each pair it uses the first intersection encountered on the source
route. Exactly 2914 of those pair diagrams have expanding slopes; every
one has a cut on that particular child's prefix with value below the
source and nonexpanding repaired slope. In 2822 cases the common endpoint
is already below the source, so selecting that endpoint is simply direct
descent. The remaining 92 have endpoints at least as large as the source
and need an earlier child cut. These are the nontrivial cases for the
finite signal; they include the 233/231 example.

This is a finite signal only. No full routes from that probe are inputs to
the join APIs or seeds of a new unresolved-source solver. It raises a
sharper question: must every such expanding pair admit a repairing cut
before the given common endpoint? No proof is supplied. The earlier
first-coalescence question only asks for some smaller witness; requiring
a witness on this particular child's prefix is stronger. A larger census
would not answer either question, so this lane stops at the exact cut law,
its constructive repair and its declared open boundary.

Connection contract: the source is a checked join with a supplied child
certificate; the target is a collection of guarded affine families.
The map cuts the child's checked prefix while keeping the original source
word. It preserves the common future and the actual grounded suffix.
It changes the child, slope, word lengths and ternary applicability guard.
Required sidecars are the cut position, exact words, source periods and
certificate provenance. Dropping those sidecars would confuse a point
proof, a decreasing dependency and a completed family.
