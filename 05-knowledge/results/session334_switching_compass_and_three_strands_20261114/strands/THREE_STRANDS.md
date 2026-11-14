# Three strands: the visual coordinate and its exact boundary

**Proved** arithmetic fibre identities, phase balance and periodicity,
and an exponential obstruction to literal Zeckendorf compilation.
No claim of full unmarked count coverage.

## 1. Provenance and arithmetic strands

The distinct-positive-parent convention, the startup holes {1,4,6}
from {2,3}, and the exact difference from the ascending complete order
(just the doubling forest) are inherited from
[THM-2422](../../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md)
and the [summand-braids repair](../arithmetic_braids_20260917_summand.md).
An apparent Fibonacci depth reset in the [early reflection](../../../07-reflections/summand-graph-fermat-zeckendorf.md)
was a cutoff artifact. Ordinary transitivity fails exactly when two
legal short arcs concatenate to a missing doubling arc: x<y<2x gives
such a failure. Generation needs the labelled hyperedge {x,y}->x+y,
not ordinary one-parent reachability.

The historical R/K/B diagonals with extra copies of 1 are in the
[reset-colors note](../reset_20260926_colours.md),
[crossroads note](../crossroads_crossing_20260926_colour.md), and
[modular Fourier/Zeckendorf routes](../modular_fourier_zeckendorf_routes_20261003.md)
section 4. The supplied 35-color word differs at 35 from the Wythoff
candidate; no uniquely established global continuation. The display
uses row r, column c, cell r(r-1)/2+c, diagonal index r-c+1.
An elementary *separate* strand chart for strict summands is

    x=1+a, y=2+a+b, n=x+y=3+2a+b, a,b>=0.

Color by gap y-x=b+1 mod3. Advancing x by one within a row changes the
gap by -2, hence the three colors cycle. Increasing n by6 adds three
strict pairs, one per color. The unordered parity lens OO,OE,EE gives
another genuine three-way chart: odd rows are OE, even rows alternate
OO and EE. These are not identifications with the old R/K/B colors.

For an exact Zeckendorf guard let alpha=phi^(-2), b(n)=floor(alpha(n+1)),
and delta(x,y)=b(x)+b(y)-b(x+y), from the inherited
[Duck/Zeckendorf carry note](../duck_zeckendorf_20260925.md).
For fixed n let B=floor(alpha(n+2))-b(n) in {0,1} and
tau=frac(alpha(n+2)). Summing floors at fixed real total gives:

    delta(x,n-x)=B-[frac(alpha(x+1))>tau].

Thus a given row contains at most two of the three possible global
carries -1,0,1. The row type determines which pair is allowed. The
geometric stripe and carry color are different quotients.

Aliquot descent is a marked projection of the same summand row:
mark left slot x iff x|n, plus the excluded equality x=n/2 if present.
These are exactly the proper divisors; read the SUM of marked labels,
not #marks. With n=3^b m, 3 not dividing m, write S_r for the sum of
proper divisors congruent to r mod3 (gap color is r+n). Then

    S_0=(sigma(3^b)-1)*sigma(m)-[b>0]*n
    S_1-S_2=product_{p^e||m} sum_{j=0}^e (chi(p)*p)^j
             -[b=0]*chi(m)*m,

where chi(1 mod3)=1, chi(2 mod3)=-1. Expand Euler factors, removing
the improper divisor. For n=2^A3^B, B>0, the signed aliquot imbalance
is (1-(-2)^(A+1))/3, unbounded, while the signed *count* imbalance is
0 if A is odd, 1 otherwise. This is a concrete demonstration that
even precise color balance can lose the absolute target. Relevant
history: the [aliquot fixed-point carrier](../../../07-reflections/perfect-aliquot-fixed-point-carrier-s645.md)
and [dyadic three-cycle exclusion](../dyadic_aliquot_three_cycle_exclusion_20261019/DYADIC_THREE_CYCLE_EXCLUSION.md).
No positivity lift from aliquot flow is assumed.

## 2. Three Hamiltonian runs and phase rank

In a two-hub word fix the hub order and lock the four adjacent
first/last endpoints or empty sentinels. Every nonanchor vertex j
then chooses independently from its eligible runs A_j subset {0,1,2}.
Take the phase to be the sum mod3 of spine-run labels. In occupation
coordinates phase is i_1+2*i_2 mod3, hence colors by diagonals. In a
path of full length t with hubs at 1-based positions p<q, the phase is
2t-p-q-1 mod3: first hub precedes t-p-1 spine vertices, second t-q.
This description is independent of hub naming.

At a primitive cube root zeta, the phase evaluation of one frame is

    zeta^(anchor charge) * product_j sum_{r in A_j} zeta^r.

A singleton contributes a character unit; two choices contribute
-zeta^(missing label), still a unit; three choices annihilate the
imbalance. Thus even one triple slot is sufficient for exact
three-way balance (cycle that vertex's choice); in a zero-triple
frame the phase counts are equal except for one signed unit, sign
(-1)^(number of binary sites), at the phase given by anchor/singleton
charges plus the missing labels. The kernel of Z^3->Z[zeta] is the
constant vector. Every global imbalance therefore comes solely from
zero-triple frames, even if most Hamiltonian mass is ternary.

Refine the phase to sum_j r_j*a_j, a_j in F3^d. Exact Fourier support
of a frame is the annihilator of V=span of the columns at its triple
sites: three-choice character sums are 0 or 3, smaller sums never
vanish. Each assignment at the nonternary sites yields one coset of
V, uniformly with fibre 3^(#triple-sites-dim V). Full 3^d-way balance
holds iff the triple columns span F3^d. This is the F3 selector
analogue, not an identification, of the F2 rank/cycle obstruction in
the [tree gauge](../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md).

Example a_j=(1,F_j,F_{j+1}) mod3. The consecutive Fibonacci pairs
traverse all 8 nonzero vectors of F3^2 from (0,1) in period8 by
(u,v)->(v,u+v). Three noncollinear triple-site residues mod8 suffice
for 27-way balance. Any four distinct such residues suffice since
a line has three points. This uses a different charge than the old
mod2 Walsh charge with three nonzero colors. With a different basis
column for each address residue j mod3^d, the Fourier support
recovers exactly the address residues with triple sites.

## 3. Exact phase periodicity, but not count periodicity

Fix any skeleton of homogeneous spine blocks, and vary positive
block lengths. In the phase enumerator let a choice of run r carry
weight zeta^r. For an occupied set S of runs in one block of length m,
the weighted surjection number is

    W_m(S)=sum_{U subset S} (-1)^(|S|-|U|)*(sum_{r in U} zeta^r)^m.

A fixed occupancy pattern is legal independently of the lengths;
block order determines adjacency across blocks, and each hub sees
only the common letter on the adjacent first or last block. Every
subset sum is zero or a sixth root. Thus phase differences are
exactly 6-periodic in each positive block length, no asymptotic or
Jordan error. For a higher-rank character the same argument works
in any homogeneous block whose charge column has nonzero pairing
with that character; trivial pairing instead leaves the ordinary
unweighted powers 1,2,3 in that block.

For the macro rail 11 0^k 2^l 311, H=1+36*2^k*3^l. Differences of
phases 0 and 1 from phase2 at (k,l)=(1,1),(1,2),(1,3) are (0,1),
(-10,-10),(19,0). An isolated excess does not follow from H=1 mod3.
The same differences persist along each of 36 residue classes while
the absolute count diverges. The replay three_strands.py compares
the phase transfer with literal run assignments on all 1365 words
of lengths 0 through 5 and gives these macro controls.

## 4. Exact analogy, then a positional-compiler lower bound

If an integer n has w canonical Zeckendorf atoms, assigning its atoms
surjectively to three labelled nonempty strands gives exactly
3^w-3*2^w+3 distinct ordered positive triples summing n, with canonical
supports disjoint and union canonical. Every subset/union is still
canonical, so all regrouping carries vanish. Uniqueness of the
canonical support makes the assignments injective. This is EXACTLY
the homogeneous Hamiltonian occupancy-surjection polynomial. The
output n is the sum of atom values, though, not the number of choices.

A direct positional compiler of Zeckendorf words into a fixed
two-hub word is exponentially expensive. Encode m canonical binary
positions by m distinct literal spine sites, allow arbitrary fixed
per-site maps from run labels to bits (including discarded choices),
and project away all auxiliary positions. Suppose the projected
language of legal paths is precisely the no-adjacent-1 words.
Every endpoint frame projects to one Cartesian product of binary
lists, possibly empty; the existence of a completion in the
auxiliary sites is coordinatewise in a locked frame. There are at
most 2(s+1)^4 frames for a spine of length s.

A rectangle in the legal bit language covers at most one maximal
independent set of the m-vertex path: the coordinatewise union of
two supports in the rectangle is again in it, and a maximal support
cannot be strictly extended. Conversely the subsets of each maximal
support are one legal rectangle, and these cover. Thus this is the
exact rectangle cover number f_m. It obeys f_m=f_(m-2)+f_(m-3) with
f_0=1,f_1=1,f_2=2: if site1 is in, remove sites1,2 and solve the
suffix; if site1 is out then site2 is in (to dominate it), forcing
site3 out and leaving the suffix from site4. Hence f_m is
Theta(rho^m), rho^3=rho+1, rho>1, and s=Omega(rho^(m/4)).

It remains true for an exact canonical lexical/value prefix of
m digits that contains all m-1 digit words with fixed leading zero
(restrict the rectangles to that zero). For ordinary Zeckendorf
weights F2 through F_(m+1), prefix [0,N) with N>=F_(m+1) does so.
This does not rule out a compiler of the same NUMBER of paths,
a nonlocal coding, or a changing core; it rules out literal local
realization of the full neighbouring-exclusion language by
polynomial padding. To make an arithmetic inverse pay, the
adjacency guard must come from outside this Cartesian architecture.

### Arbitrary words and a genuinely growing discrepancy

Since every frame contributes either zero or a unit to a phase
difference, any two phase counts differ by at most 2(s+1)^4 for a
spine of length s. Thus even beyond fixed skeletons the phase channels
track H/3 to polynomial error. The error need not be bounded. For
the word (012)^j, let the two phase differences from phase2 be (D0,D1).
The exact law is

    j=2m,   m>=1: (2m-m^2, -m(m+1))
    j=2m+1, m>=0: (-m^2+2m+2, -m^2+m+3).

Reproduction/proof: for each hub order keep the three run states
(untouched/current end disallowed/current end allowed) and phase mod3,
so an upper bound of 2*3^3*3=162 states suffices. Repeat the fixed
6-letter block 012012. Each phase difference is a linear form of
powers of its single 162-dimensional transfer matrix, starting at
the appropriate initial prefix. Subtract the stated degree-two
polynomial by adjoining a 3-dimensional unitriangular shift.
The residual has a recurrence of order at most 165 by Cayley-Hamilton;
the explicit zero checks for 165 consecutive starting terms in
three_strands.py prove identity for all m in those intervals. This is
an all-order, finite-dimensional recurrence certificate, not extrapolation
from a short tabulation. The calculation uses integer arithmetic only.

### Cost of adding more hubs to the literal compiler

For k hubs there are at most k!*(s+1)^(2k) frames: fix the hub order and
each incident run endpoint (two per hub; outside end positions require
no extra anchor). A legal consecutive-hub adjacency is just a check
in that frame. Every projected frame is still a rectangle for
coordinatewise per-site labels. Therefore the same lower bound now is
k!*(s+1)^(2k) >= f_m. If s is polynomial in m then k=Omega(m/log m):
for k>m the claim is immediate, otherwise log(k!)<=k*log m.
Increasing the hub budget by a fixed constant does not repair this
particular intended positional inverse.
