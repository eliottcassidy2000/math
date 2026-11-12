# The target budget and the apex budget control different kinds of growth

**PROVED from SCC decomposition and the cited Busch bound.** This is a
fixed-apex presentation theorem, not an enumeration of isomorphism classes.
Its ingredients and the general finite-core/rational-series mechanism are
inherited; the explicit two-budget hierarchy extends this session's two-apex
normal form. Chan–Pak remains OPEN.

Fix a labelled tournament A on k≥1 named apex vertices. Add m ordered spine
vertices forming a transitive tournament and choose all k*m cross arcs freely.
A word in the2^k-letter incidence alphabet specifies an actual tournament.
Let F_N^A(z) count those words whose Hamiltonian count is N, by spine length.
Write B(N)=floor(log_5(N^3)), computed by integer power comparison.

## 1. Exact apex-partition grammar

An ordered partition P=(B1,...,Bs) of the named apex set is compatible if every
arc of A between different blocks points from the earlier to the later block.
For a nonempty apex subset B, define K_(A[B],d)(z) to count STRONG core words
with those apexes and Hamiltonian count d. Include the singleton apex with
empty spine and d=1. No other H=1 strong core is included.

Then the exact all-order generating function is

    F_N^A(z)=sum_(compatible P) (1-z)^(-s-1)
               sum_(d1...ds=N) product_i K_(A[Bi],di)(z).       (1)

The factors di are positive integers and are ordered with the blocks.
The grammar is unambiguous, not an inclusion-exclusion count.

Proof. Every nontrivial strong component contains an apex: the spine alone
is transitive. A component's spine names form an interval, because if it
contains i<j, every intermediate spine vertex lies on i→l→j and hence is
mutually reachable with it. The apex subsets of components, in their unique
condensation order, give P. All apex arcs across blocks must point forward.
Outside these components there are exactly s+1 possible gaps containing
singleton spine vertices. Each gap has one forced incidence letter: it loses
to all preceding apex blocks and beats all later ones. Its length is an
arbitrary nonnegative integer. The block words are precisely the specified
strong cores, and path counts multiply across components. Conversely each
such list of cores and gap lengths yields that SCC decomposition, by the
same arc rules. This proves both the formula and uniqueness. QED.

A block with a single apex and no spine is essential: suppressing these
blocks would lose gap locations and falsify both the denominator and the
count-one boundary. Formula(1) is an explicit physical composition, not a
claim that arbitrary polynomials can be substituted for K.

## 2. Universal numerator degree, even as the apex count grows

The cited Busch floor implies5^(n-1)≤H(T)^3 for every nontrivial strong
n-vertex tournament. It holds also for the singleton n=1,H=1.
For a product of the s core components in(1), with total nonneutral spine
length j, multiplying these inequalities gives

    5^(j+k-s)≤N^3, hence j≤B(N)+s-k.                       (2)

Thus every core polynomial in(1) is finite, and its product degree has the
stated bound. In particular if s<k-B(N), this partition cannot contribute.

On putting(1) over the common denominator(1-z)^(k+1), its term for s blocks
is multiplied by(1-z)^(k-s). Consequently

    F_N^A(z)=P_N^A(z)/(1-z)^(k+1),   deg(P_N^A)≤B(N).       (3)

The numerator may have negative coefficients. It is not itself a physical
count bank. The degree bound is INDEPENDENT of k; the number of apexes appears
only in the maximum denominator power. Zero numerator means an empty fibre.

Equivalently, for every m>B(N), the word multiplicities obey the exact
(k+1)-st finite-difference recurrence

    sum_(j=0)^(k+1) (-1)^j binom(k+1,j) a_(m-j)=0,

with negative-index terms zero. For m≥B(N), a_m is a polynomial in m of
degree at most k. This polynomial can have smaller degree or vanish.
The target N bounds the transient numerator complexity; the apex budget k
bounds the degree of neutral repetition. These are different resources.

For N=1, only singleton blocks contribute. If A is transitive there is one
compatible singleton order and F_1^A=(1-z)^(-k-1), so a_m=binom(m+k,k).
If A is not transitive F_1^A=0. Thus(3) handles B(1)=0 exactly.
For k=0 separately, the graph is just the spine: F_1=1/(1-z), all other
fibres empty. No strong-core inequality is applied to the empty graph.

## 3. What this does and does not give a constructive proof

All words in one fixed target fibre have only B(N) nonneutral spine names,
regardless of how much neutral padding is used. Infinite word families do
not manufacture fresh target cores. Raising k can add genuinely new cores,
while increasing m after the finite-core bound only changes the reservoir
multiplicity. This distinguishes two operations often both called 'growth'.

The exact fixed-k spectrum condition is positivity of at least one core
product in(1), NOT positivity of an arbitrary numerator coefficient. A
Chan–Pak proof through a fixed apex budget would need a target-uniform
constructor for such a compatible core list. Equations(1)-(3) classify all
later repetitions but do not supply that constructor. Nor is a bounded apex
budget proved sufficient for all tournament counts.

This also explains a limitation of finite target automata. Their recurrent
states can be completely understood as neutral reservoirs while their
nonempty transient cores remain the hard arithmetic existence problem.
One should seek a recursion BETWEEN target-indexed core numerators, not try
to generate new values by pumping an already fixed target language.

## 4. Provenance and reproduction

The SCC product law, Busch's lower bound and finite-core rational enumeration
are inherited mathematical mechanisms (see the references in the wildcard
lane and Session266 multiplicity work). This note makes their fixed-apex
partition carrier and independent target/apex budgets explicit; it claims
neither a new general rational-language theorem nor a classification of all
unlabelled tournaments.

`check.py` independently constructs literal graphs, computes SCCs and path
counts, enumerates every required strong core for targets1,3,5,7,9 and up to
three named apexes, and checks the partition formula against all words
through spine length3. The general statements rest on the proof, not on this
finite check. There is no extrapolated numerical coverage claim.
