# Apex-ascent moment recurrences and deletion-marked counts

**Session250, 10 November 2026. PROVED-IN-NOTE; finite controls separately
labelled. No canon ID and no mathematical-priority claim. Chan–Pak OPEN.**

## Inheritance and purpose

Session249's `05-knowledge/results/sequence_gf_20261110/grammar/PROOF.md`
provides the positive clipping principle and a one-apex grammar which misses
13. Its synthesis explicitly asks for a richer two-apex boundary grammar.
The inherited Session191 two-apex construction already has 36 nonnegative
states and reaches 230519 in its finite census; S79 explicitly records this
in `05-knowledge/results/spine_grammars_cutwidth_apexes_hankel_seams_20261104.md`.
Sections 1–2 recover this construction self-containedly and extend its notation
to arbitrary Q. They do not claim a new two-apex transfer or a new realization
of 13. This construction preserves the actual Hamiltonian-path count, not
just a quotient statistic. The new contributions are the closed
ascent-polynomial first-moment formula in Section 7, the exact second-moment
recurrence and fluctuation limit in Section 8, and the explicit marked-census
interface. The
support construction is an application, not a new clipping theorem.

Concept board: ordered runs / positive boundary states / marked isomorphism
classes / actual-H support / the original sixteen-point carrier. Anchor:
actual-H word languages. Niche: a marked A000568-related census. Wildcard:
tensor-power moment recurrences. The original sixteen-point decomposition
does not enter the proof. Sixteen words of length two would only be a
cardinality coincidence until an operation-preserving map is supplied.

## 1. Carrier and exact run decomposition

Fix a tournament Q on named vertices 1,...,k, where k>=1. A word
w=(a_1,...,a_n) on the alphabet {0,1}^k defines T(Q,w): the new vertices
v_1,...,v_n form the transitive spine v_i -> v_j for i<j; bit a_j(i)=1
means apex i -> v_j. All remaining edges are those of Q.

Every directed Hamiltonian path determines a unique permutation pi of the
k apices and k+1, possibly empty, runs R_0,...,R_k:

    R_0, pi(1), R_1, pi(2), ..., pi(k), R_k.

Within each run the spine vertices must be in increasing order. Conversely,
an assignment of each spine vertex to one run gives exactly this candidate
path. It is valid precisely when:

* in each nonempty run its first vertex receives an edge from the left
  apex, if one exists;
* its last vertex sends an edge to the right apex, if one exists;
* an empty internal run has the required edge from its left apex to its
  right apex in Q.

These conditions exhaust all consecutive pairs. This is a bijection, not
an overcount by path presentations. Empty endpoint runs need no Q test.

## 2. Positive boundary-state theorem

For each fixed pi, scan w from left to right. For each of the first k
runs use states E,G,B: empty, or current last vertex good/bad for the
right apex. For the final run use E,N: empty/nonempty. This gives

    d_k = 2 * 3^k

states per pi and D_k=k! d_k states in total. All blocks start at their
all-empty state with weight one.

On a letter a choose the run receiving this vertex. If that run is empty
and has a left apex, require a(left)=1; otherwise disallow this choice.
Set its new state to G if it has a right apex r with a(r)=0, to B if
a(r)=1, and to N if there is no right apex. Other coordinates stay put.
Multiple choices leading to the same state add their weights. Every
transition is a nonnegative integer, and there is one matrix M_a for
each alphabet letter. Matrices have independent pi blocks.

The acceptance row f has value one exactly on states with no B, and
whose empty internal runs satisfy the corresponding Q edge. Endpoint
emptiness is accepted. If c is the initial column, then

    H(T(Q,w)) = f M_(a_n) ... M_(a_1) c.                  (1)

**Proof.** Inductively, a state's weight counts assignments of the scanned
vertices to runs with the stated last-vertex information, with every
already determined first-vertex constraint enforced. Replacing a last
vertex deliberately forgets its former right-edge condition: only the
eventual last vertex must meet that condition. At the end f enforces the
remaining conditions from Section 1. Summing the independent apex orders
counts every path once. The proof also handles n=0, yielding H(Q). QED.

The distinction E/G/B cannot be discarded casually: E triggers a first
vertex constraint, whereas B can be repaired by a later vertex in the
same run. Rejecting B early loses valid paths.

## 3. Actual-value functional equations and finite-target automata

Let C(v)=f v. Retain all D_k nonnegative coordinates. For each reachable
vector v, its word successors are M_a v. Thus the formal catalytic series

    F(t;y) = sum_w t^|w| product_i y_i^(v(w)_i)

satisfies

    F(t;y) = y^c + t sum_a F(t; phi_a(y)),
    phi_a(y)_j = product_i y_i^((M_a)_(ij)).             (2)

This follows by appending one letter, using y^(M_a v)=phi_a(y)^v.
Specialize y_i=z^(f_i) to obtain the actual-value enumerator

    sum_w t^|w| z^(H(T(Q,w))).                         (3)

Equation (2) is finite and positive; it is not a claim that (3) is
rational in both variables.

For a fixed integer target h>=0, replace every coordinate by min(h+1,v_i).
Nonnegative integer matrices satisfy

    clip(M v) = clip(M clip(v)).

Indeed a summand using an already saturated coordinate is either absent
or at least h+1; otherwise all summands agree exactly. The same argument
applies to the nonnegative output row. There are at most (h+2)^D_k states.
Consequently the words with H=h form an effectively constructible regular
language, and their word-length generating function is rational. A finite
breadth-first search either finds a witness or proves absence in this
fixed-Q family. This deliberately loose upper bound is a termination
theorem, not a practical complexity claim.

The unbounded support question, whether a target-uniform family realizes
every required value, is not discharged by separately terminating searches.

## 4. An all-orders moment recurrence engine

For r>=1 put

    B_r = sum_a M_a^(tensor r),
    m_(n,r) = sum_(|w|=n) H(T(Q,w))^r.

Tensoring (1) and summing independent letters gives

    m_(n,r) = f^(tensor r) B_r^n c^(tensor r),
    sum_(n>=0) m_(n,r)t^n
       = f^(tensor r) (I-t B_r)^(-1) c^(tensor r).     (4)

Thus every fixed raw moment over words has an explicit constant-coefficient
recurrence of order at most D_k^r, obtained from the characteristic
polynomial of B_r. The symmetric tensor subspace contains c^(tensor r)
and is B_r-invariant, so the bound improves to binom(D_k+r-1,r).
All of this follows from the specified matrices, without recurrence fits.
The zero-th moment is 2^(kn), including H=0 words by convention.

The generating index is spine length, and the universe is the fixed-apex
family. Neither superfactorial growth of unrestricted A000568 nor factorial
growth in a cyclic substitution family conflicts with (4): here
H<=k! (k+1)^n follows already from the run assignments. Growing k destroys
the fixed-dimensional bound.

## 5. A precise A000568-related quotient, especially k=2

Take k=2 and Q the edge a->b. Four letters encode the two apex incidences.
The resulting words are in bijection with isomorphism classes of pairs
(T,S), where S is an unordered two-element vertex subset and T-S is
transitive. The source and sink of T[S] canonically name a,b; the
transitive complement has a unique linear order. These structures recover
the word, proving both injectivity and surjectivity.

For an unmarked tournament isomorphism class [T], the forgetful fibre is

    d(T) = number of Aut(T)-orbits of pairs S with T-S transitive.

Thus, with N=n+2,

    sum_(|w|=n) z^H(T(w)) = sum_[T], |T|=N d(T) z^H(T). (5)

This is an orbit-resolved, deletion-marked A000568-related sequence. The
constant z=1 yields 4^n, not A000568(N). Fibre size is not the labelled
number of admissible deletion pairs and is not constant across T.
For example the transitive N-tournament has trivial automorphism group
and d(T)=binom(N,2), while a directed triangle has one pair-orbit.

The map preserves H exactly and loses the marked deletion pair. Therefore
positive coefficient support transports from (3) to actual tournament
support, but multiplicities cannot be read as the unmarked census. The
missing sidecar for unmarked enumeration is precisely d(T), and classes
with no transitive complement of size N-2 have zero weight.

## 6. Scope, reproduction, and next obligation

Run `python3 apex_words.py` in this directory. `receipts.json` records all
words through spine length five, independent subset-DP path counts through
length four, coordinatewise clipping controls, moment sums, and the first
word realizing 13. This is a proof-derived finite control, not an inherited
theorem audit or an exhaustive unmarked tournament census.

The smallest strong next step is to minimize the *reachable positive*
two-apex representation while retaining a certified H output, then find
word substitutions whose action closes on a small semigroup of boundary
vectors. This would turn fixed-target recognition into value-transforming
construction rules. A rational moment series alone cannot do that.

No identification with the original 1-5-4-5-1 decomposition is asserted
here. Such a bridge needs the original carrier and an explicit map that
commutes with one of the M_a, rather than an equality of dimensions.

## 7. Closed first moments from apex ascent polynomials

Here is a stronger reduction than the tensor dimension bound. For an apex
order pi let asc_Q(pi) count indices j with pi(j)->pi(j+1), and let
b=k-1-asc_Q(pi). Put

    P_Q(E) = (1+E)^2 / 4^k
             * sum_pi (E+3)^asc_Q(pi) (E-1)^b.

Then the first-moment exponential generating function is

    sum_(n>=0) m_(n,1) x^n/n! = P_Q(exp(2^k x)).       (6)

**Proof.** Fix an assignment of spine indices to runs. A nonempty endpoint
run imposes one independent apex-incidence bit. A nonempty internal run
imposes two, even when it is a singleton, because its left and right
apices are distinct. Distinct runs involve disjoint spine vertices, so
these constraints never overlap. All unspecified incidence bits are
free. Thus each spine vertex supplies weight 2^k; each nonempty endpoint
run adds factor 1/2, and each nonempty internal run factor 1/4. In the
exponential generating function an endpoint run therefore contributes
1+(E-1)/2=(1+E)/2. An internal run contributes (E+3)/4 when its empty
case is allowed by Q, and (E-1)/4 otherwise. Multiplication assigns the
labelled spine indices to runs with the exact multinomial coefficient;
their internal increasing order is forced. Sum over pi. QED.

The polynomial P_Q has degree k+1. Write it as sum_j p_j E^j. Then

    m_(n,1) = sum_(j=0)^(k+1) p_j (j*2^k)^n,

with 0^0=1. Hence an explicit annihilator is
product_(j=0)^(k+1) (L-j*2^k), where L is the forward shift.
This derives the recurrence before generating any coefficient. Negative
p_j can occur; positivity belongs to the combinatorial sequence and the
boundary matrices, not necessarily to this spectral basis.

For k=2 there are two orders, with ascents 1 and 0, giving

    P_Q(E)=(1+E)^3/8,
    m_n=(0^n+3*4^n+3*8^n+12^n)/8,
    sum_(n>=0) m_n t^n
      = (1+3/(1-4t)+3/(1-8t)+1/(1-12t))/8.

Thus m_0=1,m_1=6,m_2=48,m_3=432 and

    m_(n+4)=24m_(n+3)-176m_(n+2)+384m_(n+1), n>=0.

Equivalently the order-three recurrence holds starting at n=1, not n=0.
The zero-root initial term is essential. This first-moment theorem is the
new sequence compression: an arbitrary fixed apex tournament contributes
only its adjacent-ascent polynomial, with prescribed exponential bases.
It does not assert that the same ascent polynomial determines fixed-word
H values or higher moments.

## 8. An exact second-moment recurrence and fluctuation limit

For the two-apex four-letter word family let s_n=sum_(|w|=n) H(w)^2.
The proof-derived tensor automaton has 1144 reachable ordered-pair states.
The companion `moment_recurrence.py` certifies the following all-length
recurrence, with L the forward shift:

    L^3 (L-36)(L-24)^2(L-16)^2(L-12)^2(L-8)^4(L-4)^3 s = 0.  (7)

**Status: finite exact linear-algebra certificate of an all-length
identity.** We do not call a finite-prefix recurrence guess a proof.
The script first obtains a candidate via exact rational elimination of
scalar terms. It then separately checks 1144 consecutive residuals of
(7) on the full reachable tensor-state evolution. Each residual is
f2 B^n p(B)c2. Since B has dimension 1144, its characteristic polynomial
propagates these 1144 zero residuals to all n>=0. This is the finite
certificate which proves (7). The factorization is verified by integer
polynomial multiplication, without any symbolic-algebra dependency.
No claim of minimality is necessary or made.

The initial terms, full integer coefficients, all 1144 residual checks,
and OGF numerator are in `moment_certificate.json`. The denominator is

    Q(t)=(1-36t)(1-24t)^2(1-16t)^2(1-12t)^2(1-8t)^4(1-4t)^3.

The numerator has degree 16. Keeping this polynomial part matters: the
three zero roots encode an initial transient, and deleting them without
changing the recurrence's start index would be wrong.

Computing the simple-pole residue exactly gives

    s_n ~ (841/25600) * 36^n.

Every other pole has inverse modulus at most 24, with polynomial factors
from repeated poles, so its contribution is asymptotically smaller.
For a uniform word W_n of length n, there are 4^n equiprobable words and
Section 7 gives E H(W_n) ~ 3^n/8. Therefore

    E H(W_n)^2 ~ (841/25600) * 9^n,
    Var(H(W_n)) / (E H(W_n))^2 -> 441/400.              (8)

This is an exact nonconcentration constant in the deletion-pair-marked
universe, not the uniform unlabelled A000568 universe. It shows why even
complete first-moment data cannot be read as a typical H value. Equation
(7), unlike the inherited general transfer existence, is an explicit new
session output with a self-contained exact certificate.

## 9. Constructive receipt and inherited boundary

The recovered word (1,2,0) has H=13 on five vertices: apex a->b and spine
v1->v2->v3, with a->v1, v1->b, v2->a, b->v2, v3->a, v3->b.
This is an elementary constructive control that the richer grammar exits
the one-apex obstruction, not a new existence theorem for 13. All words
through length five are recorded; lexicographically first witnesses for
attained odd values below 100 and cumulative missing-value lists are in
`receipts.json`. Previous Session191 computations reach much farther;
this small census only exercises the new recurrence interfaces.

Reproduce both receipts using Python's standard library:

    python3 apex_words.py
    python3 moment_recurrence.py

The first command regenerates the word/support and first-moment controls;
the second regenerates the all-length second-moment certificate. Neither
changes inherited results or requires a large unlabelled census.

At spine length at most five, the only odd values below 100 not observed
are 7, 21, 61, 63, 87, 99; the largest attained value in that finite word
universe is 151. These are bounded-length gaps, not family exclusions.
