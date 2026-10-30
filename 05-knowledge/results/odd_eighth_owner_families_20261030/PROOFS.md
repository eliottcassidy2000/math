# Uniform owner formulas, exact continuation carriers, and signed cycles

**PROVED scoped elementary deductions / FINITE-EXACT explicitly bounded
classification / CITED and credited inherited mechanisms.** Same-agent
proof review; no independent referee or literature-priority assertion.

For a nonempty finite positive integer speed set V, write
`M(V)=max_t min_(v in V) ||v*t||`, on R/Z. Safe conditions are closed.
An isolated safe point is a witness, even when its measure is zero.

## O1. Arbitrary second-band attachments have uniform top windows

Let q>=3, let E be ANY subset of {q+1,...,2q-1}, and put

    W(q,E)={1,...,q-1} union E.
    c_u = u+q if u+q belongs to E, and u otherwise, for gcd(u,q)=1.
    s = min_(u unit mod q) c_u.

The core has maximum 1/q, with precisely the primitive times a/q.
Indeed the AP subcore bounds the maximum by 1/q and forces its equality
set to be the primitive q-grid. No added speed is divisible by q, so
each such time remains safe.

If q+1 belongs to E, then for EVERY radius

    1/(q+2) < r <= 1/q

the complete safe set is

    union_(a unit mod q)
      [a/q - (1/q-r)/c_p, a/q + (1/q-r)/c_(q-p)],
      p=a^(-1) mod q in {1,...,q-1}.                    (O1)

If q+1 does not belong to E, the same formula holds for
`1/(q+1)<r<=1/q`. Endpoints and directional owners are part of the
statement, not merely interval lengths or an unordered list of weights.

**Inheritance.** For the AP core the Farey-tent envelope is classical
and is proved directly in Session114 [GENERAL_CORE](../session114_continuation_symmetry_20261030/lrc/GENERAL_CORE.md).
For E={q+1}, that file proves the complete window description down to
1/(q+2), including the wrapped-component exclusion at a=1. We use that
proved base, then show that arbitrary remaining second-band speeds cut
only their own active residue walls. This is the extension here.

**Proof of the new assertion.** First suppose q+1 is present. Start with
the inherited windows. At a/q, an added speed v=q+j has residue
`k=v*a mod q`, with 1<=k<=q-1. On a left side whose old owner p>=2,
the uncut endpoint is `a/q-epsilon/p`, epsilon=1/q-r. Write

    v=k*p+ell*q,     ell<=1,

since v<2q. If k>=2, the left endpoint stays above its lower safe wall.
It suffices, by affinity in r, to check r=1/(q+2); there the required
inequality is

    k/q - 2v/[q(q+2)p] >= 1/(q+2)
       iff p(k-1) >= 2ell.

This holds because p>=2, k>=2 and ell<=1. If k=1, the lower wall
instead replaces p by max(p,v)=c_p. The corresponding right-side
argument uses q-k and q-p and gives exactly c_(q-p). Lower walls are
monotone along these windows, as are upper walls, so checking the
outer endpoint suffices. No extra phase branch appears: when k=1,
the endpoint's phase displacement has magnitude at most

    2(q+p)/[q(q+2)p] <= (q+2)/[q(q+2)]=1/q,

and is strictly smaller inside the open radius chamber. Thus its
positive phase cannot cross zero. The k=q-1 case is reflected.

The old p=1 side occurs only at a=1 and is already replaced by q+1.
For v=q+j, j>=2, its left endpoint at r=1/(q+2) has phase

    j/q - 2v/[q(q+2)(q+1)] >= 1/(q+2).

After clearing positive denominators the left-minus-right numerator
is `q*((j-1)*q+3*j-3)`, positive. The upper wall is handled by the
other side's owner q-1, already at least two. Reflection handles a=q-1.
Consequently nonactive walls introduce no cuts; active ones introduce
exactly the stated tighter endpoints. The intersection is (O1).

If q+1 is absent, start with the AP windows above 1/(q+1), where only
the denominator-q Farey tents survive. For p>=2 the same calculation
at r=1/(q+1) requires `p(k-1)>=ell`, which holds for k>=2.
For p=1 (a=1), v=q+j has j>=2. At the left endpoint t=r,

    v*r-1 >= r, since (v-1)/(q+1)>=1.

The right-wall argument is reflected. Active residues again replace
their old owners by c_u. These observations prove the second claim.
Every assertion includes equality at the closed endpoints of a fixed
radius window, but excludes the lower radius boundary itself. QED.

The bound `s<=q+1` is automatic: if a unit class is not raised, its
owner is <=q-1; if all are raised the minimum is q+1.

## O2. Uniform single-tail spectrum and its finite head

Append an arbitrary positive integer speed v to W(q,E).

* If q does not divide v, the maximum stays 1/q, since every primitive
  q-grid time remains safe.
* If v=qm and `2m>s`, then

      M(W(q,E) union {qm}) = m/(qm+s).                 (O2)

* If `2m<=s`, then `M<=1/(q+2)`. This inequality does not by itself
  give equality: further second-band speeds can lower the value.

In particular (O2) applies uniformly whenever `m>(q+1)/2`.

**Proof.** On any complete window from O1, every side displacement
has size at most `(1/q-r)/s`. Since qm vanishes at its center,

    r <= ||qm*t|| <= qm*|t-a/q| <= qm*(1/q-r)/s.

Thus r<=m/(qm+s). On the longest side the outer endpoint attains
this bound, provided it lies above the radius cutoff. If q+1 is present,
that is exactly 2m>s. If q+1 is absent, s=1; for m>=2 the candidate
lies above 1/(q+1), and the same argument applies. For m=1 the set
contains {1,...,q}, bounding the maximum by 1/(q+1), while t=1/(q+1)
survives all extra speeds because q+1 is absent. This gives O2 there too.
Finally 2m<=s forces q+1 to be present, and O1 excludes every r>1/(q+2).
QED.

**Maximizers and mass.** If m>1 and 2m>s, let

    a=s^(-1) mod q in {1,...,q-1}, h=(a*s-1)/q.

There are exactly two maximizing times:

    t=(a*m+h)/(qm+s), and 1-t.                         (O3)

At each, the unique binding core owner is s and the opposing tail is
qm. Equality in the displacement bound forces these two sides; distinct
unit classes have distinct integer owners, so no other side has owner s.
Inactive inequalities are strict at interior chamber radii. Hence the
last-cell slope, defined by `mu Safe(M-epsilon)=epsilon*B` for sufficiently
small positive epsilon, is

    B=2*(1/s+1/(qm)).                                  (O4)

The special m=1 case can have more than two maximizing times and must
not inherit O3-O4 without its separate analysis below.

### The complete first-tail formula, for every q and E

Let d be the least element of {q+1,...,2q} absent from E (2q is always
absent). Then

    M(W(q,E) union {q}) = 1/d,
    argmax = {a/d : gcd(a,d)=1}.                        (O5)

The speed set contains the full AP {1,...,d-1}, which proves the upper
bound and restricts all maximizers. All further speeds are <2d and
different from d, so all primitive d-grid times survive. This proves
both assertions. Its local terminal slope is

    2*sum_(u unit mod d) 1/max{v in W(q,E) union {q}: v=u mod d},

by the same active-wall rule. This is an all-size first-missing-speed
law, not an inference from a small census.

## O3. Complete classification at the four odd-eighth centers

Fix q=8 and let E run through ALL 128 subsets of {9,...,15}. Formula
O2 holds for every m>=5. The four remaining values m=1,2,3,4 are exactly
the entries of [CLASSIFICATION](CLASSIFICATION.md), certified by two
different exact finite algorithms. Together they give a complete infinite
family classification, with exactly 14 distinct functions of m.

The proof of completeness is hybrid and explicit: a universal proof
covers every m>=5, and a finite exact certificate covers 128*4 remaining
pairs. Equality of two listed profiles means equality for every m;
different profiles differ in one of the first four values or in s,
which already distinguishes them at m=5. There is no extrapolation from
the finite overlap bank.

For example E={9} recovers the q=8 instance of Session114, with
`M(m=1)=1/10` and `M=m/(8m+3)` for m>=2. In contrast E={9,10} has
`M(m=1)=1/11` and the SAME later formula. The extra speed 10 is invisible
in the upper chamber but decisive at the first tail.

## O4. Exact minimal carriers for future attachments

The state is E subset {q+1,...,2q-1}. Allowed future operations adjoin
any available speeds from the same band (set union). Queries ask for
maximum clearance after appending a named resonant tail qm.

**Upper-range consumer.** Restrict m>(q+1)/2. The smallest deterministic
carrier has exactly `2^phi(q)` states: one bit for each unit residue u,
recording whether q+u has been adjoined. Nonunit extras do not enter O2.
This bit vector updates by OR and is sufficient. For necessity, if two
states differ in the bit u, adjoin every OTHER unit extra. The first
state now has s=u, and the other has s=q+1 (or conversely). The same
large m query separates them. Every unit bit is therefore necessary.

The current scalar s alone has phi(q)+1 possible values; it is not closed
under future additions. For q=8 these counts are five current values and
16 future states.

**All-tail consumer.** Allow every m>=1. The smallest deterministic
carrier has exactly `2^(q-1)` states: the whole E. If two states differ
at a speed w, adjoin every other speed in the available band and query
m=1. By O5, the answers are 1/w and 1/(2q), distinguishing them.
For q=8 all 128 states are therefore indispensable, although only
14 complete current clearance profiles occur. This is a separation
between present profiles and future behavior, not an assertion about
arbitrary LRC sets or tournament carriers.

## O5. A weighted recurrence for every second-band attachment ensemble

List the units as u_1<...<u_h, h=phi(q). Independently include speed
q+j with probability p_j, allowing unequal probabilities and the endpoints
zero and one. Fix m>(q+1)/2 and put f(x)=m/(qm+x). Then the exact expected
clearance is Z_1, where

    Z_(h+1)=f(q+1),
    Z_j=(1-p_(u_j))*f(u_j)+p_(u_j)*Z_(j+1).            (O6)

This is the first-absent-unit recurrence: if u_j is not raised, it is
the least owner; otherwise inspect the next unit. Nonunit probabilities
cancel from this particular expectation. With fair independent bits,

    E[M]=sum_(j=1)^h 2^(-j)*f(u_j)+2^(-h)*f(q+1).

Among all 2^(q-1) sources, the owner u_j occurs 2^(q-1-j) times, and
q+1 occurs 2^(q-1-h) times. At q=8 the five multiplicities are
64,32,16,8,8. Thus a one-pass scalar recurrence computes the weighted
mean despite the larger minimal carrier required for exact future
queries. This is an exact finite-ensemble identity, not an assertion
that a mean controls every exceptional arithmetic source. Direct native
maxima verify fair, graded, empty and full laws for every q=3,...,9.

## P1. A packet's first safe arrival is the correct one-shot observable

**Concurrent-source credit.** This argument was developed here before the
incoming fetch, but the same general statement is already proved in
[L1, uniform tail banks](../odd_eighths_uniform_families_20261030/LRC.md),
commit f128982d82. It is REPRODUCED here, not claimed as a new contribution.
The explicit failure of closure under later attachments below is the
additional test supplied by this packet.

Suppose a core has complete top windows with centers in (1/q)Z,
maximum rho, and least directional owner s as in O1. For a finite
packet A of positive integers, append all speeds q*a, a in A. Define

    tau_A(r)=min{z>=0: ||a*z||>=r for every a in A},

with infinity when empty. One period [0,1] suffices. At any radius in
the core's certified top chamber,

    Safe(core union q*A,r) nonempty
       iff tau_A(r) <= q*(rho-r)/s.                    (P1)

**Proof.** Write t=center+x. All tail phases become a*(q*x). Every
window direction has |x| bounded by (rho-r)/s; one direction realizes
the entire interval of displacements from zero to that bound. Reflection
removes the sign of x. Thus a point exists exactly when the packet has
a safe phase in that initial interval. QED.

This gives uniform finite piecewise-rational recipes once the packet's
safe intervals are supplied. It does not replace those intervals under
further packet attachments. Here is an exact hostile: A={4} and B={4,28}
have identical tau(r)=r/4 for every 0<=r<=1/8. Nevertheless, with AP
core {1,...,7}, adjoining the SAME additional packet speed 32 gives

    M({1,...,7,32,256}) = 31/257,
    M({1,...,7,32,224,256}) = 3/25.

Before that addition both maxima are 4/33. The first-arrival function
is sufficient for a single composition against this core, but is not a
closed joint-attachment carrier. The complete safe-set intersection
retains the missing later phases. These four numeric instances have
independent exact maximum and interval certificates; equality of the
initial arrival functions follows from 7r<=1-r for r<=1/8.

## P2. Harmonic completion extends to every raised-first-owner core

**PROVED extension of an incoming theorem.**
[Session115, HARMONIC_COMPLETION](../session115_odd_eighth_classifications_20261030/lrc/HARMONIC_COMPLETION.md)
at 893f84da21 proves the following laws for E={q+1}. Its counting mechanism,
critical-owner correction, and affine-period recurrence are inherited.
O1 extends their hypotheses to every E containing q+1, with the complete
list C={c_u:u unit mod q} of actual owners substituted for the old list.
Assume q>=3 and q+1 in E throughout this section; put n=q+1.

For the full harmonic block

    V=W(q,E) union {qm,2qm,...,q*qm},

its maximum equals 1/n if and only if m>=s. In that preserving family,
all maximizing points are

    a_e/q - j/[q*n*m], and their reflections,
    e in C, a_e=e^(-1) mod q,
    1<=j<=floor(m/e), gcd(j,n)=1.                       (P2)

All addresses are distinct. Thus, writing mu for the arithmetic Moebius
function and ell_j=j^(-1) mod n in {1,...,n-1},

    N(m)=2 sum_(e in C) sum_(d | n) mu(d)*floor(m/(ed)),
    w_n(j)=1/ell_j+1/(n-ell_j),
    B(m)=2/(qm) sum_(e in C) sum_(j<=floor(m/e), (j,n)=1) w_n(j). (P3)

Here B is the actual physical terminal mass slope, not an unweighted
count of four abstract phases. The integer j is a displacement address;
its residue alone cannot locate the original time.

**Proof of the extension.** The q-tail AP bounds M by 1/n and forces
its base phase to be a primitive n-grid point. O1 is valid at radius
1/n because q+1 is present. On an e-owner side its displacement budget
is 1/(qne). A primitive tail phase has displacement j/(qnm); hence the
necessary and sufficient condition is ej<=m with gcd(j,n)=1. The least
possible e and j give precisely m>=s. This also proves (P2) and the
count in (P3), by disjointness of the top windows and Moebius inversion.
At such a point the fastest opposite tail owners are qm*ell_j and
qm*(n-ell_j). A tied core owner e is smaller than qm: m>=s>=2 and
e<2q<=qm. Thus the tail owners control the derivative even at a collision
ej=m. Their reciprocal sum gives (P3). This explicitly checks the only
new issue introduced by arbitrary raised owners. QED.

Put L=n*lcm(C), H_C=sum_(e in C)1/e and
H_U=sum_(a unit mod n, 1<=a<n)1/a. For every m>=s,

    N(m+L)-N(m)=2L*phi(n)*H_C/n,
    (m+L)B(m+L)-m B(m)=4L*H_C*H_U/(qn).                (P4)

Each floor bound increases by a whole number of unit periods, so these
are exact affine-period recurrences. In particular N(m)/m tends to
2phi(n)H_C/n and B(m) to 4H_C*H_U/(qn). This includes all second-band
choices satisfying the declared guard, rather than only the original
single extra speed. No assertion covers E without q+1 at this boundary:
there the available O1 chamber stops exactly at 1/n.

The native verifier checks every such E for q=3,...,8, at the failed
threshold m=s-1, equality m=s, and two larger scales. It intersects
literal safe intervals, compares every maximizing point and mass slope,
and separately checks both recurrences over an entire common period.

## T1. The Cayley transform of a local order is a literal signed cycle

For a tournament sign matrix S, put U=(I+S)(I-S)^(-1). Invertibility
follows from real skew-symmetry. For the transitive matrix S_n, with
ones above the diagonal, define

    W e_j=-e_(j+1) for 0<=j<n-1, and W e_(n-1)=e_0.

Direct column comparison gives `(I-S_n)W=I+S_n`, so U=W, exactly.
There is no unspecified conjugating basis. Consequently

    W^n=(-1)^(n-1) I,
    characteristic polynomial = z^n-(-1)^(n-1).        (T1)

For even n the phases are (2j+1)/(2n); for odd n they are j/n.
At n=4 they are precisely 1/8,3/8,5/8,7/8. At n=2^(a-1), a>=2,
they coincide with the primitive 2^a-grid maximizers of the AP with
speeds 1,...,2^a-1. This is a literal identification of marked phase
sets, not of the tournament and its speed set or their future operations.

**Complete classification.** U is a signed permutation matrix iff T
is switching-equivalent to a transitive tournament, the classical local
orders classified in [THM-468, determinant floor](../../../01-canon/theorems/THM-468-tournament-determinant-floor.md).
The old local-order classification and count are inherited. The present
operator criterion follows as follows.

Switching and relabelling conjugate U by the corresponding signed
permutation, so the forward family has this property. Conversely if U
is signed monomial, a decomposition into two or more permutation cycles
would make S=(U-I)(U+I)^(-1) block diagonal, contradicting every off-diagonal
entry of a tournament being nonzero. Therefore U is one signed cycle.
If its cycle sign product is delta, invertibility of U+I forces
`delta!=(-1)^n`, hence delta=(-1)^(n-1). Every such signed cycle is
signed-permutation-conjugate to W, and the inverse Cayley map recovers
the corresponding switching of S_n. This proves necessity.

There are `(n-1)!` directed cycles and `2^(n-1)` allowed sign patterns.
Thus the inherited labelled count L_n and its recurrence have the
literal operator explanation

    L_n=2^(n-1)(n-1)!,   L_(n+1)=2n L_n.              (T2)

The class is closed under relabelling, vertex switching and converse
(U goes to U^(-1)). Arbitrary vertex attachments need not preserve it.

**The recurrence is a complete native attachment rule.** Fix a labelled
local order on n vertices, and its actual signed cycle U. Choose any of its
n directed cycle edges j->i with sign c. Replace that edge by
j->new->i with signs a,b satisfying ab=-c. There are precisely two choices.
Then invert the Cayley map. The old tournament is retained literally,
and these 2n choices give exactly all local-order one-vertex attachments.

To see the preservation rather than infer it from the count, for a signed
n-cycle with total sign (-1)^(n-1) the finite inverse formula is

    S=sum_(k=1)^(n-1) (-1)^(k-1) U^k.                  (T2a)

For distinct old vertices its entry is the product of signs along their
unique forward cycle path, times (-1)^(path length-1). Subdividing the chosen
edge increases the path length by one exactly when it changes its sign
product by -1, so every old entry is unchanged. The whole cycle sign also
flips, as required when n grows. Conversely, deleting the new vertex from
any child cycle and replacing its two edges by sign -ab recovers the parent
cycle by the same identity. The Cayley map is injective, so there are no
duplicates. This proves the attachment rule and its completeness, including
n=1. The count itself remains the classical count credited above.

The verifier constructs every such insertion for every labelled local
order through n=5 and compares its entire child set with all 2^n native
orientations. This is a literal orientation test, not a comparison of counts.

## T2. A quartic energy certificate for the same family

Let D(T) count four-vertex vortex subsets (a directed triangle with an
external source or sink). For all tournaments,

    tr(S^4)=n^2(n-1)^2/2 - 4*binom(n,4) - 32*D(T).    (T3)

Each four-set has squared Pfaffian 1 or 9, with 9 exactly at vortices,
as already recorded and audited in THM-468. Therefore the fourth
characteristic coefficient is binom(n,4)+8D. Newton's identity with
`tr(S^2)=-n(n-1)` and second coefficient binom(n,2) gives T3. This
derivation does not claim discovery of the classical four-set formula.

The maximal possible quartic energy is attained iff D=0, hence iff
the equivalent local-order and signed-cycle conditions hold. Equivalently,
the transitive skew spectrum alone recognizes this particular class.
This exceptional sufficiency is scoped: outside this family, spectra
need not determine tournaments, Hamiltonian counts or automorphism groups.
The transitive triangle and directed triangle already share the class
and spectrum while having one and three Hamiltonian paths. The actual
81/27 control is also replayed; nothing replaces its native groups by
a cyclic-triple or spectral quotient.
