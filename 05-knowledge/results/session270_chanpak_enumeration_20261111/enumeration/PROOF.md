# Fixed-mass projective response quotient

**Status: PROVED identities and exact algorithm below; timings are finite observations.**
No polynomial-time or instantaneous unrestricted A000568 claim is made. No
literature priority claim is made. The projective quotient and its maximality are inherited in full from
Session98 `recurrence/PROJECTIVE_RESPONSE.md`. The contribution here is a
production-oriented synthesis: an odd-factorial integer realization (including
exact negative shifts), packed half-defects, a late gauge switch, and an exact
unit-tail contraction. None of those inherited GF identities is claimed new.

## Inheritance and concept board

Closest mechanisms: Session98
`05-knowledge/results/session98_tournament_growth_20261029/recurrence/PROJECTIVE_RESPONSE.md`
(including the full maximality proof repeated here for a self-contained implementation); THM-514 (`the-a000568-unit-tail-core-kernel-and-prime-tail-ladder`),
the future-boundary quotient, and the degree-local odd-factorial integral
transfer in `05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py`.
THM-5097 (`factorial-gauge-integral-cycle-response-recurrence`) rediscovers an
integer basis but does not quotient responses by fixed remaining mass.
Session249 connected diagrams give fixed Taylor jets, not fast full counts.
The hostile case is a pure unit tail: its mass is fixed, so its response is a
scalar weight, not a state coordinate. The near miss is using this quotient
unchanged for a whole prefix: the scalar depends on the requested degree.
Underused sidecar: exact target mass rather than a truncation bound.

Board: Burnside odd-cycle gas; projective weighted-automaton equivalence;
homogeneous coefficient extraction; odd-factorial gauge; analytic unit tail.
This lane counts ordinary isomorphism classes, not Hamiltonian-path support.
It has no implication for Chan–Pak.

## 1. The homogeneous generating-function identity

Fix permitted future odd lengths J containing 1. A processed core has mass d,
cycle multiset mu, and response q_j=sum_(a in mu) gcd(a,j). The future series is

 F_J(x;q)=sum_(m_j>=0) 2^[sum_j m_j(j-1)/2
                +sum_j j binom(m_j,2)+sum_(i<j) gcd(i,j)m_i m_j
                +sum_j q_j m_j] product_j x^(j m_j)/(j^m_j m_j!).

This is the existing multivariate Burnside cycle gas with external field q.
Its homogeneous coefficient obeys

 [x^R] F_J(x;q+c*j)=2^(cR) [x^R] F_J(x;q).

Indeed every contributing monomial has sum j m_j=R. Thus set

 b_j=(j q_1-q_j)/2,  b_1=0.

All b_j are nonnegative integers because all lengths are odd and gcd(a,j)<=j.
Absorb 2^(R q_1) into the processed weight and retain only b_j for j>=3.
This inherited quotient is strictly coarser than raw response equality.
The source is the full response vector; the map retains fixed-mass continuation
weights up to the stored scalar; it destroys cross-degree information; the
necessary sidecar is the target n and residual mass R=n-d.

**Maximal linear/projective equivalence.** Suppose two response fields differ
by delta, and their weights on every permitted future cycle multiset of mass R
have the same ratio. Comparing 1^R with j+1^(R-j), for every j in J with j<=R,
gives delta_j=j delta_1. Conversely that condition gives the constant ratio
2^(R delta_1). Hence this quotient is exact and maximal among quotients
identifying fields by a constant weight ratio on *every* individual future
cycle pattern. This does not assert minimality of the scalar summed answer,
where accidental cancellations/equalities could occur.

## 2. Integral recurrence with a removable unit tail

Let O_d=d!/2^v2(d!). For each state (d,b), store

 W(d,b)=O_d sum_mu [2^(e(mu)+(n-d)q_1(mu))/z(mu)],

where mu uses already processed lengths >=3, has mass d, and maps to b.
All terms are integers. The odd part of z(mu) divides O_d; its power of two
is sum_j v2(m_j!), which is at most e(mu) since each length is at least 3.

Add k cycles of odd length l>=3, put d'=d+lk and R'=n-d'. Then

 C=O_(d')/(O_d l^k O_k) is a nonnegative integer,
 D=k(l-1)/2+l binom(k,2)-v2(k!)+k R'-2k b_l,
 W' += C W 2^D,
 b'_j=b_j+k(j-gcd(l,j))/2.

The scalar correction follows by subtracting the old exponent R q_1 from
k q_l+R'(q_1+k). Since q_l=l q_1-2 b_l and R-R'=lk, the q_1 term cancels.
C is integral because (d')!/(d!l^k k!) is integral and taking its odd part
removes exactly the powers of two in the displayed factorials.

D may be negative. The right shift is nevertheless **exact on the state
weight**: each underlying term after multiplication has the integer form
O_(d') 2^(e(mu')+R'q'_1)/z(mu'). Their sum does too. This justifies exact
right shifts, not truncating division. The implementation checks divisibility.
No modular truncation is authorized by this argument: discarded bits could be
needed by a subsequent right shift.

After processing length 3, only units remain. For R=n-d the contribution to
O_n a(n) is

 W(d,empty) * O_n/(O_d O_R) * 2^[binom(R,2)-v2(R!)].

The final exponent is nonnegative for every R>=0, and final division by O_n
is exact by Burnside. This terminal contraction computes only a(n), not the
whole prefix. `gauge=False` supplies an ablation using the old responses but
the same exact unit-tail contraction.

A late switch is valid: before switching, store O_d sum 2^e/z and ordinary
responses. At any descending length l, convert each weight by 2^((n-d)q_1),
replace its response by b, and merge. This avoids unnecessarily enlarged
integers early in the computation; `switch` controls only performance.

## 3. Packing, runtime, and limits

At the moment future lengths j<=l are retained, all processed cycles have
length at least l+2. Hence j q_1<=d, and 0<=b_j<=d/2. During a transition,
retained j<=l-2 and every processed cycle has length at least l, so the same
bound holds. Using width=max(1,bit_length(n)) for each packed digit therefore
prevents carries. Coordinates above the residual budget are irrelevant and
are erased.

Let S_l be the number of retained (degree,response) states before processing l.
The recurrence uses

 O(sum_l sum_states (1+floor((n-d)/l)))

big-integer transitions, plus polynomial-size precomputed transfer tables and
one O(n S_switch) elementary-coordinate conversion. Integers have O(n^2+n log n)
bits: e and R q_1 are quadratic and factorial scales have O(n log n) bits.
Memory is the number of live states times that polynomial bit bound.
Since ordinary response states map surjectively to these projective states,
with degree preserved and transitions commuting with the map, each layer has
no more states or transitions than the same ordinary-response single-target
algorithm. Strict reduction occurs in measured instances. This alone does not
prove faster wall time, because arithmetic sizes and quotient construction
change. The partition upper bound remains exp(O(sqrt(n))) times polynomial
factors; no polynomial state bound is proved.

In particular this integrates an inherited structural compression, not a claim that
A000568 is now instant for arbitrary n. Printing a(n) already requires
Theta(n^2) bits asymptotically. Repeated coefficient queries may favor the
inherited full-prefix engine instead of rerunning this target-specific one.

## 4. Exact {1,3} terminal Laurent polynomial

The Python implementation additionally closes the last two cycle lengths in
one contraction. At a projective state (d,b_3=b), put R=n-d and r=R-3k. Define

 C_(d,k)=O_n/(O_d 3^k O_k O_r),
 E_(R,k)=binom(R,2)+3k^2+(1-2R)k-v2(k!)-v2(r!).

Then the contribution to O_n a(n) is

 W(d,b) sum_(k=0)^floor(R/3) C_(d,k) 2^[E_(R,k)-2kb].

This follows by multiplying the 3-cycle insertion and unit-tail factors from
Section2 and simplifying the exponent; no asymptotic approximation is used.
The coefficients C_(d,k) and E_(R,k) depend only on d, not the response b, so
are precomputed once per degree. Set L=min_k(E_(R,k)-2kb), evaluate the
integer polynomial P=sum_k C_(d,k) 2^(E_(R,k)-2kb-L), then contribute
(W P) 2^L. If L<0 the final right shift is exact, as each fully cleared
Burnside term is an integer. This is an exact Laurent-polynomial evaluation
at a power of four; it is a finite version of the unit/prime-tail GF ladder.

It uses one weight-times-polynomial multiplication per state instead of one
weight multiplication per permitted 3-cycle multiplicity. It still performs
that many shifts and additions, and constructing large polynomial values can
outweigh the saved multiplications. No universal runtime improvement is
claimed. `terminal3=False` replays the direct two-stage contraction.

The C++/GMP engine independently implements Sections1–3, including exact
negative-shift checks, but not this Section4 contraction. Comparing its
projective and raw switches isolates the quotient within the same language;
comparing languages does not by itself establish a mathematical speedup.
