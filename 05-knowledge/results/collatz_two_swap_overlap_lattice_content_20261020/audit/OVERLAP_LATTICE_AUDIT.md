# Independent audit: overlap roots and a uniform subexponential content bound

20 October 2026. Status: **PROVED** for the precise families below, with an
independent **FINITE-EXACT** audit. No all-length integer-cycle exclusion is
asserted. No file in the math repository or either frozen endpoint packet was
modified by this audit.

## 1. Verdict and exact scope

The proposed coprimality lemma is correct, including dependent projected
relations, nonprimitive relation rows, negative exponents, constant
polynomials, and composite moduli. Minkowski's second theorem supplies the
claimed product bound for two independent relation vectors.

The following explicit bound follows:

    g_plus  <= 38^(4 p^(2/3)),
    g_minus <= 46^(4 p^(2/3)) < exp(16 p^(2/3)).

Here p>=2, q=p+r, 1<=r<p, gcd(p,r)=1, and the word is obtained from the
zero-intercept lower mechanical source by one of the following overlapping
two-distinct-one chains. Source residue labels are a_j=(-jr mod p).

- Plus: the pair is {u,u+r}, where 0<=u<min(r,p-r); move residue u right
  first, then residue u+r right.
- Minus: the pair is {u,u+r}, where max(0,p-2r)<=u<p-r; move residue u+r
  left first, then residue u left.

The quantity bounded is the actual content g=gcd(|2^(p+r)-3^p|,B(w)), retaining
every prime-power depth. Both signs of the clock and original-cut crossings
are included. There is no restriction r>p/2.

**Audited large-range corollary:** combining this bound with the frozen
Matveev/continued-fraction certificate excludes integer outputs for every
one of these overlap chains with p>=4096. Section 10 audits that transfer.
The complete residual 3<=p<4096 is not checked here.

The pair range 0<=u<p-r by itself is insufficient for physical legality when
r<p/2. In that lower band, plus additionally requires u<r and minus requires
u>=p-2r. The algebraic lattice lemma is valid throughout 0<=u<p-r, even when
the corresponding physical chain is not legal.

## 2. Inheritance and loss ledger

Closest proved mechanism: the moving-anchor sparse formula in
`../../collatz_two_swap_norm_barrier_20261020/symbolic/TWO_SWAP_SPARSE_FORMS_AND_BARRIER.md`, especially its
Sections 2--5. The endpoint content packet
`../../collatz_two_swap_endpoint_content_20261020/ENDPOINT_CONTENT_AND_SHORT_RELATIONS.md` already
retains the common clock root, but only its two endpoint equations become
linear immediately.

Canonical hostile: the legal plus norm barrier at p=28,r=17,u=10 has norm
above 3^p but actual content 1. The repaired near miss is replacing actual
content by a norm that forgets the clock. The least-used relevant sidecar is
the full modular root together with the labelled intermediate word.

The concept board is: labelled gap legality; the common modular root; a
rank-three congruence lattice; local valuation at 5; integral resultants;
prime-power content. The useful bridge is

    actual carry -> common modular root -> two relation polynomials
                 -> nonzero integer resultant divisible by the full content.

The preserved predicate is divisibility of the actual carry and clock by g.
The resultant is only an upper bound, not an exact formula for g. The
physical residue/leader data are retained separately; lattice existence
does not confer physical legality.

## 3. All-r legality, signs, and denominator audit

The lifted source positions are

    s_j=j+ceil(jr/p)-1.

An individual right move at residue a is initially legal precisely when a<r.
An individual left move at residue a is initially legal precisely when
a>=p-r. Consecutive source ones with a gap of one have residues u+r on the
left and u on the right, with 0<=u<p-r.

For a right overlap, the rightmost one must move first; it can do so exactly
when u<r. Its move then opens the destination for the leftmost one. For a
left overlap, the leftmost one must move first; it can do so exactly when
u+r>=p-r. Its move opens the destination for the rightmost one. This proves
both displayed ranges. Distinct ones moving in opposite directions cannot
share an actual swap site without collision. Undoing one swap concerns the
same one, not this family.

For any divisor h of |D|, with D=2^(p+r)-3^p, h is coprime to 6. Choose
integers a,b with ap+br=1 and define

    z=2^a(3/2)^b modulo h.

Then z^p=2 and z^r=3/2. The common solution is unique by the same Bezout
identity. Also z-1 is a unit, because

    (z-1)(1+z+...+z^(p-1))=1.

The inherited sparse carry equation is

    1+(z-1) sum_j (2^(d_j)-1) z^(a_j)=0 modulo g.

For plus, the two summands are z^u and z^(u+r), giving

    2+5 z^u(z-1)=0.

For minus, both coefficients are -1/2, giving

    4-5 z^u(z-1)=0.

Consequently 5 is a unit modulo g in either case; otherwise reduction modulo
5 would force 2=0 or 4=0. The exact root equations are therefore

    plus:  z^u(z-1)=-2/5,
    minus: z^u(z-1)= 4/5.

These signs and constants are valid over the whole possibly nonsquarefree
modulus. The minus equation uses the cleared expression, not an unjustified
integral interpretation of the original coefficient -1/2.

The boundary p=2,r=1,u=0 is legal in both directions and has content 1. This
bound does not remove its genuine negative-cycle behavior.

## 4. Exact cleared relation polynomials

Write delta=1 for plus, delta=2 for minus, and epsilon=-1 or +1 respectively.
Thus K=epsilon*2^delta/5. Consider the lattice

    L = {(B,C,m) in Z^3 : m-Br-Cu is divisible by p}.

For a nonzero row v=(B,C,m), put

    A=(m-Br-Cu)/p,
    e=A-B+delta C.

The three modular root equations imply

    z^m(z-1)^C = epsilon^C 2^e 3^B 5^(-C).          (4.1)

For an integer t write t_+=max(t,0), t_-=max(-t,0). Define coprime integers

    n=epsilon^C 2^(e_+) 3^(B_+) 5^(C_-),
    d=             2^(e_-) 3^(B_-) 5^(C_+)>0.

The cleared polynomial is

    F_v(X)=d X^(m_+)(X-1)^(C_+)
             -n X^(m_-)(X-1)^(C_-).               (4.2)

Every denominator used to obtain (4.2) is a unit modulo g, and z,z-1 are
units. Thus F_v(z)=0 modulo the whole g.

### Nonzero, no spurious roots, and constant cases

No nonzero row produces the zero polynomial. If m or C is nonzero, a
functional identity X^m(X-1)^C=n/d is impossible, since its orders at 0 and
1 would give m=C=0. If m=C=0, the identity n/d=1 forces B=0 by the 3-adic
valuation, and then A=0 and v=0.

There is no root 0. If m!=0, exactly one term in (4.2) is nonzero at 0.
If m=0, vanishing at 0 would force n/d=(-1)^C. Its 5-adic valuation forces
C=0, and then the preceding argument forces v=0.

There is no root 1. If C!=0, exactly one term in (4.2) is nonzero at 1.
If C=0, vanishing at 1 forces n/d=1, whence B=A=m=0.

The polynomial is constant if and only if m=C=0. For opposite nonzero signs
of m and C, its two term degrees are |m| and |C|. Cancellation of the leading
term could occur only when these are equal and n/d=1, which contradicts
v_5(n/d)=-C!=0. All other nonconstant cases have plainly positive degree.
In particular,

    deg F_v <= |m|+|C|.

Two independent lattice rows cannot both yield constants, since the subspace
m=C=0 is one-dimensional.

## 5. Coprimality theorem

Let v_1,v_2 be any two linearly independent rows in L. Then F_v1 and F_v2
are coprime over Q.

Suppose instead they have a common algebraic root alpha. Section 4 shows
alpha is neither 0 nor 1, so both rational equations (4.1) hold at alpha.
Take a valuation nu at any prime over 5 in Q(alpha), normalized by nu(5)=1.
Since 2,3 and epsilon have valuation zero, the two equations give

    m_i nu(alpha)+C_i nu(alpha-1)=-C_i,   i=1,2.

If the projected pairs (m_i,C_i) are independent, the unique solution is
nu(alpha)=0 and nu(alpha-1)=-1. But the nonarchimedean triangle inequality
gives nu(alpha-1)>=min(nu(alpha),nu(1))=0, a contradiction.

If the projected pairs are dependent, choose integers s,t, not both zero,
whose combination kills m and C. The full row

    v=s v_1+t v_2=(B,0,0)

is nonzero by full-row independence. In particular B!=0. Multiplying the two
root equations to the powers s,t gives

    1=2^(A-B)3^B,

contradicting the 3-adic valuation of this rational equality. Equivalently,
because gcd(p,r)=1, write B=pj and A=-rj with j!=0, obtaining

    1=2^(-(p+r)j)3^(pj).

The reciprocal form of this last equation is equally contradictory. The
negative sign in the power of 2 is the sign supplied directly by the stated
lattice convention.

This covers rank-zero projections as well: both projections cannot be zero
under full-row independence. A nonzero constant polynomial already has no
root. No primitivity or basis assumption on the two selected rows is needed.

## 6. Minkowski: the product is genuinely p^(2/3)

The map Z^3 -> Z/pZ given by (B,C,m) -> m-Br-Cu is onto, because the
coefficient of m is 1. Hence L has index and covolume p. One explicit basis
is

    (1,0,r), (0,1,u), (0,0,p).

Let lambda_1<=lambda_2<=lambda_3 be the successive minima with respect to
the closed cube [-1,1]^3. Minkowski's second theorem gives

    lambda_1 lambda_2 lambda_3 * 8 <= 8p.

The exact form used is Theorem 3.3 in Ben Green's author-hosted
[additive-combinatorics notes](https://people.maths.ox.ac.uk/greenbj/papers/addcomb2009-3.pdf),
Section 3.7, printed page 27. Only the upper product inequality is imported.
No claim that minimizing vectors form an integral lattice basis is required.

Since lambda_3>=lambda_2>=sqrt(lambda_1 lambda_2), it follows that

    lambda_1 lambda_2 <= p^(2/3).

The cube minima are attained, and are positive integers because L is a
subgroup of Z^3. Choose independent rows v_1,v_2 with sup norms
H_1=lambda_1 and H_2=lambda_2. Thus

    H_1 H_2 <= p^(2/3).                             (6.1)

This bounds their product, not each vector separately by p^(1/3).

## 7. Explicit coefficient and resultant bounds

For any allowed u, r+u<=p-1. If H=max(|B|,|C|,|m|), then

    |A| <= (|m|+r|B|+u|C|)/p <= H.

Replacing v by -v multiplies F_v by the unit -epsilon^C, so it does not
change coefficient norms or degree. We may therefore orient every row so
C>=0. Then

    ||F_v||_1
       <= 2^(e_-+C) 3^(B_-) 5^C + 2^(e_+) 3^(B_+). (7.1)

The first term is at most 30^H:

- If B<=0, then e_-<=H+B, so it is at most
  2^H(3/2)^(-B)10^C <=30^H.
- If B>=0, it is at most
  max(10^C, 2^(H+B)(10/2^delta)^C).
  This is at most 20^H for delta=1 and 10^H for delta=2.

For the second term, e_+<=H-B+delta C, whose right-hand side is nonnegative.
If B<=0 this gives the bound 2^((delta+2)H). If B>=0, it gives

    2^(H+delta C)(3/2)^B
        <=(3*2^delta)^H <=2^((delta+2)H).

Since H is a positive integer, (7.1) implies

    plus:  ||F_v||_2 <= ||F_v||_1 <=30^H+ 8^H <=38^H,
    minus: ||F_v||_2 <= ||F_v||_1 <=30^H+16^H <=46^H.

Let R=Res(F_v1,F_v2), which is a nonzero integer by Section 5. The integral
resultant Bezout identity evaluated at the common modular root implies
g divides R. This is valid in Z/gZ without assuming it is a field or reduced
ring. In a constant case the same conclusion follows directly, or from
Res(c,F)=c^(deg F); both polynomials cannot be constant.

For nonconstant rows, Sylvester-Hadamard gives

    |R| <= ||F_v1||_2^(deg F_v2) ||F_v2||_2^(deg F_v1).

Since deg F_vi<=2H_i, this yields

    |R| <= c^(4 H_1 H_2) <= c^(4 p^(2/3)),

with c=38 or 46. The same weak inequality covers a constant row. Thus g<=|R|
proves the bounds in Section 1. Finally

    exp(4) > sum_(j=0)^6 4^j/j! >46,

so 46^(4 p^(2/3))<exp(16 p^(2/3)). No numerical approximation to a logarithm
is needed.

## 8. Exact audit and hostile controls

`audit_overlap_lattice.py` imports no repository mathematics. It independently
generates actual bitwords, follows tagged adjacent moves, computes carries,
enumerates short lattice rows, constructs integer polynomials, and computes
exact polynomial gcds and resultants with SymPy. It uses explicit checking
functions, not Python assert statements.

The completed receipt has **914,397 exact checks**:

- Every ordered actual two-edge path for 2<=p<=18, all 1<=r<p with gcd(p,r)=1:
  101 source shapes, 21,308 paths, and exactly 630 distinct-one overlapping
  paths. The proposed all-r residue and leader classifier is exact.
- Every legal same-sign overlap for 2<=p<=120:
  4,385 shapes and 175,670 paths, including 2,565 negative-clock shapes and
  4,385 original-cut-crossing paths. Of these, 6,249 have nontrivial content.
- Every algebraic triple 2<=p<=36, 1<=r<p, gcd(p,r)=1, 0<=u<p-r:
  4,759 lattices and 9,518 signed-polynomial resultant tests. The complete
  sup-norm shell enumeration checks (H_1 H_2)^3<=p^2, the coefficient bounds,
  no roots 0 or 1, coprimality, nonzero resultant, and divisibility by the
  actual full content wherever the shape is physical.
- Separate forced dependent-projection, constant, and negative-constant
  controls cover degeneracies not selected by the shortest-vector batches.
  The selected shortest pairs happened to have independent projections;
  that finite observation is not used in the proof.
- Known endpoint prime-power controls remain exact:
  plus (p,r,u)=(418,237,180) gives g=11849=17^2*41;
  minus (203,165,0) gives g=169=13^2.

There are also nonendpoint prime-power controls, including plus
(p,r,u)=(18,5,2) with g=6647=17^2*23, minus (34,19,5) with g=49, and plus
(62,5,4) with g=343. Thus the audit is not limited to squarefree content,
critical slopes, endpoint chains, or trivial gcds.

Reproduce in this directory:

    python3 audit_overlap_lattice.py > receipt.replay.json
    python3 -O audit_overlap_lattice.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json

The stored ordinary and optimized receipts are byte-identical.

## 9. Consequences and stopping boundary

The result supplies a uniform subexponential upper bound for every actual
overlapping two-distinct-one chain at k=1, not just two endpoint offsets.
For example, p>=32768 gives g/3^p<exp(-p/2), using p^(1/3)>=32 and log 3>1.

This is still not an all-length exclusion. Integrality would additionally
require |D|=g, and therefore an exceptionally small two-logarithm form.
Section 10 transfers an already certified Diophantine exclusion for p>=4096.
The complete remaining finite residual would still need verification.
Disjoint same-sign edits, mixed edits, and other two-edit families are
outside this result.

The proof mechanism generalizes to a root condition z^u(z-1)=K whenever
K has negative valuation at a prime absent from the two clock constants:
two independent projected relations force an impossible valuation pair.
Uniform coefficient constants and physical legality must be re-audited
for each such proposed extension.

## 10. Audited transfer: no integral overlap at p>=4096

This is a later, separate application of the proved content bound. It uses
the frozen arithmetic certificate in
`../../collatz_two_swap_endpoint_exclusion_20261020/large-range/` and its independent audit,
not the endpoint-specific content formula or endpoint finite residual.

The positive exponential series gives the exact certificate

    sum_(j=0)^10 (23/6)^j/j!
      =10119867931808989/219419659468800 >46.

Thus log46<23/6. Also the first two positive terms of the artanh series
give log3>13/12. If p>=4096=16^3, then

    log(g/3^p)
       <=4 log46 p^(2/3)-p log3
        <(46/3) p^(2/3)-(13/12)p
       <=((46/3)/16-13/12)p
        =-p/8<-p/9.                                (10.1)

If the word were integral, set U=2^q/3^p-1=D/3^p. Then

    0<|U|=g/3^p<exp(-p/9).

Every subsequent hypothesis of the frozen arithmetic argument holds:

- p<q=p+r<2p;
- gcd(p,q)=gcd(p,r)=1;
- U is nonzero by unique factorization;
- p>=4096 is above the thresholds 25 and 100 used for the logarithm and
  Legendre steps.

The source-dependent real Matveev estimate therefore gives p<10^12.
The bound |log(1+U)|<2 exp(-p/9), valid for either sign of U, forces q/p to
be an ordinary continued-fraction convergent of log3/log2. The frozen exact
certificate excludes all 18 such convergents with 100<=p<10^12. In
particular it excludes all 16 whose denominators are at least 4096.

The endpoint labels and r>p/2 appear only in the original source of its
content estimate. Neither is a hypothesis of the Matveev specialization,
the absolute-logarithm inequality, the Legendre criterion, or the
continued-fraction exhaustion. Equation (10.1) replaces that original
content estimate legitimately.

`audit_large_overlap_transfer.py` checks the rational threshold calculation,
replays the frozen independent large-range audit without writing to it,
compares its output byte-for-byte with the frozen receipt, and checks all
16 relevant convergent rows. Its normal and optimized receipts agree.
It recognizes both staging and canonical sibling packet names; for other
layouts, pass `--endpoint-dir /path/to/the/frozen/endpoint-packet`.

The explicit remaining universe is every legal overlap at every coprime
q=p+r with 3<=p<4096 and 1<=r<p, with both signs of D allowed. No check of
only critical clocks, only positive D, or only the two endpoint offsets
would discharge that residual.

## 11. Final producer-draft review

The complete `../OVERLAPPING_TWO_ONE_CONTENT.md` draft was inspected after the
producer's sharper-constant update. Its root-0/root-1 treatment, constant
classification, both projected-rank cases, successive-minimum product,
38/46 coefficient estimates, and separately scoped finite reduction all pass.
Its explicit partial-factor hostile is also correct: at p=5,r=3,u=1, the
stated quadratics have resultant 4225=25*13^2 while actual content is 1.
The constant-row variant gives resultant 169 without changing that content.

The mathematical statement may be promoted from proof candidate to
**PROVED**, retaining Minkowski as the core imported input and Matveev plus
the frozen exact arithmetic certificate only for the p>=4096 corollary.
The theorem must continue to distinguish that large-range exclusion from
an unperformed complete finite residual check.
