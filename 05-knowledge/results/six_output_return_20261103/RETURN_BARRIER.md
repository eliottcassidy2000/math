# Exact dyadic return gates and the missing ternary carry

**PROVED elementary statements R1--R4; FINITE-EXACT applications.** The
arithmetic grammar is closed. Its full physical tournament realization
remains **OPEN**.

## R1. Two shifted positive dyadic rays meet at most twice

Consider `U*2^m+V=A*2^n+B`, with positive integer U,A and m,n>=0. Write

    U=2^u a, A=2^v b, a,b odd,
    M=m+u, N=n+v, D=B-V.

If D=0, solutions exist only when a=b, in which case M=N gives an infinite
shifted diagonal. If D!=0, let t=v_2(D). There are only three cases:

* M=N: `(a-b)2^M=D`, giving at most one solution.
* M<N: M=t, and `b2^N=a2^t-D` determines N uniquely if it is a power.
* N<M: N=t, and `a2^M=D+b2^t` determines M uniquely if it is a power.

The last two cannot both occur: after division by 2^t their right-hand
differences are respectively strictly below and strictly above a-b.
Thus there are at most **two isolated intersections**, unless the rays
coincide up to a shift. The bound is sharp: `3*2^m=2^n+2` has solutions
`(m,n)=(0,0),(1,2)`.

`ray_gate.py` implements exactly these cases, with no exponent cutoff.
An independent 11664-case bounded grid checks the implementation. It then
compares all six output rays with all four original input rays. **Every one
of the 24 intersections is empty at all nonnegative exponents.** Therefore
none of the six outputs can be returned to those original inputs even by
matching its scalar count and arbitrarily changing clone size. This is
stronger than failure of a particular isomorphism search.

## R2. No finite dyadic-ray bank can close under a ternary offset

Let B be a finite union of rays `A_i*2^m+B_i`, A_i>0, m>=0. It cannot be
closed under `x -> 3x+c` for any fixed c. The same obstruction applies if
each point may choose one of a fixed finite collection of offsets.

**Proof.** An infinite output ray must have infinite intersection with some
input ray, because finitely many isolated intersections cannot cover it.
R1 implies coincidence up to a shift, hence

    A_j = 3 A_i * 2^z                                  (R1)

for an integer z, and `B_j=3B_i+c`. The 3-adic valuation of A_j is exactly
one more than that of A_i. Taking an input of maximal slope valuation gives
a contradiction. If the offsets vary among finitely many choices, at least
one choice occurs infinitely often and the same argument applies. ∎

Thus a finite collection of fixed cores with only one transitive-size
parameter is the wrong target for indefinite return. The obstruction is
independent of which state-dependent paths are chosen. In any b-power
version, a cycle of ray returns would force `lambda^k=b^l`; multiplicative
independence obstructs such a finite cycle as well.

Even adding a ternary height while keeping only finitely many intercepts
does not close both offsets. Uniform equality of rays forces
`B -> 3B+2` and `B -> 3B+4`. On a finite intercept set each injective map
would be a permutation; its only finite real orbit is its fixed point,
respectively -1 or -2. No nonempty finite set can satisfy both.

## R3. The exact arithmetic grammar retains height and ordered carry

For a temporal word `w=(c_1,...,c_k)`, each c_i in{2,4}, put

    C(w)=sum_(i=1..k) c_i 3^(k-i),
    D(w)=sum_(i=1..k) ((c_i-2)/2) 3^(k-i).

Then

    H_w(m)=3^k(A*2^(m-1)+B)+C(w),
    C(w)=3^k-1+2D(w).                                   (R2)

D(w) has k ternary digits in{0,1}; distinct length-k words have distinct
integer carries. Updates are `(k,C)->(k+1,3C+c)`. This is a small *formula*
with an unbounded integer state. It is not a finite bank, and a formula for
a target count does not by itself construct its source tournament.

**The length is available for free from the exact carry.** C+1 has ternary
digits `1` followed by k digits in{0,2}. Thus k is its digit length minus one;
there is no need to retain k as an independent coordinate. Replacing each
suffix digit2 by 1 gives a binary heap address b, with updates
`b->2b+(c-2)/2`. The longest common prefix of two such addresses is their
exact arithmetic LCA. Either integer encodes the complete branch word;
C gives an affine count update while b gives compact prefix operations.
Neither determines the actual state-dependent Hamilton paths used to realize
the word, which remain part of the physical source data.

For the 1457 observer, when the required physical words exist, the companion
law is

    v_w=3^k v_0 + C(w)*(1,0,1).                          (R3)

The present literal certificates realize both one-letter words and the two
mixed two-letter words. The mixed outputs differ by `(4,0,4)`, as predicted
by `C(2,4)=10` and `C(4,2)=14`. No arbitrary-length physical realization is
asserted from (R3).

In the current F49 chart, the branch law is

    (a,b,g)->(3a,3b+c,3g) mod 7.

Its formal multiplier/carry pairs `(s,C)`, s in F7*, form42 states. The
linear part is a character; the translation is a separate cocycle. At fixed
length3, the two words242 and424 have carries32 and46, so the current phase
forgets their difference14. This is the first fixed-length collision;
the words have not both been physically realized at length3 in this packet.
The collision is an arithmetic warning against replacing integer carry by
its current finite-field phase.

The ancestral F343 interval-growth coordinate is retained from the previous
packet, including its actual n=3 source and field basis. These new return
steps are not interval growth and do not advance that clock. All current
five-Q4/F49 readouts are recomputed from the actual retained root.

## R4. A focused clone-splice route has an exact parity guard

Suppose a common-frame two-complement program takes a capacity-one source
with `H(m)=A*2^(m-1)+B` to an intermediate ray `a*2^(m-1)+b`, and the last
birth simply enlarges that same transitive block by one vertex. The final
count is `2a*2^(m-1)+b`. It implements the desired ternary offset iff

    a=3A/2, b=3B+c,
    H(intermediate at m=1)=3A/2+3B+c
      =(9H(source)-3H(source at clone size 2))/2+c.        (R4)

This follows from `H_V(3)=3H_V(2)-2H_V(1)`, with
`H_V(2)=3H_T(1)+c` and `H_V(3)=3H_T(2)+c`.
Writing A=2c_2, B is odd. Hence the proposed intermediate has odd H iff
c_2 is even. The two output heads 4379 and4381 have c_2=2769, so this entire
clone-splice route is impossible at their retained ports, regardless of
which complements are chosen. This is a route-specific parity obstruction.

If the newborn is absorbed into the module and the intermediate core is
used for another such normalized return, its slope changes by A->3A/2.
It must remain even, so at most `v_2(A)-1` consecutive absorbed returns can
occur before another structural change is necessary. This is a binary
precision budget, not a claim that the allowed returns actually exist.

The initial exploratory splice search mistakenly used an arithmetic-
progression relation instead of the displayed exponential recurrence.
Its negative search output was discarded; the corrected filter and its
independent algebra tests are the only ones used here. See the correction
entry in `01-canon/MISTAKES.md` and AUDIT.
