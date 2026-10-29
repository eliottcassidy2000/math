# Nine carries decode every odd period; sparse certificates close two composite groups

2026-10-29. **PROVED elementary**: inverse-reader algorithm and sparse
primitivity certificates. **PROVED using a CITED CFSG-dependent theorem**:
G_25=A_167761 and G_55=A_312119004989. These are two finite-specific group
identifications, not a theorem for arbitrary composite Lucas degrees.

This extends `COMPATIBLE_CRT_CARRIERS.md`. The parent suggested applying
Jones after the first sparse primitivity certificate. The independent
audit lane sharpened the original 44-candidate decoder to nine candidates;
the implementation now uses the nine-candidate version.

## 1. All-odd-period inverse reader with at most nine carry candidates

Use the inherited all-odd reader bijection, with our exponent convention

    q(w)=sum_(i=0)^(n-1) w_i phi^i mod alpha,  alpha=phi^n-1=a+b phi,
    a=F_(n-1)-1, b=F_n, Norm(alpha)=-L_n.

Given **any** target t=u+v phi, not merely an integer residue, write

    -t/alpha = (A+B phi)/L_n,
    A=u(a+b)-vb, B=-ub+va.

Let c_a=floor(A/L_n), c_b=floor(B/L_n). The unique legal word can be found
by testing just

    z=a_0+b_0 phi,
    a_0 in {c_a-1,c_a,c_a+1},
    b_0 in {c_b,c_b+1,c_b+2}.

For each candidate form F=t+alpha z. Greedily subtract phi^(n-1),...,1,
recording the binary digits. Accept only if the remainder is zero and the
word has no cyclic adjacent ones. Exactly one candidate word is accepted.
The implementation uses exact integer pairs and exact signs in Q(sqrt5),
not floating-point comparisons.

**Proof of completeness.** For the true word, F is its ordinary golden
sum and z=(F-t)/alpha is integral. Put delta=F/alpha and let psi=1-phi.
The cyclic legality bound from the inherited reader proof gives
`0<=delta<=1` at the real embedding. Without even needing legality,

    -1 <= F' <= phi,

because the total of all odd negative powers of psi is at least -1 and
all even positive powers sum to at most phi. Since n is odd,
`alpha'=psi^n-1<-1`, so `-phi<=delta'<=1`.
Writing x=delta, y=delta', the coefficient offsets from -t/alpha are

    delta_b=(x-y)/sqrt5 in [-1/sqrt5, phi²/sqrt5] ⊂ (-1,2),
    delta_a=(phi*y+phi^(-1)*x)/sqrt5 in [-phi²/sqrt5,1] ⊂ (-2,1].

Integrality now gives exactly the stated candidate rectangles. For example,
if b_0<=c_b-1 then b_0-B/L_n<=-1, contradicting delta_b>-1;
if b_0>=c_b+3 then b_0-B/L_n>2, contradicting delta_b<2.
Likewise a_0<=c_a-2 contradicts delta_a>-2, while a_0>=c_a+2
contradicts delta_a<=1. These arguments include every possible fractional
part of A/L_n and B/L_n; no rounding or approximate floor is used.

The true word is therefore included. Greedy extraction is valid because,
after omitting a digit phi^i, every finite no-adjacent-one suffix is
strictly less than phi^i. Uniqueness follows from the inherited odd-reader
theorem. The algorithm uses O(n) algebraic arithmetic operations with a
fixed constant nine; this is not an O(n) *bit-complexity* assertion.
`fast_reader.py` checks 1625 source cases against enumerated banks before
producing the large certificate, and normal/-O outputs agree. The all-order
proof, not that finite comparison, establishes completeness.

## 2. Sparse period-55 primitivity certificate, including a nonreduced factor

The degree and selected root are

    N=L_55=312119004989=11²·199·331·39161,
    phi=252509579304 mod N.

The normal-form theorem makes +1 a full N-cycle. Every nontrivial block
system compatible with +1 is reduction modulo a divisor d with 1<d<N.
There are exactly 22 such divisors. `period55_certificate.json` gives,
for each d, two legal 55-bit words whose residues x,y agree modulo d but
whose reflected residues disagree. All pairs use y=x+d and x in {0,1,2}.
`verify_period55.py` independently checks all 44 word occurrences,
factorization, complete divisor coverage, root and ideal norm. It never
imports the decoder or enumerates the N words. Therefore G_55=<+1,R> is
primitive. The separate period-25 certificate uses six pairs in the same
way. Normal and optimized verifiers agree.

## 3. Jones classification now applies at these composite degrees

The primary source is Gareth A. Jones, *Primitive permutation groups
containing a cycle*, Theorem 1.2, case k=0,
[primary paper](https://arxiv.org/pdf/1209.5169), checked 2026-10-29.
For a primitive group containing a full cycle, failure to contain A_N
requires prime affine degree, projective degree
`N=1+q+...+q^(d-1)` for prime-power q and d>=2, or one of the exceptional
degrees 11 and 23. This classification uses CFSG. The theorem is used here
only after the independently certified primitivity step.

Both degrees are composite and are neither 11 nor 23. If d=2 then q=N-1
is even but is not a power of two, excluding every prime power. For d>=3,
q²<N, so it suffices to test all integers 2<=q<=floor(sqrt N) and sum
successive powers until exceeding N. `check_jones_degrees.py` finds no
solutions for either N, checking q through 409 and 558676, respectively.
This rules out even non-prime-power repunits and hence every projective
exception.

Thus both groups contain A_N. The generator +1 is even because N is odd.
For n=2m+1, reflection-fixed legal words correspond to no-adjacent-one
strings of length m: the last reflected pair must be zero, leaving
F_(m+2) choices. Therefore reflection has
`(L_n-F_((n+3)/2))/2` transpositions. These are

| n | N | reflection-fixed words | transpositions |
|---|---:|---:|---:|
|25|167761|377|83692|
|55|312119004989|514229|156059245380|

Both transposition counts are even. Hence both generators lie in A_N, and

    G_25 = A_167761,
    G_55 = A_312119004989.

The group identifications are derived from sparse exact certificates and
classification, not from generating or enumerating these enormous groups.
No identification with a tournament automorphism group is asserted. In
fact, alternating two-transitivity precludes a nonempty noncomplete
invariant loopless directed binary relation, so any tournament decoder must retain
additional symmetry-breaking data.

## Verification boundary

- All-odd reader bijection: inherited all-order theorem, source paths in the
  companion document.
- Nine-candidate inverse reader: elementary all-order proof above, sharpened
  from an initial 44-candidate implementation by the independent audit lane.
- Period25/55 sparse primitivity: direct exact word certificates.
- Alternating-group conclusions: those certificates plus Jones/CFSG and
  complete finite degree exclusions.
- No positive Collatz convergence, prime-value infinitude, or universal
  composite-Lucas group classification follows.
