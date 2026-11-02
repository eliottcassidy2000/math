# All-modulus context arithmetic: independent audit

Status: core theorems and constructions independently VERIFIED against
`/tmp/session139_research/all_moduli/PROOF.md`; finite orbit audit is complete with final status below.

## All-modulus factorial-transform theorem

For Psi=sum_(c>=1) (c-1)! p(c) (x)_c, the c-th forward difference at0
is exactly c!(c-1)!p(c). Terms of lower degree die under the difference;
terms of higher falling-factorial degree vanish at every sampled argument
0,...,c. Thus equality of the Z/MZ polynomial functions implies equality of
all normalized context coordinates. Conversely Newton interpolation at
nonnegative integers uses the integer coefficients binom(x,c), not division
by a modular nonunit. Combining with the accepted exact contextual-signature
theorem gives precisely the claimed all-M iff. This is a function identity,
not equality of formal polynomials or a claim of arbitrary-function realization.
The actual image has not been shown closed under ring operations.

At odd primes, factorial coordinates below p are units and all later
coordinates vanish; triangular evaluation at1,...,p-1 is invertible. Hence
ordinary F-functions suffice. The parity upper bound2p^((p-1)/2)-1 correctly
counts the common zero only once. Ordinary F cannot replace Psi at composite
moduli: L4 and L6 have the same mod16 context state but values3^4=1 and3^6=9.

## Exact transitive clocks

The odd-prime proof correctly recovers all j^n for1<=j<p from c<p
coordinates. Equality for n<m forces t=m-n divisible byp-1. The choice
j=p-1 then forces p^(a-1)|t: for even t,
`v_p((p-1)^t-1)=1+v_p(t)`. The c=p coordinate separates every n<a,
because its remaining difference has exact valuation n. Conversely unit
powers are periodic and nonunit powers vanish once n>=a. Thus both the
preperiod and the least period are sharp, not merely valid upper bounds.

At powers of2, q2=2^n-2 forces the exact preperioda (for a>=2), and
q3=2(3^n-3*2^n+3) forces the order of3 modulo2^(a-1). Every c>=3
coordinate has the extra factorial factor2 needed for sufficiency. The
small a=1,2 cases and the conventionM=1 are handled separately. CRT takes
the maximum transient and lcm of periods, without any hidden coprimality
requirement on the periods themselves.

`audit_moduli.py` uses the finite-difference power-sum formula directly,
not the producer's Stirling recurrence. It independently identifies the
first repeated state at all M<=200 and nine additional prime powers
(209 moduli), matching the exact formulas. The all-order theorem comes
from the preceding proof, not this finite range.

## Literal5/11 zero and projectors

An independent final-predecessor Hamiltonian DP and least-vertex partition
recursion constructs W5=C3[L2,1,1], W11=C3[L2,L2,1], and their ordinal join.
It finds orders4,5,9 and H=5,11,55; the full order9 path-cover vector is

    (0,55,1545,7495,12310,8665,2898,474,36,1).

Every normalized coordinate is zero modulo5. The exact profile functions
are x^4+4x^2 and x^5+8x^3+2x. Their product vanishes on F5, but its degree9
formal polynomial is nonzero modulo5. This hostile distinction is retained.
The supports are respectively +/-2 and +/-1, so the fourth power of W5
and square of W11 give the claimed actual positive projectors. Formal sum
of their functions is not interpreted as addition of tournaments.

Both `audit_moduli.py` normal and -O receipts agree byte-for-byte.

## Final full-function orbit audit

**Final status: PASS.** The previously pending finite-orbit audit completed.
`audit_function_states.py` uses full p-entry evaluation tuples, with neither
parity tags nor half-vector compression. Exhaustive queue closure yields
5,48,685,322100 distinct unary functions at p=3,5,7,11. Zero is reached at3,7
and absent at5,11. At5,11 the full signed point-pair support vectors at +/-1
and +/-2 are explicitly members of the reached set, so their ordinal join
realizes zero. The complete all-context class counts are therefore
5,49,685,322101, matching the independently proved upper bounds.

The same independent script checks7,776 arbitrary signed degree5 profiles
(243 per modulus at M=1,...,32) for factorial-function/signature equality
and explicit forward differences. Its normal and `-O` outputs are
byte-identical, SHA256
`db5ce7a9fe221b01e405da78d6699112b42a262dcb9f27501952d26623a6efc3`.
Script SHA256:
`06501a35b37fbb3851f299e5bd7eee6a9f4afa8ebd6553cd86194e83ef5a7b87`.

The final producer source versions, verified dependency pins, universal
P-padding corollary, and a separate permutation/cut graph enumeration are
covered in [ALL_MODULI_PADDING_AUDIT.md](ALL_MODULI_PADDING_AUDIT.md).

## Independent uncompressed prime-function orbit replay

`audit_function_states.py` stores complete p-entry functions, not parity
tags or compressed half-vectors. Exact exhaustive BFS gives5 functions at3,
48 nonzero functions at5,685 functions at7,and322100 nonzero functions at11.
At5 and11 it verifies reached disjoint-support odd function pairs, whose
ordinal product supplies the missing zero. Thus class counts agree exactly
with5,49,685,322101. Normal and -O outputs are identical.

The same independent script checks7776 signed coefficient profiles (243
degree5 profiles at each M=1..32), exact normalized-signature/function
equivalence, and recovery by forward differences. These are finite controls
for the separately proved all-modulus theorem. The latest universal padding
corollary was handed to the dedicated all-modulus auditor, not included in
this review's original scope.
