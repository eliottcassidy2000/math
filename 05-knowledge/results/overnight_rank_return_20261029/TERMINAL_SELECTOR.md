# Closing the periodic observer under terminal selectors

**2026-10-29. PROVED for the specified operation language, with FINITE-EXACT
integer certificates. Same-agent audit, no independent review or priority
claim. Arbitrary vertex re-marking and further actual complements remain OPEN.**

## S0. Inheritance, source and operations

The source is exactly [PERIODIC_COMPLEMENT E1--E5](PERIODIC_COMPLEMENT.md):
the actual complements U_(8+2r) associated with repeated two-vertex runs,
with three sources (0,1,2) and terminal vertices (n-3,n-2,n-1). Its reachable
boundary state set has size 35, bulk matrix B, initial row a, and terminal
column b. The scalar recurrence polynomial p=(X+1)f6 f11 and the finite
component-pairing transfer are inherited, not rediscovered.

We add the six possible assignments of these terminal vertices to the
three named paths. Write a permutation as its image word, for example
021. The selector pi assigns terminal vertex n-3+j to path pi(j).
Thus the endpoint of source i is n-3+pi^(-1)(i). Let b_pi be the resulting
terminal column. Sources, vertices, defect edges and the bulk clock stay
explicit. This is a change of prescribed linkage, not a tournament
automorphism or a new complement.

The closest older mechanism is the ordered cover / reversal-fixed cover
distinction in [Session24, PATHCOVERS sections 1--2](../exception_pathcover_transfers_20261021/session24-pathcovers/PATHCOVERS.md),
inheriting THM-3121, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`.
Its fixed-cover coordinate d is the least-used useful sidecar here.
The corrected near miss is inferring a free reflection pairing from an even
count; S6 supplies fixed covers. The canonical 81/27 hostile remains the
separate [full twelve-vertex calculation](DEFECT_TRANSFER.md#d6-full-three-run-matrices-of-the-8127-pair):
equal mass does not specify the joint response. No new twelve-vertex
minimality or integral recovery is inferred here.

Concept board: component pairing; selector orbit; integral observation
lattice; repeated clock factors; physical initial laws; reflection-fixed
covers. Incoming [TERNARY_EXTENSION_OBSTRUCTION](../overnight_research_20261029/TERNARY_EXTENSION_OBSTRUCTION.md)
raises the separate Paley7 clock-completion lower bound to 3^10. Its root
extension constraints and attainment question are inherited, not identified
with the terminal-selector problem below.

## S1. An exact finite operation closure

Permuting the three negative path-name tags acts on all 35 states; it
introduces no new state. Call this action P_pi. The verifier checks every
entry of the identities

    P_pi B = B P_pi,       b_pi(s)=b_012(P_(pi^-1) s).       (S1)

The second expression specifies the state action without a row/column
convention ambiguity. In particular, any interleaving of bulk insertions
and common path-name permutations reduces to a bulk power and one selector.
All six columns satisfy the SAME integer identity

    p(B)b_pi=0.                                           (S2)

Consequently form the 108-by-35 integer matrix O by stacking

    (B^j b_pi)^T,       pi in S3,       0<=j<18.             (S3)

For arbitrary boundary weights w over Z/m, the smallest exact quotient for
this language is image(O mod m), with state O w^T. Indeed equality there
implies equality of every later observation by S1--S2. Conversely every
coordinate in S3 is an allowed observation, so any exact representation
must separate different images. This is a quotient of weighted boundary
states, not a claim about the minimum among unconditioned graph prefixes.

This observation-lattice argument generalizes: for commuting clock and
finite name action, close the output lattice under both, using a monic
clock annihilator to bound powers. Its application here is exact because
S1 and S2 are checked on the full 35-state domain, not fitted on one orbit.

## S2. The new memory is multiplicity, without a new clock polynomial

**PROVED by the finite integer certificate:** rank_Q O=24; the original
single-selector rank was 18. Put

    t=sum_pi b_pi,     z=sum_pi sign(pi)b_pi,
    u_pi=6b_pi-t-sign(pi)z.

Their clock-generated spaces have dimensions 11, 1 and 12, respectively.
The exact annihilators are f11, X+1 and f6. The commuting S3 action is
trivial, sign and standard on these spaces. In rational clock-module terms,

    W_Q = Q[X]/f11  +  Q[X]/(X+1)  +  (Q[X]/f6)^2.        (S4)

The two standard copies follow from the two-dimensional standard S3
representation. Its multiplicity space has dimension six and minimal
polynomial f6: otherwise the inherited minimal degree-18 polynomial would
drop. Thus the scalar clock polynomial remains p although six additional
response directions become observable. Equal characteristic factors or
periods do not recover the relative amplitudes of repeated components.

## S3. Sharp selector counts and all-precision memory

The 24 nonzero Smith factors of O are

    1^9, 2^3, 4^6, 12^3, 24, 1344, 147676566374644608.    (S5)

Here exponents indicate repetition. The last factor is
2^7 * 3^2 * 7 * 18313066266697. Hence the exact number of weighted states
at modulus m is

    product_i m/gcd(m,d_i).                               (S6)

In particular, the closed carrier has 2^9 states modulo two and 3^18
states modulo three, versus 2^4 and 3^12 before selector closure. The
numbers of base-p digits at precision p^r are

    E2(r)=9r+3(r-1)_+ +9(r-2)_+ +(r-3)_+ +(r-6)_+ +(r-7)_+,
    E3(r)=18r+5(r-1)_+ +(r-2)_+.                           (S7)

The general precision recurrence is E_p(r+1)-E_p(r)=#{i:v_p(d_i)<=r}.
Independent unit-pivot / divide-by-p elimination checks every primary
layer of S5, including the two odd noncharacteristic primes.

A *fixed selector bank* retains all future bulk observations for its
selected terminal assignments. Its smallest size for recovering the full
six-selector observer is:

| Coefficients | Minimum fixed bank | An attaining bank |
|---|---:|---|
| Q | 2 | 012,021 |
| F2 | 4 | 012,021,102,120 |
| F3 | 2 | 012,021 |
| Z, hence every modulus simultaneously | 4 | 012,021,102,120 |

For the integral statement the row Hermite lattices agree exactly, so every
omitted row is an integer combination of the retained rows. This gives a
single reconstruction valid at every modulus. Every three-selector bank
containing 012 has binary rank at most eight, while O has rank nine.
Transport moves any bank to one containing 012, proving the lower bound
without that convention. All 32 banks containing 012 are checked. Six of
the ten four-selector banks are integrally complete; binary completeness
alone accepts eight, so even that test is insufficient over all integers.

For reference, adding one selector to 012 gives the following ranks:

| Second selector | Q | F2 | F3 |
|---|---:|---:|---:|
| 021 | 24 | 6 | 18 |
| 102 | 24 | 6 | 17 |
| 120 | 24 | 6 | 18 |
| 201 | 24 | 6 | 18 |
| 210 | 24 | 5 | 18 |

These are minima for fixed banks of full clock responses. No claim about
adaptive query cost, a single scalar read, or arbitrary vertex selectors
is made.

## S4. A positive, equal-mass separating witness

Use the inherited state order in `periodic_defect.json`. Let

    w_plus=e11+e18,       w_minus=e2+e4.                   (S8)

Both are nonnegative and have mass two. Their difference is killed by all
18 original observation rows and therefore by every B^j b_012, j>=0.
The terminal responses, before further insertion, are

| Selector | w_plus | w_minus |
|---|---:|---:|
| 012 | 1 | 1 |
| 021 | 1 | 2 |
| 102 | 2 | 1 |
| 120 | 1 | 1 |
| 201 | 1 | 1 |
| 210 | 1 | 1 |

Every supporting component state is reachable in the boundary automaton.
These are prescribed mixtures of partial-cover states. They are not
claimed to be the complete laws of two unconditioned tournament prefixes.
Thus S8 is an obstruction for arbitrary nonnegative weighted attachments;
the stronger physical-realization question is still open.

## S5. Frobenius phase and the extra field amplitude

The binary clock module has Jordan block lengths 4,2,2,1 at eigenvalue one.
The nullities of (B-I)^j on the nine-dimensional observable space are
4,7,8,9. Its clock order is four, although the actual prefix trajectory has
constant output modulo two.

Over F3 retain the inherited irreducible

    g=X^5+X^4-X^3-X^2-1.

The joint module has the exact decomposition

    W_3 = F3[X]/((X+1)^2 g^2) + F3[X]/((X+1)g).          (S9)

The kernel dimensions for powers of X+1 are 2,3; for powers of g they
are 10,15. Together with rank 18 and the full annihilator these certify
S9. Relative to the old twelve-dimensional observer, the added six
coordinates are one F3 amplitude and one F243 amplitude. They are
semisimple additions; the inherited observable square-zero layer remains.

The multiplication clock still has order 726. On each field quotient,
X has order 121 and Frobenius has order five, exactly as in E4. S1 does
not create a new Frobenius phase or identify the two independent field
amplitudes. A common field clock cannot supply their missing joint
alignment for free. Cubing still kills the old nilpotent layer and is
not an automorphism of the whole local algebra.

For the particular physical prefix a, all six scalar series obey the old
degree-twelve ternary recurrence. Their joint trajectory has period 726
with no transient. That physical trajectory is much smaller than the
arbitrary-weight carrier of size 3^18.

## S6. The precise parity boundary and the reflection near miss

For every r>=0 the six actual counts, in lexicographic selector order,
have parity

    (0,0,0,1,1,1).                                      (S10)

At n=8 their values are (6,6,6,5,3,3). Thus the old evenness theorem is
sensitive to endpoint matching. The unchanged bulk recurrence does not
force evenness for a different prescribed linkage. The verifier checks
the initial recurrence window and the constant parity solution, which
proves S10 at every size.

The natural reflection rho(v)=n-1-v is an anti-automorphism. On ordered
covers it acts by reversing each path, applying rho, and reversing the
list of paths. It preserves the identity endpoint assignment but is NOT
free. At n=8 it fixes, for example,

    ((0,4,5), (1,6), (2,3,7)).                            (S11)

Full independent cover enumeration gives:

| n | All identity-selector covers | Reflection-fixed covers |
|---|---:|---:|
| 8 | 6 | 2 |
| 10 | 38 | 6 |
| 12 | 684 | 26 |

This recovers why Session24 retained the fixed-cover coordinate d. Its
odd-order unique-fixed-vertex argument cannot be imported to these even
sources with three paths. A different path-level involution explaining
the identity-selector evenness remains an open target; S11 only refutes
the natural reflection as a free pairing.

## S7. Verification and next boundary

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/terminal_selector.py
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/terminal_selector.py
```

The checker rebuilds the reachable bulk transitions, all six terminal
columns, all common name transports, rational sector identities, integral
Smith factors, independent primary layers, all selector banks and both
weighted laws. Forty-eight native queries through n=22 agree with the
separate frontier implementation; eighteen queries through n=12 agree
with subset convolution. Direct cover enumeration checks the six n=8
matchings and reflection-fixed counts at n=8,10,12. No gate uses assertions.

The next precise questions are: which row-lattice distinctions are
realized by unconditioned prefixes; which actual complement operations
preserve a finite extension of this carrier; and whether a fixed-cover
transfer can produce a nontrivial pairing proof of S10. The general
all-size selected-path recovery conjecture and A000568 reconstruction
receive no new enumeration claim from this scoped result.
