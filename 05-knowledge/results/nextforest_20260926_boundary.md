# A carry defect with an all-height selected-descent certificate

**PROVED elementary bounds and lifting theorem / FINITE-EXACT certificate
for defect <=4096 / OPEN global Collatz.** No novelty or convergence claim.

The principal result is an actual selected-return inequality:

> For every positive odd integer n>1 whose ordered binary carry defect
> E(n), defined below, is at most4096, some j<=104 satisfies U^j(n)<n.

Here U(n)=(3n+1)/2^v2(3n+1). The quantifier has no upper bound on n.
The proof combines a local carry identity with683 finite offset
certificates and10,977 exact finite exceptions. A direct hand proof gives
j<=4 when E<=12; only n=7 needs four steps in that smaller class.
This asserts a smaller iterate, not termination: that iterate can leave
the bounded-defect class.

## 1. Inheritance, concept board and transfer

Closest mechanism: the ordered carry polynomial and exact forest gluing in
[the previous excursion note](forest_20260926_excursions.md), especially
its derivative at z=2. Canonical hostile: a growing reset can recreate
arbitrarily much positional information. Corrected near miss: near an
extremal carry pattern need not mean near a long division. Least-used
sidecar: which finite-state transitions create nonnegative slack, and the
position of the next restart after a zero-output region.

Anchor / niche / wildcard: a source-preserving selected return; the
bounded normalized jets encountered in
[THM-3000 / fixed-edge cumulant-curvature universality](../../01-canon/theorems/THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer.md);
and the exact one-bit perturbation of a direct preimage of1. The six live
concepts are state slack, output ownership, dyadic endpoint offset,
normalized jets, precision expenditure, and finite core lifts.

The older jet theorem concerns Newton ratios and root moments. Its theorem
does not apply to Collatz. The transferred operation is to normalize a
jet, separate scale from shape, and audit the resulting uniform bound.
Here the map is explicit: a binary integer goes to the probability masses
2^s/n at its occupied positions. It preserves normalized z=2 moments;
it destroys the unscaled low-bit mass as a robust coordinate. Exact low
bits are the needed sidecar. The pair differing by2 in section4 is the
decisive hostile.

## 2. A nonnegative state potential and sharp carry bounds

Write n=sum b_s 2^s, and start the odd multiplication carry at c_0=1:

    c_(s+1)=floor((3b_s+c_s)/2),
    q_s=3b_s+c_s-2c_(s+1).

Thus q_s are the bits of3n+1, and c_s belongs to{0,1,2}.
For the full polynomial C_n(z)=sum_(s>=1)c_s z^s, put

    C(n)=C_n(2),              E(n)=8n-C(n).

The independent derivative formula is

    C(n)=J(3n+1)-3J(n),       J(n)=sum_s s b_s 2^s.

Set Q(0)=Q(1)=0,Q(2)=2. Each of the six transitions has the exact slack

    e(c,b)=8b+Q(c)-2Q(c')-2c'
          =6 if(c,b)=(0,1), 2 if(c,b)=(2,1), 0 otherwise.       (1)

Multiplying by2^s and summing telescopes the Q terms: the initial state1
has Q=0 and the eventual state0 has Q=0. Consequently

    E(n)=sum_s e(c_s,b_s)2^s >=0.                              (2)

This is a genuine nonnegative boundary defect carried by a three-state
graph, not an assumed time drift. Since c_(s+1)>=b_s, and the first two
positions of an odd input contribute an additional2+4,

    2n+6 <= C(n) <= 8n.                                      (3)

Both inequalities are sharp. Upper equality means no input1 occurs from
state0 or state2; starting from1 forces the alternating word101...01
read from low to high. Thus

    E(n)=0 iff n=(4^k-1)/3 for some k>=1 iff U(n)=1.           (4)

Lower equality holds exactly when n=1 mod8 and the higher binary digits
have no adjacent1s. The first condition makes the first three inputs100
in low-to-high order, leaving state0; thereafter equality allows an
isolated1 followed by0, but not a second1 from state1.

Inverse prolongation has an exact renormalization:

    C(4n+1)=4C(n)+8,       E(4n+1)=4E(n),
    U(4n+1)=U(n),          v2(3(4n+1)+1)=v2(3n+1)+2.          (5)

The low inputs1,0 return the carry to its initial state1, proving(5)
without any asymptotic argument.

## 3. Carry defect controls an ordinary dyadic boundary offset

Let

    m=3n+1=2^t+d,      t=floor(log2 m),      0<=d<2^t.

Since m is even, d is even. Then

    E(n)>=2d.                                               (6)

**Proof with output ownership.** Output1 occurs on exactly three types:

| carry transition | slack e | role |
|---|---:|---|
| (0,1)->1 |6| restart |
| (1,0)->0 |0| exit |
| (2,1)->2 |2| continuation |

The leading output1 is necessarily an exit: after the final input1,
trailing input zeros drain any state2 to1 and then0. Every earlier exit
must be followed by a unique next restart; the state remains0 until that
restart. Assign that exit at position s to the restart at u>s. A restart
has charge6*2^u, sufficient to pay2*2^u for its own output1 plus2*2^s
for the assigned earlier exit. Restarts receive at most one assigned
exit. Continuations pay exactly2*2^s for their own output. All output1s
except the leading one are covered, proving(6).

More precisely, every restart has exactly one preceding paired exit, so

    E(n)-2d=sum_(paired exit s, restart u)(4*2^u-2*2^s).    (6a)

Each summand is strictly positive because u>s. This independent
output-pairing derivation also fixes the equality case.

If any restart occurs, the inequality is strict. In the absence of a
restart, all nonleading outputs are continuations and equality holds.
Equivalently, the finite binary expansion of n contains no pair00:

    E(n)=2d iff the binary expansion of n has no substring00. (7)

This gives the useful direction

    E(n)<=B  =>  3n+1=2^t+d with d even and 0<=d<=B/2.       (8)

The offset is nonnegative by definition; no signed shift is silently
inserted. A negative perturbation of a direct preimage generally has
large defect and is not covered by a fixed B.

## 4. The sharp hostile is repaired by a selected return

Let a_k=(4^k-1)/3 and b_k=a_k+2, k>=3. These sources have the same
dyadic height and differ in precisely the bit of weight2:

    D_(b_k)(z)=D_(a_k)(z)+z.

Nevertheless

    U(a_k)=1,       v2(3a_k+1)=2k,       E(a_k)=0,
    U(b_k)=2^(2k-1)+3, v2(3b_k+1)=1,    E(b_k)=12,
    E(U(b_k))=3*4^k+4.                                    (9)

Thus a bounded low-bit defect can generate an arbitrarily large new
defect in one growing step. Proximity to the upper extremum in(3) does
not predict a long division.

There is an all-depth version. For integers L,k>=1 set

    n_(L,k)=(2^(L+2k)+2^(L+1)-3)/3.

The numerator is divisible by3 because2^(L+2k)=2^L mod3, and its odd
quotient is positive. Its binary word is `(10)^(k-1)1^(L+1)`, so it has
no00 and

    E(n_(L,k))=2^(L+2)-4,       v2(n_(L,k)+1)=L+1.          (9a)

Its first L odd division exponents are all1, and
U^j(n)=(3/2)^j(n+1)-1>n for1<=j<=L. Yet for fixed L, E(n)/n tends
to0 as k grows. Therefore a fixed positive bound on normalized defect
E/n cannot supply any uniform return horizon. A uniform bound for the
absolute defect class E<=B, if available, must exceed L whenever
2^(L+2)-4<=B. This is an explicit logarithmic lower bound on necessary
observation depth, not a disproof of variable-depth descent.

The reset is still payable. For k>=4,

    U²(b_k)=3*2^(2k-2)+5,
    U³(b_k)=9*2^(2k-6)+1 < b_k.                            (10)

The valuation sequence is1,1,4, and the limiting multiplier is27/64.
At k=3 the last valuation is5 and U³(23)=5. For k>=5,

    E(U³(b_k))=54*2^(2k-6),

so even the descending return does not decrease the defect itself.
The successful invariant is the source-labelled affine multiplier plus
its boundary constant, not E as a scalar rank.

Equations(6) and integrality classify E<=12 completely: the possibilities
are d=0,2,6, with d=4 impossible modulo3. They yield respectively
(4^k-1)/3, (2^(2k+1)+1)/3, and (4^k+5)/3 with their valid small cases.
The d=0 family goes to1 immediately. The d=2 family descends within three
odd steps, by the lift in section5 or direct substitution; n=3 and11
are finite small controls. The d=6 family is(10); the remaining n=7 has
7->11->17->13->5. Hence E<=12 and n>1 imply descent within four steps,
and n=7 is the only member needing four.

## 5. A finite-core lifting certificate, with all quantifiers

Fix a positive even d, excluding d=1 mod3 because no exponent is then
admissible. Write d=2^a u with u odd, a>=1. Suppose a finite, checked odd
orbit of u reaches1 after J0 steps with total division exponent A0.
Append q copies of the exact1->1 transition, each of exponent2, until

    J=J0+q, A=A0+2q,     lambda=3^(J+1)/2^(a+A)<1.          (11)

This always takes finitely many additional steps because each appended
transition multiplies lambda by3/4. For every admissible source

    n=(2^t+d-1)/3,       t>a+A,

the first odd step is exactly2^(t-a)+u. The next J steps shadow the
checked core with exact affine perturbations, giving

    U^(J+1)(n)=3^J*2^(t-a-A)+1.                           (12)

For clarity, after j core steps the perturbation is
3^j*2^(t-a-A_j). The strict inequality t-a>A ensures its2-adic valuation
exceeds every next division exponent, so it cannot change that exponent.
This is the source-preserving arithmetic condition, not a population law.

Comparing(12) with n gives the exact sufficient-and-necessary inequality
within this lift:

    2^t(1-lambda)>4-d.                                    (13)

Choose an explicit cutoff K with K>a+A and(13) at K. Then(12)-(13)
prove descent for every admissible t>=K. There remain only the finite
exponents1<=t<K with d<2^t and2^t+d=1 mod3. Check first descent directly
at those exact integers. The maximum of their checked horizons and J+1
is an all-height selected-descent bound for this one offset d.

This is a **finite sufficient certificate**. It is not a proof that the
certificate succeeds for every d; core convergence and the finite
exception checks are part of its input. A shorter contractive core prefix
can replace the route to1, retaining its actual terminal constant in(13).

For B=4096, the script verifies all683 admissible positive even offsets
d<=2048. The d=0 case is(4). It verifies10,977 finite exceptions; their
sources have at most112 bits. The largest cutoff is115, the longest
large-exponent certificate has66 odd steps, and the largest first-descent
time among the finite exceptions is104. The latter is attained at
d=1406,t=53. Thus(8) proves the headline theorem for every odd n>1 with
E(n)<=4096. This is stronger than testing all n below any chosen bound,
although the covered source set has only O(log X) members below X for
each fixed B.

## 6. A constant-defect family containing27

The explicit family

    n_k=(4^k+17)/3, k>=3: 27,91,347,1371,5467,21851,...

has E(n_k)=36 for every k. Here d=18 and u=oddpart(d)=9. This is a
precise infinite arithmetic family, not an assertion that every member has
an exceptionally long orbit. The first-descent odd times are

    37 for k=3; 28 for k=4; 6 for every k>=5.              (14)

The first two and k=5 are exact finite computations. For t=2k>=12, the
first six iterates are

    2^(t-1)+9,
    3*2^(t-3)+7,
    9*2^(t-4)+11,
    27*2^(t-5)+17,
    81*2^(t-7)+13,
    243*2^(t-10)+5.

The first five exceed n_k by direct comparison, and the sixth is smaller.
The finite k=5 source347 ends at31 after six steps; its last extra
division is a small-exponent carry collision. At27 and91, the short
binary boundary causes longer transient behaviour than the eventual
large-exponent pattern. The family explains that anomaly rather than
extrapolating it to all larger members.

## 7. Why smaller cores do not finish a strong-induction proof

Indeed u<n whenever n>1: d<(3n+1)/2, and u<=d/2<(3n+1)/4<n.
Thus one might assume the smaller core converges by induction and try to
lift its proof. The missing implication is **available precision**:

    t-a > A_j.

Core convergence does not bound its necessary division budget by t-a.
If that budget is exhausted before a contractive prefix, the actual
source follows a different valuation sequence and the lift stops.

The smallest member27 of(14) is an explicit hostile. It has t=6,d=18,
a=1,u=9, hence only t-a=5 binary places. The core9 has valuations
2,1,1,2,3,..., with cumulative budgets2,3,4,6,9,... . The actual lifted
states after the initial27->41 step are31,47,71. At the next step the
core17 has division exponent2, but71 has exponent1 and maps to107.
The required5>6 has failed. The contracting core prefix would need
budget9. Merely knowing9 terminates does not certify27 through that
insufficient lift. Its first descent happens later, at the checked
37th odd step.

As B grows, the finite exception set can contain the integer one hoped
to settle. Increasing B without controlling this budget is therefore
not a strong-induction proof. The new all-height theorem is a useful
certified region, and the remaining target is now precise: an adaptive
restart must pay for the lost2-adic precision and still force a smaller
integer. It must also handle an indefinitely open forest boundary.

## 8. What the z=2 normalized jet can and cannot retain

With h=floor(log2 n), let mu(n)=J(n)/n. The occupied-bit masses2^s/n form
a probability distribution. Direct geometric sums give

    h-1 < mu(n) <= h,
    H(n)=log2 n-mu(n) in[0,2).

H is its binary Shannon entropy, since-log2(2^s/n)=log2 n-s. More
generally, every fixed centered moment

    M_r(n)=sum_s (h-s)^r b_s2^s/n

is bounded by sum_(j>=0)j^r/2^j. Thus a continuous function on the compact
closure of finitely many M_r is bounded. Adding it to positive logarithmic
height gives no eventual every-step rank, by
[THM-4507 / finite polynomial-valuation obstruction](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md),
already valid with no valuation coordinates at all. This does not exclude
singular functions on that compact boundary or exact unnormalized data.

For the pair a_k,b_k=a_k+2 in section4, the exact formula is

    M_r(b_k)=[a_k M_r(a_k)+2(h-1)^r]/(a_k+2).

Every fixed normalized moment difference tends to0, while the first
valuation is2k versus1. This explains why smooth upper-boundary shape
estimates must retain an independent low-bit obligation.

There is also an exact carry/entropy identity, not a new drift estimate:

    C(n)/(3n+1) = log2 3 + H(n)-H(U(n))
                 +log2(1+1/(3n))+mu(n)/(3n+1).             (15)

It follows by differentiating D_(3n+1) at2, and H is unchanged by removing
trailing binary zeros. Most of its main term telescopes. Rebranding(15)
as dissipated energy would supply no valuation or stopping-time control.
The nonnegative defect and finite-offset lift above use additional exact
transition/output information to produce an actual descent conclusion.

## 9. Reproduction and finite universes

Run `python -B 04-computation/experiments/nextforest_20260926_boundary.py`
and the same command with `-O`. All tests use explicit exceptions.
The independent J(3n+1)-3J(n) calculation checks the carry transducer;
131,072 odd sources below2^18 check both sharp bounds and equality
languages.398 one-bit hostile/repair pairs through800-bit scale check
the affine return identities. The27 family is independently iterated
through k=80, and320 growing-family controls cover L<=64. The finite
certificate records and their exceptions have
the SHA256 in the retained output; high representatives are directly
replayed against each symbolic affine certificate.

[Script](../../04-computation/experiments/nextforest_20260926_boundary.py)
and [retained output](nextforest_20260926_boundary.out). The analytic
proofs give all-height scope; the finite controls do not substitute for
them. The104-step constant is certificate-relative, not claimed optimal
for this entire defect class.
