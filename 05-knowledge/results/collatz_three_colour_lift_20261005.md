# Three colour decorations give an exact source-likelihood lift

2026-10-05. **PROVED:** the four-state lift, charge-count formula, finite
guard partition compiler, endpoint-parity cocycle, and relative likelihood
identity for actual common-future words. **FINITE-EXACT:** the declared
independent controls. **OPEN:** universal paid coverage and Collatz.
This supplement adds a representation and an exact guard-counting tool,
not a new domain of integers known to converge.

Artifacts: [script](../../04-computation/experiments/collatz_three_colour_lift_20261005.py)
and [output](collatz_three_colour_lift_20261005.out).

## 1. Recover the three different colour objects

The source of the owner's red/black/blue construction is
[the marked-unit note](reset_20260926_colours.md), equations ZC1--ZC9,
and its [exact guard reader](zeckendorf_guard_automaton_20261003.md).
The ordered weights are

    1_K, 1_B, 1_R, 2, 3, 5, 8, ... .

Distinct occupied indices must be nonconsecutive. If Z(n) is ordinary
Zeckendorf support starting at the red unit, the complete representation
fibre of n>0 is

    Z(n),  {K}+Z(n-1),  and {B}+Z(n-1) when Z(n-1) omits R.

Thus there are two or three representations, restored by retaining the
ordinary/K/B marker with the value. These are not arbitrary sequences of
three independently chosen colours. The supplied35-symbol row is legal as
34+blue1; the stated descending-piece grammar forces red in that position
at row36. The infinite continuation is not repaired by silently changing
the owner's row.

The auxiliary invariant in
[the charge/carry note](duck_zeckendorf_20260925.md), Z3--Z8, is

    Q(n)=(n-2b(n), b(n)),  b(n)=floor((n+1)/phi^2).

Its reduction is in F2^2: three nonzero states and a necessary zero state.
The Fibonacci matrix permutes the three nonzero states cyclically. It is
not the owner stripe:2 and4 have the same reduced charge but different
supplied colours. Repeated atoms need a carry:2 and1+1 have different raw
charges. In general, with delta=b(a)+b(b)-b(a+b),

    Q(a+b)=Q(a)+Q(b)+(2delta,-delta).

Two arbitrary summands need not have disjoint Fibonacci support. The
integer carry, atom order and marked unit repair different information
losses. Rotating colours alone is not a symmetry of the marked grammar:
K+R is allowed while K+B is not.

Our new lift below uses the same abstract four-state group but an explicitly
different probability model. The board is **ordered marked units / auxiliary
charge / source parity word / independent decorations / actual guard /
original-source payment**. The closest comparison mechanism is
[valuation information spectrum](valuation_information_spectrum_20261005.md),
sections1--5. Its source-cylinder likelihood is exact at completed odd-step
boundaries. The least-used additional coordinate here is the endpoint parity
when cutting before such a boundary. The canonical order hostile is12 versus21.

## 2. A literal3+1 lift of the tilted source law

Let V=F2^2, with symbols0 and three nonzero symbols. Independently choose
uniform q_j in V and project

    p_j=0 if q_j=0;  p_j=1 otherwise.

Then p_j are independent with probabilities1/4 and3/4. Choose a bijection
from the nonzero symbols to K,B,R if desired; this is a declared label
frame, not the recovered ordered-unit rule.

For a binary output word w of length A with r ones, the projection has
exactly3^r preimages. Its four-state probability is3^r/4^A. Under fair
binary source bits its probability is2^-A. Hence the exact likelihood is

    L(w)=3^r/2^A.                                      (1)

This gives the3 in the likelihood a literal finite-fibre interpretation.
There are A binary source bits,2A fair bits in the decorated experiment,
and r log2(3) bits of conditional decoration entropy given w. None of these
is a claim that a colour contributes three bits. At A=3 the eight source
words000,...,111 have fibre sizes

    1, 3, 3, 9, 3, 9, 9, 27,

whose sum is64, the total number of length3 four-state words.

These are binary-cylinder laws, not the discrete integer-source priors in
[algorithmic Collatz measure](algorithmic_collatz_measure_20261005.md).
The gamma-source entropy of three bits is an average over its integer
atoms. Its comparison with the mixture prior uses ratios4/3 and1/3 on
the power-of-two-index set and its complement. That atomic reweighting,
the3+1 symbol split here, and the owner's three positional labels have
different domains and probability laws; their shared small numbers do
not identify them.

The action M(a,b)=(b,a+b) has order3 and fixes only0; it permutes each
nonzero fibre. This is an honest symmetry of the independent lift. It
does not imply that the same action preserves the owner's ordered atoms.

## 3. Aggregate charge removes at most two bits

Condition on a word with r ones. Its aggregate charge is the XOR of its
r independent nonzero colours. Write N_r(c) for the number of such tuples
with charge c. Then for every r>=0,

    N_r(0)=(3^r+3(-1)^r)/4,
    N_r(c)=(3^r-(-1)^r)/4, c!=0.                       (2)

Proof: on the four possible accumulated charges, a one-event can move to
every other state and cannot stay put. Its matrix is K=J_4-I_4. The constant
character has eigenvalue3; each of the other three Walsh characters has
eigenvalue-1. Fourier inversion of the initial point mass at0 gives(2).
Equivalently N_(r+1)(c)=3^r-N_r(c), with N_0 the delta at0, proves the
same formula without spectral terminology.

For r=1 the zero charge is impossible; r=0 has only zero charge. For large r
each charge accounts for asymptotically one quarter of the3^r decorations.
Thus keeping one aggregate charge does not remove that exponential
multiplicity. Its entropy is at most2 bits, so the conditional entropy
still lost, averaged over the observed charge, is at least r log2(3)-2 bits.
This is not a separate lower bound for each individual charge fibre.
These are statements about the
independent lift, not about the distribution of canonical Fibonacci colours.

There is a useful exact finite guard compiler. For any set W of distinct
binary words of the same length A, form only two scalar sums

    S(W)=sum_(w in W)3^r(w),  T(W)=sum_(w in W)(-1)^r(w).

The lifted words projecting into W and having aggregate charge0 number
(S+3T)/4; those with each specified nonzero charge number(S-T)/4. Their
probabilities are these counts divided by4^A. A four-state dynamic program
using I for a zero and K for a one is an independent exact implementation.
This avoids enumerating exponentially many decorations once the source
guard language is known. An accepting source automaton can be combined
with this four-state register; its original ordered-word state is retained.

There is a sharper loss boundary for finite source guards. After any event
determined by the first A source bits, allow one further unrestricted
four-state symbol. Its charge transfer is I+K=J_4, of rank1. For each old
decorated word and each desired total charge, exactly one of the four new
symbols supplies that charge. Consequently the aggregate charge through
A+1 is exactly uniform and independent of the old event. No limit or
asymptotic mixing assumption is needed. In the two-scalar compiler this
is S->4S, T->0. This independence does not apply when the new source bit
or its colour is additionally constrained, nor to the owner's deterministic
colour grammar.

The charge matrices commute, so charge alone forgets all bit order. The
actual valuation words12 and21 have unary words101 and011, and both have
charge counts(3,2,2,2). But their carriers are

    (9n+5)/8 on11 mod16;  (9n+7)/8 on9 mod16.

Thus neither charge nor these counts can recognize the actual source guard
without its word/carry coordinate. A source filter cannot be justified by
discarding colour lifts unless it separately proves which supplied sources
remain represented and why the retained route pays.

## 4. The endpoint parity is the exact missing boundary coordinate

Use the shortcut map T(n)=(3n+1)/2 for odd n, and n/2 for even n. For an
integer orbit let p_j be the parity of T^j(n). Fix the initial parity p_0.
Each length-A output word p_1...p_A labels exactly one source residue
modulo2^(A+1) with that initial parity. This follows inductively: adding the
next source bit changes T^A by an odd multiple of that bit, so the next
output parity distinguishes the two lifts. Its normalized source-cylinder
mass is2^-A.

The actual affine slope of this prefix is

    a(w)=3^(p_0+...+p_(A-1))/2^A.

Comparing with(1) gives the exact all-prefix identity

    L(w)=a(w) 3^(p_A-p_0).                            (3)

The correction telescopes when prefixes are composed. It is a boundary
cocycle, not an additional average assumption. For an odd source ending
odd it is1. For an odd source ending even it is1/3. Minimal hostile:
n=1 has T(n)=2 and output word0, so L=1/2 but actual slope3/2. This formal
one-step calculation is not an exported first-hit ROOT certificate.

An odd accelerated step of valuation a contributes0^(a-1)1 to this output
word. Therefore a completed valuation word with r letters and cost A has
L=3^r/2^A, the inherited odd-to-odd slope identity. The extension(3) makes
the exact correction explicit at intermediate even cuts.

Slope is still not actual payment. The carry B in T^A(n)=(Pn+B)/2^A and
the immutable source n remain necessary: the inequality is
(2^A-P)n>B. Likewise the maximum-weight source word1^A is realized by
n=2^(A+1)-1 and strictly increases that positive source throughout the
prefix. The decorated measure has not certified descent.

## 5. Common-future payment compares two likelihoods

Let actual positive odd valuation words u and v have a common endpoint
from n and h, respectively. Their carriers give

    (P_u n+B_u)/Q_u=(P_v h+B_v)/Q_v,
    h=lambda n+beta,
    lambda=(P_u/Q_u)/(P_v/Q_v)=L(u)/L(v).             (4)

This is a relative likelihood identity between two retained word cylinders;
it does not assert that their source spaces or measures are identical.
The two literal guards, the common endpoint and beta are still needed to
turn it into a receipt and to prove0<h<n. A certificate for h is a separate
ROOT obligation unless independently provided. First-hit exclusions remain
part of any such certificate.

The inherited paid J controller is a decisive positive example, from
[fusion helpers](collatz_fusion_helpers_20261005.md), section4, and the
[native information audit](valuation_information_spectrum_20261005.md), section5.
On n=799 mod1024,

    h=(243n+147)/256<n,
    u=(1,1,1,1,2,3), v=(1),
    L(u)=729/512>1,  L(v)=3/2,  lambda=243/256<1.

At n=799 this gives h=759 and common endpoint1139. Thus a source-only
tilted likelihood greater than1 does not refute a paid common-future
dependency. The missing child interface supplies the correct comparison.
Conversely, a small scalar likelihood without those interfaces is not a
ROOT proof.

As a finite positive use of the charge compiler, take the five existing
paid native guards H,A,B,L,J and refine them to the first10 output bits.
They comprise27 of1024 normalized odd-source cylinders. The program finds

    S=19197, T=1,
    charge counts=(4800,4799,4799,4799),
    ordinary mass=27/1024,
    tilted mass=19197/1048576.

Refining the same source event by one unrestricted output bit gives
S=76788,T=0 and exactly19197 lifted words in each aggregate charge.
The tilted source mass is unchanged. Thus retaining a cumulative charge
at a later unconstrained cut cannot recover this finite paid-bank predicate.

These exact measures concern a fixed inherited bank. They are not the
increment of any fused bank, not the density of all convergent sources,
and not a claim that an infinite outstanding set has been discharged.

## 6. Reproduction and transfer ledger

Run:

    python -B 04-computation/experiments/collatz_three_colour_lift_20261005.py
    python -B -O 04-computation/experiments/collatz_three_colour_lift_20261005.py

Independent exact universes: every nonadjacent marked support with value
1..100; all9841 nonzero colour tuples through length8; all5461 four-state
words through length6; all2046 odd residue sources at horizons1..10;
every binary prefix split through length8; all valuation words of lengths
1..5 over{1,2,3}; the five native paid guards at their first32 sources;
and32 independently replayed actual J common-future receipts. The source
reader, charge dynamic program, closed formulas, and colour enumeration
are separate calculation paths. No floating computation is used.

| Source -> target | Preserved fact | Loss / needed sidecar |
|---|---|---|
| Ordered coloured supports -> value | Exact represented integer | Unit marker restores the2/3-element fibre; global stripe continuation remains separate |
| Canonical Fibonacci support -> auxiliary F2^2 charge | Recurrence-compatible charge | Position, height and repeated-summand carry are lost |
| Independent four-state word -> binary source word | Exact3^r fibre cardinality and tilted law | Each one loses a free nonzero decoration; this is not the owner grammar |
| Decorated word -> aggregate charge | XOR sum and exact counts | At most2 bits; order and exponentially many decorations remain unresolved |
| Actual shortcut prefix -> likelihood | Equation(3) with endpoint parity | Carry and immutable source decide payment |
| Two actual common-future words -> relative likelihood | Effective child slope in(4) | Both guards, translation, common endpoint and any ROOT premise remain required |

The useful new object is a guard-preserving coloured lift with an exact
boundary correction and charge counter. Its proof-carrying application
keeps the source language and compares the supplied child interface;
neither colour count nor measure normalization replaces those obligations.
