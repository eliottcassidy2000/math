# Priority-aware word acceleration, saturation and costed machine interfaces

**PROVED elementary arithmetic certificates; FINITE-EXACT implementation
controls. No general FRACTRAN halting algorithm or novelty claim.**

## 1. The exact valuation representation

Reduce every positive fraction a_i/b_i. Let P be the finite set of primes in
all numerators and denominators. Every positive input has a unique form
N=S product_(p in P) p^x_p, with S coprime to P. The factor S never changes.
Write beta_i=v(b_i) and delta_i=v(a_i)-v(b_i). Instruction i is eligible
exactly when x>=beta_i coordinatewise, and the FIRST eligible instruction
updates x to x+delta_i. Thus the vector and ordered guards are lossless for
execution on the used primes; discarding S loses literal integer identity,
not the rule sequence. Factoring arbitrary inputs is unnecessary: divide only
by the fixed program primes. The supplied implementation factors its fixed
fractions by trial division, with no polynomial factoring-cost claim.

This is the standard prime-register interpretation inherited from Sessions126,
129,132 and133. The present contribution is an executable exact repeated-word
certificate that keeps the earlier-priority guards as well as denominator
availability. Accelerating flat counter loops is an established verification
idea; see Demri, Dhar and Sangnier, [*Equivalence Between Model-Checking Flat
Counter Systems and Presburger Arithmetic*](https://lsv.ens-paris-saclay.fr/Publis/PAPERS/PDF/DDS-tcs17.pdf).
Our elementary local lemma below does not assume that an entire FRACTRAN
program is flat, and imports none of that paper's global complexity conclusions.

## 2. Exact maximum repetitions of a specified instruction word

Let w=(i_0,...,i_(L-1)) be a nonempty instruction word valid for its first
occurrence at x. Put d_0=0, d_s=sum_(j<s)delta_(i_j), and D=d_L.
At phase s of repetition t, the hypothetical state is

    x+d_s+tD.                                             (1)

For any denominator guard beta_j, the set of integers t>=0 satisfying
x+d_s+tD>=beta_j is an integer interval, possibly empty or unbounded.
For a coordinate with base z and direction d, it imposes

    t>=ceil((beta-z)/d) if d>0;
    t<=floor((z-beta)/(-d)) if d<0;
    no restriction if d=0 and z>=beta;
    impossibility if d=0 and z<beta.

Intersect all coordinates and t>=0. No numerical approximation occurs.
At phase s, the word remains correct iff its own guard is true and every
j<i_s guard is false. The initial-word assumption makes its own interval
contain0 and every earlier interval avoid0. Let R be the earliest future
integer t at which ANY phase's own guard fails or ANY earlier guard becomes
true. If no such t exists, put R=infinity.

**Theorem.** Exactly R complete consecutive copies of w can execute from x.
For finite R the state after them is x+RD and the cost is exactly RL original
FRACTRAN steps. A partial next copy may execute before its first failed phase;
we do not assert that phase0 already differs. If R=infinity, the original
program executes w forever from that state.

**Proof.** For every t<R, all phase guards and their priorities hold. Induct
through the phases and repetitions to identify the actual state with (1).
At t=R the first failing phase invalidates a complete further copy, while
all prior phases remain valid. This proves maximality. The same induction
for all t proves the infinite case. The elapsed-step count counts every
original fraction, not an accelerated unit-cost substitute. QED.

Checking one supplied word uses O(L times number_of_fractions times |P|)
integer guard operations, independently of R's magnitude. Exponent and
endpoint bit lengths still matter for arithmetic cost. No claim is made
that a literal input integer with exponentially many bits was read cheaply.

### Priority is indispensable

For the program [1/10,5/2] at input8, the actual orbit is

    8 ->20 ->2 ->5, then halt.

The second instruction can fire only once before the FIRST instruction
activates. Blindly draining all three factors2 would instead produce125,
which is a wrong halting value. Own-denominator depletion is insufficient.

### Infinite certificates are scoped

If D=0, the certificate gives an integer cycle, not necessarily a least
period unless independently checked. If D is nonzero, the infinite run need
not be periodic as integers: [2/1] doubles forever. Validity for every t also
forces D>=0 coordinatewise, so a nonzero direction grows at word boundaries.
Known PRIMEGAME inputs11 and4301 reproduce the inherited2/6 cycles. This
certificate does not establish a classification of other inputs by itself;
the sibling PRIMEGAME basin theorem uses additional global budgets.

## 3. Examples and implementation boundary

`accelerate.py` supplies the reduced program, interval and word-certificate
objects, plus a lookahead runner. Lookahead only proposes words; the exact
certificate decides whether/how often they may repeat. A finite macro budget
may return `budget_exhausted`, which is UNKNOWN about subsequent behavior,
not a halting or divergence verdict. Global discovery of a suitable repetitive
word is not guaranteed.

The guarded program [15/14,7/5,11/7] on 7*2^E executes the two-instruction
word (0,1) exactly E times, then halts at 11*3^E after2E+1 fractions. The
checker certifies E=10^18 using the factorized input vector and exact integer
arithmetic; it does not literally execute2*10^18 steps or materialize2^E.

Normal and optimized checks agree. They include93,000 guard memberships,
500 random small programs, exact finite maximality, long infinite prefixes,
first-priority and cap/residue hostiles, and fixed PRIMEGAME cycle controls.
The sibling independent audit uses a separate guard/word implementation.

## 4. Why checkpoint saturation does not preserve arbitrary microsteps

For M>=0 define cap_M(z)=min(z,M+1). For any polynomial map P whose
coefficients are nonnegative integers,

    cap_M(P(cap_M(x)))=cap_M(P(x))                       (2)

coordinatewise. Addition and multiplication obey the identity: an overflowing
positive factor cannot contribute a small positive value; a zero factor
annihilates it on both sides. Induction over polynomial expressions proves
(2), and then induction over a sequence of maps proves exact finite
checkpoint abstraction. Zero/nonzero guards also respect this quotient;
more general guards must be proved compatible separately.

Session134's positive migration matrices are the linear specialization.
They implement actual graph operations, so bounded-output acceptance can be
checked on a finite cap quotient. This is NOT a theorem about arbitrary
FRACTRAN primitive instructions, which can decrement registers.

For example, exponent values12 and18 share the threshold cap min(z,5) (M=4) and
residue0 modulo6, but
after eight decrements they have different capped values4 and5. More
generally no fixed cap and modulus can preserve every countdown step.
For [1/2], the states2^r have distinct distances r to halt. Any finite
one-step deterministic quotient preserving the current halt observation
would identify two such states r<s; after r quotient steps one would have
to be both halted and not halted, a contradiction. This does not prohibit a
finite classification of EVENTUAL behavior—all these inputs eventually halt.

Likewise PRIMEGAME's three eventual basins do not encode its full future
instruction stream or elapsed time. A decidable eventual property and a
faithful step-by-step finite abstraction are different objects.

## 5. Preserve costs when changing the time scale

A word certificate is a compressed trace with repetition count R and exact
primitive cost RL. If a separately verified implementation has per-instruction
cost c_i(z)=a_i dot z+b_i on the relevant phase domain, then its cost over
a finite number R of word repetitions is exactly

    R C_0 + R(R-1) C_1/2,
    C_0=sum_s c_(i_s)(x+d_s), C_1=sum_s a_(i_s) dot D.

This follows by summing an arithmetic progression. Affine cost annotations
are an explicit hypothesis, not assumed for every arbitrary TM implementation.
The sibling machine-translation packet proves the corresponding concrete
unary costs, including the positive-matrix round sum-of-squares formula.
Faster symbolic analysis does not mean fewer transitions of the original
FRACTRAN or literal Turing machine.

## 6. Exact finite-horizon trace caps

For each used prime p put M_p=max_i beta_(i,p). Fix an integer T>=1 and
replace x_p by min(x_p,T M_p). This preserves the first T selected-instruction
or halt observations, at times0,...,T-1. It does not preserve exact valuations
or necessarily the halt observation AFTER step T.

**Proof.** Any initially unclipped coordinate remains identical as long as
selected instructions agree. An initially clipped coordinate starts at T M_p
and can decrease by at most M_p per step. Before step j<T it is therefore
at least (T-j)M_p>=M_p, so it satisfies every denominator requirement in
both runs. Coordinates with M_p=0 never constrain eligibility. Induction
preserves every guard and hence the first eligible index, or the absence of
one. QED. For [1/2], starting exponents2 and3 agree for their first two
selected rules, but only the former is halted after those two steps.

This is a horizon-dependent trace quotient, not the fixed saturation
quotient for positive checkpoint maps in Section4. No finite universal
horizon is supplied for arbitrary programs.

## 7. Bounded-period point densities are exactly rational

Fix a finite reduced FRACTRAN list and L>=1. Membership among integer states
whose least period is at most L depends only on the capped valuations

    y_p=min(x_p,L M_p).

Indeed the first L rule observations are preserved by Section6. Return to
the initial valuation vector after t<=L means that the prefix's net update
is zero, which depends on that rule word, not the original magnitudes. A halt
before return is not a periodic orbit. The unchanged spectator factor does
not affect this test.

Each cap cell has a rational natural density among positive integer inputs.
For c=L M_p, its p-factor is

    (1-1/p)p^(-y)  if y<c,
    p^(-c)        if y=c.

The last case includes every exponent at least c; when c=0 its weight is1.
Finite-prime valuation independence, by inclusion-exclusion, multiplies
these factors across p. Summing exactly the accepted cap cells proves:

**Theorem.** For every fixed program and every specified finite period bound
L, the periodic POINTS of least period<=L have an effectively computable
rational natural density. The same holds for any prescribed subset of least
periods within1,...,L.

`periodic_points.py` is the actual finite algorithm. It reproduces the
inherited PRIMEGAME least-period-two density2654208/36409529 by all59,049
cap representatives. It does not enumerate a period-six PRIMEGAME cap box,
or claim a period bound for arbitrary FRACTRAN programs.
For [3/2,2/3], the period-two POINT density is5/12; reversing the priority
changes it to5/9. The eventual cycle basin in both cases has density2/3.
Thus periodic points, eventual basins and instruction priority are distinct.

More generally any property determined solely by the used-prime valuation
vector has a natural density: truncate each exponent at a finite upper bound,
sum the accepted exact-valuation cylinders, and bound the omitted mass by
sum_p p^(-h_p-1). This existence argument needs no decision algorithm. A total
decider makes those finite sums computable. If the fixed-support predicate
has a polynomial-time binary-integer decider, caps linear in precision bits
make the density a polynomial-time computable real. The sibling PRIMEGAME
packet gives the concrete all-input basin application and explicit bounds.

The measure here is ordinary natural density of INTEGER INPUTS. It is not
uniform measure on valuation vectors, uniform time along a trajectory, or
uniform distribution on tournament classes/rational witnesses. Neither
natural-density existence nor the bounded-period algorithm is a general
halting decider.

## 8. Prime labels are a computational gauge, not a density invariant

Let pi be a permutation of primes which moves only finitely many primes, and
extend it multiplicatively to a bijection phi of positive integers and
positive rationals. Replace each ordered fraction f_i by phi(f_i). Then

    f_i N is an integer iff phi(f_i) phi(N) is an integer,
    phi(f_i N)=phi(f_i) phi(N).

Divisibility is preserved, so the FIRST eligible index is unchanged. Thus
phi conjugates the full deterministic executions, including halting and
least cycle periods. If an output is read from a designated prime exponent,
that designation must be renamed too. An injective finite prime-support
reassignment (distinct target primes) can always be extended to such a finite
prime permutation.

Natural density of INTEGER INPUTS is not invariant under this conjugacy.
The prime swap2<->3 turns [3/2,2/3] into [2/3,3/2]. Section7 computes
period-two point densities5/12 and5/9, respectively. This is the same
abstract counter dynamics with different arithmetic weights on the carrier.
Literal integer sizes and bit-cost constants also depend on the labels.
No equality of natural measures follows from computational equivalence.

This is a useful guard on cross-domain prime-group interpretations: an
association must specify an external arithmetic observable that is retained,
not just a changeable choice of register primes. It does not deny meaningful
arithmetic interpretations when that additional structure is provided.
