# Least residue clocks for priority arithmetic

Status: **PROVED** elementary statements below; **FINITE-EXACT** replay in
`verify_clocks.py`. This is a finite-query classification, not a halting
algorithm for arbitrary FRACTRAN programs.

## The observation determines the clock

Let P=(a_1/b_1,...,a_s/b_s) be a nonempty ordered list of reduced positive
fractions. On positive integer N, execute the first i for which b_i divides
N; if none exists, halt. Instruction i is **active** precisely when no
earlier denominator divides b_i. Indeed N=b_i then selects i; conversely
an earlier divisor of b_i always preempts it. Preserve original instruction
indices when deleting inactive rules.

Put D=lcm of the active denominators and G=gcd of the active numerators.
Then gcd(D,G)=1: each prime dividing D divides an active denominator and
therefore does not divide that rule's numerator.

For m>=1, observe the pair (selected index, output modulo m), or a single
HALT symbol. The halted input itself is **not** an output observation.

**Theorem 1.** The least positive period of this observation on positive
integers is

    L_m=lcm_active_i [m b_i/gcd(m,a_i)] = Dm/gcd(m,G).       (1)

If L is any period, test N=b_i for an active i. The same rule must be
selected at N+L, so b_i divides L. Equality of output residues then forces
m to divide a_i L/b_i. Thus every displayed factor divides L. Conversely
their least common multiple preserves all active divisibility tests and
every selected output residue; inactive rules cannot affect selection.

For the second equality, at prime p with r=v_p(m), the exponent on the
left is max_i[v_p(b_i)+max(0,r-v_p(a_i))]. If p divides D, a denominator
of maximal p-valuation has numerator valuation zero, and this maximum is
r+v_p(D). Otherwise it is max(0,r-v_p(G)). These are exactly the exponents
of Dm/gcd(m,G). This proves the formula at **every modulus**, including
nonunits, prime powers, and CRT products.

The least clock for the rule index alone is D (take m=1). The numerator
gcd can remove only output precision, never a prime required by division.

**Compatible CRT charts.** For every m,n>=1,

    L_lcm(m,n)=lcm(L_m,L_n),  L_gcd(m,n)=gcd(L_m,L_n).

This follows from the primewise formula in the proof: each exponent is a
fixed increasing function of the corresponding output exponent. In
particular, if gcd(m,n)=1, the two input-residue charts must agree modulo
D; their fiber product has L_m L_n/D=L_(mn) states. They cannot be chosen
independently. For [1/2,2/3] and output moduli5,7, the charts have30 and42
residues and share the mod6 selector; exactly210 of the1260 chart pairs
are compatible. These are minimal **congruence** carriers, not a claim
that the observation itself has that many distinct values.

**Corollary 2.** A finite congruence carrier N mod M determines both the
selected rule and the next state mod M iff L_M divides M. Such an M
exists iff D=1. If D>1, any p dividing D requires v_p(M)+v_p(D) input
digits, strictly more than M contains. If D=1, the first rule is an integer
multiplier and all later rules are inactive, so every M works.

This excludes fixed congruence carriers with this exact observation. It
does not exclude finite noncongruence quotients for weaker observations.

For a horizon of h steps, keeping the branch sequence and only the final
residue mod m (with HALT absorbing and its input erased), iterating (1)
gives the sufficient input modulus

    L^(h)(m)=D^h m/gcd(m,G^h).                              (2)

The formula follows prime by prime from gcd(D,G)=1. It is generally not
minimal: priority restricts which denominator sequences actually occur.

## A complete small family with a decreasing integer rank

Consider P=(1/2,2/3). Write x=v_2(N), y=v_3(N). Its exponent transition is

    (x,y) -> (x-1,y),       x>0,          rule 1;
    (0,y) -> (1,y-1),       y>0,          rule 2;
    (0,0) -> HALT.

The natural number x+2y decreases by exactly one on every legal step.
Therefore every input halts after exactly x+2y steps, with output
N/(2^x 3^y). This is an actual well-founded certificate. Its crucial
hypothesis is the displayed transition rule; importing it into Collatz
would require a new native rank that decreases on Collatz's transitions.

**Theorem 3.** The first h rule observations, padding with HALT after
termination, have least period

    M_h=2^h 3^ceil(h/2),                                   (3)

and exactly

    C_h=1+h+floor((h+1)^2/4)                               (4)

distinct response classes, for all integers h>=0.

Sufficiency of (3): the trace is 1^x followed by (21)^y and then HALT.
At most h halvings can be observed, and at most ceil(h/2) occurrences of
rule 2 can be observed. The input residue at M_h determines the required
truncated valuations. For necessity, let L be a positive trace period.
If v_2(L)=t<h, inputs N=2^t and N+L first differ after t halvings: the
former halts and the latter executes rule 1. Thus 2^h divides L.
If v_3(L)=t<ceil(h/2), choose odd k in {1,5} such that
k+L/3^t is divisible by 3, and set N=k3^t. Both inputs are odd, and
each performs t copies of (21). Then N halts while N+L executes rule 2,
which is visible because 2t+1<=h. Hence 3^ceil(h/2) divides L.

For (4), all x>=h share the trace 1^h. For each x=0,...,h-1, there are
ceil((h-x)/2)+1 distinct tails: y=0,...,ceil((h-x)/2), with the last case
absorbing every larger y. These traces differ either at their first
non-1 symbol or their first HALT symbol. Summing proves (4).

Thus exponential-sized congruence clocks can compress to quadratically
many response classes when the query is a fixed-length instruction trace.
At h=2, (2) with m=1 gives 36 while (3) gives 12, and there are only five
response classes. No finite carrier records every complete trace: the
halting times x+2y are unbounded. A finite output is not a fixed-size
representation of all the information that future observations may ask for.

## Relation to the existing corpus and machines

The inherited [Session135 bounded-horizon exponent caps](../session135_fractran_tm_structure_20261101/SYNTHESIS.md)
retain priority and bound the precision consumed by a finite computation.
Equations (1)--(4) refine that mechanism for the declared residue and trace
queries. They are compatible with [THM-4852's selector-automaton route](../primegame_variant_selector_automaton_metagraph_20261031.md),
but a digit encoding by itself is not the operation-closed residue quotient
excluded by Corollary 2.

The **CITED** general computational boundary is Endrullis, Grabmayer and
Hendriks, [Complexity of Fractran and Productivity](https://arxiv.org/abs/0903.4366v3):
FRACTRAN is Turing complete, and uniform halting over positive integers is
Pi^0_2-complete. This does not prevent the explicit decreasing-rank proof
for the particular two-rule program above.

Transfer audit: source = ordered arithmetic program; target = residue or
trace carrier; map = input modulo (1), or its exact trace class; preserved
predicate = the stated finite query; lost information = unrequested later
divisibility/output data; required sidecar = rule order and horizon. The
cheapest hostile is h=2 for [1/2,2/3], followed by increasing the horizon.
