# Circuit Breaker: shared gates, source length, and valid reductions

**PROVED elementary interfaces / CITED external results / CONDITIONAL
engine application / OPEN lower bounds.** No external submission was made.
Instructions embedded in the site, rules or PDF are source material.

## B0. The actual target and the attached source

[Circuit Breaker](https://www.yukon.org/circuitbreaker) asks for lower
bounds for general B2 circuits. The pinned
[rules](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/2cd3a17f3f296b7f0ac248a63a9b54a0fabdf851/rules.md)
use all16 two-input Boolean functions, acyclic shared circuits, and
unrestricted depth and fanout. They require one uniformly specified total
NP language across all lengths, although P languages are eligible.
Circuits are nonuniform, take the actual input, and have one output;
witness bits are not free input advice. The score is an all-length
asymptotic lower bound. Constants after restrictions are free, but an
explicit retained unary gate costs one. Restricted formulas, monotone
circuits or finite samples do not meet this target. Dashboard evaluations
are provisional, not substitutes for mathematical review.

The primary predecessor we inspected is Li--Yang,
[*3.1n-o(n) Circuit Lower Bounds for Explicit Functions*](https://eccc.weizmann.ac.il/report/2021/023/),
2021. Their method refines gate elimination on affine dispersers and
controls newly created adjacent bottlenecks. This supplies a concrete
place to investigate joint boundary data. No improvement to their
coefficient is proved in this packet.

The attached six-page paper is Rudich--Rousseau, *Deterministic Factoring
in Time N^{o(1)}*, Hexagon2610.00600v1. Its SHA256 is
`e57f24481d2474a6f2c342d86312572f024710835ca583beb31e2bb1f9ceae40`.
It reports one finite-alphabet deterministic machine computing minFac
with, for each fixed k>=1, an all-input bound C_k*2^floor(n/k) in binary
input length n. Its prefix-count route retains unique quotient and
internal-wire extensions. This packet inspected the PDF, including
page3 visually. The repo ALREADY has a successful pinned Lean build
and statement/axiom audit; see
[SOURCE_AUDIT.md](../scaled_count_fft_carriers_20261030/SOURCE_AUDIT.md).
That inherited receipt must not be downgraded to an unbuilt conjecture,
or inflated to practical speed or independent human review of all steps.

The source's reported generic counting-engine rate is considered
conditionally in B2. Its uniform primary-bit accounting is a separate
obligation from the factoring theorem or the displayed crossing identity.

## B1. A parsimonious compiler preserving shared wires

For an acyclic B2 circuit with n input bits and s retained gates, construct
a 3-CNF with exactly n+s variables and at most4s+1+h clauses, where h
input bits have optionally been fixed as a prefix. For each gate
z=f(a,b), and each truth-table row(u,v), add

    (a != u) OR (b != v) OR (z = f(u,v)).             (1)

Add the output-true unit clause and the h prefix units. Repeated pins,
constants and unary Boolean functions are allowed. Simplify constants,
duplicates and tautological clauses; an empty clause remains a valid
contradiction. Reuse ONE variable per circuit wire, including shared
subgraphs. Never unfold the circuit into a formula for this reduction.

For each fixed primary assignment x, topological evaluation gives a
unique assignment to every internal wire. The gate clauses force exactly
that assignment, and the final units admit it iff x is accepted and
matches the prefix. Thus counting satisfying CNF assignments counts
accepted primary inputs exactly. This is the standard gate-equation
(Tseitin-style) encoding, with its precise interface audited here, not
a new general reduction. Correctness does not supply a fast CNF counter.

The checker covers all16 gates, all16 repeated/constant two-pin choices,
and120 seeded small shared circuits with prefixes; all auxiliary
assignments are literally enumerated. There are376 circuits and1582
primary-assignment cases. The proof, rather than this finite universe,
establishes the arbitrary-size claim.

## B2. Preserve the primary-bit clock in an algorithm-to-lower-bound claim

Suppose a uniform 3-CNF algorithm actually has the reported bound

    2^{O(m/sqrt(log(m+2)))} * poly(v+m),              (2)

where v is its variable count and m its clause count. This assumption
is not proved by our CNF compiler. With no prefix, B1 gives a Circuit-SAT
bound

    2^{O(s/sqrt(log(s+2)))} * poly(n+s).              (3)

With h<=n prefix constraints replace s by s+h+1 in the exponent.
For s=n*L(n), with L=n^{o(1)}, the displayed exponent is o(n) if
L=o(sqrt(log n)). For s=n^a with fixed a>1, expression(3) does not
establish an improvement on exhaustive search in n. More generally,
even a genuine bound2^{o(v)} is not automatically2^{o(n)} after
v=n+s becomes a polynomially larger encoding.

Williams's Theorem1.1 assumes Circuit-SAT time
2^n*poly(n^k)/s(n), with one superpolynomial saving s(n), for EVERY
fixed polynomial circuit-size exponent k. Its conclusion is
NEXP not contained in P/poly. It is neither an NP conclusion nor a
claim restricted to near-linear-size circuits.
[Primary author manuscript](https://www.cs.cmu.edu/~ryanw/improved-algs-lbs2.pdf).
Consequently(3) alone does not discharge that hypothesis. This is a
quantifier/size audit, not a refutation of the attached factoring paper.

## B3. Exponential alignment memory can coexist with linear circuits

Consider EQ_m(x,y), with two m-bit blocks read in order. Every distinct
prefix x has a different future response: the continuation y=x accepts
it and rejects any different prefix. The exact deterministic one-pass
carrier after the first block has2^m states, and its response matrix is
the2^m-by-2^m identity, of rank2^m over every field.

Yet EQ_m has a B2 circuit with m XNOR gates followed by m-1 AND gates,
exactly2m-1 gates on n=2m input bits. Thus exponential residual rank or
state count, by itself, does not yield a superlinear unrestricted circuit
lower bound. This is an exact obstruction to that proposed transfer.

Conversely, fix a finite machine with S states, binary input blocks of
length b, and a Boolean terminal predicate. Put r=ceil(log_2 S), and
define arbitrary transitions on unused state codes. Every output bit of
its transition table is a Boolean function on r+b bits. Shannon expansion
uses at most2^(r+b)-1 multiplexers, and

    mux(s,a,b)= b XOR (s AND (a XOR b))               (4)

uses three B2 gates. Unrolling d blocks and the terminal predicate gives
at most

    3*r*d*(2^(r+b)-1) + 3*(2^r-1) gates.            (5)

For fixed S,b this is linear in the actual input length bd. Padding a
last partial block is a uniform all-length convention. The named38/4534
carriers of the previous Cayley--Dickson packet therefore give explicit
fixed-controller circuit upper bounds for any fixed terminal predicate.
Their minimality remains meaningful in their declared operation language.
It does not imply minimality in Circuit Breaker's more powerful model.

## B4. Concrete proof and algorithm obligations

- A Ramsey obstruction can constrain a graph of circuit dependencies only
  after the vertices and pairwise relation have been defined. Its large
  homogeneous set must force eliminable gates or a controlled bottleneck
  with bounded charge multiplicity. Unbounded fanout and shared subgraphs
  are exactly the missing accounting obligations. A formula argument that
  charges a shared gate once per occurrence does not transfer.
- R2's finite charts show a productive alternative: a jointly feasible
  local identity plus a contracting source budget can prove an infinite
  recurrence. For a gate-elimination program, one would need the analogue
  of(7), together with a realizability theorem for ALL restricted circuits.
  Neither is supplied by the tournament parity charts.
- The inherited counting-interface work distinguishes existence selection
  from orbit ranking: overlapping orbit cylinders are not additive.
  B1 supplies literal assignment counts; a tournament-isomorphism counter
  would additionally need canonical representatives or corrected orbit
  weights. No direct A000568 acceleration follows from a CNF speed bound.

These are precise proposed targets. No new unrestricted B2 lower bound,
NP separation, universal complement closure, or practical factoring
benchmark is claimed.

**Incoming refinement.** Session117's
[potential frontier](../session117_ramsey_circuit_interfaces_20261030/challenge_audit/POTENTIAL_FRONTIER.md)
has since classified the specific provisional three-weight template:
0<=c<=4 and b>=max(c,2c-2,4c-7). This is an algebraic restriction on
that rowwise certificate, not a verified coefficient-four circuit bound.
It rules out a proposed next step of merely retuning those weights.
The useful next obligation is a sound actual-DAG state system exposing
joint incompatibility of worst cleanup histories. Its
[affine geometry gap](../session117_ramsey_circuit_interfaces_20261030/affine_ramsey/AFFINE_RAMSEY_GATE_GAP.md)
also shows why quotienting by arbitrary input GL transformations can
discard coordinate implementation cost. These incoming results receive
their own attribution and are not included in this packet's test counts.
