# PRIMEGAME return arithmetic contains the earlier seven-prime geometry

**PROVED elementary bridge from an inherited exact kernel.** This concerns
formal arithmetic return multisets, not PRIMEGAME's ordered executions.
Classical Fano and Hamming objects are not claimed as new discoveries.

Use Conway's fourteen fractions A,...,N in their original priority order,
as recorded in Session129 `partition/README.md`. The 10-by-14 integer matrix
V has columns equal to their prime-valuation increments, on primes
2,3,5,7,11,13,17,19,23,29. The inherited integral kernel is

    z=(g+2t+m,g+t+m,t,g,g+t+2m,g+t+2m,g,t,
       0,j,t+j,t+m,m,0).                                (1)

The nonnegative formal returns are exactly g,j,t,m>=0. Their four primitive
multisets are ABDEFG, JK, A^2BCEFHKL, ABE^2F^2LM. This rank10/kernel4
description, the user's 2+4+4+4 block height, and the complete priority-cycle
classification are inherited; they are not new in this packet.

## The seven observations and their exact arithmetic partners

Ignore the independent JK circulation parameter j, and observe the counts
of instructions D,C,E,M,A,L,B in a formal return. Their integer coefficient
columns on (g,t,m) and their reductions modulo2 are:

| instruction | integer column | binary column (LSB first) | partner prime p | square class of p+1 |
|---|---|---:|---:|---:|
|D|(1,0,0)|1|7|2|
|C|(0,1,0)|2|2|3|
|E|(1,1,2)|3|5|6|
|M|(0,0,1)|4|19|5|
|A|(1,2,1)|5|89|10|
|L|(0,1,1)|6|59|15|
|B|(1,1,1)|7|29|30|

The binary columns are every nonzero vector of F2^3, exactly once. Thus
the observed parity vectors form the binary simplex[7,3,4] code. Relations
among these seven observations form its dual, the Hamming[7,4,3] code.
These are different codes: circulation parities are simplex words; subsets
of observations giving a universal parity check are Hamming words.

The right-hand prime bank is inherited from the October28 arithmetic-code
work. Its square-class coordinates are(2,3,5). The displayed assignment
therefore gives a literal connection, not just two instances of the number7.

**Theorem.** For a subset S of these seven instruction labels, form the
squarefree integer n_S by multiplying its assigned partner primes. Then

    sigma(n_S) is a square
        iff sum_(i in S) z_i is even for EVERY formal arithmetic return z.
                                                               (2)

**Proof.** By coprimality, sigma(n_S)=product_(i in S)(p_i+1). It is square
iff the sum of the displayed binary columns is0. By(1), this is equivalent
to the selected count sum being even for every(g,t,m) in F2^3. Every such
parameter parity has a nonnegative formal-return representative, so the
universal quantifier is equivalent in both directions. QED.

There are exactly16 accepted subsets, with weight enumerator

    1+7z^3+7z^4+z^7.

The seven weight-three subsets are Fano lines. There are also seven minimal
four-point circuits; the earlier repo correction that lines do not exhaust
Fano circuits is retained. An ordinary graph cycle matroid cannot encode
this configuration: a simple rank-three graphic matroid has at most six
elements (the maximum is K4), whereas this simple rank-three configuration
has seven. The faithful carrier is binary linear incidence.

## Exact loss boundary: arithmetic return is not executable return

Session129 proves that actual periodic PRIMEGAME words are repetitions of
JK or of the six-cycle GBDEFA, with the specified source guards. In(1),
this means only j and g survive; the t and m return generators cannot be
ordered into priority-respecting periodic executions.

Consequently requiring the selected count sum to be even on **actual
periodic executions** imposes only the g-row condition, accepting64 subsets
rather than16. For example S={C} has even count on every actual periodic
orbit (C never appears), but its partner prime2 has sigma(2)=3, not a square.
The formal t-generator has z_C=1 and detects it. This is an exact hostile
to replacing formal returns by legal ones in(2).

Reduction modulo2 also loses magnitude. For example the Fano parity lines
D+M+A and D+C+E have the integer defects

    column(A)-column(D)-column(M)=2 column(C),
    column(E)-column(D)-column(C)=2 column(M).

The lost even coordinates are not zero debts. All seven line determinants
and their actual integer columns are stored in primegame_code.json.
Thus there are two separate losses: parity forgets even multiplicity;
commutative circulation forgets order and priority. Restoring only one
does not restore an executable word or its clock.

## Connection contract and relevance to the spectrum target

Source: the exact PRIMEGAME exponent kernel and its nonnegative return monoid.
Target: the October28 square-sigma bank with marked prime assignments.
Map: the seven displayed instruction-coordinate functionals, reduced mod2.
Preserved predicate: universal formal-return parity, exactly(2).
Lost data: integer multiplicities, word order, first-enabled guards and time.
Sidecar: the integral lift plus the actual instruction word and its source.

For tournament construction this supplies a precise design lesson: relations
among available operations can compress observations, but a candidate
arithmetic solution still needs an executable source word. The binary-port
classification in BINARY_PORTS.md demonstrates this boundary on actual
tournaments. Neither bridge proves an additional Hamiltonian count or the
universal spectrum theorem.

Sources: `session129_primegame_tournament_partition_20261031/partition/README.md`,
`free_bits_signed_incidence_20261028/THEORY.md`, and
`session85_arithmetic_codes_free_bits_20261028/ARITHMETIC_FREE_BITS.md`.
Bounded prior-art searches inspected the current PRIMEGAME and partition
packets for this code; no external novelty priority is asserted.

## Incoming incidence theorem makes the quotient canonical

The concurrent [entry-incidence analysis](../primegame_path_orbits_20261102/CONNECTIONS.md)
identifies the integer factorization

    Z^14 --B--> Z^19 --V--> Z^10,       Delta=VB.

Here B sends each instruction to numerator-entry minus denominator-entry,
and V takes the ten prime valuations. Its entry graph has six components,
rank(B)=13, and ker_Z(B)=Z*(JK). The exact kernel(1) of Delta therefore gives

    ker_Z(Delta)/ker_Z(B) ~= Z^3,

with coordinates(g,t,m). Equivalently this is ker(V) intersect im(B), by the
map induced by B. This identifies the three extra arithmetic circulations
created by factorization; it justifies quotienting out j intrinsically.
The entry-incidence factorization and its ranks are credited to the incoming
packet; its combination with the seven-probe code is made explicit here.

J and K individually do not descend to this quotient. Among the remaining
twelve coordinate observations, I,N are zero and the exact nonzero forms
fall into precisely seven classes:

    D=G, C=H, E=F, M, A, L, B.

Thus the seven columns in the table exhaust the distinct nonzero surviving
instruction observations; they were not selected from a larger collection
of unrelated binary vectors. Their reduction modulo2 is the Fano plane.
Four Fano-line relations lift exactly over the integers; the other three
have determinant +2 or -2. The latter become lines only in characteristic2.
The rank-three integer carrier and the binary incidence carrier are distinct.

There is also an exact Galois reading. Put K=Q(sqrt(2),sqrt(3),sqrt(5)).
A formal return with parameter parities u=(g,t,m) specifies the Galois
sign-change automorphism of K with those three signs. For every S in(2),

    sigma_u(sqrt(sigma(n_S)))
       = (-1)^(sum_(i in S) z_i) sqrt(sigma(n_S)).

Here the outer sigma_u is a field automorphism and the inner sigma is the
sum-of-divisors function; their types differ. The equality follows directly
from the displayed square classes. It transports a character evaluation,
not a PRIMEGAME execution. The incoming octonion sign-character construction
in Session143 supplies a further action-level comparison on the same binary
labels, but its nonassociative multiplication is not transported by this map.
