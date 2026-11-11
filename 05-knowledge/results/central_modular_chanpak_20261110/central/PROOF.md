# Central modular blocks as an exact Collatz edge codec

**PROVED, scoped elementary statements C1--C3 below.** C1 is inherited; C2
is an explicit decoding and relabelling of the user's observation; C3 records
the finite-word certificate needed to iterate it. No claim of priority is
made. Positive-integer Collatz convergence and the tournament Hamiltonian-path
count classification remain **OPEN**. The computation is **FINITE-EXACT** and
checks the stated finite universes; the proofs supply the universal claims.

## Inheritance and live concepts

The closest proved mechanism is the central-block identity in
[tiling_modular_atoms_20261003.md, Section 3](../../tiling_modular_atoms_20261003.md)
and **THM-4792**, [two-sheets-mod-odd-m-merged-centre-reflections-unit-shift-and-alternation-classification](../../../../01-canon/theorems/THM-4792-two-sheets-mod-odd-m-merged-centre-reflections-unit-shift-and-alternation-classification.md), part (a).
Those are fixed-modulus statements; C2 below changes the modulus along an
integer orbit. The hostile example is modulus 9: its central values 2 and 7
do not constitute the Collatz edge out of 2. The corrected near miss is
discarding the endpoint-parity bit from a finite source address. The least-used
sidecar here is the ordered affine carry, rather than the number of odd and
even steps alone.

The live board is: central unit blocks; parity and sign gauges; moving modulus;
source cylinders; ordered affine carry; prime/composite colour. The last three
are the information retained when a visually simple table pattern is iterated.

## C1. What the central four entries actually are

For odd m >= 3, index the nonzero part of the multiplication table by
1,...,m-1, and use least nonnegative residues. Put k=(m-1)/2 and let d be the
representative of 4^(-1) modulo m in {1,...,m-1}. The rows and columns k,k+1
give

    [ d      m-d ]
    [ m-d    d   ].                                      (1)

Indeed, k=-1/2 and k+1=1/2 modulo m, so their products are +/-1/4.
Consequently each of the two rows sums to m, as does each column. The
phrase “two copies of m” therefore describes the two row sums, not two
additional entries inside a four-entry block.

Equivalently, for a positive integer x:

    m=4x+1:  [ 3x+1    x   ]       m=4x-1:  [ x      3x-1 ]
             [ x       3x+1]                [ 3x-1  x    ]. (2)

This gives precisely the user's two formulas for m modulo 4. Every central
row label and every entry in (1) is a unit modulo m, even when m is composite.
For example m=9 gives [[7,2],[2,7]], while m=25 gives [[19,6],[6,19]].
Thus this patch by itself has no prime/composite test: composite zero divisors
live elsewhere in the table. The ordinary integer labels still reveal m;
the claim is about the unit/invertibility pattern of this central patch.

The nonzero-index convention matters: a full m by m table with a zero row
has odd side length and one middle entry, not this same central 2 by 2 block.
Modulus 1 and x=0 are excluded. Moduli 3 and 5 are valid base cases.

## C2. Four modulo-eight lanes encode all edges, with a decoder

Fix epsilon in {+1,-1}, and define the ordinary, unaccelerated positive-integer
map

    T_epsilon(x) = 3x+epsilon  if x is odd,
                   x/2        if x is even.              (3)

Let m=4x+epsilon and name the two central values x and y=3x+epsilon.
The following is an exact codec for the edges of (3).

| Modulus | Source x | Sign | Raw central values | Decoded edge |
|---|---|---|---|---|
| m=5 mod 8 | odd | + | x, y=3x+1 | x -> y |
| m=1 mod 8 | even | + | x, y=3x+1 | x -> (y-1)/6 = x/2 |
| m=3 mod 8 | odd | - | x, y=3x-1 | x -> y |
| m=7 mod 8 | even | - | x, y=3x-1 | x -> (y+1)/6 = x/2 |

The second row starts at m=9, not m=1. The fourth starts at m=7.
On an odd-source row the corresponding division by 6 is not an integer,
because y-epsilon=3x is odd. Thus the modulo-eight class exactly supplies
the branch guard and the sign specifies the subtraction or addition before
division by 6. The two rows per sign partition source integers by parity;
this is the precise meaning of the two “halves” of the edges.

These are decoded edges. For m=9 the raw pair 2 -> 7 is not a Collatz edge:
T_+(2)=1, and (7-1)/6=1. Literal reversal 7 -> 2 is also not a Collatz edge,
since T_+(7)=22. A reversal of the table-reading convention must therefore be
distinguished from reversing a directed edge in the Collatz graph.

There is an exact way to iterate without leaving the modulus labels. The
bijection phi_epsilon(x)=4x+epsilon maps positive integers to the domain
D_epsilon={4x+epsilon: x>=1}. Define on that domain

    M_epsilon(m) = 3m+2epsilon  if m=4+epsilon mod 8,
                   (m+epsilon)/2 if m=epsilon mod 8.      (4)

Then

    M_epsilon(phi_epsilon(x)) = phi_epsilon(T_epsilon(x)). (5)

Proof: substitute m=4x+epsilon in either branch. In the odd case the next
modulus is 4(3x+epsilon)+epsilon=3m+2epsilon; in the even case it is
2x+epsilon=(m+epsilon)/2. Both remain in D_epsilon. This proves a graph
isomorphism between the full functional graphs on their specified domains.

The plus cycle 1 -> 4 -> 2 -> 1 becomes 5 -> 17 -> 9 -> 5. It crosses from
prime moduli to a composite modulus, so restricting to primes destroys this
closed graph. The minus cycle 1 -> 2 -> 1 becomes 3 -> 7 -> 3. These are
boundary controls, not a classification of all cycles. Relabelling the
Collatz graph does not prove that every orbit reaches one of these cycles.

This moving-modulus graph is also distinct from THM-4792's alternating
fixed-modulus actions on +/-1/4. A change of representation is not an extra
dynamical contraction.

## C3. A finite automaton generates legal words; the carry certifies them

For a finite ordinary-clock word w in O,E, let a(w) and b(w) count O and E,
and let delta(w)=1 if its final letter is O, otherwise 0 (including the empty
word). An O records the odd branch of (3), an E the halving branch.

Every realized word avoids OO, because 3x+epsilon is even for odd x. Conversely,
every finite word avoiding OO is realized by infinitely many positive integers.
Its exact source set is a single residue class

    x = r_w modulo 2^(b(w)+delta(w)).                     (6)

There is also an exact affine endpoint certificate

    T_epsilon^length(w)(x)
       = [3^a(w) x + epsilon B_w] / 2^b(w),              (7)

where B_empty=0 and appending a letter updates

    O: (a,b,B) -> (a+1,b,3B+2^b),
    E: (a,b,B) -> (a,b+1,B).                             (8)

Constructive proof of (6): start with the sole class modulo 1. After O the
next E is forced, and appending it does not change the modulus in (6).
After E, or at the empty word, the two lifts r and r+2^b modulo 2^(b+1)
have opposite output parities: equation (7) changes the current output by
the odd integer 3^a. Choose the unique lift with the parity required by the
next letter. Update (8). This gives both a source-class construction and
an inductive proof that every class is nonempty. A positive representative
plus arbitrary positive multiples of its modulus supplies infinitely many
positive realizations. Formula (7) follows by composition of the two affine
branches. The construction works for both signs.

The legality automaton has just two states (last O, or not last O), and its
length-r word count is F_(r+2), with F_0=0,F_1=1. This is the inherited
golden-mean language of **THM-4528**,
[collatz-parity-in-base-phi-golden-beta-map-and-holonomy](../../../../01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md).
The two signs share the language, although their integer source classes differ.
This ordinary-clock no-OO language must not be confused with shortcut Collatz
parity words, in which every finite binary string is legal.

Order remains load-bearing. OEE and EOE both have (a,b)=(1,2), but their
carries are 1 and 2. For epsilon=+1 their maps are (3x+1)/4 and (3x+2)/4,
with source classes x=1 mod 4 and x=2 mod 4, respectively. Counts of operations
alone forget both the endpoint and the admissible source.

Likewise the final-O bit in (6) cannot be discarded. The word O has b=0,
but it is realized only by odd sources, not by all integers modulo 1.
To reconstruct a particular integer rather than an address, retain the
unbounded integer lift x=r_w+k 2^(b+delta). Neither the finite automaton nor
a finite residue observation supplies an integer descent theorem.

The related prior result
[odd_eighths_uniform_families_20261030/COLLATZ.md, E1--E2](../../odd_eighths_uniform_families_20261030/COLLATZ.md)
keeps the same original-source phase through its conditional Gram observer;
its memory grows with dyadic depth. Its clock is the shortcut clock, so its
all-binary-word polynomial is not substituted for this no-OO language. The
shared principle is that legal schedules are a small quotient while source
alignment is additional information. The analogy does not import a theorem
about one observer into another.

## Connection audit and constructive use

| Source -> target | Exact map / preserved predicate | Lost information / required sidecar |
|---|---|---|
| Odd-modulus table -> signed central pair | m=4x+epsilon -> (x,3x+epsilon); preserves equation (2) | Full table discarded; retain sign, index convention and ordinary representatives |
| Central pair -> integer edge | Four-row decoder; preserves exactly (3) | Raw pair alone lacks the even-branch decoder; retain m mod 8 |
| Integer graph -> modulus graph | phi_epsilon; all edges preserved bijectively | Nothing if full labels retained; prime-only restriction loses edges |
| Integer orbit -> finite schedule | Branch word; preserves legality | Endpoint, source and integer magnitude need (r, exponent, a,b,B, integer lift) |
| Table positions -> tournament tiling | Prior tiling_modular_atoms position map | Positions are not tournament incidence or Hamiltonian-path count; no count-preserving map is supplied here |

The reusable constructive device is a guarded affine program together with
its residue-cylinder certificate. It generates every finite legal possibility
procedurally, without tabulating enormous orbits. In a tournament construction
the analogous certificate would have to retain named ports, joint attachment
data and all phase/control coordinates used by a return rule. The table codec
does not supply that return rule. A genuine bridge must show that every claimed
affine operation is realized on the retained tournament state and preserves
the next operation's guard. Numerical similarity of an affine coefficient is
insufficient.

## Reproduction and exact finite scope

Run `python3 05-knowledge/results/central_modular_chanpak_20261110/central/check.py`
from the repository root. The checker uses independent direct products for
every odd m from 3 through 501; both signs and all x from 1 through 1000 for
the edge codec and graph conjugacy; and every ordinary parity word obtained
from all representatives 1,...,2^r for each r from 0 through 12, for both
signs. It compares direct iteration with the independent recursive cylinder
construction and checks the full source fibres, affine endpoints, Fibonacci
counts, named hostile examples and prime-colour cycle. Its ordinary and `-O`
outputs must agree. The JSON evidence hashes the checker and inherited source
files. This verifies the bounded implementation, not Collatz convergence or
the tournament spectrum conjecture.
