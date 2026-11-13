# The Collatz clock and the ordinal-triangle guard

**PROVED scoped elementary consequences, with inherited exclusions and
mutation theorem stated explicitly. FINITE-EXACT controls in check.json.
Chan–Pak and Collatz convergence remain OPEN.** No literature-priority claim.

The useful bridge is an alignment of a clock and a guard: the central-table
Collatz map lands at a modulus divisible by3 after every odd step followed
by one halving. An extreme ordinal triangle permits a count-preserving
source change at a tournament count divisible by3. The second statement
has a stronger realizability condition than the first.

## Inheritance, hostile and concept board

Closest mechanisms are [central C2–C3](../../central_modular_chanpak_20261110/central/PROOF.md),
the [four-return theorem](../../triples_ports_sieve_20261112/ports/FOUR_RETURNS.md),
and the [response-content invariant](../../odd_tables_signed_carries_20261111/tournament/CONTENT.md).
The inherited H135 reset has response content27, blocking both nonzero
ternary digits in its continuation grammar, while its H133 partner stays
primitive. That established example already shows why neutral source
improvement and arithmetic availability are distinct.

The canonical hostile here is an attained count63 with no ordinal C3
presentation at any order. The corrected near miss is replacing physical
factor availability by divisibility. The least-used sidecar is the
attainable count of the residual graph after removing a factor.
Live concepts: **clock / residue / actual factor / neutral source fibre /
ternary digit / integer realization**.

## A1. Exact clock at which the modulus regains its factor3

Use the inherited ordinary signed map T_epsilon, with odd branch
3x+epsilon and even branch x/2, for epsilon=+1 or-1. Set

    m=4x+epsilon.

At an odd source x the modulus steps to3m+2epsilon. Follow this by j legal
halvings, where1<=j<=v2(3x+epsilon). The exact endpoint modulus is

    m_j = [3m+(2^j+1)epsilon]/2^j.                    (1)

Proof: the integer endpoint is(3x+epsilon)/2^j, so its modulus is
4(3x+epsilon)/2^j+epsilon; substitute4x=m-epsilon. QED.

Because2 is a unit modulo3 and epsilon is nonzero there,

    3 divides m_j  if and only if  j is odd.           (2)

In particular the first forced halving yields

    m -> 3(m+epsilon)/2.                             (3)

Thus every post-odd one-halved endpoint has a factor3 in this modulus
gauge. Additional legal halvings alternate between a divisible and a
nondivisible modulus. The statement uses the actual number of halvings;
one cannot replace it by an unrecorded accelerated step.

Conversely, if3 divides m, either single legal modulus branch leaves this
divisibility condition:

    3m+2epsilon != 0 mod3,
    (m+epsilon)/2 != 0 mod3.                          (4)

Equations(2)–(4) explain the correct clock for comparing triangle resets
with the central tables. They do not identify a modulus with a realizable
Hamiltonian-path count or supply a graph operation for either branch.

## A2. The exact scalar guard for an ordinal triangle

Let S be the set of attained tournament Hamiltonian-path counts, including
H(empty)=1 if the empty residual graph is used. For positive integer N:

    some tournament of count N has an ordinal C3 factor
       if and only if  N=3h for some h in S.          (5)

Necessity: paths concatenate through an ordinal sum, so
H(C3 -> R)=3H(R), and likewise for an interior ordinal factor or R -> C3.
Deleting that factor leaves an actual residual tournament of count N/3.
Sufficiency: if H(R)=N/3, construct C3 -> R. QED.

For a given graph this remains only a necessary scalar condition. Even
when N/3 is attained, that particular source need not have the factor in
the required extreme position. The complete guard includes the named
induced triangle, its orientation and its uniform contacts with every
other vertex.

Even this existential guard does not identify a mutation-accessible source.
At count9, the residual count3 is attained, but the two six-vertex sources
in [matrix M3](../matrix/PROOF.md) have determinants25 and81. No ordinary
mutation word of any length connects them, and the determinant25 source
cannot reach any six-vertex ordinal-triangle presentation. This strengthens
the source-level boundary without changing the existential statement(5).

**H63 hostile, all orders.** The attained count63 is supplied by
**THM-344**, [n8-h63-complete-omega](../../../../01-canon/theorems/THM-344-n8-h63-complete-omega.md);
the literal eight-vertex witness is independently recounted below.
**THM-1370**, [h-spectrum-omits-7-21-all-n](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
excludes21 at every order. Equation(5) therefore proves that no H63
tournament has an ordinal C3 factor, despite9 dividing63.

This hostile is visible on a valid signed table trajectory:

    epsilon=-1:  x=11 ->32 ->16,
                 m=43 ->127 ->63.                   (6)

All three modulus labels43,127,63 are attained tournament counts, but the
final label cannot admit the needed triangle presentation. The plus-sign
example13 ->41 ->21 is a different failure: its final count is forbidden
altogether. In either case the missing proof is the realization guard,
not the arithmetic identity or its modulo3 calculation.

## A3. Exact triangle budget on the recovered seven-tower

More generally, some presentation with k ordinal C3 factors exists at
count N if and only if3^k divides N and N/3^k belongs to S. The factors
may be grouped at the front for sufficiency. Their deletion proves
necessity even if they were separated by other ordinal blocks.

For every r>=2, the exact maximum possible number of such factors among
all sources of count7*3^r is

    k_max(7*3^r)=r-2.                                (7)

Lower bound: use r-2 copies of C3 followed by the inherited H63 source.
Upper bound: k=r-1 would leave count21, k=r would leave count7, and
k>r would not leave an integer. Both excluded quotients are forbidden
at every graph order by THM-1370. QED.

The counts7*3^r are attained; this is not the old false gap-tower claim.
What is unavailable is the last two powers of3 as *separate ordinal
triangles*. The invariant distinguishes arithmetic valuation from
physical factorization. At count189 one triangle can be exposed, while
at count63 none can.

## A4. Why this reset alone cannot close the nonzero ternary induction

Let R denote any finite composition of tournament-return blocks of at
most four ordinary Fomin–Zelevinsky mutations. The inherited four-return
theorem gives H(R(T))=H(T). It follows by composition, without any
all-length mutation-invariance assumption.

Either desired digit map sends every integer H to a nonmultiple of3:

    3H+2 = 2 mod3,       3H+4 = 1 mod3.               (8)

Therefore after a realized nonzero digit, no composition of these
neutral return blocks can produce an ordinal C3 factor at the same
count. In particular an extreme-triangle transport cannot be the sole
rearming mechanism between successive nonzero digits. This is an
all-input obstruction; it does not rely on a finite failed search.

For the prime and7p lanes it is already restrictive at the input. A
prime count can have this factor only at H=3. If p is prime and H=7p,
divisibility by3 would force p=3 and H=21, which is forbidden.
Thus no attainable7p source is in this particular triangle-reset domain.

The strongest surviving use is within the compatible3-divisible count
fibres, or as one stage of a larger construction with a different next
guard. A new operation can change the count back to a multiple of3; a
new reset family can instead work at counts1 or2 mod3. Both require an
actual construction. Multiplying a formal response vector by3 and then
dividing it again does not establish a graph operation.

## A5. A precise composition problem in place of scalar analogy

The neutral return relation stays inside each set of actual sources with
one fixed H. An arithmetic count-changing step may be attempted only
from the members with its full attachment guard. The correct target is:

    given an admissible source, find a neutral source change,
    realize one required digit from the changed source,
    and prove the resulting source has a usable next guard.

The source includes the graph and named ports, not merely H. A finite
controller for residues can select candidates; an integer/factor/port
certificate must establish each chosen edge. The central-table clock
and the neutral mutation relation preserve different predicates:

| Map | Preserved predicate | Lost data / required sidecar | Cheap hostile |
|---|---|---|---|
| Integer orbit -> central modulus | Signed Collatz edges under4x+epsilon | Ordinary clock, sign and positive integer lift |13->41->21 |
| Odd step plus j halvings -> modulo3 guard | Factor3 iff j odd | Attainable cofactor and actual factor |43->127->63 |
| Short tournament return -> H fibre | Exact total H and fixed vertex set | Root readiness, factor position, joint attachments |H63 has no C3 factor |
| Neutral return then nonzero digit | Count equation when physically realized | A new reset guard after loss of3-divisibility |Equation(8) |

The nonzero-digit obstruction is compatible with the existence of useful
repeatable neutral source changes. It identifies the additional family
of guards a universal prime/7p induction must supply.

## Reproduction and bounded controls

Run `python3 -B check.py` and its `-O` version in this directory.
The checker verifies both signs at every odd x<=511 and every legal
halving depth; recounts the inherited eight-vertex H63 graph by subset
DP and all8! permutations; recounts C3 -> that graph at count189; checks
the displayed clock hostiles, the scalar factor-budget consequences
through r=16, and the ternary residue obstruction for every odd H<=501.
The all-order exclusion and short-return theorem are inherited proved
dependencies, not extrapolations from these finite controls.
