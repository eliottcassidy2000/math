# A genuine role for 105: singular golden bridge words

**PROVED, scoped to the seven-letter single-bridge language:** its full
three-coordinate observer can be made divisible by a positive integer m,
uniformly over all input states, exactly when every prime factor of m divides
105. Explicit mixed words supply primitive contractions at 3, 5 and 7.
**FINITE-EXACT:** the complete row-subspace graphs at these primes, the
displayed integer products, and the independent literal-atom audit.
**OPEN:** conversion into a reusable source switch, the named two-port
17-coordinate induction, and the unrestricted Chan–Pak spectrum conjecture.

## 1. Inheritance, operation and actual ports

The closest proved mechanism is [golden Hamilton RECURRENCE G1–G4](../../golden_hamilton_transfer_20261103/RECURRENCE.md),
the all-order one-reversed-bridge join. Its observer is

    v(T,a)=(A,B,C)=(H(T)-E(T,a), E(T,a), C(T,a)),

where C counts paths starting at the named port a plus spanning two-path
covers with a path starting there. The new block B has first vertex b and
last vertex c. All old-to-new arcs are forward except b->a; the new active
port is c. Counts are ordinary integer path counts on actual tournaments.

The literal word and all physical vertex names remain sidecars. The three
coordinates answer every future operation of this particular join language;
they do not answer arbitrary cuts, insertions, both-end contexts, or queries
at two retained old ports. In particular, retaining actual old names through
the join is not a proof that the 17-coordinate two-port observer closes.

The hostile is the inherited 81/27 **automorphism-order** pair, not the rings
Z/81 or Z/27. The corrected near miss is treating a scalar 55->111 arrow as
executable from every H55 source. The least-used sidecar here is the
normalized p-adic phase after a visible zero residue. The board is:

| Concept | Precise role | Boundary |
|---|---|---|
| Actual bridge word | Legal Boolean graph construction | Does not preserve H |
| Three-coordinate observer | Exact future bridge responses | Does not replace joint ports |
| Singular directions | Products can remove every visible residue | Lifted phase can survive |
| 55 and source choice | Two golden H55 sources fail the inherited 111 branch | Scalar equality does not identify sources |
| 105 | Exact prime support of uniform observer divisibility | No exact spectrum coverage follows |
| Golden packing clock | A quotient can hide displacement or phase | No packing-to-tournament conjugacy |

The seven matrices in RECURRENCE G4 have determinants

    s:1, c:1, q0:1, q1:9, q2:5, q4:9, q5:21.

We use row states and right multiplication, in **processing word order**.
The code independently recompiles these matrices from literal graphs and
endpoint/cover subset DP, rather than importing their numerical table as
its only evidence.

## 2. Three explicit primitive contractions

Let M(w) be the product of the matrices of w. The following are actual words:

    W3 = q5 s q2 q5 s s q5,
    W5 = q2 c c q2 s s s s q2,
    W7 = q5 q0 q2 q5 q2 q5.

Their vertex costs are respectively 19, 22 and 24. Direct multiplication gives

    M(W3) = 3 U3,
    U3 = [[51188,12543,88569],
          [41096,10070,71107],
          [72061,17658,124685]],

    M(W5) = 5 U5,
    U5 = [[160526,32615,346665],
          [134642,27356,290767],
          [295163,59970,637421]],

    M(W7) = 7 U7,
    U7 = [[1309872,349995,2246175],
          [1028111,274709,1763010],
          [1730050,462265,2966700]].

Their exact Smith data are

| Word | det M | SNF(M) | det U | SNF(U) |
|---|---:|---|---:|---|
| W3 | 46305 = 3^3·5·7^3 | (3,3,5145) | 1715 = 5·7^3 | (1,1,1715) |
| W5 | 125 = 5^3 | (5,5,5) | 1 | (1,1,1) |
| W7 | 231525 = 3^3·5^2·7^3 | (7,21,1575) | 675 = 3^3·5^2 | (1,3,225) |

In each row U_p is invertible over Z/p^k for every k>=1. In particular the
W5 normalization is an integer unimodular change of coordinates. This is a
stronger exact statement than observing that one scalar output is divisible
by p. In particular, for every nonzero integer state v, the gcd of the three
coordinates of v M(W5) is exactly five times the gcd of the coordinates of
v. The shorter word q4 q1 has cost 8 and matrix

    [[117,36,186],[90,27,144],[144,42,231]],

divisible by 3, but its normalized determinant is 3. It is a useful hostile
to the unsupported inference that every visible contraction is primitive.

**Proposition G105 (uniform divisibility).** For an integer m>=1, there
exists a golden bridge word w with every entry of M(w) divisible by m iff
all prime factors of m belong to {3,5,7}.

**Proof.** Necessity follows from the determinant product. If a prime p
divides every matrix entry then p^3 divides its nonzero determinant, whose
only possible prime factors are 3,5,7. Conversely, write m=3^a5^b7^c. The
literal word W3^a W5^b W7^c has matrix

    m U3^a U5^b U7^c,

which is integral. The empty word handles m=1. This supplies a construction
at every exponent, not a bounded residue census. QED.

Thus 105 is the squarefree support of the singular primes in this exact
operation language. The connection is an operator identity, not a claim
that a count105 tournament belongs to a special golden isomorphism class.
The modulus55 is excluded: its factor11 is invertible in every letter.
At every other prime outside {3,5,7}, every word matrix is invertible.

The same necessary obstruction applies even to demanding that H be divisible
by such a prime for **every actual nonempty input**. The three actual
singleton prefixes have state rows (0,1,1),(1,0,2),(2,1,3), determinant2.
They span over every odd prime. If all three outputs had H=0 mod p, then
M(w)(1,1,0)^t=0 mod p, impossible for an invertible matrix. At p=2 each
tournament H is odd. This assertion concerns a uniform operation, not the
existence of an individual tournament with H divisible by 11 or 55.

## 3. Why mixed singular choices work, and what they do not reset

The procedural search operates on row **subspaces**, not millions of
Hamilton-count targets. Starting from F_p^3, apply S -> S M_l and put the
result in reduced row-echelon form. There are exactly

    2p^2+2p+4

subspaces in dimension three. The seven-letter graph reaches every one at
p=3,5,7: respectively 28,64,116 states, including zero. This reachability
is FINITE-EXACT. The graph at p=11 contains only the full space. Breadth-first
search proves shortest **letter counts** to zero of 2,7,6 at 3,5,7;
W7 is one shortest word. W3 and W5 were chosen for primitive normalization,
not asserted to be shortest. No global vertex-cost minimum is claimed.

The mechanism is rotation of the surviving subspace between singular
letters. Repeating a single q5 need not kill its surviving unit directions;
interleaving q0 and q2 reaches the next kernel. This is a concrete model
for paying a valuation cost by changing the active singular direction.

**Proposition Gphase (exact remaining memory).** For p in {3,5,7}, every
nonzero integer row v and k>=0 satisfy

    min_i v_p((v M(Wp)^k)_i) = k + min_i v_p(v_i).

Moreover, for every precision r>=1, the normalized output

    p^(-k) v M(Wp)^k = v Up^k  (mod p^r)

determines v mod p^r uniquely.

**Proof.** M(Wp)^k=p^k Up^k, and Up^k and its inverse are integral over
Z_p because det Up is a p-adic unit. Such an automorphism preserves the
largest power of p dividing all coordinates. The second assertion is the
same invertibility reduced modulo p^r. QED.

W7^2 and W7^3 therefore give explicit zero residues modulo49 and343, but
they have not erased the phase: the quotient by49 or343 retains it
bijectively. These are controlled contractions, **not H-neutral native
resets**, and they do not solve the older mod49 row-lattice obstruction
for a different quotient/operation language. Any proposed transfer between
the two must supply an actual intertwining map and carry the joint ports.

The inherited golden55 sources make this loss visible. Their states are
(34,21,61) and (32,23,63). Both have H55 and next-singleton H95. After W7,
both raw states vanish mod7, while their normalized phases are respectively
(3,0,2) and (2,3,5) mod7. A reset controller remembering only the raw zero
would identify different surviving data. The full rows and states are
retained in CHECK.json. The concurrent [target audit](../targets/receipt.json)
strengthens their source-compatibility obstruction; its all-order statements
are separate evidence, not a consequence of these modular calculations.

## 4. Recovered 55, 105, golden and packing connections

The earlier [golden seam recovery R1–R4](../../golden_seam_recurrence_20261103/RECOVERY.md)
already made the following exact distinctions:

- At order6 the fixed-path complement/class observer has visible rank56,
  omitted rank55 and total111. It observes a different operation and output
  from the three-coordinate bridge observer.
- The square-sum graph Q56 has111 edges, a supplied Hamilton path with55
  edges, and cycle-space dimension56. Its boundary exact sequence has
  ranks56,111,55; no canonical intertwiner with class observation was supplied.
- The golden source codec gives F10=55=40+15 at source order9, with a
  terminal flag. The corresponding marked source count is184. Its Fibonacci
  clock is not the singleton bridge's tribonacci clock.

The older [PRIMEGAME source clarification](../../session126_primegame_symmetry_20261031/primegame/SOURCE_ATTRIBUTION_CLARIFICATION.md)
gives a different operational 105: at input14 one ordered fraction rack
chooses15/2 and reaches105, while the alternate rack chooses15/14 and
reaches15. The ordered guard is essential. Neither rack is identified with
the present matrix language. The new G105 proposition gives 105 a separate,
proved operator role; the numeral alone supplies no map between them.

The packing connection is through an **explicit loss of lift data**, not a
claim about a physical optimal packing. In the inherited
[golden metric note, sections2–3](../../session91_packing_golden_carriers_20261029/packing/GOLDEN_METRIC_AND_PALEY_TORI.md),
Q=[[0,1],[1,1]], C=4I-Q, det C=11, and Q^5-I=Q^3 C. Therefore Q^5 is the
identity on Z^2/C Z^2 while the integer lift is not the identity. The
missing coordinate is the lattice displacement. The same note proves
Q^t N Q=-N for N=[[2,1],[1,-2]]; this is a split arithmetic norm, not a
positive-definite packing metric. No current physical-packing literature
claim is made or needed here.

The [norm-controlled observer theorem, sections3–5](../../session91_packing_golden_carriers_20261029/zeckendorf/NORM_CONTROLLED_OBSERVERS.md)
uses the same N: a primitive observation (a,b) has observation determinant
a^2+ab-b^2. At11 the scalar golden eigenchannel can compress two registers,
but the ordinary integer channel cannot; the second eigenchannel repairs
the loss. This is not the same prime mechanism as golden Hamilton bridges:
**every bridge matrix is invertible at11**, while the Fibonacci observation
map can be singular there. An operator determinant and an observation
determinant must not be substituted for one another.

Finally, the [guarded flip leakage theorem](../../session91_packing_golden_carriers_20261029/operators/GUARDED_FLIP_LEAKAGE.md)
proves rank(PB(1-P)BP)=I(G)-2^b for independent-set legality, with b the
number of bipartite components. On the five-cycle,11 feasible words sit
inside an unguarded complete55-edge increment graph, but only15 insertion
arcs are legal. The guard sidecar, like a tournament's actual ports, cannot
be inferred from a scalar count. This is a precise source/target comparison:
both compress a fully specified transition system; each needs its own
observer and legal-operation proof. No map from a square packing to every
tournament source is asserted.

The constructive next ingredient is consequently sharper: find a legal
source-changing operation that transports **normalized phase plus joint
ports** back to the input predicate for both ternary branches. A positive
phase contraction is now available in this restricted language; closure of
that stronger return predicate remains the unproved part.

## 5. Audit, provenance and reproduction

The prior exact matrices and arbitrary-block closure are inherited, with
their proofs in RECURRENCE G1–G4. Existing mixed-return work in
[RETURN_BANKS sections3–5](../../adaptive_state_induction_20261103/RETURN_BANKS.md)
already identified q5 as a 3-singular escape channel and proved an infinite
source-pair identity. The present result is a different full three-dimensional
modular semigroup statement; it does not claim priority for singular
channels. Targeted searches for the displayed words, 105-smooth/congruence
eraser language, and golden reset/Smith claims found no prior version of
G105 in the inspected repository. No external novelty claim is made.

`check.py` reconstructs every atom from Boolean graph data, regenerates all
three matrices, computes Smith invariants from gcds of minors, exhausts the
declared small subspace universes, and checks the normalized phase bijections.
An eight-vertex composite is counted directly by subset DP. The large
19/22/24-vertex blocks use the inherited all-order bridge theorem, rather
than a disguised giant Hamilton-path enumeration. Normal and optimized
receipts are identical; all 59 explicit checks survive Python -O.

    python3 -I -B 05-knowledge/results/chanpak_55_105_golden_20261110/golden/check.py --output CHECK.json
    python3 -I -B -O 05-knowledge/results/chanpak_55_105_golden_20261110/golden/check.py --output CHECK-optimized.json

The negative controls are prime11 invertibility, the nonprimitive q4 q1
contraction, and the two retained H55 normalized phases. There is no claim
of full spectrum coverage, a way to construct H7 or H21, or a replacement
for the 81/27 automorphism hostile.
