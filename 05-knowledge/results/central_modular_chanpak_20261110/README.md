# Central multiplication, exact phase, and constructive address arithmetic

**PROVED in the scoped notes below / FINITE-EXACT controls. Chan–Pak
completeness and Collatz convergence remain OPEN.** This session turns the
user's central-table observation into an exact graph codec and obtains a
new classification inside the actual golden tournament bridge grammar.
It also identifies a restriction that a recursive seed-switching proof
must escape. No external priority claim is made for the in-session results.

## 1. The strongest positive result: an exact content monoid

In the inherited seven-letter bridge grammar, there exists a word whose
matrix is rU, with U an integer matrix of determinant1, **if and only if**

    r = 5^b 9^k 21^c,     b,k,c >= 0.

Equivalently, the attainable universal gcd multipliers are exactly that
monoid. This is stronger than the preceding session's divisibility result
for arbitrary products of3,5,7. The new words are

    W9  = q1 s s q1 s s q1,
    W21 = q5 s q0 c q0 q5 q0 c s s s q5.

Together with the inherited W5, they give the constructive word
W5^b W9^k W21^c at every exponent. The determinant valuations prove that
no other universal exact content multiplier is possible in this grammar.
In particular3 and7 separately fail this stronger requirement, while21
succeeds. The legal singular operation q5 couples the two primes.

For the bridge state v=(H-E,E,C), division of vM by r leaves vU, from
which v can be recovered integrally and modulo **every** modulus. Thus
W5 W21 gives a globally primitive105 operation: its normalized phase
survives simultaneously modulo49,343 and all other precisions.
The literal word remains an actual graph construction; its old vertex
names survive. The normalized vector need not itself be an actual graph,
and this is not H multiplication or an H-neutral reset.

See [information G1–G4](information/PROOF.md) for the theorem, explicit
matrices, simultaneous subspace alignment and information budget. The
search used small row-subspace orbits; it did not enumerate large
tournaments or scan a new target-count interval.

## 2. The user's four central cells encode the two Collatz signs

In the nonzero-residue table modulo odd m, write m=4x+epsilon,
epsilon in{+1,-1}. The two central values are x and3x+epsilon, each twice.
Each row sums to m. These are the two copies of m in the observation;
there are four entries altogether.

| Modulus class | Source | Decoded destination |
|---|---|---|
| 5 mod8 | odd x, plus sign | 3x+1 |
| 1 mod8, m>=9 | even x, plus sign | ((3x+1)-1)/6=x/2 |
| 3 mod8 | odd x, minus sign | 3x-1 |
| 7 mod8 | even x, minus sign | ((3x-1)+1)/6=x/2 |

This is an exact decoding of all edges, with parity separating the two
source halves for each sign. The raw pair is not itself an even-source
Collatz edge. On moving modulus labels the entire graph is conjugate to

    m -> 3m+2epsilon   on odd sources,
    m -> (m+epsilon)/2 on even sources.

The central identity is inherited from **THM-4792**,
[two-sheets-mod-odd-m-merged-centre-reflections-unit-shift-and-alternation-classification](../../../01-canon/theorems/THM-4792-two-sheets-mod-odd-m-merged-centre-reflections-unit-shift-and-alternation-classification.md)
and the earlier modular-tiling note. [Central C1–C3](central/PROOF.md)
supplies the decoder, exact domains and recursive finite-word cylinders.
The ordinary-clock legal words form the inherited no-OO language; ordered
carry and the final-odd bit recover source alignment. An unbounded integer
lift remains necessary for a particular orbit.

## 3. The prime/composite split depends on the map being observed

Multiplication by a modulo m has gcd(a,m) points in each nonempty fibre.
Every nonzero row is reversible exactly when m is prime. A composite
modulus still has reversible unit rows, including both central row labels.
Indeed the central block has the same unit pattern for every odd modulus.

For uniform inputs, the row's lost information is exactly log2 gcd(a,m).
The central two-coordinate block is a different observer: it retains only
the difference x-y and loses the diagonal coordinate, even at prime m.
Smith invariants extend this exact fibre calculation to all matrix
moduli at once. This makes the information-transfer intuition precise
without asserting that primes universally minimize every representation.

Prime-only table dynamics are not closed: the ordinary Collatz cycle
becomes5 ->17 ->9 ->5. Likewise, the golden exact-content monoid admits
the composite21 while excluding the individual factors3 and7. The
information question must specify both the observation and the legal
operation grammar.

## 4. A recursive proof needs more than bounded suffix switching

The actual two-flank compiler has the form

    Phi(A,B)=d(g(A)+g(B)-2)+qAB+eA+sB+t+u(A mod2)+v(B mod2),

where d is the deletion-core Hamilton count and g is the inherited
binary recurrence. A translated four-point stencil extracts exactly -2d.
After any fixed binary suffix the coefficient of g remains1.
Consequently a uniform suffix/cylinder rule Phi_new=lambda Phi_old+kappa
forces d_new=lambda d_old. This also holds for finite policies controlled
only by bounded low bits.

Therefore a bounded core bank cannot support unlimited uniform3H+2/3H+4
returns by that policy class. This does not exclude unbounded carry
guards, height-dependent choices, non-dyadic rebases or actual core growth.

The positive construction is an explicit nonnegative matrix lift with
m+3 integer coordinates for all g(mn+r), 0<=r<m, for **every fixed odd m**.
It shares the height n and two parity indicators. This generalizes the
inherited multiplier-three H13 rebase. It makes an infinite arithmetic
family manageable with two digit matrices and actual graph reconstruction.
These are finite coordinates with unbounded values, not a finite-state
solution of the integer problem. See [compiler C1–C4](compiler/PROOF.md).

The [affine-gauge lemma](gauge/PROOF.md) explains the boundary of the
most direct scalar bridge. The common even law (H+1)/2 and odd law3H+delta
force H=(delta+2)x+1. For delta=0,2,4 these domains contain7 or21, and
even updates13->7 or41->21 leave the attained spectrum. Arithmetic can
therefore act safely on a **constructor address**; its equality to a
Hamilton count needs a separate proof.

## 5. The remaining constructive gate

The full range conjecture is [Chan–Pak, Conjecture6.2, section6.5](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf),
which credits the2016 formulation to bof and Gordon Royle. We retain its
OPEN status. For the prime and7p lanes the previous
[one-active-padding reduction](../chanpak_55_105_golden_20261110/h5/TERNARY_PADDING.md)
still applies: prepare one singleton on each side, then the possible
low-ratio readings reduce to three named-cut expressions. This session
supplies no universal source switch between those branch-ready cuts.

The next precise target is an actual positive chart-switching rule with
an exact count equation3H+2 or3H+4, whose output is usable again. Its
certificate must contain the named ports and their joint response,
the named F49/F343 charts, normalized integer phase, and the physical source.
The residue rings Z/49 and Z/343 are not the fields F49 and F343. Our
unimodular theorem concerns the three-coordinate integer bridge observer;
it does not identify it with all the
[named-port field codecs](../session180_170807_returns_20261104/phases/candidate_phases.py).
The [219-key seed packet](../recursive_seed_switching_20261103/README.md)
already warns that finite scalar registers modulo49/343 do not recover
the full port phases. The81/27 hostile
concerns automorphism orders, not arithmetic modulo81 or27. Neither the
three-coordinate golden state nor the carry coordinates silently replace
the two-port17-coordinate attachment observer.

The incoming H13 carry rebase demonstrates one exact zero fibre; the new
odd-multiplier lift makes more such fibres explicit. The new simultaneous
primitive word removes one genuine phase-loss defect of previous mixed
singular words. What remains is to prove exact target reachability and
reusability in the positive, physically realizable state set.

Incoming Session255 supplies a finite exact-equality compiler for comparable
nonnegative states. The [height/cancellation bridge](compiler/INCOMING_HEIGHT.md)
tests its applicability here. Two actual eleven-vertex tournaments at right
addresses37 and39 both have H=611, but their carry difference is
(-10,2,0,0): the height rises while g falls. In general any same-parity
fixed-chart collision with a positive linear address coefficient requires
that mixed sign. Thus the incoming support theorem applies where
comparability is proved; these exact repair collisions require signed
integer cancellation data or a different proved common frame. This
identifies another specific missing coordinate, rather than another
unexplained failure of a finite search.

The final incoming
[iterable ready chart](../uniform_return_metagraph_20261111/reset/ITERABLE_RETURN.md)
is also compatible with this distinction: its actual loop reaches3H+2
and then returns at4H-33, losing one old vertex and retaining the others.
It is indefinitely repeatable, but its final target changes. Our obstruction
concerns exact uniform ternary closure in the specified suffix language;
it neither retracts this physical loop nor turns its target drift into a
target-preserving induction.

## Validation and files

- [Central proof](central/PROOF.md), [information proof](information/PROOF.md),
  [compiler proof](compiler/PROOF.md), [gauge proof](gauge/PROOF.md).
- [Session inheritance and decisions](SESSION.md), audit notes in the
  component directories, and source pins in SOURCE_PINS.json.
- Run `python3 -B 05-knowledge/results/central_modular_chanpak_20261110/reproduce.py --rerun`
  from the repository root. The normal and optimized checker receipts,
  source hashes, local links and artifact manifest are checked.

Universal conclusions come from the displayed proofs and explicit words.
The finite controls include composite/prime rows, the raw-edge hostile,
complete small parity fibres, actual small tournament recounts and exact
matrix products. They supply no new full-spectrum coverage record.
