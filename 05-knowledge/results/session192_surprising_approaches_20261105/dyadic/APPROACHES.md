# Arithmetic/dynamics connection board and decisive probes

Status labels distinguish new deductions, inherited machinery, conditional reductions, and failed transfers. The main proved packet is `PROOF.md`; no new Collatz theorem or generalized undecidability assertion is made.

## 1. Retained past digits versus fresh future digits: successful transfer

- Sources: THM-2555 / `natural-extension-sheet-charge-and-future-digit-boundary`; THM-4473 / `collatz-digit-chains-rotation-repunits`.
- Target/map: one-apex word -> left binary fraction and reverse-complement right fraction.
- Preserved: actual H via its exact dyadic correction. Lost in the continuum: individual integers, marked continuation data, multiplicities when passing to support.
- Decisive test: uniformly sampled strong words have exactly the S26 product-of-uniforms limit. This was recovered, then correctly credited to S26, not claimed new.
- New survivor: independently selectable end cylinders control magnitude uniformly over an arbitrary middle. The middle can consequently carry a separately proved residue program. The simultaneous theorem in PROOF sections 4–5 is the useful transfer.

## 2. Finite-range arithmetic interactions: successful transfer

- Source: the exact one-apex slot expansion, read as a finite-memory digital dynamical system modulo 2^k.
- Map: truncate pair distances to at most k, because all farther weights vanish modulo 2^k.
- Needed sidecar: the k-bit guard on each side of a switch. Without it, the exact 01->10 increment depends on the full prefix and suffix.
- Decisive test: the padded integer 2^j-2 before a switch and k ones after it give increment 2^j modulo 2^k in *every* exterior.
- Proved output: a length-(2k^2-2) bank with k-1 independent switches; every odd residue occurs inside any fixed exterior. The bank was exhaustively checked for k<=8 and many short exteriors.
- Limit: the size error grows with k. The available inequality is never tight enough to force an exact integer from its residue.

## 3. Collatz-style inverse ternary descent: exact conditional reduction

The arithmetic side of the ternary target grammar has a particularly small complete entry bank. Let S be a set of positive odd integers containing

    1, 3, 5, 23, 63, 67,

closed under multiplication, and satisfying

    h in S, h>=3  =>  3h+2 and 3h+4 are in S.          (* )

Then S contains every positive odd integer except possibly 7 and 21.

**Proof by strong induction.** For an odd N, select its unique odd inverse parent:

    N=0 mod 3: q=N/3;
    N=2 mod 3: q=(N-2)/3;
    N=1 mod 3: q=(N-4)/3.

Except at the listed small entries, q is a smaller positive odd integer. A parent q=7 creates the exceptional targets 21,23,25; a parent q=21 creates 63,65,67. The targets 21 and 7 are allowed exceptions; 23,63,67 are seeds; 25=5^2 and 65=5*13, while 13=3*3+4. The q=1 boundary is covered by seeds 3,5 and exceptions 1,7. Every other parent is neither 7 nor 21, so induction and multiplication by 3 or (*) applies. QED.

The seeds have actual tournament witnesses, checked in `seed_check.py`: one-apex words 10,100,11000,1001100 give 3,5,23,63, and the two-apex word 12301 gives 67. A singleton gives 1. Multiplication of H is physically available through ordered joins. Thus the sole unproved step in this reduction is the **universal existential closure** (*), not a number-theoretic convergence statement. No fixed macros, preservation of a marked state, or strong-child assertion is included. S191's 13->41->{125,127} end-growth obstruction defeats one proposed physical implementation without refuting this existential condition.

This is a Collatz-inspired termination proof template, but its descent is elementary and strictly decreases. It must not be advertised as resolving the physical closure, or as new finite experimental coverage.

## 4. Christoffel free-cycle coprimality: failed transfer with a minimal short hostile

- Source: THM-4716 / `balanced-collatz-cycles-are-three-smooth-ladders-and-are-free`.
- Attempt: use the same balanced binary necklaces as apex words and transfer the coprimality to the clock D=2^q-3^p.
- Preserved: the binary word and its adjacent-transposition order. Lost: the Collatz affine carry and its linear recurrence modulo primes dividing D.
- Decisive hostile: q=3,p=1, the balanced word 100 is strong and has H=5, while D=8-3=5. Thus gcd(H,D)=5, not 1. This exact hostile is in verify.py.
- Strongest survivor: the ordered binary word can still organize monotone H ladders through the inherited adjacent-switch law, but the increment is a sum of left and right binary responses, not a single 2^i 3^j carry. No prime-lane conclusion transfers.

## 5. FRACTRAN return states: scoped stopping reason

- Source: `05-knowledge/results/collatz_procgen_20260922_barrier_atlas.md`, Conway/register-machine row.
- Candidate: encode construction-state coordinates by prime exponents, with a special return marker for the scalar H consumer.
- Required map: actual tournament operation -> legal register update -> return to the same usable marked interface. A code of the number is not such a map.
- Cheapest hostile already available: words 00 and 01 both have H=1 but their appended-zero counts are 1 and 3. In the two-apex grammar S191 proves exact rational linear continuation dimension 17. Scalar return states discard necessary registers.
- Stopped result: no tournament compiler or persistent return invariant was found. No undecidability consequence for this fixed tournament problem follows from universality of a broad family of generalized Collatz maps. Retaining the full physical continuation carrier is the useful methodological instruction, not a theorem transfer.
