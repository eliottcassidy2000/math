# Tournament Theory as Coding Theory

**Session:** kind-pasteur-2026-03-23-S20dp
**Status:** DEEP REFLECTION

---

## The Dictionary

A tournament on n vertices is a binary word of length m = C(n,2). The entire theory translates:

| Coding theory | Tournament theory |
|---------------|-------------------|
| Binary alphabet F_2 | Forward/backward arc |
| Codeword of length m = C(n,2) | A labeled tournament |
| The full code F_2^m | All 2^m tournaments |
| Group action S_n on coordinates | Vertex relabeling |
| Orbit = equivalence class | Isomorphism class |
| Quotient code F_2^m / S_n | The metagraph G_n |
| Hamming distance d(x,y) | Number of arc flips between T and T' |
| Single-bit error | Single arc reversal |
| Hamming graph Q_m | The complete arc-flip graph on all tournaments |
| Quotient graph Q_m / S_n | The metagraph G_n |
| Weight of codeword | Usually Hamming weight; here H(T) = Hamiltonian paths |
| Minimum distance of code | Minimum Hamming distance between orbits |
| Covering radius | Max distance from any tournament to nearest class representative |
| Undetectable error | Neutral arc (arc flip that preserves the iso class) |
| Error-detecting capability | d_min - 1 of the orbit code |
| Dual code | The complement tournament space |
| MacWilliams identity | Complement symmetry: H(T) + H(T^c) structure |

---

## What Each Object Represents

### A tournament IS a codeword

Every binary string of length C(n,2) encodes a tournament. There are no invalid codewords. The "code" is the entire space. But the INTERESTING structure is the quotient by S_n: which codewords are equivalent under coordinate permutation.

### An isomorphism class IS a coset

The S_n-orbit of a codeword = all codewords equivalent to it under coordinate permutation. The orbit size = n!/|Aut(T)|. The number of orbits = A000568(n). This is the "code" in the quotient sense.

### The metagraph IS the channel

The metagraph G_n = Q_m / S_n is the communication channel where:
- Input: an iso class (the "message")
- Noise: a single arc flip (a "bit error")
- Output: the resulting iso class

The channel's transition graph IS the metagraph. Its properties (diameter, connectivity, spectral gap) determine the channel's capacity and error behavior.

### H(T) IS the weight — but a NONLINEAR one

In standard coding theory, the weight = number of 1s. For tournaments, the natural weight is H(T) = number of Hamiltonian paths. This is a HIGHLY nonlinear function of the binary representation.

The OCF H(T) = I(Omega(T), 2) says: the weight function IS the partition function of a hard-core lattice gas on the conflict graph at fugacity 2. This is a statistical mechanical weight, not a combinatorial one.

### Neutral arcs ARE undetectable errors

A neutral arc = a bit position where flipping doesn't change the orbit. This is exactly an undetectable error: the error goes through the channel without being noticed.

SL_orbits = the number of (message, error-position) pairs where the error is undetectable. The twin_SL formula counts the dominant mechanism: two coordinates that are interchangeable (transposition symmetry in the code).

### Self-complementary IS self-dual

A tournament T is SC iff T ~ T^c (complement = bitwise NOT). In coding theory, a code is self-dual if C = C^perp. The SC tournaments are the "self-dual codewords" — fixed points of the complement involution.

---

## The Three Codes

### Code 1: The Full Tournament Code (rate 1)

All 2^m binary words of length m. Every word is valid. The S_n action partitions this into A000568(n) orbits. This is a "trivial" code in the coding sense (no redundancy) but the GROUP ACTION gives it rich structure.

### Code 2: The Regular Tournament Code (rate ~ 0)

The 24 regular tournaments at n=5 form a (10, 24, 3) code. This is a VERY GOOD code:
- Length 10, size 24 = 2^{4.58}, minimum distance 3
- Detects 2 errors, corrects 1
- Gilbert-Varshamov bound: 2^10 / V(10,2) ~ 1024/56 ~ 18.3. We have 24 > 18.3. Exceeds the bound!
- The 24 codewords are equidistributed (constant weight = 5, since score = 2 everywhere)
- The complement of any codeword is also a codeword (the code is self-complementary)

At general n (odd): the regular tournaments form a constant-weight code of weight (n-1)/2 * (n-2)/something... actually weight = number of "1" bits varies depending on encoding. But all regular tournaments have the same Hamming weight (same number of forward arcs).

### Code 3: The SC Backbone Code

The SC classes form a code within the orbit space. At n=8 it splits into 7 components = 7 SUB-CODES. The 7 sub-codes correspond to the 7 lines of the Fano plane.

The Fano plane IS the parity check matrix of the [7,4,3] Hamming code. So the 7 SC backbone components at n=8 correspond to the 7 PARITY CHECK EQUATIONS of the Hamming code.

---

## The Hamming Code at n=7

The [7,4,3] Hamming code has:
- 7 positions (= 7 vertices of QR_7)
- 4 information bits
- 3 parity checks (= 7 - 4)
- Parity check matrix = Fano plane incidence matrix
- Minimum distance 3

The Paley tournament QR_7 encodes this code as a tournament:
- 7 vertices = 7 code positions
- Connection set {1,2,4} mod 7 = quadratic residues = parity check support
- The 7 directed 3-cycles = the 7 lines of the Fano plane = parity check equations
- |Aut(QR_7)| = 21 = the automorphism group of the Hamming code restricted to cyclic shifts

The tournament STRUCTURE (directions of arcs) encodes the SIGNS in the parity checks. The Hamming code's syndrome decoding IS the cycle structure of QR_7.

---

## Error Correction in the Metagraph

### The Channel: Arc Flip Noise

A "message" is an iso class. The channel randomly flips one arc. The receiver observes the new iso class.

- If the flip is a SELF-LOOP (neutral arc): the message is unchanged. Error undetected.
- If the flip produces a NEIGHBOR: the receiver detects an error (different class received).

The self-loop probability = SL_orbits / T_n ~ 1 - 2E/T ~ 1/(2m) at large n.

So for large n: a random arc flip almost always changes the class. The "error rate" (probability of detectable error) approaches 1. This is a NOISY channel — almost every transmission corrupts the message.

### Decoding: Finding the Original Class

Given a corrupted tournament T' (one arc flipped from T), can we recover the original class?

If T and T' are in DIFFERENT classes: we know an error occurred. But we don't know WHICH arc was flipped (there could be multiple arcs whose flip connects these two classes).

The META-GRAPH DISTANCE between two classes = minimum number of arc flips to transform one into the other. This is the Hamming distance in the quotient space.

### The Weight = Hamiltonian Path Count as Decoding Metric

H(T) measures how "ordered" the tournament is:
- H = 1: perfectly ordered (transitive, unique path)
- H = max: maximally cyclic (regular/Paley, many paths)

In the coding interpretation: H is a RELIABILITY METRIC. High H = high redundancy (many Hamiltonian paths = many valid orderings). Low H = low redundancy (unique ordering = fragile).

The transitive tournament (H=1) is the most FRAGILE codeword: it has exactly one decoding (ordering). Any error destroys the unique structure.

The Paley tournament (H=max) is the most ROBUST: it has many valid orderings. Errors are less likely to destroy ALL of them.

---

## The Spectral Gap as Channel Capacity

The metagraph's spectral gap lambda_2 controls:
1. **Mixing time** of random walks = how quickly random arc flips explore all classes
2. **Expansion** = how well-connected the code is
3. **Shannon capacity** of the arc-flip channel

The spectral gap ~ 2/n (from the project's data). This means:
- Mixing time ~ n/2 (logarithmic in V_n but linear in n)
- The channel capacity ~ 1 - O(1/n) bits per arc flip
- At large n: the channel is almost perfect (each flip gives almost 1 bit of information about the class)

---

## The Deepest Interpretation

The tournament metagraph G_n = Q_m / S_n is the **moduli space of binary codes on the complete graph K_n**. Each tournament is a "direction assignment" to edges — a binary code on the graph. Isomorphic tournaments give equivalent codes.

The metagraph structure reveals:
- **How many distinct codes exist** (V_n = A000568)
- **How codes relate to each other** (edge structure)
- **Which codes are robust** (high H = many decodings)
- **Which codes are fragile** (low H = few decodings)
- **Which codes are self-dual** (SC tournaments)
- **How error correction works** (neutral arcs = undetectable errors)

And our main formula E(G_n) ~ (T_n - twin_SL)/2 tells us: **the number of distinct single-bit-error channels between codes** is controlled by the Burnside representation of S_n, with the twin mechanism (interchangeable coordinates) as the dominant self-correction mechanism.

Tournament theory IS coding theory on complete graphs, with Hamiltonian paths as the decoding metric.
