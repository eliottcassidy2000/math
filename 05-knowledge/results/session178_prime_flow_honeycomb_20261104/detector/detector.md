# Exact local Collatz detection on a space-filling triangular lattice

Use the unaccelerated map C(n)=n/2 for even n and C(n)=3n+1 for odd n, on positive integers. The triangular lattice (hexagonally arranged vertices, degree six) is not the degree-three honeycomb graph. All claims below use the former unless explicitly described as a digit strip.

## 1. A fully specified space-filling labeling

Use axial coordinates (q,s) in Z², neighbors given by the six displacements
(0,-1), (-1,0), (-1,1), (0,1), (1,0), (1,-1).
The radius is r=max(|q|,|s|,|q+s|). Label the origin 1. Radius r≥1 starts with label B(r)=2+3r(r-1) at (r,0). Traverse the six directions above for r steps each, including each side's starting corner but excluding its ending corner. A site j steps along side i, with 0≤i≤5 and 0≤j<r, has label B(r)+ir+j. This is a bijection N⁺→Z²: shell r has exactly 6r vertices and the labels through that shell end at 1+3r(r+1). The script supplies both exact directions of this bijection.

Given two sites, decode a,b. The complete directed-edge oracle is
- a→b iff a=2b, or a odd and b=3a+1;
- b→a iff b=2a, or b odd and a=3b+1.
Testing both orientations determines the direction without assuming the larger label is the source. For this unaccelerated map, a positive distinct pair cannot satisfy both orientations. No arrow is inferred merely from proximity.

For this particular spiral, the entire nearest-neighbor Collatz edge set is exactly
1→4, 2→1, 3→10, 28→14, 30→15, 32→16, 34→17, 36→18.
This is an exhaustive finite calculation, not an extrapolation: if the larger-labeled endpoint is in shell R≥2, its neighbor has label at least A=2+3(R-1)(R-2), while the larger label is at most 1+3R(R+1). Thus the label gap is at most 12R-7. Every Collatz edge has larger-minus-smaller ≥ smaller, so A≤12R-7, or 3R²-21R+15≤0. Hence R≤6. The implementation exhausts all sites through shell 6 (and their neighbors); shells 30 and 60 provide redundant checks.

This shows why an actual space-filling arithmetic visualization should expose long edges or wire routes rather than quietly delete them. It does not disprove every possible labeling: this finite-neighbor conclusion is specific to the stated spiral.

## 2. Constructive constant-state local wire checker

The useful replacement for an impossible fixed-size numerical lookup is a finite-state arithmetic strip. Let x_i,y_i be the least-significant-first binary digits of positive integers x,y, padded with zeros to a common finite length L at least their bit lengths. To certify y=A x+B, impose at each site

A x_i+c_i = y_i+2c_(i+1),    c_0=B,    c_L=0.

For y=2x, take A=2, B=0 and carry alphabet {0,1}. For y=3x+1, take A=3, B=1 and carry alphabet {0,1,2}. The odd branch additionally requires x_0=1. Each transition depends only on one pair of bits and its incoming/outgoing carry; all alphabets and tile types are fixed independently of x,y.

Soundness is exact: multiplying by 2^i and summing telescopes the carries, giving A x+B=y. Completeness follows from the standard recurrence c_(i+1)=floor((A x_i+c_i)/2); acceptance checks the corresponding output bit and final zero carry. Bounded carry sets follow directly from A≤3 and bit inputs. This is also a completely specified deterministic finite-state transducer, so no appeal to an unspecified numerical approximation is involved.

To test directed n→m, run either the doubling checker with x=m,y=n, or the 3x+1 checker with x=n,y=m and n odd. Run these with n,m exchanged for the reverse direction. These are edge-presence and direction tests, not Collatz convergence tests.

Example: 13→40 uses x bits (1,0,1,1,0,0), y bits (0,0,0,1,0,1), and carries (1,2,1,2,2,1,0). Every cell satisfies 3x_i+c_i=y_i+2c_(i+1).

Embed the arithmetic strip in triangular coordinates by using sites (i,0), i=0,…,L, which are adjacent. A tile stores its two bit symbols, incoming/outgoing carries, and a branch tag. The local tile constraint checks arithmetic and neighboring tiles check that both carry symbols and branch tags match, keeping A uniform across the entire strip. Start/end marker tiles impose the boundary conditions, including x_0=1 at the odd-branch start. Thus a valid complete strip certifies an edge using O(log max(n,m)) cells and a constant finite alphabet. A reader moving along the strip needs constant carry state and O(log max(n,m)) steps. Parallel local verification of a supplied certificate does not mean one site can collect the verdict in constant time: end-to-end acceptance still needs communication.

Crucial interface condition: the bits must actually encode the endpoint labels. Copying/checking endpoint data requires wires or a larger certificate. An arbitrary supplied strip certifies its own encoded integers, not automatically whatever remote lattice vertices someone points to. The construction can attach a checked binary representation to each selected endpoint, but this is not a one-cell-per-integer finite-alphabet bijection storing all digits at zero cost. The original hex spiral remains the injective global labeling; digit strips are expanded checking gadgets.

A fixed modulus alone is only a rejection filter. For any modulus M, the true pair (n,2n) and false pair (n,2n+M) have the same ordered residue pair mod M. Therefore residues without magnitude or longer digit information cannot be an exact doubling oracle over all integer pairs. This does not preclude computing exact arithmetic from unbounded coordinate labels.

## 3. The precise connection to prime-factor flow

On even steps, prime exponent coordinates have an exact unit move: v₂(C(n))=v₂(n)-1, and v_p(C(n))=v_p(n) for all odd primes p. This is the natural factor-space transport for division by two.

For odd n, set k=v₂(3n+1). The 3x+1 transducer produces k initial zero output bits followed by a one. Reading those bits gives the entire immediate sequence of halving steps, and the odd-to-odd map U(n)=(3n+1)/2^k is obtained by dropping that zero prefix. The finite transducer can emit/count k in unary along a wire; an unbounded numeric k cannot be kept in a fixed finite internal state.

The valuation is governed by exact congruences: k≥j iff n≡-3⁻¹ (mod 2^j), and k=j iff n≡(2^j-1)3⁻¹ (mod 2^(j+1)). More generally, for p≠3, v_p(3n+1)≥j iff n≡-3⁻¹ (mod p^j), while v₃(3n+1)=0. These are valid residue detectors at each chosen depth. They do not produce a uniformly bounded-depth detector because valuations are unbounded, and the additive +1 means odd-prime exponents are not independently transported from the factorization of n.

## Verification

Run `python detector.py`. The dependency-free tests verify 100,000 label→coordinate→label cases, 40,401 coordinate→label→coordinate cases, 250,000 ordered candidate pairs against direct arithmetic, 100,000 actual Collatz transitions and carry bounds, and the listed spiral-neighbor edge set. The finite-shell argument above supplies the global justification beyond tested ranges. These are elementary constructive results and computational checks; no novelty or resolution of Collatz is claimed.

## Prior-work boundary and cycle convention

The session survey identifies earlier repository work in session 140's `principles/PRINCIPLES.md` on triangular-site conjugacy, shell-radius growth, and geometric jump bounds, and in `collatz_hex/PROOF.md` on binary supply and residue nonclosure after inverse division. Those broad themes are prior work, not new claims here. The additions in this note are the explicitly specified spiral convention with its complete eight-edge list, the finite-state carry certificate with its exact interface and communication costs, and the executable paired verification. The precise repository links should be resolved by the integrating session, which supplied these source identifiers; this note does not claim an independent literature novelty search.

The standard positive Collatz cycle is 1→4→2→1. In this particular spiral, 1→4 and 2→1 are nearest-neighbor arrows, but 4→2 is not. Thus the eight-edge list neither omits a reverse arrow under a mutual-edge convention nor claims to geometrically realize the whole cycle. It uses directed edges of the unaccelerated C only. Different rotations, origins, shell orderings, or shortened Collatz conventions can change the list.
