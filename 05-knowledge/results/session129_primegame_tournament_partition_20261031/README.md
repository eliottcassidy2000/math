# Session129: exact PRIMEGAME–tournament bridges

**Status: PROVED scoped algebra and dynamical classification; FINITE-EXACT tournament census; explicitly chosen encodings.** Three verified publication checkpoints. No global open conjecture is claimed solved, and no literature-priority claim is made.

## 1. The proposed partition really has structure

The four blocks are

- X0={11/13,13/11};
- X1={55,1/7,1/17,17/91};
- X2={19/51,77/19,77/29,29/33};
- X3={78/85,23/38,95/23,15/2}.

Retain Conway's original first-eligible instruction order. With e_p=nu_p(N), define

`T=e_2+e_5+e_23`, `U=e_11+e_13+e_17`, `V=e_19+e_23+e_29`,
`Phi=2e_7+2U+3V-4T`.

[The height theorem](partition/README.md) proves that every X1 instruction lowers Phi by2, every X2 instruction raises it by1, and X0/X3 leave it unchanged. This is the unique primitive character supported on the ten program primes with these block-constancy and horizontal conditions, up to sign. X0 consists exactly of the horizontal moves that change no data counter. Thus the user's split can be recovered by an additive grading plus the control/data distinction.

For every executed word, `Delta Phi=#X2-2#X1`. The nonnegative physical height `Psi=Phi+4T` gives the global prefix restriction

`2#X1-#X2-4#(55 instruction) <= Psi(start)`.

This is a genuine elimination/transition budget. It does not discard unbounded register values or assert a finite-state simulation.

## 2. A weighted triangle and a lossless rational encoding

[Multiplying the blocks](bridge_audit/README.md) gives1 and

`P=55/637, Q=539/153, R=585/34`.

Control-prime cancellation turns P,Q,R into the arcs13->11,17->11,17->13: an actual transitive tournament, with orientations determined by signed valuations. Its data holonomy `PR/Q=2025/4802` is nontrivial. An integral chain map explains exactly which boundaries and formal cycles are preserved.

The three products are multiplicatively independent, and their valuation image is a primitive direct summand. For any integers(a,b,c), put z=P^a Q^b R^c. Then

`a=nu_2(z)+nu_5(z)+nu_23(z)`,
`b=nu_11(z)-a`,
`c=nu_13(z)+a`.

This is an explicit integer decoder. However, product integrality is not executability:34*(585/34)=585, while none of the24 permutations of X3 executes from34. The control graph and its gains need guards and first-eligible priority to describe computation.

## 3. How the factorization can encode all ten nodes

[The tournament packet](tournaments/REPORT.md) independently reconstructs merged-node counts **2,3,10,34** at orders3,4,5,6. At order3 the reciprocal pair can model the single **complement line with two endpoints**, not a single merged node.

Order4 has three merged nodes, canonically ordered by Hamiltonian-path count1,3,5. Endpoint deletion from every marked-path tiling of an order5 node gives a three-coordinate count vector. Divide by its coordinate gcd. The ten resulting primitive vectors are

`(1,0,0), (1,1,0), (2,3,5), (0,1,0), (2,0,7),`
`(1,4,4), (1,0,10), (0,2,3), (0,5,8), (0,0,1)`.

They are distinct, so give a complete order5 classifier. This mechanism is inherited from **THM-796, three-sorted-recursive-tiling-line-node-incidence**, and independently reconstructed here rather than claimed as new.

Assign the three ordered four-fraction blocks to these three coordinates. The rational code P^a Q^b R^c is then a faithful encoding of **all ten nodes using precisely the user's products**. The assignment is declared, not forced by PRIMEGAME. The independently derived incidence table is what gives ten; factoring the fractions alone does not predict that count.

The alternative selected integer code `7^a17^b55^c` gives an exact repair of the earlier symmetry triple. The three pure-blue order5 nodes have precisely the coordinate-axis profiles, and therefore receive7,17,55. Decoding uses valuations at7,17,5,11;55 need not be prime.

These are **arithmetic node labels**, not the earlier incorrect bit-mask interpretation. Also55 remains unbalanced for the separate F=U+S divisor statistic. At order6, `1 -> strong_4 -> 1` is mixed but has an axis profile. The report proves the general self-converse sandwich mechanism. Even a complete static fingerprint is not a sufficient recursive state: two lower tilings in the same node can have different lift counts.

## 4. Complete periodic families and a global cycle restriction

The [circulation classification](partition/README.md) solves the entire nonnegative integer kernel of the ten-by-fourteen increment matrix. Its free generators are the instruction multisets

`JK`, `ABDEFG`, `A^2BCEFHKL`, `ABE^2F^2LM`.

Letters follow Conway's original A..N order. These are formal multiplicative relations until priority is imposed.

Priority greatly sharpens them. Every genuine periodic orbit is either an explicitly classified JK two-cycle, or uses only A,B,D,E,F,G, each equally often. Every non-JK period is a multiple of6 and needs at least three control tokens. In particular all cycles with at most two control tokens are completely classified.

The first six-cycle exhibited is

`4301 -> 17765 -> 16302 -> 9867 -> 8671 -> 23023 -> 4301`,

with word GBDEFA. **All least-period-six cycles** are exactly the constant rescalings of this cycle by

`X=S*2^a*11^r*13^v*19^w*23^s`,

where a,r,v,w,s are nonnegative, a*w=0, and S uses only primes outside the ten program primes. The unique G-state is4301X. A symbolic priority argument proves that every non-JK periodic orbit reaches exactly these guards and immediately returns in six steps. Thus **all finite cycles have least period2 or6**; the displayed families are complete, not merely a short-period census. A saturated3^10-state check and an independent proof audit supplement the symbolic argument.

[Uniform cycle formulas and density](periodic_density/README.md) reduce the scale restrictions to coprimality: two-cycles11M↔13M have gcd(M,58815393)=1; six-cycles use gcd(X,51765)=1 and38∤X. Exact inclusion-exclusion gives the total periodic-point density144064009693952/1969552535775825. This counts points already on cycles, not their basins.

A four- or ten-state deterministic machine cannot preserve the future block output:2^a emits exactly a consecutive X3 instructions before leaving X3, giving infinitely many continuation-distinguishable states. The inherited divisor-loop construction gives the same obstruction on seed2's reachable orbit. Finite control is useful only with its counter/phase sidecars.

## 5. Triangular and Fibonacci bridges with exact scope

[Figurate transport](figurate/THEORY.md) maps every integer stateN to the harmonic triangle `(N,N+1,N(N+1))`, whose third coordinate is2T_N. Its Fibonacci-recurrence window is `(1,N,N+1,2N+1)`. The decoder recoversN, and enabled FRACTRAN steps transport exactly. This is a generalized recurrence window; onlyN=1,2 lie on its unit-Cassini golden slice.

For an instructiona/b, the canonical seed(b,a) has quadratic norm `Q=b^2+ab-a^2`. The uniform eraser law is `Q(q,1)=2T_q-1`. Hence1/7 gives55, and1/17 gives305. The latter is not an identified bridge to the unrelated cubic continued-fraction305 in Session127.

The fourteen instructions meet exactly thirteen orbits under the stated Fibonacci-unit/Galois action. The unique paired instructions are11/13 and15/2, via `(x,y)->(x-y,2x-y)`. The reciprocal pair11/13,13/11 has different absolute norms and is not that golden pairing.

A separate explicit observer assigns the ten prime registers to the ten edges of K5 and reads exponent parity as tournament orientation. Factoring each instruction gives its exact edge-toggle mask; the observer reaches all ten merged classes. It is not autonomous:1 and4 have the same observer, but their successors55 and30 have different merged nodes. The edge assignment is a declared coordinate gauge.

[Prime-renaming covariance](encoding/README.md) explains the boundary uniformly: prime names can be changed while preserving every instruction trace, but55 then need not remain triangular or Fibonacci. Numerical identities are retained as facts about Conway's chosen encoding, not promoted to label-free dynamical theorems.

## 6. Evidence and next precise boundary

See [verification](VERIFICATION.md), [independent audit](independent/AUDIT.md), and the portable `verify_session.py`. Checks cover exact guards, integer kernels, all fixed-path cubes through order6, independent full-labelled orbit enumeration, numerical identities, negative controls, optimized execution, relocation, and source hashes.

The strongest reusable result is the separation of three levels: an additive block grading, a lossless arithmetic encoding of static tournament incidence, and an independently checked executable language. Preserving one does not automatically preserve the others. Remaining targets include classifying eventual cycle basins, all-order completeness of deletion profiles, and a genuinely natural operation-preserving map between tournament growth and PRIMEGAME beyond the explicitly chosen codecs.
