# One-apex modular synthesis and the exact binary-word interface

Status: all-order elementary theorem, with independent direct Hamiltonian-path DP controls. This is a candidate for canon insertion, not a claim of literature priority. It does not settle global spectrum completeness.

## Main theorem

For every positive odd integer q and every residue r modulo q, there is a strongly connected tournament T, obtained by adjoining one apex to a transitive chain, such that

    H(T) = r (mod q),             |V(T)| <= 3q^2.

Consequently this thin construction class meets every odd residue modulo 2q. For any prescribed finite list of odd moduli and compatible residues, one strong one-apex tournament meets them simultaneously, by applying the theorem to their least common multiple. No collection of odd-modulus congruence restrictions alone can account for missing exact odd Hamiltonian-path counts.

This does NOT assert arbitrary prescribed 2-adic residues, exact realization of every permitted integer, or that the witness value is near its target residue.

## 1. Exact interface, without the odd-cycle formula

For a binary word w=b_1...b_m, let T(w) have apex a, chain 1->2->...->m, and a->i iff b_i=1. Define s(w) as the binary integer represented by w, p(w) as its final bit (zero for the empty word), and h(w)=H(T(w)). The empty state is (0,0,1). Appending bit b gives exactly

    (s,p,h) -> (2s+b, b, h+(1-b)(s+p)).                 (1)

Proof. A newly appended chain vertex z dominates no old chain vertex. If b=1, it is a global sink, so every Hamiltonian path ends at z, and h is unchanged.

If b=0, partition new paths by whether z is final. Those ending at z correspond to old paths not ending at a. There is exactly 1-p old path ending at a: it must traverse the entire transitive chain first. (For the empty chain this is also one.) Thus this part contributes h-(1-p).

Every other new path must use z->a. Its old-chain vertices preceding z form an increasing sequence, and those following a form another increasing sequence. Choose the latter subset R. The empty R gives one path. If R is nonempty, its minimum i must satisfy b_i=1, and every vertex greater than i can independently belong to R or its complement. Therefore there are

    1 + sum_{i:b_i=1} 2^(m-i) = 1+s

such paths. The total is h-(1-p)+1+s=h+s+p. This proves (1) for every word, without OCF or an unproved enumeration identity.

The familiar one-apex odd-cycle formula is recovered independently:

    h(w)=1+2 sum_{i<j, b_i=1,b_j=0} weight(j-i),
    weight(1)=1, weight(d)=2^(d-2) for d>=2.             (2)

It follows either by summing (1), separating adjacent pairs, or by counting the odd cycles through the apex. Our tests compare (1), (2), and a completely separate subset DP.

Important interface boundary: this is three scalar registers, but s carries m bits. It is constant-dimensional, not a bounded-information exact state and not a compression of arbitrary tournaments. Reduction modulo fixed q is genuinely finite-state.

## 2. Three append gates

Start with the word 1. Thereafter append only blocks 1, 01, and 001, so the current final bit remains p=1. Set x=s+1. On (x,h), these blocks act as

    B (append 1):     (x,h) -> (2x,h),
    A (append 01):    (x,h) -> (4x-2,h+x),
    D (append 001):   (x,h) -> (8x-6,h+3x-2).           (3)

These are exact integer affine gates, derived by composing (1), not fitted recurrences. For odd q all three are permutations of (Z/qZ)^2. Therefore the semigroup they generate is a group: each inverse is a positive power of the corresponding permutation. Negative powers used below are proof devices only; a forward breadth-first search constructs actual positive append words.

Use functional composition from right to left. In the generated group put

    C = A B^(-2):     (x,h) -> (x-2,h+x/4),
    F = D B^(-3):     (x,h) -> (x-6,h+3x/8-2).         (4)

Fractions mean modular inverses, which exist because q is odd. Two combinations give

    Z = F C F^(-1) C^(-1): (x,h) -> (x,h+3/4),
    E = C^3 F^(-1):        (x,h) -> (x,h+(3x+22)/8).   (5)

For clarity, C^3(x,h)=(x-6,h+3x/4-3/2), and
F^(-1)(x,h)=(x+6,h-3x/8-1/4). These yield the constant 22 in (5).

At a fixed x, both E and Z leave x unchanged. Hence E^2 Z^(-x) translates h by 22/4. The Bezout identity

    15*3 - 2*22 = 1

shows that Z^15 (E^2 Z^(-x))^(-2) translates h by 1/4 at that fixed x. Since 4 is a unit modulo odd q, every h-coordinate is attainable while holding x fixed. C changes x by -2, which runs through all x modulo q. The group therefore acts transitively on (Z/qZ)^2. This argument works also when 3 divides q; retaining the 001 gate and the constant term in E is essential for that case.

## 3. Strong tournament realization and size bound

The initial word 1 has state (x,h)=(2,1). By transitivity, a positive block word reaches the state

    (x,h)=(2,r-2) modulo q.

Append a final bit 0. Formula (1) adds s+p=x=2 to h, so the resulting count is r modulo q.

The full word begins with 1 and ends with 0. Thus

    a -> 1 -> 2 -> ... -> m -> a

is a directed Hamiltonian cycle, making the tournament strongly connected.

There are q^2 states in the directed graph with the three append gates as edges. A shortest route between reachable states has at most q^2-1 edges. Each block has at most three bits; the initial and final bits add two. Thus m<=3(q^2-1)+2 and |V(T)|=m+1<=3q^2. For q=1 the word 10 is the three-cycle and meets the same bound. This finishes the proof.

A constructive algorithm is BFS in q^2 states, with three outgoing edges each. It takes O(q^2) arithmetic operations and O(q^2) storage in the standard modular-arithmetic model; reconstructing one witness costs its word length. Exact H may have O(q^2) bits, so modular synthesis does not claim constant bit complexity.

## 4. Decisive finite tests

Companion script: modular_apex.py. Frozen transcript: modular_apex.out.

* Every binary word of lengths 0 through 9: 1,023 comparisons of recurrence, independent cycle sum, and direct Hamiltonian-path subset DP, up to tournament order ten.
* Every odd modulus 1 through 101: BFS reaches exactly q^2 gate states, constructs each target residue, checks the exact integer H against (2), and checks the 3q^2 order bound.
* Every state modulo each odd q=3,...,25: direct verification of both affine identities (5).
* All checks raise explicit exceptions, so optimized Python retains them.

For q=7, one strong residue-zero witness is word 1001101010, with H=413=7*59, order 11. This is an actual 7-divisible strong construction, not a construction with H=7. Other BFS examples are recorded in the transcript. The theorem's value is all-modulus coverage, not efficiency of these examples.

## 5. Targeted archaeology and what it changes

Sources below are precise repo paths, read directly. 01-canon/MISTAKES.md was not read or scanned.

1. `01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`: recent all-order Fibonacci fan gives H=2^k-1 for k>=4 with three back arcs. Its known one-apex pair formula supplied the construction direction, but the direct path-partition proof of (1) removes its 1-not-in-A convention and all boundary ambiguity. The new modular synthesis uses arbitrary apex words, not a claim that three back arcs suffice for arbitrary residues.
2. `01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`: exact deletion interfaces and two-prime-lane reduction. Spectrum completeness remains equivalent to realizing all ordinary prime atoms p!=7 and all exceptional atoms 7p for odd prime p!=3, with the finite carry atoms. Modular coverage is not either exact prime-lane theorem.
3. `01-canon/theorems/THM-4104-selected-order-eleven-strong-ear-solid-interval.md`: exact (Start,End,Q) cut interface, finite prefix through 80,405, and explicitly open interval propagation. Our interface is much smaller only because the chain-plus-apex family is highly restricted. It does not replace Q for general ears.
4. `01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`: full endpoint jet is needed for general cyclic substitution. This reinforces the operation-specific interface distinction: one-apex append operations close on (s,p,h), but arbitrary substitution does not.
5. `07-reflections/the-source-sink-arc-is-the-pivot-s530.md`: a single reversed source-sink arc gives 1+2^(n-2), a strongness-preserving construction seed. Here first-bit 1 / last-bit 0 is the explicit Hamilton-cycle guard.
6. `07-reflections/the-apex-recursion-modular-tree-recursive-H-and-recursive-carving-s531.md`: disjoint reverse-arc blocks multiply but nested ones couple. The exact affine gates make one such nested coupling computable and invertible modulo odd q. Its broad LRC prose is not used as a theorem.
7. `07-reflections/tournaments-as-proof-engines-a-generative-catalog.md`: the SCC-factorization sieve and realizability-hole toolkit are useful, but this new theorem marks a limit of a residue-only sieve. The source's old finite claims about 21 are superseded by the global exclusion source below.
8. `07-reflections/H-impossibility-spectrum-and-its-propagation-s599.md`: historical candidate pattern 7*3^k is obsolete already at 63, as the later Fibonacci fan demonstrates. The reflection itself says its evidence was finite/sampled. Do not promote its early wording to a global obstruction.
9. Exact permanent-gap authority: `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`, as explicitly quoted and relied upon by THM-4094 and THM-4748. Global gaps 7 and 21 are known. Other missing sampled values are not automatically global gaps.

## 6. Honest research boundary

The theorem establishes that any proof of a new global forbidden odd value cannot consist merely of showing its residue fails to occur modulo some odd integer: every such residue already occurs among strong one-apex tournaments. Exact equality remains a size-sensitive, non-modular issue. The next construction question is whether suitably bounded lifts of these residue targets can be controlled, not whether the finite modular state graph is connected. No such bounded-lift theorem is claimed here.
