# Chan–Pak and Rule 30: how far the campaign has come (evaluation, 2026-11-14)

mac-mini-2026-11-13-cuts. **Both problems are OPEN.** This is a map of what is proved, what is only verified, which
proof architectures are alive, and which are closed, with the evidence for each verdict. Detailed records:
- [Chan–Pak waves 1–17](../chanpak_waves_20261113/README.md);
- [Busch staircase and Fermat-block cuts](../busch_staircase_fermat_block_cuts_20261113/README.md);
- [Rule 30 continued fraction and 2-adic profile](../rule30_continued_fraction_2adic_20261113/README.md).

## Chan–Pak: every odd N except 7 and 21 is H(T) for some tournament T

### What is proved

| result | where |
|---|---|
| H(T) is odd (Rédei); H is multiplicative over strong components | classical |
| Fermat blocks: H(R_j) = 2^j + 1; the gluing law N = H(A) + 2^j·L(A, v) | Session 324; busch §2 |
| The port sets Λ(1) = {2^t} and Λ(3) = {3·2^t}, and the exact reason 7 and 21 fail | busch §2 |
| Frame identity H = Σ 2^α 3^β for two-hub tournaments; the nested-frame lemma | waves 6, 8 |
| Every n is a sum of 3-smooth numbers none dividing another (the decoder lemma) | two_three_problems note |
| One hub: the hub port equals the binary value, so H(x1^j0) = H(x) + 2^j(1 + X(x)) | busch §3 |
| Base-b gaps: m = 3 is a gap in every base, and m = 10 in every base except 4 and 8 | waves 2, 3 |
| Terrace corridors: det K_q = −T(−q−3), rank one exactly at q = 0 and 13 | wave 16 |
| Conditional: the Staircase-Run Conjecture implies the exceptions are exactly {7, 21} | busch §5 |

### What is only verified

| family | coverage |
|---|---|
| one-apex joins (THM-4875 (D)) | **every odd N < 2^32** except 7 and 21 |
| two-hub words, length ≤ 18 | every odd N < **226,079,703** except 7 and 21 (new today; ratio X(L)/3^L = 0.584) |
| interval-hub words (u2's out-set an interval), m ≤ 23 | every odd N < **26,610,535** except 7 and 21 (new today; binary-scale horizon 3–5·2^m) |
| O/Z/J root-change words (KNOWN: THM-4904, HYP-9494 (D)) | every odd N in (2647, 10^9]; the wave-19 decoder adds random samples near 10^12 and 10^15 |
| bounded runs ≤ 5, the 63 matrices of HYP-9677 | every odd N in (10,391, 2^24] (wave 20; the record was 2^19); runs ≤ 4 (31 matrices): every odd N in (130,247, 2^25] |
| Fermat-block cuts of order-≤9 sources | every odd N in [69, 398,129] |
| span-2 ribbons | every odd N ≤ 100,000 except 7, 21 and 47 |

### Proof architectures

| route | status | decisive evidence |
|---|---|---|
| local ternary or binary recursion (one vertex added gives 3n+{0,2,4} or 2n±1) | **closed** | P18 per configuration; wave 17 per value (323/1231) |
| affine numeration in one register (a common eigenvector) | **closed** | no common left eigenvector for λ = 2, 3 (wave 17) |
| decoder heredity (closure, shortlex, bounded backtracking, residue steering) | **closed** | P25, P30, wave 16, P33 |
| multiplicative closure alone | **closed** | Session 324 countermodels: closed sets with extra holes |
| fewest-triangle terraces alone | **closed** | density O(X^0.98), s* = 0.9636 (wave 14) |
| two-hub horizon law X(L) ≥ c·3^L | **alive** | 0.487, 0.515, 0.569, 0.584 for L = 14, 16, 17, 18 |
| Shallow-Cut Conjecture and port-wedge induction | **alive** | order-9 sources reach 398,129; the §4 transfer closes linearly on (H, D, h, E, S) |
| rational base 3/2 frame digits (P34) | **new, untested** | nested frames weigh 2^(A−g)·3^(g−1) |
| interval-hub numeration (P35/P36) | alive, superseded | complete through 26,610,535; but the O/Z/J root-change grammar (KNOWN) is smaller: 4-dimensional state, 3 letters |
| **cofiniteness of the O/Z/J semigroup orbit (HYP-9640, HYP-9677)** | **the repository's main formulation** | one-chain equations are near-products (wave 19); a proof must use chains of length ~log N; 14–16 matrices already appear to suffice beyond 130,247, and the run cap has a trichotomy with critical cap 2–3 (wave 20) |
| dichotomy: a third exception must belong to an arithmetic family | **alive (heuristic)** | every candidate family checked so far is complete |

### Assessment

Every reduction found so far trades the conjecture for an equally hard statement: the range of an N-rational series
(the 17-dimensional two-hub automaton) or a two-dimensional numeration (ports against counts). The computations are
overwhelming (2^32), and the obstruction at 7 and 21 is completely understood, as a port-set phenomenon. The
missing ingredient is a **self-similar covering lemma**: a family whose state space is small enough that
"intervals of attained values at size X map onto intervals at size cX" can be proved by induction with a finite
base. The one-hub family is such a system, with state (H, Y = 1 + X) and block moves (H, Y) ↦ (H + 2^j Y, 2^(j+1)Y − 1) for j ≥ 1, but
it reaches only density 1/2 (HYP-9469). Two hubs reach density 1 but have no low-dimensional invariant.

**Update (wave 18).** The smallest complete system has now been found. In the interval-hub family the second hub
beats exactly one interval of the spine. It covers every odd N < 26,610,535 except 7 and 21, at binary scale (m ≈
log2 N − 2). For each interval its count is a quadratic pseudo-Boolean function of the first hub's bits. That
follows from the bi-quadratic frame law, which is proved for all two-hub words. The part of the word before the
interval acts through a 6-dimensional state. So the next target is concrete: closed forms for the scaled one-hub
numerals that each interval produces, then a carry decoder proving that their union is complete (P36).

**Update (wave 19).** The repository search caught up with the construction. The "caterpillar" family behind the
new decoder is the O/Z/J root-change grammar (THM-4904, THM-5148, HYP-9494): 3 letters and a 4-dimensional root
state. It is smaller than the interval-hub family, and complete through 10^9 apart from 29 values up to 2647. The
wave-19 decoder carries the reader's factors into the source equation and ends in one exact check. Random targets
near 10^12 and 10^15 all decode. One chain's equation is a near-product of the source's high bits with its reversed
low bits, so each chain covers density about 1/(A + B). The open step is the union over chains whose length grows
with log N, as in HYP-9677.

## Rule 30: the three prize problems

### Repository program (36 theorems, THM-3456 … THM-4210)

- **Reductions.**
  - The trace–boundary bijection for left-permutive rules (THM-3456).
  - Periodicity from a given onset and period is exactly the vanishing of one inverse boundary word (THM-4048).
  - Balance is equivalent to discrepancy on dyadic shells (THM-4048).
  - The dyadic block-current Cartier tree (THM-4210).
- **Tame pieces.**
  - Fixed left-front columns are eventually dyadic-periodic through an affine monodromy clock (THM-4047), with
    infinitely many eventual-zero columns (THM-4065).
  - The right edge is a 2-adic odometer (THM-3458).
- **Barrier.** The center is the moving diagonal c_t = ℓ_t(t), and no fixed-column onset bound reaches it.

### This session

- Prize 1 is irrationality of α = Σ c_t X^t ∈ F_2[[X]], equivalently of β = Σ c_t 2^t ∈ Z_2.
- Through 2^22 bits the linear complexity is exactly n/2 with |2L − n| ≤ 25, and the partial quotients are
  geometric. So any eventual period has onset + period ≥ 2,097,152, and no LFSR shorter than n/2 − 12.5 produces
  any prefix.
- Through 2^16 bits the 2-adic complexity is n/2 ± 10.5. The same holds for the 3x+1 conjugate x_30, whose Collatz
  parity record is the center column; it would be rational if the column were eventually periodic.
- The shortest register with any nonlinear feedback (maximum-order complexity) has length 44 = 2 log2 N at 2^22. No
  44-bit window repeats, so onset + period ≥ 4,194,261.
- All 32,767 kernel elements with k ≤ 14 are distinct.
- Controls: Python's Mersenne Twister is caught at linear complexity 19,938, and a periodic control stops at 10.

### Assessment

There is no route through the moving diagonal yet. The irrationality form puts Prize 1 inside transcendence theory
in positive characteristic. A proof that the center's 2-kernel is infinite would settle it, and the kernel has the
recursive structure that THM-4210's tree describes. That is the most concrete proof target. Prize 3 is now
quantified in the LFSR and FCSR models; the quadratic-register model is next.

## Ranked next moves

1. Prove cofiniteness, or density one, for a sub-semigroup of the O/Z/J grammar whose chains grow like log N
   (HYP-9677); the wave-19 decoder supplies exact instances (Chan–Pak).
2. Test P34: realize rational base-3/2 digit vectors with nested two-hub frames (Chan–Pak).
3. Look for a provably infinite 2-kernel of the Rule 30 center through THM-4210's tree (Rule 30, Prize 1).
4. Measure the quadratic span of the Rule 30 center (Rule 30, Prize 3).
5. Extend the two-hub horizon to L = 19 or 20 with a sharded bitmap (P21).
