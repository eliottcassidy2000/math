# Placing prescribed primes with an exact binary inverse: what dimensions 8 and 24 do and do not buy

**Proved:** an absolute inverse for the established fixed-root binary continuation family, its address-volume constant, an obstruction to stable independent digits in that family, and a contrasting interleaved scalar repair. **Witnesses:** two 6-vertex sources followed by literal address words reach the two predetermined primes 1,000,000,007 and 1,000,000,009, on 32 and 31 vertices respectively. **Open:** selection of a physical source for every odd prime other than 7; no claim about new coverage bounds for the unrestricted conjecture.

## 1. Exact inversion, not a congruence search

At the retained root let (g=H(T-r)), and let (E,O) be the integer weighted even/odd path-start coordinates used in the [fixed-63 continuation](../fixed63_missing_lane_20261116/README.md). A weight is the Hamiltonian count of the complement of the root-start path. Start with (O=E+g), (E>0). Write (A=2E) and (J=H-A-g). Read the bits of an address (z>0) from low to high: 0 adds an ordinary sink (S), 1 a near-sink beating only the root (N). The last bit is 1. Address 0 is the empty word. Trailing S's, if present in a continuation, do not affect its final count. The established transfer is
[
S:(H,g,E,O)mapsto(H,g,E+O,E+O),qquad
N:(H,g,E,O)mapsto(H+2E,g,E+O-g,E+O).
]
As in the [all-moduli construction](../forty_two_7p_plane_20261116/ALL_MODULI.md), the exact output is (A(z+1)+g f(z)+J), with (f(0)=1, f(2u+1)=f(u), f(2u)=f(u)+u+(umod2)) for positive (u). The general (J) is an additive constant: the operations on (E,O,g) do not read (H).

Let (q=2^k), (kge1), (z=qb+r), (0le r<q), and let \(\operatorname{rev}_k\) reverse the (k) bits including leading zeros. Then
[
f(qb+r)= f(b)+f(r)-1+b(q-1-\operatorname{rev}_k(r))
                 +[r<q/2](bmod2).
]
Proof: at every zero low bit (i) the contribution from the quotient (b) is (2^{k-i-1}b); only at (i=k-1) does the parity of that high quotient also enter. All other low contributions are exactly (f(r)-1); stripping all (k) bits exposes (f(b)).

Define (K_r=Aq+g(q-1-\operatorname{rev}_k(r))) and (C_r=A(r+1)+g(f(r)-1)), and let (X=P-J) for prescribed count (P). The physical high digits must lie in
[
\max{0,\lceil(X-C_r-2g)/(K_r+g)\rceil}
\ \le b\le\ \lfloor(X-C_r-g)/K_r\rfloor .
]
Indeed (1le f(b)le b+1); the upper bound follows by induction, using (f(b)le b) for odd (b). The final parity term is at most 1. Check the exact identity on every surviving (b). At (q) of order (sqrt X), each interval has (O_{A,g}(1)) candidates. Thus this is a complete inverse in (O_{A,g}(sqrt Xlog X)) elementary arithmetic steps, versus checking (O(X)) addresses. It returns every address or a certificate of absence from this *particular* source. It is not polynomial in bit length.

For example these inputs were **prescribed**, not obtained by testing output for primality:

| (P) | Outgoing bit rows of source | Root | Address |
|---:|---|---:|---:|
| 30011 | 2,4,9,19,39,15 | 2 | 1561 |
| 30013 | 2,4,9,19,39,15 | 1 | 2363 |
| 1000000007 | 2,4,9,19,39,15 | 2 | 47827936 |
| 1000000009 | 2,12,25,49,35,7 | 1 | 30268243 |

The fixed H63 source with (A=52,g=11,J=0) misses all four values. The larger pair has different terminal profiles, not a common continuation. Run `python decoder.py` and `python -c 'from seed_sweep import literal_checks; literal_checks()'`. Run `python seed_sweep.py` to reproduce the target-first search through degree-ordered isomorphism representatives up to order seven. It found six distinct eligible source triples for each large prime; the code verifies literal sources by independent subset DP and replays the full transfer words exactly.

## 2. A volume calculation and a covering calculation are different

For fixed (A,g,J), the number of binary **addresses** (with multiplicity, not distinct values) giving count at most (X) is
[
X\frac{\log(1+g/A)}{g}+O_{A,g,J}(\sqrt X).
]
For a residue (r) the exact formula sandwiches the eligible (b)'s between the linear bounds of §1, now interpreting the lower bound as a sufficient upper cutoff for membership in ({Hle X}). It gives (X/K_r+O(1+X/q^2)) addresses. Summing, (K_r/q) runs over the evenly spaced values (A+g j/q). Their reciprocal average is the Riemann sum for (int_0^1 dt/(A+gt)), with error (O(1/q)). Choose (q\asymp\sqrt X). The three banks ( (A,g)=(28,13),(30,13),(52,11)) together have address mass approximately (0.0744731X), far below the (X/2) odd integers. This is an unconditional obstruction to covering **all odds** with those banks using these words; a linear count bound alone does not imply any prime obstruction.

Conversely, for every large (P) one fixed physical source comes within (O_{A,g,J}(\sqrt P)). Let (z=q b-1) with (b\ge1), so (f(z)=f(b-1)); then (H=Aq b+g f(b-1)+J). Choose (q\asymp\sqrt P), (b=\lfloor(P-J)/(Aq)\rfloor). Rounding the center costs at most (Aq), while the deposit is at most (gb). Closeness and address volume thus coexist with holes. A genuine coverage certificate must include interval overlap and a decoder that pays the final error without breaking the prescribed profile.

## 3. 8 and 24: real mathematical bridges, with the loss of information marked

The unit-ball volume is (V_d=\pi^{d/2}/\Gamma(d/2+1)); the ratio (V_{d+2}/V_d=2\pi/(d+2)) peaks at dimension 6 along the even subsequence, not at 8 or 24. In all integer dimensions the maximum is at 5 ((V_5/V_6=16/(5\pi)>1)). Optimal packing fractions, by contrast, are (\pi^4/384) for (E_8) and (\pi^{12}/12!) for Leech, from their covolumes 1 and squared minima 2 and 4. Those fractions are both below one; a packing optimum is not a no-holes result. See [Viazovska](https://arxiv.org/abs/1603.04246) and [Cohn–Kumar–Miller–Radchenko–Viazovska](https://annals.math.princeton.edu/2017/185-3/p08).

* [THM-481](../../../01-canon/theorems/THM-481-paley-gauge-generates-gleason.md) gives the actual Paley gauge route from the 7- and 23-tournaments to the length-8 and length-24 extended codes. The [archival reflection and errata](../../../07-reflections/sphere-packing-eight-twentyfour.md) explicitly retract a mechanism from the gaps 7/21 to the packing dimensions. The perfect codes for correction are the **punctured** [7,4,3] Hamming and [23,12,7] Golay codes; the extended codes do not supply a perfect packing of Hamming balls.
* A literal transfer worth pursuing is an absolute chart (H=C+2\sum_{i=0}^{d-1}2^i b_i). A physical center for each Hamming codeword and any one bit repair would cover 128 consecutive odds; the Golay analogue with up to three repairs would cover (2^{23}). The ball sizes are 8 and 2048 respectively, and minimum distance gives disjointness; cardinality gives **coverage** of the entire finite cube. To iterate, place the next chart and keep terminal data under control. Neither physical hypothesis follows from Paley or packing.
* A different criterion uses quadratic shells. (E_8) has every even nonnegative squared norm: repeat each coordinate of a four-square representation twice to get a vector of its (D_8) sublattice. Leech has no squared norm 2, but all even squared norms at least 4. Indeed its [theta coefficient at half-norm (m)](https://uva.theopenscholar.com/files/ken-ono/files/006_8.pdf) is proportional to (sigma_{11}(m)-\tau(m)); for (m>1), Deligne's bound yields (|\tau(m)|\le d(m)m^{11/2}\le m^{13/2}<m^{11}\le\sigma_{11}(m)). A physical, continuation-preserving norm decoder would therefore suffice; no such decoder is asserted. This reformulates the main target as four independent square increments rather than an inventory of high counts.
* [The Leech deep-hole swap](../leech_deep_hole_swap_theorem_20261022.md) leaves a controlled cusp-form correction; it inspires keeping the interaction sidecar instead of equating the two spectra. [The Bott/Collatz marked periodicity](../collatz_bott_marked_periodicity_20261007d.md) likewise retains the exact arithmetic corner rather than assuming an eight-step physical return; equal dimension or a stable class is not a decoder. The old [Lonely Runner relation-lattice analysis](../../../07-reflections/the-density-floor-is-cohn-elkies-poisson-summation-over-the-relation-lattice-macmini-S24.md) distinguishes an exact Fourier reformulation from a pointwise positive certificate. Here that distinction is made explicit by the inverse and its volume bound.

## 4. The simplest exact code chart is obstructed, but nonlocal edits change the interaction

For literal address bits (b_i), the recurrence gives
[
f=1+\sum_{i<j}(1-b_i)b_j(2^{j-i-1}+[j=i+1]).
]
Every pair of literal bit toggles (0\to1) has nonzero mixed difference (-g(2^{j-i-1}+[j=i+1])). More generally for two nonzero edits supported in separated, **nonadjacent** blocks (I<J), the mixed response is
[
-g\Big(\sum_{i\in I}\delta_i2^{-i}\Big)\Big(\sum_{j\in J}\epsilon_j2^{j-1}\Big).
]
Neither signed binary sum can vanish, as its largest weight exceeds the sum of smaller weights. Hence independent scalar digits on separated blocks of a single fixed-root address cannot simply be sewn together.

Actual *interleaving* can remove the mixed term. With the cyclic triangle and root 0, four addresses 289, 297, 278, 286 realize (721,723,725,727); the masks of the two edits are 8 and 55, disjoint, with the latter straddling the former. In the pair ((H,E)) their increments are ((2,-16)) and ((4,28)) throughout the square. After a common (k)-N tail their signed half-steps become (1-16(2^k-1)) and (2+28(2^k-1)). These are coprime for even (k): seven times the first plus four times the second is 15, and neither 3 nor 5 divides their gcd at an even height. Seven of the first and four of the second would cancel (E) and leave 30 in (H), **if** those independent copies could be realized. The square itself supplies only one of each.

There is a categorical obstruction to a stable code chart on one fixed source: terminal words of length (L\ge1) give (E'=2^L E+g x) with (0\le x<2^{L-1}) and each (x) read from the lower bits. Distinct lengths have disjoint bands if (2E\ge g), as they do at all four displayed sources. At a fixed length the (E') register is injective. Thus no two of those distinct addresses share the same continuation; the interleaved patch is a final scalar selector, not yet a repeatable digit.

The next constructive target is a source-changing (or nonuniform) edit with a coordinated ((\Delta H,\Delta E,\Delta g)) budget. The inverse can prescribe the exact center, but needs a bank large and structured enough to place it. Positive density, modular steering, perfect codes on abstract labels, and fixed-core scalar repairs are each weaker than that goal.
