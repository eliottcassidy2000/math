# A fixed palette cannot pump primitive positive Collatz cycles

Use the ordinary shortcut map T: even x→x/2, odd x→(3x+1)/2.
A finite nonempty parity motif u of length l(u) and s(u) odd letters has the formal affine map
F_u(x)=r_u x+c_u, where r_u=3^s/2^l and c_u≥0. Since r_u≠1, write z_u=c_u/(1−r_u) for its rational fixed point. A "cycle word" below is the actual word of an integer orbit, not simply a formal composition. Its primitive length N is its least period.

**Theorem.** Fix any finite palette P of finite nonempty parity motifs. If a primitive positive cycle of length N>max_P l admits a cyclic decomposition into m nonempty blocks, each a repetition of a motif in P, then, writing M=max_P l,

    N ≤ C (log₂(3))^m (1+M+log N+log m),

where C is effective and absolute. In particular the statement also allows the motifs to vary with N: if M≤N^θ for fixed 0<θ<1, then m≥[(1−θ)log N−O(1)]/log(log₂ 3).

In particular, m=Ω_P(log N). For a fixed bound on m there are only finitely many such primitive positive cycles; this holds for independent and unbounded repetition exponents, not just for affine rays of exponents. The finite upper bound on N is effective (though far too large for a practical cycle exclusion).

## 1. Cut the circular word to bound the actual value at one block boundary

Let a_i=r_(u_i)^n_i be the multiplier of block i, and a=∏a_i=3^S/2^N. For a positive cycle a<1: the affine period identity reads (1−a)x=c with c>0, since an all-even word cannot be a positive cycle.

Rotate at a block boundary just after a minimum of the cumulative sums of log a_i, considering both endpoints. Every rotated prefix logarithm is at least log a: before wrapping it is nonnegative, and after wrapping it is log a plus a nonnegative difference. Hence all suffix multipliers A_i=a_i⋯a_m belong to (0,1], including A_(m+1)=1. The repeated motif's affine formula is a_i x+(1−a_i)z_(u_i), and thus the exact cycle value at the chosen cut is

    x₀ = Σ_i z_(u_i)(A_(i+1)−A_i)/(1−a).

Each difference is at most one in magnitude. With Z=max_P |z_u| this gives x₀≤m Z/(1−a). The carry estimate can also be centered: the differences sum to 1−a, so |x₀−z_*|≤m (max z−min z)/(2(1−a)) for the midpoint z_* of the palette. This is a block-level cyclic shift argument, akin to cutting a ballot path at its low point, but here it bounds the affine **carry**, not just the slope.

By the fixed-base case of [Matveev's theorem](https://www.mathnet.ru/eng/im314), 0<N log 2−S log 3≥c N^(−K), with effective absolute constants c,K: the logs are those of the fixed positive rationals 2,3, and the two integer coefficients have magnitude at most N. (When S=0 there is no positive cycle.) For 0<δ≤1, 1−e^(−δ)≥δ/2, and for δ>1 the same quantity is bounded below by 1−e^(−1). It follows that x₀≤C m N^K with an effective C depending on the palette. There is also a palette-free polynomial anchor at the minimum **vertex** from the standard product formula (repository note [Collatz log gap](../collatz_log_gap_20261014a.md)); the block cut is useful because its anchor is already at a permitted boundary.

## 2. A block is bounded by the 2-adic precision of its entry

Write z_u=p_u/q_u in lowest terms, q_u odd since the raw fixed-point denominator is 2^l−3^s. If a block u^n with length L=n l takes integer x to integer y, its affine fixed-point identity is

    2^L(q_u y−p_u)=3^(ns)(q_u x−p_u).

Thus 2^L divides q_u x−p_u. If x=z_u, one copy of u already returns to x (the actual orbit does follow u), so the primitive cycle length is at most l, contrary to the theorem's hypothesis. Otherwise the nonzero multiple gives 2^L≤|q_u x−p_u|≤C'_P(x+1).

For any positive x, T(x)+1≤(3/2)(x+1). Let L_i be the total number of letters from the cut through block i, with L₀=0. The entry to block i has x+1≤(3/2)^(L_(i−1))(x₀+1), so

    L_i ≤ log₂(3) L_(i−1) + log₂(C'_P(x₀+1)).

For a motif of length at most M, its raw affine carry numerator is at most 3^M−2^M, by bounding the expression by the all-odd motif term by term in the affine composition; the raw fixed-point denominator has absolute value at most 3^M. In particular Z and C′_P can each be bounded by 2·3^M. Hence log₂(C′_P(x₀+1))=O(1+M+log N+log m), with absolute effective constant. Iterating the geometric series proves the theorem. Since m≤N, if M≤N^θ for fixed θ>0 the parenthesis is O(N^θ), which gives the asserted seam tradeoff.

## 3. Stronger special case away from resonance

If a≤1−ε for some fixed ε>0, the block-cut value satisfies x₀≤m Z/ε even without a linear-forms theorem. The same argument now gives N≤C_(P,ε) (log₂ 3)^m (1+log m). For a fixed list with exponents n_i(k)=d_i k+e_i, d_i nonnegative integers not all zero, the leading drift coefficient Σ d_i log r_(u_i) is nonzero by unique factorization. If positive, a(k)>1 eventually, excluding all such cycle words. If negative, a(k)→0, and the previous bound excludes all large primitive lengths. This recovers the fixed-ray obstruction without Matveev.

## 4. Exact eventual dichotomy on a specified affine ray

In the negative-drift affine case, the cutting boundary can be fixed for all sufficiently large k: compare the finitely many affine functions of k representing the boundary prefix logarithms, first by slope and then by intercept. Each comparison is a comparison of products of fixed positive rational powers. The cut expression for x₀ is a rational function of finitely many positive rational exponential sequences in k. Each suffix base is at most 1 by the cut and the total multiplier base is less than 1, so x₀(k) converges exponentially to a rational L. If L is not a positive integer, eventually x₀(k) cannot be a positive integer. If L is a positive integer, eventually the only possible integral value is L. After multiplication by 1−a(k), the equality x₀(k)=L is a finite exponential-polynomial equality. Either all coefficients vanish after grouping equal bases, or the largest surviving base fixes a nonzero sign for all sufficiently large k. In the former case x₀(k)=L identically for every sufficiently large k.

An integral affine fixed point really follows the parity word: a word of length N with S odd letters and affine numerator 3^S x+C has a unique actual parity class modulo 2^N. This follows by lifting one bit at a time; the two lifts differing by 2^j have values after the first j steps differing by the odd number 3^(number of odd letters in that prefix). The actual class satisfies 3^S x+C≡0 (mod 2^N), which has a unique solution because 3^S is odd. Any integral fixed point belongs to this same class. Hence the identity case above is sufficient: eventually every word in the ray is a presentation of the same genuine orbit. Otherwise only finitely many parameters give positive integer cycles. This is a coefficient-checkable classification for a specified affine ray, not an exclusion for general exponents near the resonant hyperplane.

## Provenance and remaining gap

The affine map and rotation/cocycle law are inherited (repository note "Collatz pentagon operator coherence"; THM-4634). THM-4716 uses balanced Christoffel rotations for a different smooth-gap analysis. The fixed-run theory for actual Collatz cycles is older: the palette {0,1} specializes to the logarithmic lower bound on the number of parity runs, in the spirit of the classical [m-cycle bounds of Simons and de Weger](https://deweger.net/papers/%5B35%5DSidW-3n%2B1-ActaArith%5B2005%5D.pdf). The contribution here is the common statement for **arbitrary fixed mixed-parity motifs**, with a block-boundary carry estimate and a simple forward precision recurrence. Nothing here rules out a primitive cycle with enough switches or a growing motif palette, let alone all nontrivial positive cycles. The same issue is visible in the tournament construction: a fixed palette under increasing size may preserve the scalar relation while running out of allowable changes of state.

## 5. Two-sided sharpening when every motif has an odd letter

Suppose the odd fraction of every motif is at least \(\rho>0\). Let \(\alpha=\log_2 3\) and \(\beta=1+(1-\rho)/(\rho\log_2 3)\). The forward bound after j blocks from the cut is \(O(D\alpha^j)\), where \(D=1+M+\log N+\log m\). In reverse, the block equation also gives \(3^{ns}\mid q_u y-p_u\) at its exit y. The zero case again would make a genuine cycle return after one copy of u, since the equation then gives the entry \(x=z_u\). Inverting an even step grows \(x+1\) by at most 2 and inverting an odd step does not increase it on a positive orbit. If the reverse prefix has S odd letters then it has at most \((1-\rho)S/\rho\) even letters, as it comprises complete blocks. The next block's odd count is at most \(\log_3(C'_P(x_0+1)) + (1-\rho)S/(\rho\log_2 3)\). Therefore after j reverse blocks its odd count is \(O_\rho(D\beta^j)\) when \(\beta>1\), and its full length is at most \(1/\rho\) times this. If \(\rho=1\), there are no positive cycles because every motif expands.

Split the blocks into a forward prefix of \(j\) and a reverse suffix of \(m-j\), each anchored at the same cut, then sum their length bounds. Balancing the exponents gives

    N ≤ O_ρ(D exp(γm)),  γ = log(α) log(β)/(log(α)+log(β)).

This improves the one-sided exponential base for any fixed positive \(\rho\) and does not apply if purely even motifs are permitted. The example palette \(\{01,10\}\) has \(\rho=1/2\) and the estimate holds, but it has a simpler total-slope obstruction if its motif structure is known; the meaningful cases can mix multiple odd fractions.

## 6. A matching obstruction for divergence

The forward precision recurrence in Section 2 does not require a cycle or Matveev. On any nonrepeating positive trajectory starting at x, any actual finite word of length N admitting a decomposition into m repetitions-blocks with motif length at most M satisfies

    N ≤ C α^m(1+M+log(x+1)).

If an entry equalled the motif's fixed point, the very next copy would repeat an orbit value, so the nonzero-divisibility argument still applies. Thus an escaping trajectory cannot have all prefixes encoded by a bounded number of repetitions from a fixed palette; such a compression must grow at least logarithmically in its observed length. In particular an eventually periodic parity template cannot describe a divergent positive integer orbit. This is an obstruction to compressed itineraries, not a proof against divergence.
