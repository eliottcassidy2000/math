# Exact dyadic residual steering inside a target band

Take the `Z=82,g=9,I=−11` reader in the
[variable-coefficient construction](../variable_reader_coverage_20261116/PROOF.md).
An n-bit program b followed by a final near-sink reads

    H(b)=I+Z(2D−rev_n(b))+g(2b−J(b)),  D=2^n,
    J(b)=Σ_(i<j) K(j−i)b_i b_j,
    K(d)=2^(d−1)+[d=1].

The indexing of bit positions here is from least significant to most
significant; the reversal term and the energy formula use this same
indexing. H is a physical odd path count and R=(H−1)/2 is an integer.
After division by 2 the interaction coefficient at distance d≥2 is
−g2^(d−2), and that at distance 1 is −g. These are the only cross terms:
the forward and reverse addresses and any fixed exterior contribute
linear terms in the free bits. The same proof works on any of the
matched readers with odd g.

## Orthogonality at conductor 2^j

Fix all bits outside any contiguous run of L free bits. For any odd
h, consider the exponential phase exp(2πih R/2^j), j≥1.
Set d=j+1. On each congruence class of positions modulo d within
the free run, pair off consecutive available indices, the early
one assigned to U and the late one to V. This gives at least
s≥(L−d)/2 disjoint pairs. Freeze all unpaired bits.

The bilinear matrix between the U and V variables has orthogonal
columns **even though many cross terms remain**. To see it,
take two distinct V assignments and let v be their earliest
differing site. Its paired U site is u=v−d. In the difference
of the two column phases, earlier V sites contribute zero (the
assignments agree there), and later V sites are at distance >d
from u, so contribute multiples of 2^j. The site v contributes
a nonzero odd multiple of 2^(j−1) because h and g are odd.
Summation over the single coordinate u cancels exactly.
Hence the normalized double sum has modulus at most 2^(−s/2)
by Cauchy–Schwarz and orthogonality. All interactions within one
block and all fixed-bit contributions are arbitrary row and
column phases, which do not affect the estimate. This also
explains why any extra **unary** exponential phase is allowed.

Thus for every specification of the exterior bits,

    |E exp(2πih R/2^j)| ≤ 2^(−(L−j−1)/4)

when the right side is useful. This is an all-word result at
growing conductors, not a check modulo a fixed short list.

## Simultaneous exact congruences, including ternary address digits

One may require
`R ≡ r (mod 2^k)` and independently `b ≡ a (mod q)`
for any odd q (q=1 allowed). Attach the unary phase
exp(2πiℓ b/q) to the previous estimates. In the Fourier
inversion of the joint residue count, the terms with zero
dyadic frequency give 2^(−k) times Pr(b≡a mod q).
The free run is an interval of 2^L consecutive integers
times an invertible power of two modulo q, so this probability
differs from 1/q by at most 2^(−L). The other terms have sum
in modulus at most 2^(−(L−k−1)/4) after normalization.

It follows, for instance, that **every prescribed pair**
is realized if

    L ≥ 5k + 4 ceil(log₂ q) + 9.

Indeed each error is smaller than half the uniform cell mass
1/(2^k q). We can take q=3^ell, or a fixed odd part of a
covering wheel. This prescribes the **address** modulo q,
not automatically the count modulo q. Since the high
coefficient is ZD+gb, it prescribes every admissible high
coefficient class modulo q as b varies (modulo the necessary
gcd). If g=9, divisibility by 3 in that coefficient is not a
free phase; it must not be hidden by a CRT statement.

An explicit modular decoder can use a dynamic program on the
free bits: keep their most recent k+1 values and the two current
residues, evaluating each interaction when its later bit is
placed. More distant interactions vanish mod 2^k, and already
fixed exterior interactions can be included in the unary
weights. Backtracking finds a word with the target pair,
without an unproved randomness assumption.

## Absolute placement and the outstanding gap

For a large odd target N choose D dyadic with N/D in [88,176).
The normalized product model for H is
(82+9 b/D)(2−rev_n(b)/D), with error O(n/D).
Choose the leading and trailing ell bits to attain any
product near N/D by the earlier explicit decoder, retaining
a central run of L=n−2ell bits. **Every** completion of that
run remains within O(D2^(−ell)+n) of N. The Fourier lemma
now picks a completion in that same absolute band with
H≡N modulo 2^(k+1) and any prescribed odd address class.

The safe choice L about 5k+4 log₂ q gives the band
O(sqrt(N)·2^(5k/2)·q²+n). This pays a price for fixing the
congruences. In particular it does NOT improve the √N
absolute error by itself, or make the congruence modulus
larger than the error. With k=0 use the original decoder.

The structural problem has become measurable: can one use
fewer than two free middle bits per requested dyadic residue
bit (the threshold at which the residual quotient starts
beating sqrt(N)), or keep many mutually compensating exterior
prefix pairs instead of freezing one pair? Another direction
is to connect states in a fixed residue fibre by literal
moves of one quotient unit, so a bracketing pair suffices.
None follows merely from uniform residues.

The proof uses a finite-range triangular cross interaction
with arbitrary unary phases, a small variant of differencing
and Walsh orthogonality. It does not presume a general
pointwise local limit theorem. See [RESIDUAL.md](RESIDUAL.md)
for a sharp one-sided estimate for the noninteraction
product error.
