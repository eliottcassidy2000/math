# A finite-state transfer that keeps the exact repair window

The older [October 17 quasi-Riemann/ROOT bridge](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/preprint_connections_solfeggio_digits_20261017.md)
distinguished average-energy amplification from a transfer that contracts
the remaining scale. The [October 5 quasi-Riemann manuscript](https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex),
§2, embeds a target sum in a larger family, passes to Fourier variables
to estimate its mean square, and then seeks a power contraction. Here the
transfer is an exact finite action on address digits, and the same
distinction matters: modular transfer by itself is not absolute descent.
The claims below are proved directly and do not assume a zero-free region.

Use the verified three-center family
`H_j = 432xy+129x+258y+15[x odd]+30[y odd]+45(p(x)+p(y))+135+2j`
with j=0,1,2. The checksum obeys p(2x+1)=p(x) and
p(2x)=p(x)+x+[x odd]. State variables are over Z/m, m odd.

## 1. Exact completion from any prescribed odd prefix

Put w=p(x), keep x odd, let A append 1 and B append 01. Then
A:(x,w)→(2x+1,w) and B:(x,w)→(4x+1,w+x+1).
Using right-to-left composition,

    C=A^(-2) B:         x→x−1/2,    w→w+x+1,
    E=C^(-2) A C A^(-1): x→x,       w→w−3x/2−1.

First move x to 0 with powers of C; E then translates w by −1.
Set w in anticipation of the final C steps to the desired x.
This programs any pair (x,p(x)) modulo *every odd m*. A^ord_m(2)
and C^m are identities, so all inverse steps are realizable as
positive appends. No negative or formal address counts are used.

Take y=1 while keeping x odd. Then H_j=561x+45p(x)+303+135+2j.
The first two coefficients have gcd 3, and j covers all classes
modulo 3. For a target N choose j to leave N−(303+135+2j)
divisible by 3, and use 187·13−15·162=1 to specify a compatible
state modulo m. Thus every odd modulus admits exact count
residues, independently of the absolute size of the output.

## 2. Complete while fixing the *actual sum* of the address pair

Append paired complementary bits to a prefix with x odd and y even.
Write w=p(x)+p(y), and use only the macro boundaries at which x stays
odd and y stays even. A appends 1/0, B appends 01/10, respectively.
With S=x+y and d=x−y, their exact maps are

    A: (x,y,w) → (2x+1,2y,w+y),
    B: (x,y,w) → (4x+1,4y+2,w+x+2y+3).

As before let C=A^(-2) B and E=C^(-2) A C A^(-1). Evaluation gives

    C: (x,y,w) → (x−1/2,y+1/2,w+d+3/2),
    E: (x,y,w) → (x,y,w−3d/2−1/2).

Both preserve S modulo m. First move to d=0 with C (d decreases
by 1); the E shift is then the unit −1/2. Lastly move x to the
desired state. This selects all pairs (x,w) with the same S.
For C^m=Id, sum the affine increment along the m translations;
the triangular term vanishes since m is odd. The order of A
divides ell=ord_m(2) also on the w coordinate, because on that
coordinate A^ell adds (2^ell−1)y. Every word representing C or E
has length divisible by ell when A^-1 is replaced by A^(ell−1)
and C^-2 by C^(m−2). Pad with A^ell to a single length K
independent of the selected residue, K=O(m² ell). Complementary
appends of a fixed length K send the **integer** sum S to
2^K(S+1)−1, not only to a fixed modular residue.

## 3. Joint near-target placement and residue selection

Start with the companion complementary-phase bracket at scale L=√(N/432).
It gives x,y=Theta(L), phase g=2^s=Theta(L^(1/3)), magnitude
|x−y|=O(√(Lg)), and x+y≡g−1 mod g. All these are literal addresses;
the existing three counts are within O(N^(1/3)).

Fix any m≥3 coprime to 30 (odd and coprime to 15). For sufficiently
large N, g≥2^(K+1). Drop the final K+1 low bits of each address, append
the first pair of bits 1/0 to establish boundary parity, and reserve
the remaining K bits for the transfer. All replacement pairs have the
same actual sum as the old pair because the dropped bits are
complementary. Each address moves by at most 2^(K+1). Pick j first
by N modulo 3, choose any x residue modulo m, infer y from the fixed
sum, and solve for w modulo m using invertibility of 45. The resulting
designated branch has exactly the residue N modulo 3m. Its two other
physical variants retain their exact +2 spacing.

The absolute error stays O_m(N^(1/3)). For the bilinear term, a change
(x,y)→(x+δ,y−δ) costs O(|δ|·|x−y|+δ²). For the checksum, let u,v
be the unchanged high quotients at scale g. A low word b contributes
`p(gu+b)=p(u)+A_b u+p(b)+z_b[u odd]`, with A_b the sum of
2^j over zero bits of b counted from the most significant end, and
z_b the indicator that the leading low bit is zero. For the
complementary word A_{g−1−b}=g−1−A_b, so the total variation is at
most (g−1)|u−v|+2g+2, hence O(√(Lg)+g). Parities and center offsets
are bounded. Since m is fixed, 2^(K+1) is constant, proving the claim.
To preserve a specified side of N, begin with the companion construction
for N±C_m N^(1/3), for a constant strictly larger than the error
bound, while still programming the residue of N itself. This yields
two exact windows, one on each side, with a branch in each precisely
aligned modulo 3m to the original target.

The two sides can in fact use the **same integer sum** if one
allows scanning short blocks. In the companion construction, let
u be its selected diagonal, w_0=gu+g/2, and take
x=w_0+t, y=w_0−1−t on the same sum. Partition t≥0 into blocks of
length B_0=2^(K+1), assuming g/2 is divisible by B_0. Each block
has a steered representative, using its reserved B_0 suffixes.
In the first block, both high quotients are u and the discarded
low-bit checksum terms are nonnegative. Moreover
p(g/2+t)≥p(g/(2B_0))=g/(2B_0) by repeated appending and the
checksum of a power of two. Thus this whole block is above N for
large N, since its quadratic deduction is O(B_0²).
Near t=C√(Lg), the smooth quadratic deduction exceeds the
O(Lg) diagonal surplus for a sufficiently large fixed C, while
the companion's complementary-phase error is O(N^(1/3))
uniformly on this range (the high quotients differ by
O(t/g+1)). So the last block is entirely below N.
Consecutive steered representatives have t-distance ≤2B_0
and their smooth readouts differ by O_m(√(Lg)); adding the
uniform phase errors gives O_m(N^(1/3)) gaps between them.
The first crossing therefore gives both sides with the same sum
and the designated residue. This same-sum proof does not claim
the logarithmic search time of the unsteered bracket.

This does not correct the final integer residual. It replaces it by
an integer multiple of 3m of size O_m(N^(1/3)); the theorem uses
**fixed** m. Making m grow increases the suffix budget and the
largest local address perturbation, so there is not yet a demonstrated
power contraction. Joint bracketing for 3m with gcd(m,30)=1 is narrower
than the unrestricted odd-modulus assertion of §1.

## 4. An unconditional spectral transfer and its scale barrier

On a complete p×p address square modulo a prime p not dividing 432,
the discrete Fourier transform in the output residue of the bilinear
phase plus any unary weights, even completely arbitrary checksum
phases, has |S(a)|≤p^(3/2) for a≠0. Indeed the matrix of phases
e(a·432xy/p) has orthogonal rows of norm sqrt(p), and the unary
factors have Euclidean norms sqrt(p). Pair this matrix bound with the
cyclic triangular window W_h(s)=#{0≤i,j<h:i−j=s mod p}. It has support
within distance h−1 and nonnegative spectrum, with sum of Fourier
coefficients p h. Thus its count at any residue has main term p h²
and error at most p^(3/2) h. Taking h>sqrt(p) ensures an actual covered
residue in its support. This is an all-phase Fourier estimate, not a
claim of exact surjectivity or a claim on short address intervals.
For arbitrary unary phases, −216x² and −216y² turn the bilinear phase
into a multiple of (x−y)², an explicit reason exact residues need
more structure than a unitary estimate.

The analogous obstruction for absolute windows is quantified by
the gap energy. Normalize the odd grid to integers. A gap of size
d between consecutive represented sites contributes exactly
(d³+2d)/12 to squared distance if d is even and (d³−d)/12 if d
is odd. Sum the squares of the two mirrored ramps to see this.
Consequently a square-energy bound of O(X^(1+epsilon)) on an octave
would imply a largest gap O(X^(1/3+epsilon/3)), with boundary sites
included or controlled separately. At a positive density of missing
sites, square energy is already Omega(X), so this amplifier alone
cannot beat the newly proved one-third exponent. A higher moment or
a scale-preserving transfer would have to supply a genuinely
stronger conclusion.

Replay: `python verify.py` checks literal and streamed single and
paired states, all residues for selected odd moduli, and exact
sum retention after the joint correction. These replays check the
transfer, not any conjecture about absolute coverage.
