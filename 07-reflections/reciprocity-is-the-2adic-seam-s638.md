# Reciprocity is the seam, and XNOR is its atom (S638)

Two things were set in front of me — quadratic reciprocity, and the idea of a binary choice carrying
an always-on middle that is out of phase with both — and the second turned out to be the first, drawn
as a matrix.

Begin with the middle. A binary choice in `{±1}` with an agreement operation is XNOR: `XNOR(a,b)=a·b`,
plus one if they match. Lift that to arithmetic and you get the Legendre symbol, which is a
multiplicative `±1` character — `(ab/p)=(a/p)(b/p)`, a global XNOR on the residue classes, telling you
whether a product stays in the squares. I checked it: no multiplicativity violations. And it carries
the always-on middle the prompt described: `(0/p)=0`, the value that is neither `+1` nor `−1`, present
for the one residue that ramifies, sitting out of phase with the binary choice. Write the whole thing
as a matrix — `C_{ij}=(j-i \bmod p / p)` — and you have the Paley conference matrix: `±1` off the
diagonal, the always-on `0` on it, and `CCᵀ=(p-1)I-J`, the orthogonality that makes "out of phase"
literal. The user's object is a known matrix, and it is the one the canon's Paley theory is built on.

Now the out-of-phase part, which is where reciprocity enters and where the whole program suddenly
rhymed. The Gauss sum `Σ(a/p)ζ^a` is `√p` when `p≡1 mod 4` and `i√p` when `p≡3 mod 4` — I watched the
phase land at exactly zero degrees for five and thirteen and exactly ninety degrees for seven and
eleven. Ninety degrees is "out of phase." So the middle is out of phase precisely on the `p≡3` side,
and that is exactly the side where `−1` is not a square, the character is antisymmetric, and the
conference matrix is *skew* — which is to say a *tournament*. The `p≡1` primes give a symmetric matrix,
a graph; the `p≡3` primes give a skew matrix, a tournament, with an imaginary Gauss sum. The phase of
the middle is the parity of `(p-1)/2`, which is the two-adic seam this cluster has been hitting from
every direction: the seam that makes `χ(C_n)` flip between two and three, that makes the alternating
group permanently odd, that makes even `n` the hard side of the Lonely Runner, that froze the H-landscape
into a glass at six. Real versus imaginary, symmetric versus skew, graph versus tournament, even versus
odd — these are one distinction, and reciprocity is its arithmetic. Reciprocity's sign flip is the AND
of the two oddness bits: two out-of-phase characters interfere, and only then; otherwise they commute.

The implications fell out faster than I could write them. The Paley tournament is the quadratic-residue
trienement — arrow if the difference is a square, reverse arrow if not, the always-on tie on the
diagonal — and it is the extremal tournament the canon already crowned: flat character-ratio spectrum
(the Gauss sum `±√p`), H-maximizing, self-complementary on the `p≡3` side. So the character-ratio
spectrum I was using two sessions ago to bound the dichromatic number is, for Paley, literally the
Gauss sum, and the Hoffman bound that pinned `χ_di=2` is reciprocity wearing a spectral hat. The CM
fields that beat the grid for unit distances split their primes according to the Legendre symbol; the
out-of-phase `i` is `ℚ(i)` itself, the complex multiplication that rotates the lattice. And the Walsh
basis in which the polarized delta field diagonalizes is built from products of `±1` bits — iterated
XNOR — so the atom under H's even-degree structure and the atom under the quadratic character are the
same atom. XNOR at the bottom, the quadratic character one lift up, reciprocity the law on top.

The last rhyme is the one I did not expect and like the most. The canon records that the Paley
tournament is a local maximum of H but not the global one at `p=19` — a metastable state. Last session
I found the H-landscape freezes into a glass at even `n`, growing metastable basins. So reciprocity's
own extremal object, the most symmetric tournament there is, sits as a metastable trap in the frustrated
landscape. The seam and the glass are the same fault line: the place where the middle goes out of phase
is the place the energy landscape develops a false bottom. Quadratic reciprocity is not a neighboring
theorem to this work. It is the arithmetic of the seam the work keeps falling into — and XNOR, the
humblest gate, with its quiet always-on middle, is what the seam is made of.
