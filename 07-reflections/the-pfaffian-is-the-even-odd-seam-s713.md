# The Pfaffian is the even/odd seam, made algebraic (S713)

"Understand even odd and the Pfaffian." The two are the same thing. The even/odd seam the project keeps
hitting — `chi(C_n)` even vs odd, the glass transition at even `n`, the parity defect, the `(p-1)/2`
reciprocity bit, `v_2(T(n))=(n-1)/2` — is not a list of coincidences that happen to be governed by
parity. It is one structure, and the Pfaffian is its cleanest algebraic name.

Start with the definition itself. The Pfaffian of a skew-symmetric matrix exists *only* in even
dimension; in odd dimension a skew matrix is singular and `Pf` is zero. A tournament's skew-adjacency
`S = A - A^T` is the canonical skew matrix with `±1` off the diagonal, so the question "even or odd `n`?"
is literally the question "does this tournament have a Pfaffian?" The seam is the Pfaffian's domain of
definition. Everything else follows from looking at what survives on each side.

On the even side there is a single integer, `Pf(S_T)`, and it is always **odd**. The reason is the
seam again, one level down: mod 2 the skew matrix forgets its orientation entirely — every `±1` becomes
`1` — and `S` collapses to `J - I`, the adjacency matrix of the complete graph. The Pfaffian mod 2 is
then just the number of perfect matchings of `K_n`, which is `(n-1)!!`, a product of odd numbers. So
`Pf(S_T) ≡ 1 (mod 2)` for every tournament, no matter how the arcs are oriented. This is the exact
even-`n` twin of Redei's theorem that the number of Hamiltonian paths is odd. Redei's count is a
near-Hamiltonian (spanning *path*) object alive at every `n`; the Pfaffian is a perfect-*matching*
object alive only at even `n`; both are forced odd because mod 2 the orientation washes out and you are
left counting a structure on `K_n` whose count is intrinsically odd. For even `n` a tournament carries
both odd invariants at once. I checked whether they are more tightly linked — whether `|Pf| ≡ H` mod 4
— and at `n=4` it holds perfectly, but it is a small-number mirage: it breaks at `n=6`. The honest
shared content is precisely the oddness, and the oddness is precisely the mod-2 collapse of orientation.

On the odd side the Pfaffian vanishes, but it does not leave nothing behind. The skew matrix is now
singular with a one-dimensional kernel, and that kernel is spanned by the vector of *cofactor
Pfaffians*: `w_i = (-1)^i Pf(S_{T-i})`, the signed Pfaffian of the tournament with vertex `i` deleted.
Deleting a vertex from an odd-`n` tournament lands you at even `n-1`, where the Pfaffian is defined and
odd — so the kernel of an odd-`n` tournament is a vector all of whose entries are odd integers, each one
the Pfaffian of an even subtournament. The even and odd sides are not two separate phenomena; they are
two rungs of a ladder, and the rung-to-rung map is vertex deletion, `n <-> n-1`. This is the same
`n-1`-things-to-`n`-things correspondence the owner keeps circling (the degree-`n` polynomial with `n+1`
coefficients and `n` roots; the `n` runners with `n-1` speeds). Here it is exact and algebraic: the
odd-dimensional object is assembled from the even-dimensional objects one vertex down, glued by the
cofactor identity `S w = 0`.

Putting both sides together gives one clean invariant that does not care about parity at all:
`det(I + 2A)` is always an *odd perfect square*, hence `≡ 1 (mod 8)`. On the even side it is `Pf(S)^2`
with `Pf` odd; on the odd side it is `(1^T w)^2` with `1^T w` an odd sum of odd kernel entries. The
perfect-square half is the old THM-174; the oddness is the trivial observation that `I + 2A ≡ I` mod 2;
together they pin it to `1 (mod 8)`, the residue every odd square occupies. The determinant of the
skew part, `det(S) = Pf^2`, then reads the seam off as a 2-adic valuation: even valuation (a square)
when `n` is even, infinite valuation (zero) when `n` is odd. That is the same ledger as THM-305's
`v_2(T(n)) = (n-1)/2`; the Pfaffian is where the valuation bookkeeping becomes a single matrix
invariant.

Finally the seam has a finer, mod-4 grain, and the Pfaffian carries that too. Complementation reverses
every arc, `S -> -S`, and a `2m`-dimensional Pfaffian scales by `(-1)^m`, so `Pf(T^op) =
(-1)^{n/2} Pf(T)`. For `n ≡ 0 (mod 4)` the Pfaffian is complement-invariant; for `n ≡ 2 (mod 4)` it
flips sign. This is the project's recurring mod-4 phase — the `√p` versus `i√p` of the Gauss sum, the
"middle element out of phase," the reciprocity sign — appearing now as the action of complementation on
the Pfaffian. The even/odd seam splits even `n` again into `0` and `2` mod 4, and the Pfaffian's sign
is the indicator.

So the answer to the prompt is a single sentence with four faces: the Pfaffian is defined exactly on
the even side of the seam (existence), is forced odd there by the mod-2 collapse of orientation
(Redei's even twin), reconstitutes the odd side as a ladder of odd cofactors one vertex down (the
`n <-> n-1` map), and records both the 2-adic valuation (`det = Pf^2`, `det(I+2A) ≡ 1 mod 8`) and the
mod-4 reciprocity phase (`(-1)^{n/2}` under complement). Even/odd is not a property the Pfaffian
happens to have. The Pfaffian is what even/odd looks like when you write it as algebra.
