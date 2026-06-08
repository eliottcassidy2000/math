# LRC dual conformal symmetry: the observer is spacetime, the obstruction is dual (S717)

The prompt borrowed a word from scattering amplitudes — dual conformal symmetry — and asked what the
lonely runner conjecture can and cannot say in spacetime. The word turns out to fit exactly, and the
answer is a clean division of labor: the observer lives in spacetime, the obstruction lives in the dual,
and the quantity that decides the conjecture needs both.

There are two spaces, just as in amplitudes. The spacetime is the cylinder, time crossed with the circle,
on which the runners draw worldlines `x_i(t) = frac(v_i t)` and the observer is the single vertical strand
at zero — the picture the cluster already built as a pure braid. Loneliness is a local, geometric event
there: a time at which every strand is far from the observer's strand. The dual space is arithmetic: the
speeds as residues, the multiplier group `(Z/m)*` that the cluster calls perspective and that is really
the Galois group of a cyclotomic field, and the autocorrelation / homometry data. Ordinary conformal
symmetry acts on the first; dual conformal symmetry acts on the second; and the second is hidden from the
first.

The ordinary symmetries are the ones you can watch on the worldlines. Rescaling all speeds is rescaling
time, and the gap `M` does not change — conformal weight zero, verified on every test. Reversing all
speeds is reversing time, and `M` does not change. Relabeling runners is permuting strands. These are the
visible conformal group, and they are unremarkable precisely because they are visible.

The dual symmetry is the interesting one, because it is invisible. Multiply every speed by a unit `a`
modulo the shell `m`, and on the witness orbit the gap is exactly preserved — three hundred out of three
hundred, and the proof is one line: multiplication by `a` just permutes the witness times `b -> ab`, so
the maximum over witnesses is unmoved. But this multiplier scrambles the worldlines beyond recognition;
there is no spacetime motion it corresponds to. It is a symmetry of loneliness that lives only in the
dual. And it has a special-conformal generator: inversion, `v -> v^{-1}` mod `m`. Inversion is exactly the
dangerous move — `a = v_i^{-1}` sends runner `i` onto the central band, which is the bad multiplier of
THM-420. So the inversions are the special conformal transformations, and the reason the hard core (the
transversal, Paley-like configurations) is hard is that there the inversions `{+-v_i^{-1}}` exhaust the
whole multiplier group: every rotation is killed by some runner's inversion, and no dodge survives. The
hardness is a statement purely about the dual symmetry being maximal.

Then the prompt's real question — what cannot be expressed in spacetime, and what cannot be expressed in
the dual — answers itself in two complementary facts. The observer cannot be expressed in the dual. The
gap `M` is not translation invariant (move all the speeds and it changes), but the distance multiset, the
autocorrelation, the entire dual data, is translation invariant. So the choice of origin — which strand is
the observer — is a spacetime fact the dual literally cannot see. `{1,2}` and `{3,4}` have the same dual
data and different gaps. The observer is irreducibly spacetime. And conversely the gap cannot be expressed
in the dual either: homometric configurations, genuinely non-congruent sets with identical autocorrelation,
have different `M` — `{0,1,2,6,8,11}` and `{0,1,6,7,9,11}` share every distance and give gaps `1/5` and
`2/9`. So `M` is dual-symmetric (the multiplier preserves it) but not dual-determined (the autocorrelation
does not fix it). Loneliness carries spacetime information beyond anything the dual records.

That is the whole division of labor, and it is the amplitudes lesson transposed. In the amplitude the hard
structure lives in the dual, where the hidden symmetry is largest, and locality is a partial, emergent
description. Here the obstruction to loneliness lives in the dual — it is maximal dual symmetry, the flat
autocorrelation of the Paley core, the inversions covering every multiplier — while the observer, the
origin, the thing that makes the question well posed at all, lives in spacetime and nowhere else. The gap
`M`, the answer to the conjecture, is the pairing of the two: a dual-symmetric quantity anchored to a
spacetime origin. You cannot compute it from the worldlines without knowing the arithmetic, and you cannot
compute it from the arithmetic without knowing where the observer stands. The conjecture is hard for the
configurations where the dual symmetry is so large that the spacetime side has no handle left — exactly
the place the amplitude is hardest for the same reason. Ordinary conformal symmetry on the worldlines plus
dual conformal symmetry on the residues is the LRC's Yangian, and the perspective key the owner keeps
returning to is its dual half.
