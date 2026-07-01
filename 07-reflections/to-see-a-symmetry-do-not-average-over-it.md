# To see a symmetry, do not average over it

*mac-mini-2026-07-01-S90. Reflection on HYP-3817.*

The owner named the next instrument before I built it. After the skew-spectrum turned out blind to the
covering-hardness — blind, as it happens, to symmetry itself — the tempting move was to reach for a finer
transform: the adjacency spectrum, which at six vertices did see the automorphism group, or the Seidel
matrix, or the determinant lens. The owner's instinct was to stop. The next instrument, they said, has to be
built to be sensitive to fixed points, not blind to them by symmetry — a covering invariant, or the variance
of H, rather than another clever transform. It is exactly right, and the reason is a small tautology worth
saying slowly.

A spectral transform is built to be invariant under relabeling. That is its whole virtue: it does not care
which vertex you call one. But the thing the covering problem cares about — how symmetric a tournament is,
how large its automorphism group, whether it is one of the few-copy needles a thin subcube cannot cover — is
*precisely* a statement about how relabeling acts on it. An invariant of the group action has averaged over
the group action. It has integrated out the very coordinate you are trying to read. So a transform can be
blind to symmetry not by weakness but by design, and the adjacency spectrum's apparent sight of the
automorphism group was a small-case accident: it leaked the symmetry through the score degeneracies, and at
seven vertices the leak closes and a single spectrum holds tournaments of different symmetry again. No
transform is a reliable instrument for a property of the group, because a transform is what you build to not
see the group.

The instruments that see fixed points are the ones built from the action rather than invariant under it. A
covering number is defined by how the orbits pack the cube; it feels a high-symmetry class directly, as a
needle with too few copies to be caught by a coarse net. A moment — the variance of H over the labeled
ensemble — is a sum weighted by orbit size, so the automorphism group is written into it as a weight, not
integrated away. And the automorphism group itself is the most honest instrument of all: it is the count of
fixed points, the trace of the fold, the Burnside number the whole involution atlas turned on. These are not
cleverer transforms. They are the opposite kind of object — engaged with the group instead of quotiented by
it — and that is why they can measure what the transforms average over. There is a complementarity here, not
a ranking: the transform separates the asymmetric twins the covering number cannot; the covering number
flags the symmetric needles the transform cannot. You do not want the finer instrument. You want the one
whose blind spot is not the thing you are hunting.

And this, finally, is the lonely runner's lesson wearing tournament clothes. The covering minimum is a fixed-
point problem — its extremal loneliness is an atom, a single point, the deep symmetric well no perturbation
reaches. For years the analytic side reached for transforms: the Fourier expansion, the Delsarte bound, the
Fejér kernel, every one a PSD, averaging, symmetry-respecting instrument. And every one hit the same wall,
the spectral gap, blind to the pointwise atom exactly as the skew-spectrum is blind to the automorphism
group. The instruments that made progress were the other kind: the second-moment floor, a moment; the lazy-
cut, a covering. The variance of H on the tournament side and the danger-relation's second moment on the
runner side are the same instrument, and they are the right one for the same reason — they are built from
the group and the geometry, not invariant over them.

The pattern that transcends the theorem: **to measure a symmetry, use an instrument built from its action,
never one invariant under it — invariance is averaging, and you cannot see what you have averaged away.** When
a fine structure keeps escaping your transforms, ask whether the structure is a property of the group itself,
and if it is, put the transforms down. The right next instrument is not more resolution in the symmetric
direction; it is a covering, or a moment — something that engages the fixed points instead of politely
declining to notice them.
