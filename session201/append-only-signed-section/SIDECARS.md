# First-growth gap and a parity obstruction to boundary collisions

**PROVED, independently audited in the stated scope.** These two short consequences are
separate from the frozen main proof. They do not establish target coverage.

Throughout this sidecar, a physical source means a member of the nontrivial
orbit S defined in PROOF.md; the singleton is excluded.

## 1. The first two possible append-only section targets

For every physical source s, the two smallest distinct values in Dtail(s)
are exactly

    D_s and F(s)=D_s-E_s+U_s.

Indeed, a word which has not used Z or J keeps D_s. If its first such letter
follows r O letters, its new path count is

    F(s) if r=0, and D_s+2^r U_s if r>=1.

The latter is strictly greater than F(s). Subsequent D values do not
decrease. Conversely the empty word and the one-letter word Z attain D_s
and F(s). Also F(s)>D_s because U_s>=2d_s>=2E_s and U_s>0.

By the main theorem, if f_s>0 these are also the two smallest distinct
section targets attainable by appending letters to s. In particular,
merely arranging that the signed target exceeds D_s is insufficient.

An exact hostile is

    s=t.ZOOZ=(15,1,23,1), f_s=8, L(s)=33,
    D_s=15 < 33 < F(s)=37.

It has no target-preserving append-only repair into T, although target33
is already represented: t.OOOOZ=(33,1,63,1), whose W extension is in T.

This has an infinite version. For lambda=3^r,

    s_r=t.Z.(OJ)^r.OOZ=lambda*(15,1,23,1)

has target33lambda strictly between its first two reachable section
targets15lambda and37lambda. The number33lambda is represented by the
balanced source t.OOOOJ.(OJ)^r. Its W extension supplies a section witness.

## 2. A proposed one-letter signed-boundary collision is impossible

Let z=(-1,-7,-1,-7) be the formal signed boundary. Its entire O/Z/J orbit
has

    D odd, d odd, U congruent to E modulo2.

These facts hold at z and are preserved by each displayed raw operation.
In particular no point of that formal orbit can have d=0.

For an arbitrary formal raw vector w, exact equality of distinct
one-letter outputs has the following complete descriptions:

    O(w)=J(w) iff w=lambda*(1,1,1,1);
    Z(w)=J(w) iff w=lambda*(1,0,2,1);
    O(w)=Z(w) iff w=(D,U,U,0).

The last collision locus is therefore unreachable from z, at every word
length. This excludes one proposed source of globally equal signed
readouts before entering a positive boundary. It says nothing about
longer-word collisions or the usefulness of the other two ray collisions.
The vector z and its descendants are formal; no physicality is assigned to
them merely from these parity calculations.
