---
source: claudebox-2026-06-03-S602
status: REFLECTION — the worry-set's two block-diagonalizations; the apex is the 2-block of the
  even modulus n; the additive face mod 2n-1 is odd, apex-free, and rigid.
tags: [block-diagonalization, two-faces, additive, dynamical, doubling, apex, 2-block, mod-2n-1,
  64-classes, transversal, primitive-root, 3-adic, LRC]
---

# Two faces of one seam: the apex is the 2-block of the wrong modulus

**Prompt (human):** see everything from block-diagonalization. The worry-set's static modular
rigidity (the 64 classes mod 27 = 3³) is the additive-face shadow and survives, while the dynamical
doubling face fragments at 2q — the same 2·7 seam from both faces.

The instruction and the insight together resolve something I had been circling. The worry-set is one
object, but it block-diagonalizes under two different group actions, and the Lonely-Runner
obstruction is a block that exists in only one of the two bases.

## The dynamical face, mod n

The phase dynamics is the doubling map `x↦2x mod n`. For even `n` it 2-adically collapses, and the
gcd-divisor decomposition of the danger (S599) carries a rank-1 **2-block** — the apex, the
zero-divisor, the single point where the single-modulus corrector dies (HYP-2063). Everything hard
about `n=14` lives in this one block of the mod-14 picture.

## The additive face, mod 2n−1

But the worry-set has another skeleton: the `64 = 2⁶` self-converse classes, realized as antipodal
transversals **mod `2n−1 = 27`** (oracle-S570: 8191 of them, all lonely). This is the *additive*
face — the residue structure of the speed set in the doubled-and-shifted modulus. And `2n−1` is
always **odd**. So the doubling is a clean permutation, the divisor decomposition has no even part,
and there is **no 2-block**. The apex has no counterpart here. The static rigidity — the Boolean
flip-lattice of the 64 classes, with the AP at the tight bottom — sits in this face and is
untouched by the seam.

The verification made the duality vivid: at `n=14`, the dynamical face mod 14 is dead, while the
additive face mod `27 = 3³` is *maximally* alive — 2 is a primitive root mod `3³`, so the doubling
there is a single 18-cycle, the most ergodic possible. The static modular rigidity is not merely
surviving; it is at its healthiest exactly where the dynamical rigidity has collapsed. Same prime
`2`, two opposite verdicts: it kills the even modulus and is absent from the odd one.

## The resolution: transport the obstruction to the odd face

This reframes the whole apex story. The apex is not intrinsic to loneliness; **it is the 2-block of
the wrong modulus.** Compute in the even modulus `n` and you meet a rank-1 obstruction; compute in
the odd modulus `2n−1` and it is not there. The repo already discovered the cure empirically —
HYP-2075, "multi-sieving has no apex" — and the two-faces picture says *why*: the pair-sum sieve
takes sums `v+w`, pushing the structure into the additive (odd, `2n−1`-flavored) face, where there
is no 2-block to obstruct. The single corrector fails because it stays in the dynamical face; the
multi-sieve succeeds because it is, secretly, the additive face. Resolving the conjecture is
transporting the obstruction from the face where it is a block to the face where it is nothing —
and the static rigidity (every transversal lonely) is the guarantee that nothing survives the
transport.

## A symmetry of heights

There is even a duality of depths. The dynamical face has the 2-adic height `v₂(n)` (H6) — `n=14` is
height 1. The additive face, when `2n−1` is a prime power, has its own `p`-adic height: at `n=14`,
`2n−1 = 3³` is a height-3 tower `{3,9,27}` (the repo's curvature-cubed Cassini atom). The hardness
lives in the small 2-tower of the dynamical face; the rigidity lives in the tall 3-tower of the
additive face. The seam is the place where a single factor of 2, fatal on one side, is simply
absent on the other.

## The transcending pattern

Block-diagonalization is the right altitude because the obstruction is basis-dependent: it is a
genuine block in one symmetry-adapted basis and identically zero in another. The Lonely Runner is
hard only in the dynamical (multiplicative, mod-`n`) eyes, where an even modulus manufactures a
2-block; in the additive (mod-`2n−1`) eyes the worry-set is a rigid Boolean lattice with every
transversal lonely. To prove it is to look with the additive eye — and the human had already named
that eye: the 64 classes mod `3³`, the shadow that survives.

**Artifacts:** `04-computation/lrc_two_faces_additive_dynamical_s602.py` (+`.out`); new **HYP-2150**.
Builds on HYP-2145 (divisor blocks), HYP-2140, HYP-2063 (apex/2-block), HYP-2075 (pair sieve),
HYP-2097 (64 classes), HYP-2117 (doubling).
