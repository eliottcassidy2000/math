# Typed analogy atlas: Cauchy–Schwarz, central limit versus Tychonoff, reciprocity as a tournament, and the two fundamental theorems

2026-10-25, klein. Each entry follows the connection contract (source, target, map, preserved
predicate, lost information, sidecar, decisive test) of `00-navigation/RESEARCH-PROTOCOL.md` §5 and
the card "Type every analogy and every implication". Statuses: **EXACT** (a theorem, proved or
classical), **ANALOGY** (shared mechanism, no identification), **COINCIDENCE** (a shared numeral
with no map). No historical novelty is claimed for the classical parts.

## 1. Cauchy–Schwarz is the triangle inequality — EXACT, in inner-product spaces only

**Statement.** In a real or complex inner-product space, `|<x,y>| <= |x||y|` for all `x,y` is
equivalent to `|x+y| <= |x| + |y|` for all `x,y`.
*Proof.* `|x+y|^2 = |x|^2 + 2 Re<x,y> + |y|^2`, so the triangle inequality for all pairs is
`Re<x,y> <= |x||y|` for all pairs; replacing `y` by `e^{i theta} y` gives `|<x,y>| <= |x||y|`; the
converse is the usual expansion. QED.
**Map.** Identity on the space; the predicate is the same inequality read on `|x+y|^2`.
**Lost.** The inner product itself: a norm satisfies the triangle inequality without any
Cauchy–Schwarz unless it comes from an inner product (Jordan–von Neumann: the parallelogram law).
**Decisive test.** `l^1` on `R^2`: `|(1,0)+(0,1)|_1 = 2 = 1+1`, triangle inequality holds, no bilinear
form reproduces the norm, so no Cauchy–Schwarz statement exists. The analogy is exactly as strong as
the parallelogram law.

## 2. Central limit theorem and Tychonoff — ANALOGY with a precise shared step

**Shared mechanism.** Tychonoff: a product of compact spaces is compact. The weak-convergence
compactness used in the central limit theorem (Prokhorov/Helly: a tight family of laws has a weakly
convergent subsequence) is proved by a diagonal argument on distribution functions, i.e. by the
compactness of `[0,1]^N` — Tychonoff for a countable product of intervals, which needs no choice.
So the existence of subsequential limits of the normalized sums `S_n/sqrt n` is Tychonoff's theorem
in the product-of-intervals form, after Chebyshev's inequality supplies tightness.
**What Tychonoff does not supply.** Identification of the limit. That the limit is Gaussian needs the
characteristic-function (Lévy continuity) or Lindeberg replacement argument: the Gaussian is the
unique fixed point of the renormalization map `mu -> law of (X+X')/sqrt 2` among finite-variance
centred laws. Tychonoff gives a cluster point; uniqueness of the fixed point makes it a limit.
**Lost.** Rates (Berry–Esseen), dependence structure, and the finite-variance hypothesis: without it
tightness can hold at another scaling and the fixed points are the stable laws — same compactness,
different fixed point, which is the decisive test that the analogy is about compactness, not about
the Gaussian.
**Map.** Laws on `R` -> points of a compact product space via distribution functions; preserved:
subsequential limits exist; sidecar: the fixed-point equation.

## 3. Quadratic reciprocity as a tournament — EXACT for one half, COINCIDENCE for the numeral four

**Exact statement.** Let `P_3` be the primes `≡ 3 (mod 4)` and `P_1` those `≡ 1 (mod 4)`. Quadratic
reciprocity says `(p/q)(q/p) = (-1)^{(p-1)(q-1)/4}`. Hence on `P_3` the Legendre relation
"`p` is a square modulo `q`" is antisymmetric: exactly one of `(p/q), (q/p)` equals `+1`. **The
Legendre relation is a tournament on `P_3`**, and a symmetric graph on every pair involving a prime
of `P_1`. (Compare the Paley tournament, which is the same construction inside one field `F_p`.)
**The numeral four.** `A000568(4) = 4` counts isomorphism classes of 4-vertex tournaments (`TT4`,
source over a 3-cycle, 3-cycle over a sink, the strong one). Reciprocity's "four cases" are the
residues `(p mod 4, q mod 4) in {1,3}^2`. There is no map between a count of orbits and a table of
residues; MISTAKE-214's warning applies. What is true and checkable: among the first 120 primes of
`P_3` every one of the four classes occurs as an induced subtournament, with frequencies
`0.368, 0.125, 0.126, 0.381` against the random values `3/8, 1/8, 1/8, 3/8`, and triangle density
`0.2534` against `1/4` (`qr_tournament.py`; sampled at 2,000 and 20,000 primes the frequencies are
within `0.002` of random). See [QR_TOURNAMENT.md](QR_TOURNAMENT.md) and HYP-9375.
**Lost by the tournament reading.** The `P_1` half (symmetric), the supplementary laws for `-1` and
`2`, and the field structure that makes Paley tournaments doubly regular; the prime tournament is
only random-like, not regular (scores on the first 120 primes range from 47 to 68 around 59.5).

## 4. Area bands, Möbius strips, the two fundamental theorems, and the spiral — ANALOGY, typed

**The signed area band.** For a real function the region between the graph and the axis is a
2-chain whose orientation is the sign of `f`; the fundamental theorem of calculus is Stokes'
theorem `int_a^b dF = F(b) - F(a)` for the 1-form `dF = f dx`. The band is orientable; it changes
sign, not side, at a zero. The honest Möbius object is elsewhere: the **square root** `s^2 = f` over a
loop around a simple zero returns with `s -> -s`, i.e. the real line bundle with monodromy `-1`, the
Möbius band. In this session's own objects this is exactly the sign deck `(mu_2)^7` of the hourglass
root surface: each root coordinate `s_i = sqrt(d_i L_i)` is a Möbius band over a loop around the
branch line `L_i = 0`, and BOUNDARY §5's residual torsor is the statement that over a non-closed
field the band may fail to have a rational section.
**The spiral.** The fundamental theorem of algebra by the argument principle is Stokes' theorem for
the closed 1-form `d arg p = Im(dp/p)` on a disk minus the roots: `(1/2 pi) oint d arg p = deg p`
on a large circle. `arg p` is single-valued only on the universal cover of `C^x`, the infinite
spiral surface of `log`; a loop on the base lifts to a path on the cover displaced by
`2 pi i · deg p`. "The strip returns" with a **shift** (`Z` monodromy, log) where the Möbius band
returns with a **flip** (`Z/2` monodromy, square root). Both are the same mechanism, monodromy of a
covering, with different groups.
**Hostile to "the only topology".** Algebraic closure is not a topological property: `C_p` is
algebraically closed and totally disconnected, with no winding numbers at all. The winding-number
proof needs `pi_1(k^x) != 0` (true for `C`, false for `R`, where real polynomials indeed need not
have roots); other proofs of the fundamental theorem of algebra (Galois-theoretic, Liouville) use
different sidecars. The repaired statement: the argument-principle route to roots is available
precisely over fields whose multiplicative group has a nontrivial universal cover, and the spiral is
that cover.
**Pairing FTC with FTA.** Both are degree statements: `int_{partial} = int d` with degree one for
`[a,b]` and degree `deg p` for the circle. The decisive difference: FTC is exact over any ordered
field with the intermediate value property (it fails over `Q`), FTA over `R` fails (no winding
obstruction) — so the two theorems have different minimal fields, which is why they cannot be
identified.

## 5. Method note

Three of the four analogies became exact statements only after naming the object the numeral or
shape was standing for (the parallelogram law, the product-of-intervals compactness, the
antisymmetry of the Legendre relation, the monodromy group). The fourth produced a repaired
statement and a hostile. That is the repository's standard: a resemblance earns one cheap probe, and
the probe either yields a map or a stopping reason.
