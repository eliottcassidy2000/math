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

**Shared mechanism.** For iid real input of mean zero and variance one,
Chebyshev gives tightness of the normalized-sum laws. Compactness of
`[0,1]^Q` extracts limits of their distribution functions at rational
arguments; monotonicity and tightness reconstruct a probability law.
This is the compact-product step related to Tychonoff. Bare product
compactness alone does not establish weak convergence of probability laws.

**Identification.** With characteristic function phi,
`phi(t)=1-t^2/2+o(t^2)` gives `phi(t/sqrt(n))^n -> exp(-t^2/2)`.
Characteristic-function continuity identifies the standard Gaussian.
It is also the unique normalized-convolution fixed point **at mean zero
and variance one**. The former sentence saying that a unique fixed point
alone forces convergence was false: on the compact discrete set {0,1,2},
`0->0, 1->2, 2->1` has a unique fixed point and a two-cycle. Every cluster
point must be identified; uniqueness cannot substitute for that step.

**Map.** Laws on R -> distribution-function coordinates in a compact
product. Preserved: subsequence existence after tightness and reconstruction.
Required sidecar: identification of every possible limit. Rates, dependence
conditions and finite variance are additional data; other stable laws can
occur with other normalizations and assumptions.

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

**Stronger recovered arithmetic result:** the [prime Legendre tournament is
the countable random tournament](../root27_machine_bridges_20261026/anchor/LEGENDRE_SIDECAR.md),
a known CRT/Dirichlet construction. Every finite tournament occurs; increasing
prime-prefix frequencies remain a separate question. The continuation also
gives a limited four-object dictionary preserving the converse involution,
without asserting an arithmetic or native-composition equivalence.

## 4. Area bands, Möbius bundles, the fundamental theorems, and the spiral

**Signed area.** For C1 real f, the parametrization
`Phi(x,t)=(x,t f(x))`, 0<=t<=1, satisfies
`Phi*(dx wedge dy)=f(x) dx wedge dt`. Its Jacobian changes sign when f does;
at a zero a vertical fibre collapses. This is a signed region with a
degenerate parametrization, not nonorientable gluing. If F'=f on the real
interval, FTC gives `integral f = F(b)-F(a)`.

**Möbius bundle versus square-root cover.** Around a simple complex zero,
continuation of a square root changes its sign. The associated real sign
line bundle over that loop has a Möbius band as a bounded-fibre subbundle.
The covering itself is not that band: over the loop it is a circle, and
its complex covering surface is orientable. The hourglass root surface
has this sign monodromy; arithmetic squareclasses additionally decide
rational lifting. Trivialization of a complex topological cover does not
by itself decide a rational section. The former wording conflated these
three objects.

**FTC--FTA connection.** Along a parametrized loop avoiding the zeros of p,
choose a continuous logarithm lift L. Ordinary FTC on the parameter interval
says `integral p'/p = L(1)-L(0) = 2pi*i*w`, where w is the winding number.
For a large circle, dominance by the leading term of a degree-n polynomial
gives w=n. If p had no zero in the disk, p/|p| would extend to the disk and
force winding zero. This proves a root exists without first assuming a
factorization into roots. The infinite logarithm cover stores integer
winding; an m-th-root cover retains its residue modulo m. The sign bundle
retains only parity. Logarithm branches also exist on simply connected
subdomains avoiding zero, not only on the full universal cover.

**No uniqueness-of-topology claim.** The complex plane itself satisfies FTA.
The algebraic numbers inside C are algebraically closed and totally
disconnected in their subspace topology; C_p is another algebraically closed,
totally disconnected example. The old assertion that the argument-principle
route is available *precisely* over fields with a nontrivial multiplicative
universal cover was unsupported and is withdrawn. The proved argument above
is over C, with its analytic and topological hypotheses retained.

Real FTC likewise needs explicit regularity and a specified field. A
polynomial intermediate-value property over a real-closed field alone does
not supply arbitrary analytic FTC: over the real algebraic numbers, the
step at the missing cut pi is locally constant, with derivative zero, but
changes between 3 and 4. See the [2026-10-26 continuation](../orientation_monodromy_rule30_20261026/README.md)
for the complete connection maps and hostiles.

## 5. Method note

Three of the four analogies became exact statements only after naming the object the numeral or
shape was standing for (the parallelogram law, the product-of-intervals compactness, the
antisymmetry of the Legendre relation, the monodromy group). The fourth produced a repaired
statement and a hostile. That is the repository's standard: a resemblance earns one cheap probe, and
the probe either yields a map or a stopping reason.
