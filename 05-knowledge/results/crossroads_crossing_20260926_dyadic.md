# Dyadic potentials, the Bernoulli boundary, and poset projections retaining one integer

**Status: PROVED elementary identities and obstructions; FINITE-EXACT controls;
CITED literature scope; OPEN arithmetic coupling. No Collatz proof.**

## 1. Inheritance and the question we are actually testing

The closest proved mechanism is the bounded-correction obstruction in
[the blueprint energy audit, sections 3--4](collatz_blueprint_20260921_energy.md).
The Bernoulli identities and jump-to-valuation bridge were already proved
in [the B1 boundary note](bernoulli_boundary_20260925.md); section 3 below
applies that inherited boundary discipline to the user's printed series.
The canonical hostile is an arbitrarily long all-odd prefix. The corrected
near miss is replacing an orbit increment by a spatial derivative. The
least-used sidecar is the nonzero Fourier modes of the arithmetic source
filter: uniform extension counting keeps only its zero mode.

Our six concepts are dyadic phase, Bernoulli boundary values, actual crossing
events, valuation depth, poset extension laws, and positive atomic weights.
The anchor is a potential along one fixed orbit; the niche is an exact
poset-to-residue projection; the wildcard is why B_1 remembers a boundary.

## 2. The two displayed functions have distinct, exact behaviours

Interpret the user's first function, for x>0, as

    F(x)=2^[floor(log_2 x)+cos^2(pi {log_2 x}/2)].

Write t=log_2 x, k=floor(t), u={t}, and g(u)=cos^2(pi u/2).
Then F(2x)=2F(x), and F decreases within each dyadic interval. At a seam
x=2^m its left limit is 2^(m-1), but its value is 2^(m+1): a factor-four
jump. Thus dyadic homogeneity does not mean continuity or monotonicity.
The logarithm is

    log_2 F(x)=log_2 x+g(u)-u,

a bounded correction to logarithmic height. The inherited obstruction
therefore excludes one-step descent everywhere, and even excludes a
globally bounded lookahead for descent on odd accelerated orbits. It does
not exclude an input-dependent stopping time or an unbounded extra state.

The explicit verbal rendering of the second formula gives

    G(x)=(pi/x)[3/2+(1/2)cos(2 pi log_2 x)].

Here G(2x)=G(x)/2 and pi/x<=G(x)<=2pi/x. An antiderivative is

    E(x)=(3pi/2) ln x+(ln 2/4) sin(2 pi log_2 x).

Normalizing E by (3pi/2)ln2 gives

    V(x)=log_2 x+sin(2 pi log_2 x)/(6pi).

Its derivative with respect to log_2 x is
1+(1/3)cos(2pi log_2 x), between 2/3 and 4/3. Consequently V is strictly
increasing: V(y)<V(x) iff y<x. It is a smooth height coordinate, but it
cannot by itself discover descent that ordinary height misses.

Writing omega=2pi/ln2 gives the exact complex-power expansion

    G(x)=(3pi/2)x^-1+(pi/4)(x^(-1+i omega)+x^(-1-i omega)).

This is the legitimate Mellin connection to log-periodicity. It supplies
three frequencies, not an Euler product or a prime-distribution theorem.

## 3. The printed fractional-part series is a floor series

For u not an integer, the correct identity is

    {u}=1/2-(1/pi) sum_(k>=1) sin(2pi k u)/k.

The user's expression u-1/2+(1/pi)sum sin(2pi k u)/k therefore equals
floor(u), not {u}. At u=1/4 its value is 0, whereas the fractional part is
1/4. At integer u the sine series is zero: the Fourier expression takes
the midpoint value 1/2 for the fractional-part sawtooth. Exact integer
boundary values must be supplied separately.

This is a precise connection to the earlier B_1 observation. In the usual
generating-function convention B_1(x)=x-1/2 and B_1=-1/2. Its periodically
extended polynomial has a jump; the Fourier series retains the midpoint,
not either chosen endpoint. Higher odd Bernoulli numbers vanish because
t/(e^t-1)+t/2 is even, as a formal power series. The surviving linear term
is an exact symmetry/boundary correction, not evidence by itself for a
Collatz drift. The Fourier identity and its open-endpoint restriction are
[DLMF 24.8.2, n=0](https://dlmf.nist.gov/24.8.E2).

The analogous Collatz correction is explicit: an odd multiplication
contributes log_2(1+1/(3n)). Its accumulated value can be bounded on a
hypothetical divergent orbit by THM-4476. Its sign is positive, so simply
retaining this boundary term does not create contraction.

## 4. A poset construction that preserves the source exactly

Use the shortcut T(n)=n/2 for even n and (3n+1)/2 for odd n. For a word
w of length L with a odd steps, its affine map is

    T_w(n)=(3^a n+C_w)/2^L,
    C_w=sum_(j:w_j=1) 2^j 3^(number of later ones).

Its unique starting residue is r_w=-3^(-a) C_w modulo N=2^L.
By [THM-4503, parity posets and height selection](../../01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md),
words with every prefix multiplier strictly above one are the linear
extensions of an explicit two-chain poset, whenever this word set is
nonempty. This is a sufficient class of no-descent words. Carry-supported
no-descent words outside it must not be discarded in a Collatz reduction.

For a fixed such poset P with a odd labels, introduce the cyclic polynomial

    F_P(z)=sum_(w in LE(P)) z^(C_w)   modulo z^N-1.

The exact number of extensions compatible with an integer n is

    [z^(-3^a n mod N)] F_P(z)
      = (1/N) sum_(ell=0)^(N-1) exp(2pi i ell 3^a n/N)
                    F_P(exp(2pi i ell/N)).

This equals 0 or 1, because each residue has a unique parity word. The
same formula with a comparison indicator inserted selects its exact
arithmetic value. In contrast, F_P(1)=|LE(P)| gives the uniform extension
law: this is only the ell=0 mode of the source projection.

This construction identifies the lost information rather than hiding it
inside a heuristic probability. Source: labelled linear extensions.
Target: a cyclic group-algebra polynomial with exact carry exponents.
Preserved predicate: compatibility with n modulo 2^L. Lost by augmentation
z->1: every nonzero residue Fourier mode. Necessary sidecar: those modes,
or equivalently the coefficient vector. Cheapest hostile: the four words
of length six with five ones and positive prefix multipliers:

| word | C_w | r_w | A_4 precedes B_1 |
|---|---:|---:|---|
| 110111 | 287 | 27 | no |
| 111011 | 251 | 39 | no |
| 111101 | 227 | 47 | yes |
| 111110 | 211 | 31 | yes |

The uniform comparison probability is 1/2. For source 27 it is 0; for
source 31 it is 1. Randomizing the reveal order of facts about the fixed
source does not randomize its actual parity word.

The attached Kahn--Saks and balance-barrier arguments concern a specified
uniform law on linear extensions; their detailed scope and unresolved
external dependencies are recorded in the
[previous paper audit](crossroads_poset_20260926_barrier.md) and
[width audit](crossroads_poset_20260926_width.md). No result there automatically
controls the nonzero Fourier modes above. A new weighted inequality would
need its law and source filter in the statement. Even a proof of the full
1/3--2/3 conjecture would not, by this construction alone, prove Collatz.

## 5. Two constructive routes that survive the audit

**OPEN: a decorated crossing forest.** Actual upward and downward dyadic
crossings match in last-in-first-out order, giving nested excursion
intervals and a containment poset. Unmatched upward crossings form its
boundary. Label each event by its exact integer, valuation, and carry.
This preserves the original orbit while allowing alternative orders for
verifying commuting facts about disjoint excursions. A uniform linear
extension theorem can then organize a proof search, but a sign inequality
must still hold for the labelled arithmetic costs. Chains of all-odd
growth show that width cannot be assumed large. The excursion
27->41->62->31 ends in the same dyadic band at a larger value, so
returning across a boundary alone is insufficient.

**OPEN: positive theta weights with an exact residue filter.** The supplied
theta square has coefficient r_2(n), which vanishes at n=3. Thus it cannot
give positive mass to every possible exceptional starting integer.
The four-square theta series repairs full support: r_4(n)>0 for all n>=1.
For 0<q<1 weight each actual no-descent survivor n>=2 by r_4(n)q^n.
The resulting mass tends to zero iff every n>1 eventually falls below
itself, hence iff Collatz holds. This follows from positivity and a
summable divisor bound, not from modularity. The arithmetic lane proves
the explicit divisor formula. For example r_4(n)<=4n(n+1), so the
geometric-polynomial tail is explicitly summable. A surviving singleton has positive
mass and cannot disappear through a density normalization.

The actionable question is whether a positive theta/divisor weight and
the decorated crossing order produce a new decay estimate with controlled
endpoints. A modular transformation for the unmasked theta function does
not commute with the orbit-dependent survivor mask. That commutator,
alongside the residue Fourier modes, is the specific quantity to study.

## 6. Literature boundary and reproduction

The rotation connection is not a novelty claim. Kontorovich--Miller's
[2005 paper, section 5.2](https://arxiv.org/abs/math/0412003) explicitly
distinguishes source-distribution limits from one fixed trajectory.
The January 2026 preprint
[An Explicit Near-Conjugacy Between the Collatz Map and a Circle Rotation](https://arxiv.org/abs/2601.04289)
studies a shifted base-six logarithmic coordinate for the standard map.
Our crossing note uses an elementary unshifted identity and the repository's
THM-4476 reciprocal summability; it does not certify that preprint's
numerical claims or its full argument.

Run `python -B 04-computation/experiments/crossroads_crossing_20260926_dyadic.py`.
The [retained output](crossroads_crossing_20260926_dyadic.out) enumerates
every parity word at lengths 1--14, independently checks its actual source,
and verifies the exact singleton source projection and balanced-pair hostile.
All-odd identities are checked through length 300. Floating evaluations
of the two potentials are labelled diagnostics; the proofs use derivatives,
exact dyadic limits, and integer affine identities. No floating computation
decides a mathematical claim. All checks remain enabled under `python -O`.
