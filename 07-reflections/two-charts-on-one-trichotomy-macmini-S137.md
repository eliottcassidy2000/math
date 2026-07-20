# Two charts on one trichotomy: `{−1,0,1}` and `{1,½,0}` are the same object

*mac-mini-2026-07-20-S137. The owner: "look back at the summand and multiplicand graph, and
… think about how `{-1,0,1}` and `{1,1/2,0}` are functionally equivalent but have shown up in
this repo many times each. this is what i mean by even/odd vs positive/negative."*

*They are functionally equivalent because `x ↦ 2x − 1` carries one onto the other. What makes
the observation worth a file is that the repo has been computing in **both charts for months
without naming the change of variables** — and once named, four separate results this session
turn out to be the same statement.*

---

## The dictionary

| | `{0, ½, 1}` | `{−1, 0, +1}` |
|---|---|---|
| chart | **additive** / probability | **multiplicative** / sign |
| natural operation | averaging | `(−1)^k` |
| the involution | `p ↦ 1 − p` | `S ↦ −S` |
| the owner's name for it | **positive/negative** | **even/odd** |

`x ↦ 2x − 1` and `x ↦ (1+x)/2` are inverse. This is the Fourier–Walsh ↔ probability
dictionary: a `±1` variable with `P(X=1) = p` has `E[X] = 2p − 1`. And it is the repo's own
"two arithmetics" (additive tiling hypercube vs multiplicative `H`-norm) one level down, on
the value set rather than on the algebra.

**The whole content is the fixed point.** Both charts have exactly one, and it is the same
point: `p = ½ ↔ S = 0`. That is the **tie** — the diagonal, the draw, the balance.

## Four instances, all verified this session

**1. A tournament is a probability matrix, and its skew-Seidel matrix is `2p − 1`.**
Read `p_ij ∈ {0, ½, 1}` — `1` if `i` beats `j`, `0` if it loses, `½` on the diagonal. Then
`S = 2p − 1` is exactly the skew-Seidel matrix of THM-1440, landing in `{−1,0,1}`, skew, with
zero diagonal. Verified on 200 random tournaments, and **complementation `p ↦ 1−p` IS
`S ↦ −S`**, with the tie as the unique fixed point.

Paley makes it literal: the skew-Seidel matrix of the Paley tournament **is the Legendre
symbol** `χ(j−i)`, whose value set is exactly `{−1,0,+1}`; its probability avatar `(1+χ)/2`
is the quadratic-residue indicator with `½` sitting at `0`. Checked at `q = 7, 11, 19`.

**2. The fiber fraction is a tie probability.** CLAUDE.md records
`f(n) = (½)_{n−2}/(n−2)!` with GF `(1−x)^{−1/2}`. That equals

> `C(2k,k)/4^k` — **the probability of an exact tie in `2k` fair coin flips**

and also `E[U^k]/k!` for `U = T²/2`. All three agree exactly for `k = 0..8`. So the repo's
`½` is a **fair coin**, its "two-sheeted branched cover" is that coin's square root, and the
tie it computes is the `0` of the sign chart.

**3. THM-1500's admissible ladder *is* the trichotomy.** The master equation
`Φ(−s·g(s)) = 1/(1+s)` with `U = χ²_d/2` gives `Φ(x) = (1−x)^{−d/2}`, admissible exactly when
`2/d ∈ ℤ₊`:

| `d` | exponent `d/2` | `2x−1` | `n = 2+d` | `g(s)` | status |
|---|---|---|---|---|---|
| 0 | `0` | `−1` | 2 | `log(1+s)/s` | **obstructed** — GMC(2) open |
| 1 | `½` | `0` | **3** | `2 + s` | **minimal counterexample** |
| 2 | `1` | `+1` | 4 | `1` | the four-term example |

The exponent ladder `{0, ½, 1}` maps under `x ↦ 2x−1` onto `{−1, 0, +1}` — and **the minimal
counterexample sits at the fixed point.** The same `½` that is the fiber fraction's fair coin
is the `½` that makes `(1+s)^{2/d}` a *perfect square*. Square root and square are one
one-half.

**4. The telescoping principle is a statement about the sign chart.** THM-1520/1540: the
nullcone is the strictly one-sided charge support. In sign coordinates, take the trichotomy of
charge signs `{−1, 0, +1}`:

> the nullcone is the two **pure corners** — support inside `{+1}`, or inside `{−1}`. Anything
> touching `0`, or touching both signs, is outside.

So "forbidding one variable telescopes" (S135) and "the nullcone is one-sided" (S136) are the
statement that **the nullcone is the boundary of the trichotomy and the hard case is its
interior**, with the balanced point the deepest interior. That is instance 1's fixed point,
one level up.

## Census — the repo has been in both charts all along

| marker | files | chart |
|---|---|---|
| Legendre | 84 | sign |
| `sgn` | 20 | sign |
| skew-Seidel | 7 | sign |
| odd function | 6 | sign |
| merged metagraph | 69 | probability |
| fiber fraction | 29 | probability |
| `(½)_k` | 15 | probability |
| `(1−x)^{−1/2}` | 9 | probability |

## What to do with it

The useful form is a **habit**, not a theorem:

> When a repo object takes three values, or an involution has a fixed point, ask which chart
> you are in — and write the other one down. Sign coordinates make parity, characters and
> determinants easy; probability coordinates make averaging, measure and generating functions
> easy. The translation is free.

Two concrete payoffs already visible. First, it explains why the `½` keeps recurring in
unrelated-looking places: the fiber fraction, the `d=1` MGF exponent, the perfect-square
substitution and the merged metagraph's `Z₂`-quotient are all the fixed point of the same
involution. Second, it suggests a question I have not answered: **the `n=3` GMC counterexample
sits at the fixed point — is that forced?** If the minimal case of a construction always lands
on the self-dual point, that is a principle, not an accident, and it would predict where the
minimal case sits in other families.

## The honest limit

This is a **change of coordinates**, not a theorem, and it proves nothing on its own. It
should not be cited as evidence for anything. Two specific cautions:

- The trichotomy in instance 4 is the sign of a `ℤ`-grading, not a three-valued quantity; the
  "balance" statistic used to illustrate it is crude and mishandles charge-0 monomials. The
  precise statement is THM-1540's, and it should be quoted from there, not from here.
- The recurrence of `½` is partly structural (instances 2 and 3 share a genuine generating
  function, `(1−x)^{−1/2}`) and partly just that `½` is the commonest number in mathematics.
  I have checked the first three instances share an actual identity; I make **no claim** that
  every `½` in the repo is this one.

---

*Cross-links: THM-1440 (skew-Seidel, the sign chart), THM-1500 (the exponent ladder),
THM-1520 / THM-1540 (charge telescoping, the nullcone as pure corners),
`07-reflections/the-telescoping-principle-macmini-S135.md`, CLAUDE.md (fiber fraction, the two
arithmetics). Artifacts: `04-computation/two_coordinate_systems_macmini_S137.py` (+out).*
