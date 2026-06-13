---
source: claude-2026-06-03-S605
status: synthesis (carrier-level unification) + verified concrete observations
tags: [resonance-lattice, cancellation-family, nowhere-zero-flow, redei, OCF, involution, permanent, determinant, collatz, baker, formal-group, LRC, riemann]
---

# Resonance lattices everywhere

After poking around the repo I stopped seeing separate problems and started
seeing one object wearing different clothes. Two existing threads were already
most of the way there; this note welds them to the lattice picture and adds the
concrete reason LRC and Rédei sit on opposite sides of the line.

## The one object

Given quantities `x_1,…,x_m` (speeds, prime-logs, tournament arcs, `log2/log3`),
form the **resonance lattice**

```
L = { c ∈ ℤ^m : Σ_i c_i x_i = 0 }.
```

Its short vectors are the **resonances**. Every member of the repo's
"cancellation family" (S599) asks the same thing: **control the sign /
non-vanishing of a weighted sum over `L`**,

```
Σ = Σ_{c ∈ L} w(c)     — the PERMANENT face (all-orders cancellation, Vitali wall).
```

The clean identification: my LRC Poisson formula (HYP-2154)
`p_0 = Σ_{c∈L(V)} ∏_i κ(c_i)` **is** S599's `(★)` inclusion–exclusion sum
`p_0 = Σ_S(−1)^{|S|}meas(⋂_S D_i)` re-indexed by the lattice, and (S537o) its
**full-support terms are nowhere-zero flows** on the speed dipole. Three threads,
one carrier.

## The line: SOLVED ⇔ a symmetry of `L` collapses permanent → determinant

S599's resolution template says every solved member turns a permanent-shaped sum
into a determinant-shaped (certifiable) one. The lattice picture says *how*: by a
**symmetry of `L`** — a sign-reversing involution, or a unimodular/flow structure
that factors the sum.

**Rédei is the solved face [verified].** `#HamPaths(T)` is a permanent-shaped sum
over the cycle lattice; it varies wildly (`{1,3,5,9,11,13,15}` at `n=5`) but its
**parity is forced to 1** for every tournament — a sign-reversing involution over
GF(2). This is the Rosetta stone: the project's own theorem is the cancellation
family's *solved* instance, and it is solved by a lattice symmetry.

**LRC resists the same trick — and now I can say exactly why [verified].** The
obvious symmetry of `L(V)` is negation `c ↔ −c`. But the LRC kernel `κ` is
**even** (`κ(−c)=κ(c)`), so negation is sign-*preserving*: it doubles terms,
never cancels them. The trivial involution fails. LRC needs a *non-obvious*
sign-reversing symmetry of `L(V)`, and none is known — that gap *is* the open
problem. (A pretty side-fact: `κ((n+1)/2·k)=0`, arc-width Fourier zeros, e.g.
`κ(3)=0` at `δ=1/6`.)

**Collatz is the near-lattice face [proved/analogy].** `log2, log3` are
ℚ-independent, so the exact lattice is trivial; the controlling object is the
**near-lattice** — the convergents of `log₂3`, i.e. Baker's linear forms
`|K log2 − L log3|` (the cycle gap `|2^E − 3^k|`). My HYP-2147 identity
`K log2 − L log3 = log n + D(n)` lives here, and the convergent gaps
`0.405, 0.288, 0.118, 0.052, 0.0136, …` shrink but stay bounded below — Collatz's
Vitali wall. And the **formal-group logarithm `arctanh = log_F`** (the repo spine,
my first Collatz session) is precisely the *linearizer* carrying the
multiplicative resonances `2^K` vs `3^L` into this additive relation lattice. The
arc closes.

## The atlas

| problem | eval `x_i` | resonance lattice `L` | the `(★)` sum | certifying symmetry | status |
|---|---|---|---|---|---|
| LRC | speeds `v_i` | `{c: c·V=0}` | `p_0=Σ∏κ(c_i)` | sign-rev involution? (`κ` even ⇒ trivial fails) | open |
| LRC / flows | `v_i` mod `n*` | full-support `c` | flow polynomial | Tutte/Seymour | partial |
| Collatz | `log2, log3` | near-lattice (convergents) | `K log2 − L log3` | Baker linear forms | open |
| Rédei / OCF | tournament arcs | cycle lattice | `#HamPaths` (permanent) | GF(2) sign-rev involution | **solved** |
| Riemann | `log p` | zeros `ρ` (dual) | `ψ(x)−x=−Σ x^ρ/ρ` | zero-free region | open |
| Goldbach | primes | additive `a+b=n` | `r(n)=ΣΛΛ` | circle method | open (3-prime done) |
| P vs NP | — | — | permanent vs determinant | barrier theory | open |

## Why this is more than a slogan

It is **actionable**: it says where to push each open member — *find the symmetry
of its resonance lattice that pairs terms of opposite sign*. For LRC the concrete
sub-goal becomes sharp: `κ`-negation is sign-preserving, so look for a *twisted*
involution on `L(V)` (e.g. a flow-reversal combined with a residue shift) that is
sign-reversing on the full-support shell — the shell that S537o already
identifies with nowhere-zero flows. That is a finite, well-posed search on the
flow lattice, not a vague hope.

## Open / next

1. **Twisted involution on `L(V)`.** Search the full-support (flow) shell for a
   sign-reversing involution combining `c↦−c` with a coordinate permutation /
   residue shift; even a partial one would bound `p_0` structurally.
2. **Flow-polynomial ↔ `p_0`.** Make the S537o flow count quantitative inside the
   Poisson sum: express the full-support contribution to `p_0` via the dipole
   flow polynomial evaluated at a root of unity.
3. **Unimodularity test.** For which speed sets does `L(V)` have a basis turning
   the `(★)` sum into a determinant (a TU-like structure)? Those would be the
   *certifiable* (poly-decidable) LRC instances — the easy side of the
   structure↔complexity dichotomy (HYP-2154).
4. **Carry the lens back to OEIS/partition functions:** linear recurrences are
   resonance lattices (kernel relations); the repo's extended OEIS sequences are
   sums over those lattices — is there a sign-reversing involution explaining the
   integrality / sign patterns there too?
