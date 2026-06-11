# Pentagonal Lyapunov and the 72-code gate

**codex-2026-06-11-P1.** Prompt: extend the random-sign Lyapunov constant theme through
Euler's pentagonal theorem, the partition function, and the self-dual `[72,36,16]`
code problem.

## The useful reframing

Euler's pentagonal theorem is not just a recurrence for partitions. It is the
statement that the sparse sign denominator

```text
E(q) = 1 - q - q^2 + q^5 + q^7 - q^12 - q^15 + ...
```

has the hidden product factorization `prod(1-q^n)`. That product keeps all zeros
on the unit-circle boundary and gives the reciprocal `1/E(q)` the partition
numbers, whose growth is subexponential:

```text
log p(n) ~ pi sqrt(2n/3).
```

So Euler's signs are a zero ordinary-Lyapunov point in a much larger sparse-sign
space. If one randomizes the signs on the same generalized-pentagonal support, the
reciprocal recurrence usually starts behaving exponentially. The script
`04-computation/pentagonal_lyapunov_code72_codex.py` tested `160` paired-random
and `160` term-random samples through `n=650`; every sample had positive
finite-window slope. This is not a theorem, but it gives a precise theorem to
prove:

> A random pentagonal sign denominator almost surely has a zero inside `|q|<1`.

That is the whole Lyapunov constant in one sentence: nearest interior zero gives
the exponential rate. Euler's exception is not numerology; it is the product.

## The [72,36,16] gate

The exact Type II Gleason computation is similarly crisp. A putative extremal
doubly-even self-dual `[72,36,16]` code must have

```text
W = A^9 - 126 A^6 B + 3015 A^3 B^2 - 4398 B^3
```

where `A=x^8+14x^4y^4+y^8` and `B=x^4y^4(x^4-y^4)^4`. The forced vanishing at
weights `4,8,12` is the finite coding-theory analogue of the pentagonal
cancellation gate. The scalar gate is open:

```text
A_16 = 249849
lambda_5 for the weight-16 design = 78
W(1,1) = 2^36
all coefficients nonnegative
```

So the famous problem is not "does the formal modular enumerator exist?" It does.
The problem is support realization: can one realize that scalar partition
function as a self-dual binary support system whose minimum words form the
required `5-(72,16,78)` design?

This matches the current literature status checked during the session. Error
Correction Zoo still lists existence of `[72,36,16]` as unknown, Hannusch-Major
2023 gives newer necessary conditions via neighborhoods, and Borello 2013
restricts automorphism groups. The pressure is all on support and neighborhoods,
not on the scalar enumerator.

Sources checked:

- Error Correction Zoo, self-dual code page:
  https://errorcorrectionzoo.org/c/self_dual
- Hannusch and Major, "Neighborhoods of binary self-dual codes":
  https://arxiv.org/pdf/2206.05588
- Borello, automorphism restrictions for a self-dual `[72,36,16]` code:
  https://arxiv.org/abs/1303.4899
- Bell, "Euler and the pentagonal number theorem":
  https://arxiv.org/abs/math/0510054

## The common proof shape

The session's synthesis is HYP-2411: eta, Gleason, and OCF are all cancellation
gates.

- Eta: an infinite product hides behind sparse pentagonal signs.
- Gleason: a finite invariant ring hides behind low-weight vanishing.
- OCF / Reed-Muller digit ladder: odd-cycle collections hide behind low 2-adic
  digits of `H`.

In all three, a generic signed object wants positive Lyapunov/exponential growth,
while the mathematical object has a gate that cancels the dangerous low layers.
The scalar gate is necessary but lossy. It preserves the vanishing layers and the
growth signature; it destroys the support geometry that a proof must recover.

## Tournament Analysis

I used sign laws as tournament vertices, explicitly not runners/arcs. Alternate
vertex sets considered: individual pentagonal exponents, code weight layers,
Gleason monomials, automorphism-prime cases, root-angle packets, and proof
obligations. The chosen quotient preserves cancellation strength and destroys
root location and support design.

Pilot observable: lower finite-window reciprocal Lyapunov slope wins. Tie path:
the listed order of sign laws. The pilot tournament is transitive, with score
histogram `{0:1,...,8:1}`, no directed 3-cycles, singleton SCCs, and one
Hamiltonian path. That is useful mostly as a warning: scalar cancellation ranking
is too clean. The next observable should be two-dimensional, e.g.

```text
(nearest-zero radius, root-angle rigidity)
```

or

```text
(low-weight suppression depth, support-design compatibility).
```

The `[72,36,16]` problem should appear exactly where these gauges disagree.

## Proof targets

1. Prove the random pentagonal interior-zero theorem. A Jensen/Rouche/small-ball
   route should be possible because the support is sparse but not too sparse
   (`g_k ~ 3k^2/2`).
2. Classify zero-Lyapunov deterministic pentagonal sign laws. Euler is one; the
   all-plus control has low finite-window slope and must be audited before making
   any uniqueness claim.
3. Build the length-72 support problem from the weight-16 design upward:
   `5-(72,16,78)` first, then self-dual closure and MacWilliams suppression of
   weights `4,8,12`.
4. Define a code-facing Lyapunov diagnostic: random signed supports with the same
   first layer should typically leak low dual weights unless the self-dual gate is
   real.

The emotional shape of the result is pleasant: Euler's old pentagonal signs are
not just a formula for partitions; they are the prototype of a cancellation gate.
The length-72 code asks for the finite, binary, support-level version of the same
miracle.
