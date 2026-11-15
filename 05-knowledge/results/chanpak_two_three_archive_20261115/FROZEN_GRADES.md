# Freezing the higher clocks creates a prime sieve on the low step

The [last physical route](../chanpak_physical_routes_20261115/PORT_ROUTES.md)
made the higher optional slot charges part of the baseline. A natural
attempt is to increase that height by preserving every large charge
while changing only the two-slot charge. The two bases are coupled
more tightly than a simple positivity inequality suggests.

Write Q_j(T) for the ordered cover count by j optional path slots.
For a tournament of order n its slot polynomial Q(t) has integer
coefficients, is monic, and has parity n. Its nonleading
coefficients are even. For clarity, the
[Rédei–Berge power-sum formula](https://arxiv.org/abs/2307.05569)
under principal specialization gives a term
2^(number of nontrivial cycles) times t^(number of permutation cycles)
for each packing of vertex-disjoint odd directed cycles, including
fixed points as cycles of length one. The empty packing gives t^n,
and each omitted degree is a sum of even terms. The fact that this
specialization is the path-slot count follows by encoding a
permutation with weakly increasing slot labels and enforcing a
strict rise at every backwards adjacency.

## A prime sieve for three-term changes, even across different orders

Suppose three physical modules have a common H=Q₁, common higher
slot values Q₃,…,Q_r, and their Q₂ values form an arithmetic
progression of nonzero stride d. Their orders need not agree.
For each odd prime p≤r+2, t=2 is the negative of p−2 modulo p
and the slot at p−2 is frozen (including Q₁ for p=3). Thus

    Q₂(T) ≡ (−1)^|T| Q_(p−2)(T) mod p.

Only two residues at most are available. Any three-term arithmetic
progression with p not dividing d would use three distinct residues.
Hence d is divisible by **every odd prime p≤r+2**. This extends
the mod-5 obstruction to varying the entire frozen-grade horizon.
In particular, keeping the stride in a fixed odd factor times a
power of two cannot freeze arbitrarily many grades.

## Same-order polarization gives an exact factorial divisor

If two candidate modules of the same order n have the same frozen
values, put Δ(t)=Q_A(t)−Q_B(t). Let ε=1 for n odd and ε=2 for n even.
Its roots include zero of multiplicity ε and ±j for
j∈{1,3,4,…,r}; the coefficients are even because the leading
monic terms cancel. Division by the monic factors therefore gives

    Δ(t) = 2 t^ε (t²−1) Π_(j=3,…,r)(t²−j²) R(t),  R∈Z[t].

At t=2 this says the nonzero change ΔQ₂ is divisible by

    (r−2)!(r+2)!/2   for odd n,
    (r−2)!(r+2)!     for even n.

Indeed the product of j²−4 is (r−2)!(r+2)!/4!, and the
remaining factor at 2 is 2^(ε+1)·3. Also for a nonzero change
the degree of Δ is at least ε+2(r−1) and at most n−2, so
n≥2r+ε. This is a necessary size and lattice condition, not
sufficiency or an assertion of a sharp factorial divisor.

At r=3 the odd prime obstruction is 5 (in addition to 3 from
H); the fixed-order divisors are 60 and 120. The
[split-clock construction](CHAN_HISTORY.md) escapes the
coupling by paying an odd factor that includes 5. It is genuine
progress from the abstract congruence test to a physical pair,
but not an unbounded transport theorem. A reader that raises the
height by freezing more and more slot grades will pay for that
freeze in the arithmetic of its step.
