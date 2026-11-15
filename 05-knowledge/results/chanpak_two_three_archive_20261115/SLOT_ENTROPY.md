# The physical 2-slot and 3-slot clocks lie in an exact growth cone

A literal (Q_2) step in the previous [port routes](../chanpak_physical_routes_20261115/PORT_ROUTES.md) was easy to make
high-dyadic, but it was too large compared with height. Is that just
a defect of constructing the module from singletons? No, if the
reader's step costs a fixed fraction of its two-slot count.

**Slot comparison.** For every nonempty tournament D and integers
b≥a≥2,

    Q_b(D) ≤ Q_a(D)^(log_a b).

For b>a, equality holds precisely when D is transitive.

**Proof.** By the [Rédei–Berge power-sum formula](https://arxiv.org/abs/2307.05569),
as used in the [prior Cayley readers](../chanpak_bott_cayley_periods_20261115/CAYLEY.md),
the integer polynomial Q(t) is a sum of monomials t^d with nonnegative
integer multiplicities. These monomials count odd directed cycle
packings together with independent signs on each nontrivial cycle;
the empty packing is the single term t^n. Index all terms with
multiplicity. At a the term weights are the positive numbers x_i=a^d.
At b they are x_i^λ, λ=log_a b≥1. Thus

    Q_b(D) = Σ_i x_i^λ ≤ (Σ_i x_i)^λ = Q_a(D)^λ.

For λ>1 equality needs only one term; any directed triangle of a
nontransitive tournament adds terms. A transitive D has Q_j=j^n
by labeling each vertex by its slot, so equality is attained.
The ratio Q_b/Q_a^λ is exactly Σ_i (x_i/Q_a)^λ, an explicit
concentration of the packing weights, not an estimate of a
random tournament.

The same calculation gives a quantitative penalty. The largest
packing term at a is the empty packing a^n. If D has τ directed
triangles then Q_a≥a^n+2τ a^(n−2). Since a probability p obeys
p^λ≤p·p_max^(λ−1), the sum of powers gives

    Q_b / Q_a^λ ≤ (a^n/Q_a)^(λ−1)
                ≤ (1+2τ/a²)^(1−λ).

So the transitive comparison is not only sharp; one can bound
its loss by a genuine local combinatorial invariant. In the
2-to-3 case a triangle already imposes a uniform strict loss.
This does not make Q₂ a substitute for the full higher profile.

## Consequence for an absolute return at a small modulus

Suppose a fixed marked reader sees one module through slot grades
at most r, so that its count is a fixed integer linear combination
of Q₁,…,Q_r, possibly with a fixed constant. There is a constant
C depending only on the reader for which

    N ≤ C Q₂(D)^(log₂ r)       (r≥2).

Use the comparison term by term with absolute coefficients; Q₁≤Q₂.
If this device's physical step M is at least c Q₂(D) for fixed
c>0, it follows that

    M ≥ c' N^(log_r 2).

In particular at rank three the best possible exponent in this
architecture is at least log₃2≈0.631, not the earlier precision
exponent 1/160. Achieving the latter with a step proportional to
Q₂ would require r≥2^160, or a change of architecture. This applies
also to a split-clock pair with common transitive padding where its
Q₂ difference is a fixed nonzero fraction of Q₂ itself. Replacing
the singletons by more sophisticated modules cannot evade the
inequality within these assumptions.

An important escape remains: use two physical modules with a step
o(Q₂), such as a **small signed difference** despite large carrier
coordinates, or change the reader rather than scale its same charge.
Independently programming congruences does not break this absolute
growth cone. The [frozen-grade sieve](FROZEN_GRADES.md) describes a
separate arithmetic restriction on that signed-difference alternative.
