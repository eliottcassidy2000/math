# Keeping the logarithms: a sharper transfer budget and an exact prime-power envelope

2026-10-19. **PROVED elementary identities and conditional transfer sharpening.** The analytic supply and ambient-cancellation inputs (S),(A) are exactly those stated in the attached Howe–Huffman–Kohlmeyer–Kruer paper. This note does not independently prove those analytic inputs.

## 1. The modulus is a join of predecessor valuations

For a selected prime bank P, the paper chooses L=lcm_{p in P}(p−1). At every prime q,

    v_q(L)=max_(p in P) v_q(p−1).

This is an actual join in the coordinatewise divisibility lattice, not a sum of predecessor sizes. It is the same retained coordinate used by the repo's Carmichael interface (N,lambda,e), where the slack is

    v_q(e)=v_q(N−1)−max_(p|N)v_q(p−1).

The prime bank and the final selected subset have different lcms; the compiler deliberately uses one common multiple to satisfy every selected predecessor obligation. The converse is not implied: many Carmichael subsets may satisfy their own smaller lcm while failing product-one modulo the whole bank's L.

For real X≥2 and 2≤Y≤X define the maximal smooth predecessor envelope

    E(X,Y)=product_(q≤Y prime) q^floor(log X/log q).

Every Y-smooth positive integer at most X divides E(X,Y). Hence L divides E(5y,y^δ) under (S). The exact logarithm is

    log E(X,Y)
       = sum_(q≤Y prime) floor(log X/log q) log q
       = sum_(j≥1) theta(min(Y,X^(1/j))),                    (1)

where theta(z)=sum_(q≤z prime)log q, and terms vanish once X^(1/j)<2. This is a two-cutoff prime-power carrier: one cutoff bounds the prime support, the other bounds the power. It retains information that the radical alone loses.

At Y=X, E(X,X)=lcm(1,…,floor X) and log E=psi(X). Thus the previous session's prime-to-six norm-torsion order is obtained by taking this full envelope at odd degree and deleting its 2- and 3-primary factors. Equality of the numerical orders does not identify their groups: the supported norm lattice has elementary p-torsion multiplicities, whereas an lcm is an integer carrying maximum valuations.

## 2. The paper's exponent loss is partly a bookkeeping choice

The paper uses

    log L≤y^δ log(5y),  log|G|≤log L,
    D=max(1,ceil(log_2 |G|)),  R=ceil(64Q/N).

It then upper bounds both log|G| and D by U=y^(δ+σ), and R by O(y^σ). This is a convenient sufficient envelope. Keeping the original logarithms gives instead

    log|G| = O(y^δ log y),
    D = O(y^δ log y),
    R^2 = O(y^(2σ)).                                       (2)

No prime number theorem or improved smooth-prime input is needed: N≥y^(1−σ) and Q≤5y suffice. The original ambient condition

    β+2δ+2σ<1                                             (3)

still ensures (R−1)A≤Q and hence the same character-energy lemma.

For expected size t=y^η, the nontrivial-character error is at most

    exp(−c t/(D R²)).

To make this negligible relative to 1/|H|, it suffices that

    [t/(D R²)] / log|G| → infinity.

By (2), the ratio is bounded below up to a fixed positive factor by

    y^(η−2δ−2σ)/(log y)^2.                                (4)

Therefore the sharper sufficient power condition is

    η>2δ+2σ,                                              (5)

in place of the paper's η>2δ+4σ. Bernoulli size concentration also dominates 1/|H|, since (5) implies η>δ. Retained size n≥N/2 and ρ=t/n≤1/2 require η<1−σ as before.

The original proof now runs unchanged with a group-size envelope C y^δ log y. There is no need to replace the loss in log(n/t) by 2σ:

    log(n/t)≥(1−σ−η)log y−log 2.                          (6)

Consequently the Bernoulli-collar transfer yields the lower counting exponent

    ((1−ν)/(1+ν))(1−σ−η),                               (7)

under (3), (5), η<1−σ and 0<ν<1. This is a sharper conditional parameter theorem; letting all small parameters tend to zero still gives the same exponent-one conclusion as the paper, not a finer Pomerance asymptotic.

A fixed auxiliary finite group or cardinality modulus changes only constants in these estimates. The exact-cardinality transfer developed alongside this note removes the collar separately; that is an additional theorem, not a consequence of (1) alone.

## 3. Boundary and next questions

- At η=2δ+2σ, the lower estimate in (4) tends to zero. This proof gives no conclusion at the boundary; it does not prove the theorem false there.
- Improvements to constants 64,100,200 do not change these power thresholds. Improvements to the ratio Q/N, the pruning depth D, or the energy mechanism can.
- The envelope (1) can be evaluated exactly to compare a concrete bank's actual lcm with the worst-case budget. Prime support size, maximum valuation and log-group size must be kept distinct.
- A bank with a common predecessor modulus is sufficient for Korselt selection. Passing to an arbitrary graph or Collatz word set requires a new certificate map; the numerical similarity of an lcm or clock modulus is not enough.

## Attribution

This is an elementary refinement of the supplied paper's Sections 3–6. Its estimates intentionally trade logarithms for small powers for convenience; retaining them is not a new analytic input. The lattice/prime-power interpretation links the paper to the repo's existing Carmichael successor interface and the preceding integral norm-lattice synthesis.
