# Finite-surplus prime-composition jets

**PROVED:** the finite-horizon reduction and recursive count formula below.
**FINITE-EXACT:** 81 projection/substitution checks with aggregate prime data
through order 15, plus independent direct cyclic-palette orbits. No novelty,
minimality, or global complexity improvement claim. This classifies
module-preserving inflations, not arbitrary one-vertex attachments.

## Recovery and actual delta

Read-only source pin: `5862c94555eafea4e9a0638a17f7c661f8e45d37`.
The closest mechanisms are Session79 `PRIME_CYCLE_INDEX.md` (aggregate
plethystic recovery), Session86 `A000568_COMPOSITION.md` (recomposition and
universal-palette separation), and the THM-4803/4814/4816 prime-family tower.
All paths are under `05-knowledge/results/` except the named theorem entries.
Session79 already supersedes needing to enumerate nonrigid prime representatives
for the aggregate ordinary index. Session86 already warns that this recovery,
when fed the tournament Burnside index, is a structural round trip, not an
independent discovery of A000568.

The delta is a **surplus-adapted finite-horizon carrier**: a root with k slots
being inflated only up to k+h vertices needs at most D(h) rational coordinates,
independent of k, rather than its entire action index. Its state is refined
when the allowed surplus increases. This is a different truncation from total
root arity, and from the future-gcd-response moment carrier for Burnside sums.
The current note does not establish whether this precise refinement already
occurs elsewhere in the literature or the full repository.

Canonical decomposition uses Gallai's theorem, stated as Theorem 16 of
[Boussaïri–Chergui–Ille–Zaidi](https://arxiv.org/html/1805.04607v1), consulted
2026-10-29. The algebraic reduction below is proved directly.

## 1. The carrier and its theorem

Let Q be a prime tournament on k>=3 vertices, G=Aut(Q), and

    Z_Q(p) = |G|^(-1) sum_(g in G) product_(d odd) p_d^c_d(g).

Every cycle of a tournament automorphism has odd length. Let

    J_h(Q;z) = [weighted degree <=h] Z_Q(1+z_1,1+z_3,...),
    weight(z_d)=d.

Equivalently its coordinate at z^alpha is the mixed binomial cycle moment

    M_alpha(Q) = |G|^(-1) sum_g product_d binom(c_d(g),alpha_d),
    sum_d d alpha_d <=h.                                      (1)

For a family of distinct prime classes, sum their J_h; tagged multiplicities
are also allowed, but must not be silently deduplicated. Let the child palette
contain exactly one singleton and arbitrary other nonempty tournament
isotypes, with ordinary series

    C(x)=x(1+b_1 x+b_2 x^2+...),  E(x)=C(x)/x-1.

**Theorem.** The number generating function I_Q of Q-rooted inflations obeys

    I_Q(x) = x^k J_h(Q; E(x),E(x^3),E(x^5),...)
             modulo x^(k+h+1).                               (2)

**Proof.** Canonical prime-root decomposition identifies inflations precisely
under G's slot action. Burnside therefore gives I_Q=Z_Q(C(x),C(x^3),...).
For every term of Z_Q, sum_d d c_d(g)=k. Factor x^k to obtain exactly
Z_Q(1+E(x),1+E(x^3),...). Since E(x^d) has valuation at least d, every
monomial of weighted degree greater than h contributes only beyond k+h.
Deleting those monomials proves (2). The binomial theorem proves (1). QED.

The singleton condition matters: it selects one canonical baseline color.
For all ordinary tournament children it holds. A palette with multiple tagged
size-one colors needs a Taylor expansion around that different baseline.

There are at most

    D(h)=sum_(r=0)^h p_odd(r),
    sum_(h>=0) D(h)t^h = 1/((1-t) product_(d odd)(1-t^d)),       (3)

coordinates. For h=0..8 this is 1,2,3,5,7,10,14,19,25, regardless of k.
Additional action-specific vanishing can reduce the number. In particular
J_h is not a claim of a coarsest carrier for the single fixed tournament
child series: its h+1 output coefficients could be a smaller observation.

## 2. What changes when the size increases

Projection J_(h+1) -> J_h simply discards weight h+1. No old coordinate needs
changing. The new coordinates are mixed cycle moments indexed by partitions
of h+1 into odd parts. Thus this classification has a concrete refinement law,
not a conjectural scalar recurrence on prime counts.

Write M_j=M_(alpha_1=j), N=M_(alpha_3=1). With all tournament children,
E(x)=x+2x^2+4x^3+12x^4+..., formula (2) gives surplus coefficients

    h=0: M_0 = 1 (or the number of roots in an aggregate family),
    h=1: M_1,
    h=2: 2M_1+M_2,
    h=3: 4M_1+4M_2+M_3+N,
    h=4: 12M_1+12M_2+6M_3+M_4+M_(alpha_1=1,alpha_3=1).       (4)

In particular only fixed-point moments occur before surplus 3; the explicit
three-cycle variable becomes visible at surplus 3, five-cycle at surplus 5,
and a mixed 1/3 term at surplus 4. This does NOT mean longer cycles have no
influence at lower surplus: the normalized fixed-point moments themselves
reflect the full group action. It means no additional long-cycle coordinate
is needed at that horizon.

For C3 the surplus counts h=0..8 are
1,1,3,9,26,104,664,8248,208568. Removing its z3 term at h=3 gives 25/3,
not 9; the missing contribution is 2/3. Individual moments can be fractional,
although the final orbit counts are integers.

For the three prime roots of order 5 aggregated together, the surplus counts
are 3,11,44,154,539,2115,11374,113828,2528900. Their scalar count 3 alone
cannot supply these coefficients: the slot action is still necessary.

## 3. An exact triangular generating-function recursion

Let a_n be A000568, a_0=1; let s_n count nonempty strong classes, with the
singleton strong. Set E(x)=sum_(r>=1) a_(r+1)x^r. Let J_(k,h) be the summed
h-jet of the prime roots of arity k. Then

    s_n = [n=1] + sum_(k=3)^n [x^(n-k)]
              J_(k,n-k)(E(x),E(x^3),...),
    a_n = sum_(j=1)^n s_j a_(n-j).                            (5)

The second line is ordered SCC composition, t(x)=1/(1-s(x)).
In the first line a child has size at most n-k+1<=n-2, so only earlier
ordinary counts are used once the required prime jets are supplied. Root
arity n contributes its zeroth moment, the number of n-vertex primes.
Thus the missing new-prime information has NOT been wished away.
The all-order generating identity remains

    t(x)=1/(1-x-P(t(x)-1,t(x^3)-1,...)),

but (5) shows exactly which finite piece of P each coefficient consumes.
Prime jets can be obtained from the inherited plethystic inversion, hence
(5) is fully executable without a prime-class census; when so obtained from
Burnside input it retains the inherited noncircularity warning.

## 4. Metagraph boundary and useful next coordinate

The observer is the count of module-preserving inflations through a specified
surplus, not the adjacency observer. A new vertex with nonuniform incidence
inside an existing child module can destroy that module and change the prime
quotient. Such an attachment does not stay in the fixed-Q inflation family.
Neither J_h nor the full Z_Q retains the aligned masks needed to decide this.
Consequently this is a smaller exact **composition branch** carrier, not an
all-attachment or root-switch carrier. The constructive next bridge would
pair these jets with ports recording which old modules a new attachment
splits; no sufficiency theorem for that paired state is claimed here.

## 5. Reproduction

`inherited_recreate.py` is copied unchanged from the pinned Session86 script.
It is an explicit inherited dependency, not a new independent prime recovery.
Run in this directory:

    python3 verify_surplus_jets.py > receipt.json
    python3 -O verify_surplus_jets.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

Checks use explicit exceptions. The 81 horizon checks compare jet substitution
against direct full-index substitution for every root arity 3..15 and surplus
0..min(8,15-k), including empty order-four prime bank. Separate direct orbit
enumeration on TT1/TT2/TT3 palettes under cyclic groups C3,C5,C7 agrees at every
possible surplus, without calling the cycle-index routine to find the orbits.
A separate implementation of the triangular jet recurrence (5) reconstructs
all a_n through 15 and agrees with the full-index reconstruction. Every
reported inflation coefficient is checked integral before conversion.
These are checks of cyclic slot actions; no claim is made that every chosen
cyclic action is the full automorphism group of a specified prime root.
The C3 missed-three-cycle hostile is checked exactly. Normal and optimized
receipts are byte-identical. No files in the source repository were modified.
