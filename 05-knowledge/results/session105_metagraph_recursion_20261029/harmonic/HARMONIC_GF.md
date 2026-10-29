# Deletion harmonics: exact generating transform and doubled pole ladder

**PROVED relative to the attachment/deletion algebra and the inherited Session98 global pole theorem.** This is an analytic refinement for an explicitly defined class-space statistic. It is not a new enumeration of modular-prime tournament classes.

The source-free primitive count and complete integral deletion classification are **inherited** from [October19 Hasse/deletion work, §5](https://github.com/eliottcassidy2000/math/blob/c49370cdbf6c331a735a17c690a8e865beb3ac2b/05-knowledge/results/tournament_hasse_frobenius_20261019/PADDING_FROBENIUS_PROOF.md), blob `9f035d2324589187ee3fb4bc3ff1b3131eb4a917`. This packet's new assembly combines that count with the normalized uniform-refresh spectrum and the newer Session98 global pole theorem. The old source-free induced-count basis is not being claimed as a new discovery.

## 1. A recursively defined dimension

Let a_n=A000568(n), a_0=1, a_(-1)=0. The [attachment/deletion theorem](../updown/UPDOWN_HARMONIC_AUDIT.md) defines H_n=ker D_n on the rational vector space of tournament isomorphism classes. Normalized attachment E and deletion D obey [D,E]=I. Hence

    h_n=dim H_n=a_n-a_(n-1),
    V_n=direct-sum_(j=0)^n E^(n-j) H_j.

The first h_n are 1,0,0,1,2,8,44,400,6424,184656,9541520. In particular the 44 at order six is a deletion-harmonic dimension, not the number of prime tournaments. Nonzero harmonics in positive degree have both positive and negative class coefficients. Their dimension nevertheless has a separate positive interpretation: it counts source-free tournament classes, as proved below. This does not identify those class basis vectors themselves with elements of ker D.

As formal ordinary series, H_raw(z)=(1-z)A_raw(z). Both have radius zero. This identity by itself is not a smaller closed scalar recurrence for a_n: one still needs the h_n.

There is nevertheless an exact metagraph return identity. For the vertex-refresh operator K_n=E_(n-1)D_n/n,

    n trace(K_n)=sum_(j=0)^(n-1) a_j,
    a_n=(n+1)trace(K_(n+1))-n trace(K_n),               (1)

where 0 trace(K_0) is defined as zero. Proof: its eigenvalue (n-j)/n has multiplicity h_j; summing and telescoping gives the first formula, then subtract. This is a class-space trace with the exact multiplicity-weighted operator, not an average over distinct metagraph neighbours. It supplies a global consistency test and a different counting interface, not a claimed speedup for computing the traces.

### Positive class interpretation: strip the unique universal source

A tournament has at most one universal source. Adjoining a universal source to any (n−1)-vertex tournament and removing that source are inverse operations on isomorphism classes. Therefore h_n is also the number of n-vertex tournaments **without a universal source**, with the empty tournament included at n=0. Repeated source stripping gives a unique initial transitive string followed by a source-free core. This is an elementary structural interpretation, not the unrestricted modular-prime atlas and not a claim that source-free basis vectors already lie in H_n.

## 2. Normalized generating transform

Write

    beta_n=n! a_n / 2^binomial(n,2),
    eta_n=n! h_n / 2^binomial(n,2),
    B(z)=sum beta_n z^n,   J(z)=sum eta_n z^n.

For n>=1,

    eta_n=beta_n-n 2^(1-n) beta_(n-1).

Consequently, initially as power series for |z|<1,

    J(z)=B(z)-z B(z/2)-(z²/2) B'(z/2).                (2)

The coefficient derivation is exact: set m=n-1 in the subtracted sum and use sum(m+1)beta_m w^m=B(w)+wB'(w). Eta_0=1. The first positive-degree coefficients eta_1 and eta_2 are zero, while eta_3=eta_4=3/4 and eta_5=15/16, consistent with the first intransitivity degree three.

## 3. Every power of two is now a pole, with a precise order

The inherited theorem [Session98 global pole expansion](../../session98_tournament_growth_20261029/limits/ALL_POLES_CORE_EXPANSION.md) proves that B extends meromorphically to C with poles exactly 4^k, order 3k+1, k>=0. Applying (2) gives a meromorphic continuation of J.

At z=4^k, only B(z) is singular: z/2 cannot be a power of four. Thus the pole remains of exact order 3k+1.

At z=2·4^k, B(z) is analytic. If B(w) has order p=3k+1 at w=4^k, B'(w) has order p+1 with nonzero leading coefficient. The factor z²/2 is nonzero at the pole. The other shifted term has only order p, so cannot cancel that leading pole. Thus J has exact order 3k+2 there.

There are no other poles. Equivalently, at every z=2^r, r>=0, the exact pole order is

    floor(3r/2)+1.                                    (3)

This is an interlacing of an inherited square-scale ladder and its differentiated doubling translate. It does not identify these analytic scales with Collatz time or prove an arithmetic convergence theorem. Since there are infinitely many finite poles, J is not D-finite and eta_n is not P-recursive, by the usual finite-singularity property of a solution of a polynomial-coefficient linear differential equation.

## 4. Each symmetry core gives an explicit paired term

For an inherited odd nontrivial-cycle core mu, let m=sum mu, c=number of parts,

    E_mu=sum_i (mu_i-1)/2+sum_(i<j)gcd(mu_i,mu_j),
    z_mu=product_p p^(m_p)m_p!,
    A_mu=(m!/z_mu) 2^(E_mu-binom(m,2)),
    lambda_mu=2^(-(m-c)).

Its term in B is A_mu z^m/(1-lambda_mu z)^(m+1). Substitution into (2) yields the exact pair

    A_mu z^m/(1-lambda_mu z)^(m+1)
      -(m+1) A_mu 2^(-m) z^(m+1)
          /(1-lambda_mu z/2)^(m+2).                    (4)

This makes the additional pole explicit for each core, rather than only at aggregate level. The empty core contributes

    1/(1-z)-z/(1-z/2)^2.

Its coefficient for n>=1 is 1-n 2^(1-n), explaining the first cancellation at n=1,2. The nontrivial cores begin at n=3. The negative sign of the second rational summand is an algebraic subtraction, not a probability. In fact every full paired core has nonnegative Taylor coefficients: for n=m+f with f>=1, its subtracted/original coefficient ratio is f*2^(1-c-f)<=1. Independently, h_n is a dimension by the Weyl algebra argument.

The nonnegativity has an exact Burnside meaning. For an odd permutation type with nontrivial core mu and f fixed vertices, every universal source must be fixed by the permutation, because it is unique. A specified fixed vertex is a source when c+f−1 independent edge-orbit bits point outward. These f source events are disjoint. Thus the proportion of fixed tournaments with a source is f·2^(1−c−f), and the source-free count multiplies that cycle-type contribution by

    1−f·2^(1−c−f),

with factor one for f=0 (including the empty object). This is precisely the ratio in the paired-core Taylor coefficient. It gives a nonnegative Burnside refinement for h_n, while individual weighted cycle-type contributions need not be integers. Positivity does not turn the signed harmonic projection into a positive class basis.

The inherited normally convergent core expansion may be differentiated on compact sets away from its poles. Alternatively the meromorphic identity (2) already establishes (3), without termwise differentiation. No new global convergence assumption is being smuggled into (4).

## 5. What the recurrence improves

The structural improvement is a canonical all-degree deletion decomposition and a complete refresh spectrum, with an analytic signature computed from the established A000568 series. Formula (1) recovers counts from return traces; (2) transports the existing meromorphic theory to the new harmonic dimension. Neither bypasses the full class count by declaring h_n to be a pre-enumerated prime bank.

The separate enumeration lane improves arithmetic realization of the projective Burnside recurrence and honestly compares it with the faster existing prefix engine. The native Hamilton-path metagraph lane retains path-word data because its class-pair quotient fails extension lumpability. The three operators share a subject, not interchangeable transition rules.
