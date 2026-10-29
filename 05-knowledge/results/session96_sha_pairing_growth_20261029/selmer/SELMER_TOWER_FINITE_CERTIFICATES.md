# Selmer towers: exact finite stopping certificates and the pairing-limit obstruction

**Status:** elementary group-theoretic statements proved below; arithmetic use is a criterion requiring certified descent bounds, not a proof that such bounds always become available. No claim of arithmetic realization is made for the abstract hostile groups.

## 1. Arithmetic input and a concrete target

For an abelian variety A over a number field K, the standard descent sequence is

    0 -> A(K)/p^n A(K) -> Sel_(p^n)(A/K) -> Sha(A/K)[p^n] -> 0.

The Selmer groups are finite. The Cassels–Tate pairing has divisible kernels; passing to the quotient by the maximal divisible subgroup is essential. See [Milne, Arithmetic Duality, §3](https://www.jmilne.org/math/xnotes/Tate100.pdf), and [Poonen–Stoll, introduction and §2](https://math.mit.edu/~poonen/papers/sha.pdf). For elliptic curves, [Poonen's ICM survey, §2.3](https://math.mit.edu/~poonen/papers/icm2018.pdf) also gives the finite and p-infinity exact sequences. These are cited arithmetic inputs, not new results.

A tractable target is therefore: certify a single exact plateau in the *Sha torsion layers*, after subtracting the actual Mordell–Weil contribution. The square-or-twice-square theorem does not supply that plateau.

## 2. One plateau proves fixed-prime finiteness

Let G be a p-primary torsion abelian group with G[p] finite. Put G_n=G[p^n], G_0=0, and

    b_n=log_p |G_n|,        a_n=b_n-b_(n-1), n>=1.

All G_n are finite: multiplication by p on G_n has kernel G_1 and image in G_(n-1), giving |G_n|<=|G_1|^n.

Multiplication by p induces an injection

    G_(n+1)/G_n -> G_n/G_(n-1).

If px lies in G_(n-1), then p^n x=0, so x lies in G_n; this proves injectivity. The quotients are elementary abelian p-groups. Hence

    a_1 >= a_2 >= ... >= 0.                              (1)

**Plateau criterion.** If b_(n+1)=b_n for any n>=0, then G=G_n is finite. Indeed equality of the nested finite subgroups gives G_(n+1)=G_n. An element of any larger p-power order would yield, by multiplication by a suitable power of p, an element of exact order p^(n+1), a contradiction. Conversely, a finite G eventually has such a plateau.

This is a finite positive certificate. It does not give a universal finite cutoff at which every still-growing sequence can be classified as finite or infinite.

## 3. The exact recurrence profile and its missing stopping clock

Since (1) is a nonincreasing sequence of nonnegative integers, it stabilizes at some r>=0. Put

    m_e=a_e-a_(e+1), e>=1.

Only finitely many m_e are nonzero, and direct summation yields

    b_n = r n + sum_(e>=1) m_e min(n,e).                 (2)

No classification theorem for infinite abelian groups is needed for this cardinality statement. In the usual cofinite p-primary decomposition, r is the divisible corank and m_e the multiplicity of cyclic factors of order p^e; equation (2) itself follows just from the injections above.

Consequences:

* G is finite iff r=0, equivalently iff some a_n is zero.
* The jump m_e is recovered exactly from a second finite difference.
* With E=max{e:m_e!=0}, b_n is affine for n>=E.
* The count sequence C_n=|G_n| then obeys C_(n+1)=p^r C_n for n>=E.
* Its generating series is rational, but the transient endpoint E is part of the data.

For the logarithmic counts,

    sum_(n>=0) b_n z^n
      = [r z + sum_e m_e z(1-z^e)]/(1-z)^2.             (3)

An eventual one-step recurrence does not bound E. Taking G=(Z/p^E)^2 makes b_n=2 min(n,E), with arbitrarily late stopping despite square order at every layer. Taking G=(Q_p/Z_p)^2 makes b_n=2n forever. The same initial count prefix can belong to either case.

## 4. A certified Selmer-ratio stopping test

Write A(K)=Z^R plus its finite torsion subgroup T, and define

    S_n=|Sel_(p^n)(A/K)|,
    t_n=log_p |T[p^n]|.

The finite descent sequence gives the exact identity

    log_p S_n = nR+t_n+b_n,                             (4)

where now G=Sha(A/K)[p^infinity]. Thus

    S_(n+1)/S_n
       =p^R * |T[p^(n+1)]|/|T[p^n]| * |G_(n+1)|/|G_n|. (5)

If R and T are certified, equality of (5) to its first two factors proves G=G_n. A guessed rank or an unproved Selmer estimate cannot be inserted into this equality.

More usefully, a rank *lower bound* can suffice. Suppose:

1. r0 is a certified lower bound for R, obtained for example from independently verified rational points;
2. the torsion ratio tau_n=|T[p^(n+1)]|/|T[p^n]| is certified;
3. L_n<=S_n and U_(n+1)>=S_(n+1) are certified positive bounds.

Then the finite inequality

    U_(n+1) <= L_n p^r0 tau_n                           (6)

forces **R=r0 and Sha(A/K)[p^infinity]=Sha(A/K)[p^n] finite**.

Proof: the actual ratio S_(n+1)/S_n is at most U_(n+1)/L_n, while (5), R>=r0 and G_n subset G_(n+1) make it at least p^r0 tau_n. Equality throughout is forced. In particular the rank difference is zero and the Sha layer ratio is one, so §2 applies. This proof also shows why both bound directions matter.

Equation (6) is a concrete target for a verified higher-descent computation. The session does not produce those arithmetic input certificates for a new curve, nor prove that a search for them terminates.

At n=0 take S_0=L_0=1 and T[p^0]=0 as a group of order one. Then the familiar case of a tight p-Selmer bound proves G[p]=0, hence G=0. At larger n the criterion can certify nontrivial finite p-primary Sha. The [classical two-descent calibration](../example/EXACT_TWO_DESCENT_CONTROL.md) for y²=x³-x supplies an actual arithmetic n=0 case: rank lower bound zero, four rational two-torsion points and exact two-Selmer order four give Sha[2]=0 and therefore Sha[2^infinity]=0. This calibration is classical, not a new finiteness result.

### Wrong-rank hostile

A formal tower with R=0, T=0, and G=Q_p/Z_p has S_n=p^n. Subtracting the unsupported rank guess one would produce a fake zero Sha growth rate. A certified rank lower bound cannot make this mistake, since one is not <=R. This is an abstract exact-sequence model, not a claim that such a Sha occurs for an abelian variety.

## 5. Perfect finite pairings do not pass isometrically to the divisible tower

Let H_n=(Z/p^n)^2 and give it the perfect alternating pairing

    beta_n((x1,x2),(y1,y2))=(x1 y2-x2 y1)/p^n mod Z.

Its order is the square p^(2n), with no bound as n grows. The natural injections

    i_n:H_n -> H_(n+1),       x -> p x

have direct limit (Q_p/Z_p)^2. Their pairing transport is

    beta_(n+1)(i_n x,i_n y)=p beta_n(x,y),              (7)

not beta_n(x,y). A change of unit cannot repair this: any perfect alternating pairing on H_n has coefficient u_n/p^n with p not dividing u_n. Isometric compatibility would require

    p u_(n+1) = u_n mod p^n,

impossible after reduction modulo p. Thus one cannot silently glue these finite perfect pairings into a perfect pairing on the direct limit with unchanged target Q/Z.

In fact **every bilinear pairing D x D -> Q/Z is zero when D is a divisible torsion group**. Given y of order m and x, choose z with mz=x. Then B(x,y)=B(z,my)=0. This elementary observation matches the need to quotient by divisible kernels in the arithmetic theorem; it does not assert that the divisor is trivial.

The stage clock also controls how much of a finite ambient pairing is visible. Embed H_n in H_E for E>=n by multiplying coordinates by p^(E-n). The restricted pairing is

    p^(E-2n)(x1 y2-x2 y1) mod Z.                        (8)

It vanishes identically if E>=2n. In general its radical has order

    p^(2 min(n,E-n)).                                  (9)

For E<2n, the coefficient form shows that a radical vector has both coordinates divisible by p^(2n-E), giving p^(E-n) choices per coordinate. For E>=2n every vector is radical. At n=E the ambient perfect pairing is recovered.

## 6. A finite-depth hostile that retains group maps and pairings

Fix a cutoff N and choose E>=2N. Compare:

* finite H_E with its perfect alternating pairing;
* infinite D=(Q_p/Z_p)^2 with its necessarily zero bilinear pairing.

Through levels n<=N, their p^n-torsion groups, inclusion maps and restricted pairings are identical: the groups are H_n, the maps are x->px, and (8) makes every observed pairing zero. One ambient group is finite and one infinite. The finite group has square order; the infinite group's quotient by its divisible radical is trivial with its nondegenerate pairing.

One may orthogonally append the same fixed finite paired group to both models. For p=2, appending C2 with self-pairing 1/2 gives the twice-square finite version. The hostile still agrees through every prescribed finite cutoff. These are group-theoretic models satisfying the stated pairing patterns, not constructed abelian varieties or counterexamples to Sha finiteness.

This does not conflict with the plateau certificate: neither hostile has reached a plateau in the observed window. It proves that a bounded growing prefix, even with all its restricted pairings, gives no uniform negative answer about eventual stopping.

## 7. Fixed-prime versus all-prime debt

For a torsion abelian group, global finiteness requires both:

* every p-primary component is finite;
* only finitely many primes have nonzero components.

The second condition cannot be dropped. The abstract direct sum over all primes of (Z/p)^2 has finite square-order p-primary components and an alternating pairing formed by summing its finite-support hyperbolic pairings, yet it is infinite. Its pairing has zero radical but its algebraic dual is much larger; nondegeneracy is not finite perfect duality. This is again a structural warning, not an asserted Sha example.

Thus a proof of general Sha finiteness needs arithmetic control of both late p-power layers and prime support. The square/twice-square classification constrains the finite nonsingular pairing; it does not pay either unboundedness debt by itself.

## 8. Relation to the repository's recurrence/carrier boundary

[Session95](../../session95_niche_recurrences_20261029/README.md) distinguishes a short scalar recurrence from a faithful carrier for all future operations. Here the precise source is a torsion filtration with actual inclusion maps and restricted pairings. The map to b_n retains torsion-layer cardinalities but loses pairing transport and embedding depth E. The finite-depth hostile shows that even restoring those data at bounded depth still leaves a future height obligation. This is a common logical mechanism, not an isomorphism of Sha with a tournament grammar.

`check_selmer_towers.py` verifies finite torsion layers, monotone increments, exact moment/recurrence formulas, the certified-ratio implication on explicit model towers, and pairing scaling/radicals. All checks use exact integers and remain active under Python optimization.
