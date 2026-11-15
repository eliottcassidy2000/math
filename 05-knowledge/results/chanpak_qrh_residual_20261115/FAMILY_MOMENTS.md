# What low-phase averaging can and cannot correct

The earlier [quasi-Riemann / ROOT session](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/preprint_connections_solfeggio_digits_20261017.md) separated *family mean square*, *replication*, and *scale contraction*. Here the analogous object is the residual on a family of **literal complementary phases**. Its moments are exactly computable, but the family mean is not itself a tournament. The statements below do not assume the quasi-Riemann hypothesis or import a bound for L-functions.

Fix g=2^s ≥2 and high parts u,v≥1. For b uniform in [0,g), write bar b=g−1−b, x=gu+b, y=gv+bar b. Remove the exact product and linear portion P(x,y) and the center constant k_j from the three-window count. Normalize the remaining correction by 45:

    T_b=p(x)+p(y)+(parity(x)+2 parity(y))/3.

Write X_j=1−2b_j for the independent Rademacher bits of b, numbered j=0,…,s−1 from the most significant end. Write w_d=2^max(0,d−2). From the one-apex checksum,

    p(b)+p(bar b)=sum_(i<j) w_(j−i) (1−X_i X_j).

This is because adding the two weighted inversion sums counts each unequal digit pair exactly once; retaining the factor two turns the XOR into 1−X_i X_j. The weights sum to g/2−1 (group by distance). From the [complement transfer](COMPLEMENT.md),

    A_b=(g−1)/2 + (1/2) sum_j 2^j X_j,
    z_b=(1+X_0)/2.

Finally x,y are of opposite parity, and the outer parity term of T is 1/2+X_(s−1)/6. Thus one has the **exact** mean

    E T = p(u)+p(v)+(g−1)(u+v+1)/2
                        +(parity(u)+parity(v))/2

and centered expansion

    T_b−E T = sum_j c_j X_j − sum_(i<j) w_(j−i) X_i X_j,
    c_j = (u−v) 2^j / 2
          +(parity(u)−parity(v)) [j=0]/2
          +[j=s−1]/6.

The degree-one and degree-two Walsh characters are orthogonal, so the **exact variance** is

    Var T = sum_j c_j² + sum_(i<j) w_(j−i)²
          = sum_j c_j² + (4^(s−1)+6s−7)/9.

The weight identity follows by summing (s−d)4^max(0,d−2) for d=1,…,s−1 and holds also for s=1. In particular as g grows the standard deviation is already Ω(g) when the high halves agree. If u≠v its linear part also has standard deviation Ω(g |u−v|), by the top coefficient, up to the bounded parity term. Accordingly full-orbit second moments alone do not supply a power contraction of the residual scale O(g+g|u−v|), even after the exact smooth part P has been removed and the exact mean has been charged to the baseline.

Under b↔bar b the degree-one part flips sign and the degree-two part stays fixed. This cancels the asymmetric term exactly in the two-element mean, but neither that mean nor a cancellation between distinct graph counts is necessarily physical. The needed new input is *target-conditioned phase selection* or a new literal gadget; a smaller average is not an exact selector. Compare the earlier [guard-composition session](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_complement_routing_20261004.md), which likewise kept the final legality check separate from identities in the auxiliary composed family.

This obstruction concerns the **full complementary phase orbit** only, not a specially chosen subset of phases, the whole three-window model, or Chan–Pak itself.

Replay: `python moments.py` verifies the expansions and both moments with exact rational arithmetic on 72 phase families.
