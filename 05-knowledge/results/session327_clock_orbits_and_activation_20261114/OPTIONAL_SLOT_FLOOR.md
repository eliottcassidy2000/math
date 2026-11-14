# Busch's shared vertex has a slot cost: a new source-grade floor

Let Q_T(z) count ordered covers of all vertices of a tournament by z directed paths, with empty slots allowed. In particular H=Q(1). This note gives an all-order strong-tournament lower bound for every integer z≥2 and exact small-order activation thresholds for the slope-99 source machinery.

## 1. The optional-slot merge inequality

Suppose A,B induce subtournaments of T, meet in exactly v, and together cover T. Then

    Q_T(z) ≥ Q_A(z) Q_B(z) / z.

To see this, take ordered covers of A and B in which v occupies the same slot. Slot permutation symmetry makes their number Q_A(z)Q_B(z)/z. In each other slot use the path-merging insertion lemma from [Busch](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/) to merge the two paths, preserving the order on each factor. In the common slot split each path at v and merge the two prefixes and the two suffixes separately. In each case the directed join to v is legal since each candidate neighboring endpoint was legal in its original path. The resulting cover restricts to the original two slot covers, so this is an injection. The divisor z is the precise alignment constraint in this construction; setting z=1 recovers the usual merge lower bound for H.

## 2. A uniform grade bound for the two-piece clock and higher slots

For integer z≥2, put b=z²+2 and d=z(z²+4). For any strong tournament T of order n≥3,

    Q_T(z) ≥ z b^((n−1)/2)                 if n is odd,
    Q_T(z) ≥ z d b^((n−4)/2)               if n is even.

At orders three and four the exact normalized counts Q/z are b and d. For n≥5, Busch's structural lemma gives either a near-transitive tournament (one reversed extreme edge of a transitive order), or two induced strong subtournaments meeting in one vertex and spanning T. For the overlap, normalized Q/z is supermultiplicative by the previous section. Put l=n−1. If both factor lengths l_i are even, or one even and one odd, multiplication gives the claimed bound with equality in the algebra. If both are odd, use d²≥b³, since their difference is 2z⁴+4z²−8>0. This handles the overlap induction.

For the near-transitive case, an odd cycle uses the reversed edge and an odd subset of the interior vertices. Since all cycles share that edge, the odd-cycle packing formula gives

    Q_R(z)/z = z^(n−1) + (z+1)^(n−2) − (z−1)^(n−2).

Let m=n−2≥1. The difference of this expression at n+2 and b times its value at n is

    −2z^(m+1) + (2z−1)(z+1)^m + (2z+1)(z−1)^m > 0,

since already (2z−1)(1+1/z)^m≥2z+1−1/z>2z. Induction in steps of two from n=3,4 finishes the proof. Only the lower bound is asserted for general order, not equality or classification of extremizers.

In particular Q(2) has normalized growth base √6 per unit of glued length, with exact lower bound 2·6^((n−1)/2) at odd n, and 32·6^((n−4)/2) at even n. At n=7 this gives Q≥432 rather than the one-triangle or three-triangle floor used in the previous atlas. Thus no source with target N<99·432=42768 can contain a strong SCC of order at least seven, without any finite tournament census. This is why a chamber around 30,000 can be computed from small recursive pieces.

## 3. The first six exact activation levels

At orders 3 through 8, the small literal catalogue yields the following exact minima for Q(2). At each order the same displayed minimizing pair also attains Busch's exact minimum for H. The replay enumerates every isomorphism class at these orders by extending all classes at the preceding order and canonically choosing within outdegree classes.

| new SCC order n | minimum (H,Q2) simultaneously | first target N allowing this or a larger source SCC, actual 99 digit bank |
|---:|---:|---:|
| 3 | (3,12) | 1191 |
| 4 | (5,32) | 3173 |
| 5 | (9,84) | 8325 |
| 6 | (15,216) | 21399 |
| 7 | (25,548) | 54277 |
| 8 | (45,1416) | 140229 |

For the lower direction, any source containing a strong SCC C has H≥H(C), Q≥Q(C) by ordinal multiplication. Each digit is positive, and the actual bank phase is nonnegative. For n=3,…,8 the analytic lower bound for all larger orders exceeds the tabulated Q minimum at n (and the Busch H floor increases), so no larger order can enter earlier. For the upper direction, take the minimum core alone in the c=1, phase-zero chart (the literal S325 bank contains it). Its exact output is 99Q+H. These are activation levels of a *source grade* in this chart, not thresholds of full Hamiltonian-spectrum coverage. In particular the change near 30,000 is not a new n=7 core: the whole chamber from the first n=6 entry at 21,399 up to the first n=7 entry at 54,277 allows at most six-vertex strong components. Changing the target scale changes these bounds in a controlled way.

This is one exact meaning of least resistance: before optimizing coefficients, respect the combinatorial alignment cost at the recursive seam. It is a different scale law from the primorial's wheel period, even where both guide the same examples.
