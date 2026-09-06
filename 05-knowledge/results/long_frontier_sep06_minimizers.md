# All minimizing sequences for the sharp signed-root stability quotient

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
Promoted as [THM-4455 / three-atom minimizing-sequence rigidity](../../01-canon/theorems/THM-4455-three-atom-minimizing-sequence-rigidity.md).
This classifies
the limiting geometry of every minimizing sequence for
[THM-4454 / sharp global signed-root duplication stability](../../01-canon/theorems/THM-4454-sharp-global-signed-root-duplication-stability.md).
It does not transport the quotient into actual Laurent coefficient rows.

## 1. Statement

For each n let r^(n) be a finite real list satisfying

    p1 = sum r_i = 1, p2 = sum r_i^2 = 1,
    E = (1-p4)/2 > 0.

Use the definitions from THM-4454:

    D = (5-8p3+3p4)/6, J = D/E,
    c_* = (13-8sqrt(2))/3,
    d2 = 2-sqrt(2)(a+b), R = (J-c_*)/d2,
    K3 = 4sqrt(3)/[3(1+sqrt(2))(1+sqrt(3))],

where a>=b>=0 are the two largest positive entries, padded with zero.
Let c>=0 be the third largest positive entry. Define the squared distance
to the permutations of three equal positive atoms by

    Delta3 = inf_perm ||r-perm(1/sqrt(3),1/sqrt(3),1/sqrt(3),0,...)||_2^2
           = 2-(2/sqrt(3))(a+b+c).

Zero padding allows the infimum formula even when fewer than three
positive entries occur. Then the following are equivalent:

1. R(r^(n)) tends to K3.
2. Delta3(r^(n)) tends to zero.
3. a,b,c tend to 1/sqrt(3), and the square mass outside these three
   positive coordinates tends to zero.

The equivalence holds for arbitrary varying finite lengths and arbitrary
sign patterns. In particular there is no other asymptotic minimizing
shape, including at the zero-energy and zero-distance boundaries.

The signed sum of the remaining coordinates tends to 1-sqrt(3). Their
positive and negative sums separately need not converge: both may diverge.
The number of negative coordinates necessarily tends to infinity.

## 2. Inheritance and the missing boundary coordinate

The closest proved mechanism is THM-4454's complementary signed secant
and convex tail packing. Its actual sharpness family proves that the
three-positive-atom shape is attainable as a limit. The hostile examples
are the one-atom zero-energy boundary and the two-equal-positive-atom
zero-distance boundary: the unnormalized objective vanishes at both.
Thus merely classifying zeros of that objective cannot classify R.

The least-used sidecar is the denominator E*d2. We retain it in the two
boundary expansions below. The compactness argument uses only the two
largest positive coordinates and their total remaining square mass;
it does not silently identify weak convergence with norm convergence.
The final tail estimate recovers all lost norm explicitly.

The corrected near miss is to infer that the *negative* first moment
must tend to 1-sqrt(3) from square-norm convergence. Only the signed
remaining first moment has this forced limit. Section 7 supplies an
actual normalized hostile with divergent positive and negative dust.
This distinction was caught before promotion, not inherited as a theorem.

Targeted searches at inherited HEAD f0521b87281f for minimizing sequences,
near minimizers, concentration compactness and three-atom rigidity found
no existing classification of this quotient. THM-4454 explicitly leaves
the question open. No external priority claim is made.

## 3. The two degenerate boundaries cannot minimize

Write u=sqrt(2), h=1/u, A=2-u and B=u-1. The exact identity

    g = B-p3+A*p4 = (3E/4)(J-c_*)

will also control R.

### One positive atom

Suppose a tends to 1. Set delta=1-a^2. Every remaining coordinate has
absolute value at most sqrt(delta), so its signed cubic sum is
O(delta^(3/2)) and its fourth-power sum is O(delta^2), uniformly in list
length. Therefore

    p3 = 1-(3/2)delta+O(delta^(3/2)),
    p4 = 1-2delta+O(delta^2),
    D = delta+O(delta^(3/2)), E = delta+O(delta^2).

Eligible lists have delta>0. Hence J tends to 1, b tends to zero, and

    R tends to (1-c_*)/(2-u) = u-2/3 > 1/2 > K3.       (1)

### Two equal positive atoms

Suppose a,b tend to h. Put

    delta=1-a^2-b^2>=0, w=a-b>=0, eta=delta+w^2.

We have eta>0 on actual eligible lists: eta=0 would leave only two
positive roots h and contradict p1=1. The exact relation
a+b=sqrt(2-2delta-w^2) and the same tail moment bounds give, uniformly
over all ratios of delta to w^2,

    p3 = h-(3h/2)delta+(3h/4)w^2+O(eta^(3/2)),
    p4 = 1/2-delta+w^2+O(eta^2),
    d2 = delta+w^2/2+O(eta^2), E tends to 1/4,
    g = (3h/2-A)delta+(A-3h/4)w^2+O(eta^(3/2)).

Since delta+w^2/2>=eta/2, division is uniform. Set

    Ldust = (28u-32)/3,
    Lsplit = (64-44u)/3.

Then

    R = [Ldust*delta+Lsplit*w^2/2]/[delta+w^2/2]+o(1).

Here Ldust-Lsplit=24u-32>0, and Lsplit>1/2 (square the
equivalent positive comparison 125>88u). Consequently

    liminf R >= Lsplit > 1/2 > K3.                    (2)

THM-4454 proves K3<4/9, so the displayed separation is exact.
The formal two-atom family is used only for its boundary expansion;
it is not asserted to satisfy the actual first-moment constraint.

## 4. The two leading coordinates must tend to 1/sqrt(3)

Use z=1/sqrt(3), gamma=(2-u)(3-sqrt(3))/4 and

    C0=(1+u+sqrt(3)-u*sqrt(3))/2,
    C=C0-gamma(a+b),
    F_actual=1-C-p3+C*p4=(3E/4)*d2*(R-K3).

The quantities E and d2 are bounded, so R tending to K3 implies
F_actual tending to zero. Pass to any subsequence on which (a,b)
converges and which lies entirely on one side of b=z.

For b<=z, THM-4454 Section 3 gives F_actual>=F_env(a,b)>=0.
The continuous envelope has precisely the zeros (z,z) and (1,0) in
its closed domain. The second is excluded by (1).

For b>=z, put v=sqrt(1-a^2-b^2). THM-4454 Section 4 gives
F_actual>=F_pack(a,b,v)>=0. Its only zeros in
a>=b>=z>=v>=0, a^2+b^2+v^2=1 are

    (a,b,v)=(z,z,z), (h,h,0).                         (3)

For completeness, the strict concavity in a+b at fixed v forces any
zero to an interval endpoint. On b=z the strictly positive secant
boundary factor excludes a>z. On a=b=x the factorization is
v^2(x-z) times a strictly positive brace, so x=z or v=0. These are
exactly (3); this uses the strictness already proved in THM-4454,
not an equality assumption in the actual signed comparison.
The second triple is excluded by (2).

Every convergent subsequence of (a,b) therefore has limit (z,z).
Compactness in these two real coordinates proves a,b tend to z.

## 5. The third positive atom captures all remaining norm

Let m=1-a^2-b^2, so m tends to 1/3. Let c be the third largest
positive coordinate, padded with zero. Every positive tail coordinate
is at most c<=sqrt(m). For all sufficiently large n, the function

    f_C(t)=t-C*t^2

has derivative at least a fixed positive epsilon on [0,sqrt(m)],
because C tends to (3-sqrt(3))/2 and its limiting minimum derivative
on [0,z] is 2-sqrt(3)>0. Negative tail coordinates have f_C(t)<0,
while f_C(c)>=0. Thus

    sum_tail (r_i^3-C*r_i^4) <= m*f_C(c).

Define the continuous three-positive-coordinate expression

    F3=1-C-a^3-b^3+C(a^4+b^4)-m*f_C(sqrt(m)).

It tends to zero at (a,b,sqrt(m))=(z,z,z). Hence

    0 <= epsilon*m*(sqrt(m)-c)
       <= F_actual-F3 -> 0.

It follows that c tends to z. The three positive square masses now
sum to one in the limit, so the entire remaining square mass tends
to zero. This proves (1) implies (3) in the theorem statement without
an assumption on bounded support size or signed dust mass.

Conversely, square-norm convergence to a permutation of three z atoms
implies convergence of p3 and p4, and of the top-two sum a+b. For
example p_k is k-Lipschitz in the square norm on the unit ball for
k=3,4 by integrating its gradient; the top-two positive sum is
sqrt(2)-Lipschitz. Consequently

    E -> 1/3, d2 -> 2-2sqrt(2/3)>0, R -> K3.

The formula for Delta3 shows that (2) and (3) are equivalent. This
completes the classification.

## 6. Negative dimension must diverge; a quantitative count consequence

Let N_minus be the number of strictly negative entries and let
Delta3 be as above. The total negative magnitude is at least
max(a+b+c-1,0), while the total negative square mass is at most
Delta3. Cauchy--Schwarz therefore gives the explicit bound

    N_minus*Delta3 >=
      max(sqrt(3)-1-(sqrt(3)/2)*Delta3,0)^2.             (4)

Thus every minimizing sequence has N_minus tending to infinity.
This is a consequence of its actual p1 normalization, rather than
an assumption that all small entries have the same sign.

## 7. An actual hostile to separate first-moment rigidity

For each integer n>=4 set N=n^4, and let a_n be the unique root in
(1/2,z) of

    3a^2+[n^2+(n+3a-1)^2]/N=1.                       (5)

Existence and uniqueness follow because the left side is increasing
for a>=0, is less than one at a=1/2 and greater than one at a=z.
The former follows from
n^2+(n+1/2)^2 < (n+1)^2*2 <= (25/8)n^2 < n^4/4
for n>=4. The last inequality is strict already at n=4.

Take the actual finite list of three entries a_n, N entries n/N,
and N entries -(n+3a_n-1)/N. Equation (5) gives p2=1; its signed
sum is exactly 1. Both dust magnitudes are below a_n, since
(n+3a_n-1)/N<(n+1)/n^4<=5/256<1/2. The three leading coordinates
tend to z and all remaining square mass tends to zero.

By the proved equivalence, R tends to K3. Yet the positive dust
first moment is n, and the negative magnitude is n+3a_n-1. Both
diverge. Their signed difference tends to 1-sqrt(3), as required.
Thus neither nonnegative-dust absence nor convergence of the negative
first moment is an additional consequence of optimality.

## 8. Verification and scope

The proof's universal inputs are the audited analytic inequalities of
THM-4454, uniform tail moment bounds, elementary compactness in two
coordinates, and the exact identities and expansions displayed here.
The [independent proof/source audit](long_frontier_sep06_minimizers_audit.md)
passes without correction. The [exact source](../../04-computation/long_frontier_sep06_minimizers.py)
checks91 always-active gates: symbolic quotient and uniform boundary
identities, exact radical signs, and seven actual algebraic mixed-dust
lists indexed by n=4,5,6,8,12,20,50, with no exclusions. Normal and
optimized execution match the [frozen output](long_frontier_sep06_minimizers.out)
byte-for-byte. No numerical experiment is used to infer the universal
classification.

```bash
python3 -B 04-computation/long_frontier_sep06_minimizers.py
python3 -B -O 04-computation/long_frontier_sep06_minimizers.py
```

Raw LF hashes are source
6d3be5780138c24f58312733270981d3375eae93b9281b26710529bf955dc7cd
and output
19d6c54228a78db3b7231e30e9d6897bdc9a71ba904dd5159b9c861ab2040b20.

The new source-to-target map retains the actual moment normalization,
the quotient denominator and the ordered leading roots. Square-norm
convergence preserves p3,p4 and the quotient at the three-atom limit,
but destroys separate positive and negative first-moment information;
the signed sum and estimate (4) are the necessary first-moment sidecars.
