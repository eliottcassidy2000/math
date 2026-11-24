# 7, 9, 15: the Pythagorean parent functional and a Collatz boundary

2026-11-24. Continuation of the repo's Lemoine same-pair projection,
[weighted diagonal and finite census](../lemoine_weighted_diagonal_20261124/),
and its split/inert and 61/63/65 themes. Order triples by odd leg a,
even leg b, and hypotenuse c. Let T(x)=oddpart(3x−1) on positive odds.

## 1. The 3-4-5 identification is rigid, and extends to the entire tree

Write (a,b,c)=(m²−n²,2mn,m²+n²), where m>n>0 are coprime and
opposite parity. On every primitive triple,

    b+c=(m+n)²,      c−b=(m−n)²,
    k := 3c−2(a+b) = (m−2n)²+n².

Thus b+c=a² exactly on the consecutive-parameter spine m=n+1.
Its first triple gives a+b=7 and b+c=a²=9.

The more surprising invariant is **k = the parent's hypotenuse**.
Here is a direct derivation of the Berggren descent. Apart from
the root (m,n)=(2,1), form the parent by ordering (n,|m−2n|)
decreasingly. Depending on m/n the inverse of that parent map is
one of the three branches, written on the parent pair (u,v):

    A: (2u−v,u),  1<m/n<2;
    B: (2u+v,u),  2<m/n<3;
    C: (u+2v,v),  3<m/n.

The two reduced boundary ratios 2 and 3 give the root or violate
opposite parity. Each branch preserves primitivity and parity;
the parent hypotenuse decreases strictly since |m−2n|<m.
So every pair has a unique path to the root and the formula for k
is exact. The root can be regarded as having a degenerate parent
of hypotenuse 1. On the child triple we have

    (3c−k)/2 = a+b,

and that division by 2 is exact in the odd Collatz sense:
3c−k=2(a+b) with a+b odd. Thus the same affine 3x-minus-intercept
branch applies on each generation; its intercept is **not** fixed
to 1 except at the root. On the root, T(5)=7 and T(7)=5.

Indeed, fixing k=1 forces (m−2n)²+n²=1 and yields only the root.
Equivalently, put U=a+b and V=a−b; the right-triangle cone
U²+V²=2c² and the linear constraint 2U=3c−k give

    (c−3k)² + (2V)² = 8k².

At k=1 the only positive nondegenerate solution up to leg swap
is c=5, |V|=1. More generally, if (3U−k)/4=c as well, elimination
forces c=5k and U=7k, giving only (3k,4k,5k); primitivity forces k=1.

### Split/inert and multiplicity

The intercepts are exactly 1 and the positive integers whose
prime divisors all split in Z[i] (are 1 mod 4). More precisely:
for k>1 supported on split primes there are
3·2^(ω(k)−1) primitive triples with intercept k, otherwise zero.
The descent identifies the fiber with the three children of
each primitive triple of hypotenuse k. In Gaussian factorization
of a primitive norm k, each split prime power must be allocated
entirely to one of its two conjugate factors; allocation to both
would divide both legs of the Euclid pair. Units and simultaneous
conjugation identify the sign/leg order, leaving 2^(ω(k)−1).
An inert odd divisor of k would force both Euclid parameters to
be divisible by it, and even k is incompatible with their opposite
parity. This is an exact use of the repo's split/inert distinction,
not a claim that splitting predicts Collatz basins.

## 2. The first three children and the Lemoine transition at 21

The ordered children of (3,4,5) are

    A: (5,12,13),  B: (21,20,29),  C: (15,8,17).

In general on the consecutive-parameter spine, the child odd legs
are a+2, a(2a+1), a(2a−1). This follows directly from the three
parameter branches: A preserves m−n=1 and adds 2 to m+n, while
B and C multiply the new difference and sum to give the two products.
Therefore the old Lemoine diagonal seeds 9=3·3, 15=3·5 and
21=3·7 use the old triangle's odd leg, hypotenuse, and leg sum,
respectively. Exactly at a=3, the B/C leg multipliers 2a±1
agree with a+b and c respectively: insert b=(a²−1)/2 and
c=(a²+1)/2: the B equality gives (a+1)(a−3)=0, and
the C equality gives (a−1)(a−3)=0. The positive a=1 case is degenerate.
So the two odd legs 15 and 21 are not arbitrary matches.

There is a direct Lemoine reason for the small exceptions.
Remove q=2 and diagonal p=q:
* 7 is below 3+2·3; the only unrestricted pair is (3,2).
* 9 is the corner (3,3), plus (5,2) on the boundary.
* In a multiple of 3, p+2q=0 forces p≡q mod 3. If one is 3,
  both are 3; otherwise both are ≥5. At 15, equality forces
  (5,5); the remaining unrestricted pair is (11,2).
* For distinct odd primes in a common nonzero class mod 3, the
  first two in residue 2 are 5,11 and in residue 1 are 7,13.
  The minimum weighted sum uses the smaller in the doubled
  slot. Hence the first recovered multiple of 3 is
  21=11+2·5 (the residue-1 lane cannot open before 27).
  This coexists with the old diagonal 7+2·7.

Thus 21 is exactly the first off-diagonal recovery in the mod-3
lane even though it is a tournament count hole. This distinguishes
the two spectra rather than explaining one by the other. The preceding
repo census found only {7,9,15} through 10^8; it does not show that
there are no further removed-diagonal exceptions.

One more rigidity at 15: it is the unique positive odd N>4 with
T(N)=N−4, so the accelerated Collatz target literally equals the
prime left after removing the doubled 2. If 3N−1=2^j(N−4),
then (2^j−3)(N−4)=11. For j≥2 the only possible exponent is 2;
j=1 would give N=−7. This gives N=15 and T(15)=11.

## 3. What actually enters the known Collatz cycles

The three known positive 3x−1 cycles, in the odd convention, are

    (1), (5,7), (17,25,37,55,41,61,91).

This is the negative 3x+1 problem under sign reversal; the list
is not claimed exhaustive. See [Lagarias, §2.6](https://web.williams.edu/Mathematics/sjmiller/public_html/372Fa15/addcomments/Lagarias_3x%2B1AndItsGeneralizations.pdf).
There is nevertheless an exact interpretation of all three labels
as cycle statistics, provided one says which clock is used. For
each odd-to-odd edge write j=v2(3x−1). It contributes one step to
the odd clock, j steps to the shortcut clock (combine the multiply
and first divide), and j+1 to the raw clock (separate every divide).
The valuation words from the displayed representatives are
(1), (1,2), (1,1,1,2,1,1,4). Thus:

| cycle | odd steps | shortcut steps | raw steps |
|---|---:|---:|---:|
| through 1 | 1 | 1 | 2 |
| through 5 | 2 | 3 | 5 |
| through 17 | 7 | 11 | 18 |
| total | 10 | 15 | 25 |

In particular **7** is the long odd period, **9** counts the odd
sites of both nontrivial cycles together, and **15** counts all
sites of the three cycles in the shortcut map. These are exact
identities, not the same statistic applied to three cycles.
There is also a precise 9+14=23 split for the two nontrivial cycles
(odd+shortcut=raw); this gives a clean meaning to those values
without making them tournament counts.

Our three inputs are **not** representatives of three basins:

    7 → 5 → 7,
    9 → 13 → 19 → 7,
    15 → 11 → 1.

The absent basin enters as 17, paired with 15 in the C-child
(15,8,17). The other child hypotenuses follow 13→19→7 and
29→43→1. The hypotenuse labels of the three children do thus
meet all three known basins, one each.

Here is an infinite exact prefix law. For every e≥1,

    T^e(2^e+1) = oddpart(3^e+1).

For the first e−1 steps the state is 3^i·2^(e−i)+1 and the
raw numerator has v2=1; the last gives the stated odd part.
It explains 9→7 in exactly three steps, since 9=2³+1 and
3³+1=4·7. On the even-parameter ray of Euclid pairs (2^j,1),
the first three primitive triples are (3,4,5), (15,8,17),
and (63,16,65), and the two endpoints have the exact prefix laws

    T^j(4^j−1)  = oddpart(3^j−1),
    T^(2j)(4^j+1) = (3^(2j)+1)/2.

On the minus side each of the first j−1 steps divides by 4;
the last takes the remaining odd part. At j=2 this sends 15 to 1
and 17 after four steps to 41 on its own cycle. At j=3, 63
reaches 13 then the 5-cycle, and 65 reaches 365 then eventually
1. The repo's neighboring 61 already lies on the 17-cycle.
So {61,63,65}, unlike {7,9,15}, actually samples the three known
basins once each. This is a finite exact observation with an
infinite prefix rule, not a proposed rule for the final basin.

## 4. What the two linked papers contribute

[The Weil-ground-state paper, §8.6](https://arxiv.org/html/2611.13129v1#S8.SS6)
keeps the transition and pole modes when eliminating the
coercive block by a Schur complement. As a literal finite
analogue, form the directed Collatz graph U from the three known
cycles and the routes of 7,9,15 to them. Give each directed edge
unit energy (x_v−x_T(v))². Fix measurements x7,x9,x15 and
minimize over the other coordinates. The reduced form is exactly

    E_eff = (x9−x7)²/3.

The three-edge chain from 9 to 7 sets the denominator; every
other measured-basin vertex can be constant at its terminal
value. The 17 component has an undetermined constant, so an
unqualified inverse on the eliminated block would be illegal.
Ground it, quotient its constant, or keep 17 as a fourth
terminal. The full form has three zero modes whereas these
three terminal observations see only two. For the composition
operator P(f)(v)=f(T(v)), the same finite set gives
det(I−zP)=(1−z)(1−z²)(1−z⁷): order transients before the
cycle blocks and triangularize. These are graph identities,
not a claim about the Weil eigenfunctions or RH.

[The local-smoothing note, §§1.2, 2.2.2, 3.3–3.5](https://arxiv.org/html/2611.30274v1)
keeps interference in a positive matrix before taking the
quaternionic determinant, and retains a normal cofactor
under the Schur decomposition. The relevant small exact
model here is a positive, rank-one light-cone matrix:

    [[c+b, a], [a, c−b]] = (m+n,m−n)^t (m+n,m−n).

Its two diagonal readings are squares but they do not replace
the marked off-diagonal reading a. For the parent functional,
the separate positive matrix [[1,−2],[−2,5]] has determinant
1 and Schur complement 1; completing its square exposes the
coordinate m−2n whose sign/order chamber chooses the branch.
The transferable idea is to keep the joint block until after
elimination, not that the PDE or its critical estimate applies
to Collatz or primes.

## Replay

python3 replay.py checks all primitive Euclid pairs with m≤360,
and compares the fiber sizes up to intercept 10,000 against the
independent Gaussian count. It also checks the finite routes,
the Collatz formulas, and the Lemoine boundary. No arbitrary
completion of the functional graph into a tournament is used;
the missing arcs are essential to the reduced energy.
