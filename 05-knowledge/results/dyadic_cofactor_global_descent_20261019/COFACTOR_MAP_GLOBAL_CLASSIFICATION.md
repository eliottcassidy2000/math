# A global classification of the dyadic cofactor map

19 October 2026. **PROVED, with independent audit PASS.** This map agrees with the aliquot map only on its explicitly retained prime domain. The global classification below is not a claim that general aliquot or Collatz orbits converge.

## 1. The map and the strengthened aliquot consequence

For each positive even integer n, write n=cp with c=2^a, a>=1, and p>=1 odd. Define

    F(n)=(c-1)p+2c-1 = n-p+2c-1.

F maps positive even integers to positive even integers. When p is an odd prime, F(n)=sigma(n)-n exactly. When p=1 or p is composite, that equality is not assumed.

**Theorem.** Every orbit of F eventually reaches one of:

- a fixed point c(2c-1), c a power of two at least two;
- the four-cycle 8 -> 22 -> 14 -> 10 -> 8;
- the three-cycle 16 -> 46 -> 26 -> 16.

These are all finite cycles, and there are no divergent F-orbits.

**Aliquot corollary.** No nontrivial aliquot cycle of any length can consist entirely of numbers 2^a p with a>=1 and p an odd prime. An aliquot trajectory that remains in this class forever must eventually reach an even perfect number.

Indeed the genuine trajectory equals its F trajectory as long as the prime guard holds. Both nontrivial F cycles contain a power of two, whose odd part is 1 and is outside the prime domain. The remaining possible fixed point is c(2c-1), with 2c-1 prime, which is exactly the usual even-perfect form.

This strengthens the separately proved all-exponent three-cycle exclusion, whose size/congruence proof remains valid. The new argument uses short, exact blocks of the actual cofactor-map orbit, rather than a fixed cycle-length template.

## 2. A short-block descent statement

Except for fixed points and the five values

    {2,4,8,12,16},

every positive even n has F^j(n)<n for some 1<=j<=6. The exception set is a set of starting values for this descent statement, not a list of cycle members.

We prove this by cases. All comparisons concern exact integers.

Write n=cp as above and put q=p+1. Then

    F(n)-n=2c-q.

If q>2c, there is immediate descent. If q=2c, n is fixed. We therefore only need the abundant case q<2c.

## 3. Non-power-of-two q gives descent in two steps

Suppose q is not a power of two. Write q=dk, where d=2^{v_2(q)}>=2 and k>=3 is odd. Since q<2c and c is a power of two, d<c and in fact 2d divides c.

The first image is

    n_1=c(p+2)-p-1=c(q+1)-q.

Its exact dyadic cofactor is d, since

    n_1/d=(c/d)(q+1)-k

is odd. A direct second step gives

    d(F²(n)-n)=(c+d)(2d-p-2)+p+1.

Because p+1=dk with k>=3, we have 2d-p-2<=-d-1. Thus

    d(F²(n)-n)
    <=-(c+d)(d+1)+p+1
    <-(c+d)(d+1)+2c
     =c(1-d)-d(d+1)<0.

The strict inequality uses q<2c. Hence F²(n)<n.

## 4. The power-of-two starting values

If p=1, n=c is a power of two. For every C=2^a>=8,

    F(C)=3C-2,
    F²(C)=3C/2+2,
    F³(C)=3C/4+4.                                  (1)

The first two intermediate values have exact dyadic cofactor 2: for the second, C>=8 is used. Hence the formulas are genuine consecutive F steps. For C>=32, (1) is strictly below C. The remaining powers 2,4,8,16 belong to the declared exception set.

From now on p>=3 and q=d=2^t>=4. Abundance q<2c forces c>=d, since both are powers of two.

## 5. Power-of-two q with d<c<d²

Suppose d<c<d² and put e=c/d. Then e is a power of two with 2<=e<=d/2. Directly,

    n_1=d[e(d+1)-1],
    n_2=e(d²-1)+d.

The exact dyadic cofactors are d at n_1 and e at n_2. The second follows because d/e is an even integer. The next step gives

    n_3=e d²-d²+e+d-d/e,
    n_3-n=-d(d-e)+e+d-d/e.

Since e<=d/2 and d/e>=2,

    n_3-n <= -d²/2+3d/2-2<0   (d>=4).

Thus there is descent within three steps.

## 6. Power-of-two q with c>d²

Put h=c/d², a power of two with h>=2. The same first two formulas now give

    n_1=d[(hd)(d+1)-1],
    n_2=d[h(d²-1)+1].

Both have exact dyadic cofactor d. A third step yields

    n_3-n=-h(d-1)+3d-2.                            (2)

For h>=4 this is at most -d+2<0.

For h=2, equation (2) instead says n_3=n+d. Since

    n_3=d[2d(d-1)+1],

its exact dyadic cofactor is again d. A fourth step gives

    n_4-n=-2d²+5d-2=-(2d-1)(d-2)<0.

This proves descent within four steps throughout c>d².

## 7. The two resonant equalities c=d and c=d²

If c=d, then n=d(d-1) and its first image is d², a power of two. For d>=8, apply (1) with C=d² to get

    F⁴(n)-n = -d²/4+d+4<0.

The remaining d=4 gives n=12, an exception.

If c=d², then n=d²(d-1), and the first two images are

    n_1=d[d(d+1)-1],
    n_2=d³.

For d>=8, apply (1) with C=d³:

    F⁵(n)-n = -d³/4+d²+4<0.

For the remaining d=4, n=48 has the exact six-step chain

    48 -> 76 -> 64 -> 190 -> 98 -> 52 -> 46,

so it too has a smaller iterate within six steps. These exhaust the power-of-two q case with p>=3.

## 8. Global conclusion

Sections 2–7 prove the short-block descent statement. The exceptional trajectories are explicit:

    2 -> 4 -> 10 -> 8 -> 22 -> 14 -> 10,
    12 -> 16 -> 46 -> 26 -> 16.

The same displayed cycles cover starting values 8 and 16. The fixed points are exactly q=2c, equivalently p=2c-1.

Starting anywhere else, pass to a smaller iterate within six steps. Repeating this operation cannot create an infinite strictly decreasing sequence of positive integers. It must reach a fixed point or one of the exceptional starting values, whose fate has just been listed. This proves global eventual periodicity and excludes additional cycles.

For the native aliquot corollary, retain the odd-prime condition at every step. The formulas sigma(p)=p+1 used to identify F with the aliquot sum fail at p=1 and at composite p; no conclusion about the aliquot orbit after leaving the class is made. In particular the formal fixed point 120=8·15 is not an aliquot fixed point, and the two displayed nontrivial F cycles are not aliquot cycles.

## 9. Relation to the proof program

The source formula is inherited from the repository's aliquot core compiler. The earlier three-cycle argument combines size, congruence and valuation to close its A3/D3 frontier. Here the same retained dyadic cofactor supplies a uniformly bounded descent block, with an explicit finite exceptional set. That extra source-specific mechanism is what turns local affine contraction into a global statement for F.

This does not imply such a descent-block theorem for Collatz. The step formulas, residue transitions, exceptions and native prime domain above are all specific and proved. Historical priority for this auxiliary-map classification or its aliquot consequence has not been established.
