# Quartic AP-square walls force nonregular Berggren ancestry languages

**Status: PROVED, with FINITE-EXACT implementation controls.** No theorem number or external priority is claimed.

Companions: [checker](../../04-computation/experiments/berggren_ap_square_nonregular_20261019.py), [output](berggren_ap_square_nonregular_20261019.out).

This extends [the AP-square/Pythagorean classification](ap_square_pythagorean_20261019.md).

## Inheritance and scope

The closest established mechanism is THM-3536,
`01-canon/theorems/THM-3536-berggren-angle-languages-signed-c4-and-harmonic-density.md`:
finite inverse wall orbits, with orientation retained, give regular ancestry
languages. This note proves a converse criterion for a single irrational cut,
and uses it to obstruct regularity for the AP-square cone. It does not claim
that every quadratic wall has a finite inverse orbit, does not import a
continued-fraction periodicity theorem, and proves no arithmetic-height or
word-density formula. The canonical hostile control is the quadratic fixed
wall `1+sqrt(2)`, whose middle-branch orientation really reverses.

The active concept board is: (1) the AP-square cone polynomial, (2) outer-first
word residuals, (3) opposite-parity Euclid slopes, (4) finite wall-orbit automata,
and (5) the orientation sign. The source-to-target connection sends a PPT to
its Euclid slope and its Berggren word, preserving cone membership. Forgetting
word order or orientation destroys the residual recurrence; both are retained.

## 1. Coding and density

Let `J=(1,infinity)`, root `t0=2`, and

    L(t)=2-1/t,   M(t)=2+1/t,   R(t)=t+2.

These are bijections from `J` onto `(1,2)`, `(2,3)`, and `(3,infinity)`,
respectively. Only `M` reverses orientation. For a root-to-leaf word
`w=w1...wk`, its slope is

    E(w)=f_wk o ... o f_w1 (2).

For an outer-to-inner word put

    H(w)=f_w1 o ... o f_wk (2)=E(reverse(w)).

The root orbit is exactly

    S={n/m>1: gcd(n,m)=1 and n,m have opposite parity}.

Indeed all three forward maps preserve these properties. Conversely, if
`t=n/m` belongs to `S` and is not `2`, the unique inverse branch is

    T(t)=1/(2-t)  for 1<t<2,
         1/(t-2)  for 2<t<3,
         t-2      for 3<t.

The resulting numerator/denominator pairs are, respectively,

    (m,2m-n),   (m,n-2m),   (n-2m,m).

They remain positive, coprime, opposite-parity pairs with first coordinate
larger than the second, and their sum strictly decreases. The seam `3=3/1`
is excluded by parity. The only remaining seam is the root `2=2/1`, so descent
terminates there. Finally, reduced fractions with denominator `2^k` and odd
numerator are in `S`; their spacing `2^(1-k)` tends to zero. Thus `S` is dense
in `J`.

## 2. Exact finite-orbit criterion

For irrational `eta in J`, define outer-first cut languages

    D_eta^-={w:H(w)<eta},   D_eta^+={w:H(w)>eta}.

**Proposition.** `D_eta^-` is regular if and only if the actual inverse orbit
`{T^j(eta):j>=0}` is finite. The same is true of `D_eta^+`.

All inverse iterates remain irrational: the branch maps are invertible
rational Möbius transformations. Consequently no iterate is a seam `2` or
`3`, so its branch itinerary is unambiguous.

Write `eta_j=T^j(eta)`, and let `d_j` be the unique branch with
`eta_j in f_d_j(J)`. Let `p_n=d_0...d_(n-1)`, and let `sigma_n` be `-` when
this prefix contains an even number of `M`s and `+` otherwise. Direct inversion
gives the exact left-quotient identity

    p_n^(-1) D_eta^- = D_eta_n^(sigma_n).                 (1)

Here `p^(-1)D={v:pv in D}` is language quotient notation. The word is read
outermost letter first; no root-to-leaf/outer-to-inner identification is made.

Different pairs `(eta,sign)` give different cut languages. For two different
thresholds with the same sign, density of `S` supplies a slope strictly between
them. For opposite signs, density supplies a slope in `(1,min(eta,eta'))`,
which belongs to the lower-cut language and not the upper-cut language. Hence
all residuals in (1) have exactly their indicated threshold and orientation.

If `D_eta^-` is recognized by a finite automaton, its left quotients are finite
in number. Equation (1) and the preceding distinguishability imply that its
inverse orbit is finite.

Conversely, if the inverse orbit is finite, use undecided states
`(eta_j,sign)` and two sinks. The unique branch containing the current wall
updates the wall by `T` and flips the sign precisely for `M`. Either other
branch lies wholly on one side of the wall and goes to the appropriate sink.
At the end of a word, an undecided state accepts exactly when the rational
root `2` satisfies its indicated inequality. This finite automaton recognizes
the cut. (It uses at most twice the number of inverse-orbit points plus two
states.) This proves the proposition.

## 3. Algebraic obstruction, without a coding-periodicity assumption

**Corollary.** If a wall `eta` is neither rational nor quadratic over `Q`,
its cut language is nonregular. This includes every algebraic wall of degree
greater than two and every transcendental wall.

A finite inverse orbit has `eta_m=eta_n` for some `m<n`. By reconstruction,
`eta_m` is fixed by the nonempty branch composition

    F=f_d_m o ... o f_d_(n-1).

Every branch is represented by an integer matrix of determinant `+1` or `-1`:

    L=[[2,-1],[1,0]],  M=[[2,1],[1,0]],  R=[[1,2],[0,1]].

The composition `F` cannot be the identity, because
`F(J)` is contained in the proper branch interval `f_d_m(J)`.
If `F(t)=(a t+b)/(c t+d)`, its fixed point satisfies the nonzero polynomial

    c t^2+(d-a)t-b=0.

Thus `eta_m` has degree at most two. Since `eta` is another rational invertible
Möbius transform of `eta_m`, it has the same degree, a contradiction.

This proves what is needed directly from repeated actual thresholds. It does
not infer a repeated threshold merely from repetition of symbolic letters,
and does not require contraction of parabolic cylinders or any Lagrange
continued-fraction theorem.

## 4. Application to the AP-square cone

For the Euclid PPT

    (a,b,C)=(n^2-m^2,2nm,n^2+m^2),   t=n/m>1,

the cone condition `C^2>4ab` is exactly

    P(t)=t^4-8t^3+2t^2+8t+1>0.                         (2)

Writing `t=cot(theta)` with `0<theta<pi/4`, the ratio `4ab/C^2` is
`2 sin(4theta)`. Consequently (2) is equivalent to

    t in (1,alpha) union (beta,infinity),
    alpha=cot(5pi/24),   beta=cot(pi/24).

In particular, `1<alpha<2<3<beta`. Both walls are roots of `P`.

The polynomial `P` is irreducible over `Q`. Its only possible integral roots
are `+1,-1`, and `P(1)=P(-1)=4`, so there is no linear factor. By Gauss's lemma,
a factorization must therefore be into monic integral quadratics. Their
constant terms multiply to one, so they have the form

    (t^2+u t+epsilon)(t^2+v t+epsilon),  epsilon=+1 or -1.

Comparing coefficients gives `u+v=-8`, `epsilon(u+v)=8`, hence
`epsilon=-1`, and then `uv=4`. This would make `u,v` integral roots of
`z^2+8z+4`, whose discriminant is `48`, a contradiction. Thus both walls have
degree four. Equality in (2) never occurs at a rational tree slope.

Let the root-to-leaf admissibility language be

    A={w:E(w) in (1,alpha) union (beta,infinity)}.

Suppose `A` were regular. Regular languages are closed under reversal, so

    K=reverse(A)={w:H(w) in (1,alpha) union (beta,infinity)}

would be regular. Since `L(J)=(1,2)` excludes the entire upper cone,

    L^(-1)K={v:L(H(v))<alpha}
           =D_(1/(2-alpha))^- .                        (3)

A left quotient of a regular language is regular. But the wall
`1/(2-alpha)` has degree four, since `alpha=2-1/(1/(2-alpha))` recovers it by
a rational Möbius transformation. Section 3 contradicts (3).

**Conclusion.** The finite root-to-leaf Berggren words whose terminal PPT
satisfies `C^2>4ab` form a nonregular language. Its complement is also
nonregular. A DFA, NFA, regular expression, or finite-state wall machine
cannot exactly classify this cone at all depths.

## 5. What this result does and does not transfer

- THM-3536's angle-wall languages are regular because that theorem verifies
  finite four-point inverse orbits, including orientation signs.
- The AP-square walls have degree four, so no repeated actual inverse state
  can occur. This is an exact all-depth obstruction, not a conclusion from
  seeing many distinct numeric iterates.
- The union of two cuts could in principle hide a nonregular cut; equation (3)
  prevents that cancellation here by isolating one boundary in one quotient.
- Nonregularity does not prohibit an exact algorithm using unbounded memory,
  a polynomial-time numerical/exact slope test, a uniform short constructive
  suffix into the cone, a rational limiting frequency, or some other counting
  recurrence. None of those stronger obstructions is claimed.
- This note does not assert the converse “quadratic wall implies regular”.
  The proved equivalence is with finite actual inverse orbit.

## 6. Reproduction

Run the companion with ordinary and optimized Python:

    python3 -B 04-computation/experiments/berggren_ap_square_nonregular_20261019.py
    python3 -B -O 04-computation/experiments/berggren_ap_square_nonregular_20261019.py

The companion uses exact rational arithmetic, explicit failure checks rather
than `assert`, rational bisection for the quartic wall, direct root/outer word
checks, inverse residual checks, a quadratic positive control, and a wrong-
orientation hostile control. Its finite transcript is evidence for the
implementation, not the proof of infinite nonregularity.
