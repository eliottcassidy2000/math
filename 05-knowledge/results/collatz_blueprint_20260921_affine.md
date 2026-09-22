# Blueprint audit: the exact guarded affine structure

**Status: PROVED elementary algebra and domain results; FINITE-EXACT replay.**
The proposed triangle-group and odd-perfect-number arguments fail. No Collatz
convergence theorem or complete cycle classification is established here.

Inheritance: the [summand note](arithmetic_braids_20260917_summand.md) supplies
the dyadic forest and shortcut map; the [signed-cycle note](arithmetic_braids2_20260917_signed_cycles.md)
supplies the ordered carry and warns that sign conjugacy is not arrow reversal.
The closest hostile is loss of the ordered carry despite equal exponent counts.
The least-used sidecar here is the integer-domain guard on an affine map.
This note rederives these mechanisms for the attachment's specific matrices;
it makes no novelty claim. Parent integration owns correction-ledger entries.

The live board has four objects: partial integer inverse maps (legality),
projective real matrices (scale-invariant spectrum), ordered words (carry),
and divisor sums (multiplicative prime-power data). The anchor is faithful
inverse Collatz; the niche is affine group versus semigroup structure; the
wildcard is the proposed divisor-sum transfer. Dropping guards or prime-power
data is the explicit loss in the latter two proposed bridges.

## 1. Domains come before residue arrows

Write `C(n)=n/2` for even positive `n`, and `(3n+1)/2` for odd positive `n`.
Its inverse branches are

```text
D(n)=2n                         (all positive n),
E(n)=(2n-1)/3                   (n=2 mod3).
M3(n)=(4n-1)/3 = E(D(n))        (n=1 mod3).
```

Thus the attachment's third generator is redundant. D produces even
predecessors; E produces odd predecessors. D does not act on the odd states
alone. On all six residue classes D is not a permutation: `D(1)=D(4)=2 mod6`.

Neither division branch descends to a function modulo six:

```text
E(5)=3, E(11)=7;       5=11 mod6 but 3!=7 mod6.
M3(1)=1, M3(7)=9;      1=7 mod6 but 1!=9 mod6.
```

Input modulo 18 is sufficient to determine output modulo six. More
precisely, replacing n by n+6 changes E(n) by 4 and M3(n) by 8; replacing
n by n+18 changes them by 12 and 24. The three lifts of a legal mod-six
input therefore give three distinct odd output classes. Modulo nine is
already sufficient for each branch's output modulo six; modulo 18 is the
natural lift that also retains input parity. Keeping only mod-six labels
loses one ternary digit, and iteration requires further digits.

For an odd target u, the accelerated predecessors are `(2^k u-1)/3`.
If u=2 mod3 then k is odd; if u=1 mod3 then k is even; if 3 divides u
there is no odd predecessor. These are the guarded words `E o D^(k-1)`.
Increasing k by two changes the predecessor by `n -> 4n+1`, recovering
the previously proved inverse-fibre braid. None of this forces a basin.

## 2. Projective normalization rules out the claimed triangle group

Use real projective matrices

```text
A1 = [[2,0],[0,1]], A2 = [[2,-1],[0,3]], A3 = [[4,-1],[0,3]].
```

For positive determinant the scale-invariant classifier is
`J(A)=tr(A)^2/det(A)`. Determinant-one normalization gives elliptic for
`J<4`, parabolic for `J=4` except the identity, and hyperbolic for `J>4`.
Here

```text
J(A1)=9/2, J(A2*A1)=49/12, J(A3*A1)=121/24.
```

All three tested elements are hyperbolic. Bare traces are not projective
invariants: multiplying a representative by a nonzero scalar leaves its
map unchanged. In particular `det(A3)=12`, whereas `det(A3/2)=3`.

All generators fix infinity and are orientation-preserving affine maps.
In fact their unrestricted generated group is exactly

```text
{x -> a*x+t : a=2^i*3^j, i,j in Z, t in Z[1/6]}.
```

Proof: the displayed class is a group containing D and E. Conversely
`D E D^-1 E^-1(x)=x-1/3`. Its inverse cubed gives translation by one.
Composing E on the left with translation by 1/3 gives pure scaling by
2/3; together with D this gives all slopes `2^i 3^j`. Conjugating
translation by one by these pure scalings gives `1/(2^a 3^b)` for all
a,b>=0; integer sums produce all translations in `Z[1/6]`.
This is a metabelian affine group, and it has no nonidentity finite-order
elements: a finite-order positive slope must be one, and a nonzero
translation has infinite order.

It is also **not discrete**: `D^-j [D,E] D^j` translates by
`-1/(3*2^j)` and tends to the identity. Consequently it is not the
Fuchsian triangle group Delta(3,4,infinity). The latter's order-three
and order-four elliptic generators cannot occur here. The commutator
uses inverses in the unrestricted real group; it is not asserted to be
a legal positive-integer inverse-Collatz path.

No action or quotient map to S6 is supplied. An actual six-dimensional
sphere cannot be a space in bijection with Z, since it is uncountable.
The standard Fuchsian action is on a real two-dimensional hyperbolic
plane; its existence would not identify that plane's orbifold with S6.

## 3. The carry is a complete legality and cycle coordinate

Let `w=g1...gr` be a chronological word in D,E: g1 acts first. Let m
count its E letters. Define B initially zero, updating in chronological
order by

```text
D: B <- 2B;
E: B <- 2B+3^m, then m <- m+1.
```

Direct matrix multiplication proves

```text
W(n)=(2^r n-B)/3^m.                                    (1)
```

**PROVED iff:** for positive integer n, all intermediate guards are
valid iff `2^r n = B mod3^m`. In particular there is exactly one legal
source class modulo `3^m` (the whole class modulo one if m=0).

Necessity follows from integer evaluation. For sufficiency, let Bj,mj
be the carry and E-count after j steps. Later E contributions are
divisible by `3^mj`, so `B=2^(r-j)Bj mod3^mj`. Reducing the final
congruence modulo `3^mj` and cancelling the unit `2^(r-j)` proves
`3^mj | 2^j n-Bj`. Thus every prefix is integral. At an E step this
forces the previous integer to be 2 mod3; a positive such integer is
at least two, so the output is positive. Induction proves all guards.

For a nonempty word, `2^r != 3^m`. If m>0 then B>0. Therefore

```text
W has a positive integer fixed point
iff m>0, 2^r>3^m, and (2^r-3^m) divides B;
then the fixed point is n=B/(2^r-3^m).                  (2)
```

The congruence in (1) is automatic at that integer fixed point, so (2)
certifies every intermediate step, not just a rational return. For m=0
the only affine fixed point is zero, outside the positive domain.
Repeating a word repeats a cycle; word length need not be minimal period.

**Spectral hostile.** `w=DEDE` gives `(16n-7)/9`, fixing 1; `w=DDEE`
gives `(16n-5)/9`, fixing only 5/7. Both matrices have trace 25,
determinant 144, and `J=625/144`. Even the full eigenvalue pair forgets
the carry that decides integer closure. Baker-type exponent estimates
alone cannot replace the missing divisibility condition.

In fact 5/7 belongs to the nontrivial positive rational accelerated cycle
`5/7 -> 11/7 -> 5/7`, with exact halving exponents 1 and 3 (valuation
extends to rationals with odd denominator). Clearing the denominator
gives the integer `3n+7` cycle `5 -> 11 -> 5`. This refutes the blueprint's
literal rational uniqueness claim, not the integer Collatz conjecture.

**Free-semigroup survivor.** Distinct D/E words define distinct affine
maps. Equality first forces equal r,m by unique prime factorization of
the slopes, and then equal B. Choose a positive integer in their common
legal source class. Starting at its common endpoint, parity uniquely
recovers the last letter: even means D, odd means E. Applying C reverses
that letter. Repeating for r steps recovers the complete word. Thus the
semigroup is free although its containing group is solvable. This is
freeness of transformations, not a free action on integers: E o D fixes 1.

**Exact ternary/dyadic progression correspondence.** If s is the least
positive legal source of w, then `1<=s<=3^m` and positivity plus (1) give
`1<=W(s)<=2^r`. Thus w is a bijection between the positive progression
`s+3^m*j` and `W(s)+2^r*j`, for j>=0. Its endpoint class is exactly

```text
y = -3^(-m)*B mod2^r.
```

For fixed r these endpoint classes, one per D/E word, partition all
`2^r` residue classes: distinct words cannot share an endpoint, because
its r backward parity decisions recover the word. This is a precise
bridge between the ternary source guard and the dyadic parity address.
Uniform counting of endpoint residues gives every length-r word mass
`2^(-r)` and E-count distribution `binomial(r,m)/2^r`. It does not assert
any limiting frequency along a fixed positive integer trajectory.

## 4. What survives from the odd-perfect discussion

The Euler and Touchard restrictions are valid necessary conditions:
`N=pi^alpha*m^2`, with pi prime, `gcd(pi,m)=1`, `pi=alpha=1 mod4`, and
`N=1 mod12` or `N=9 mod36`. A primary corroborating source is
[Starni, Some Extensions to Touchard's Theorem](https://arxiv.org/abs/1709.05286).
For completeness, the residue deduction is elementary: if 3 does not
divide N, then pi=2 mod3 would make `sigma(pi^alpha)=0 mod3`, impossible
because `sigma(N)=2N`. Thus N=1 mod3 and N=1 mod4. If 3 divides N,
pi is not 3, so its exponent in the square factor is even and 9 divides N;
combine this with N=1 mod4.

The next implication fails. A projective determinant can be rescaled
arbitrarily and has no established relation to divisibility of N.
Moreover, M3 is not integral at any N=3 mod6, the very state the
argument attempts to force. A missing accelerated Collatz predecessor
does not constrain the separate divisor equation `sigma(N)=2N`.
Its multiplicative prime-power data have no specified transport into
these affine maps. No odd-perfect-number nonexistence result follows.

## 5. Replay and remaining proof obligation

Run `python 04-computation/experiments/collatz_blueprint_20260921_affine.py`.
The [stdlib script](../../04-computation/experiments/collatz_blueprint_20260921_affine.py)
checks all D/E words through length ten, source/guard equivalence,
complete source residue classes through length six, every endpoint
residue partition through length ten, the arithmetic-progression
bijections, independent rational and matrix evaluation, parity decoding,
hyperbolic invariants, fixed points,
and the explicit hostile controls. It uses exact integers/Fractions and
explicit checks that remain active under `python -O`. Its
[matching JSON](../../04-computation/experiments/collatz_blueprint_20260921_affine.json)
is a finite replay, not an all-cycle census. Normal and optimized runs
agree byte-for-byte: 2,047 words, 5,461 complete-residue checks, and all
listed hostile and positive controls pass.

Connection contract: ordered guarded inverse words map to triangular
matrices with their carry and legal source class; this preserves exact
integer return. Projecting further to slope/trace loses B and legality;
the four-letter hostile is decisive. The forward global obligation
remains proving that every positive orbit reaches the known cycle,
including excluding infinite escape and all additional positive cycles.
