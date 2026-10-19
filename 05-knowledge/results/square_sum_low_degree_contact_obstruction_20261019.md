# Low-degree rational contact maps: degree one is trivial; degree two must split over Q(i)

**Status: PROVED scoped, with exact formula controls.** This note develops the initial lower-degree obstruction and exhaustive pole reductions. The complete [degree-minimality theorem](square_sum_cubic_contact_degree_minimality_20261019.md) now closes the remaining cases and proves that useful cubic maps are minimal.

## 1. Statement

Put

\[
 f(t)=\frac{4t(t^2-1)}{(t^2+1)^2},\qquad
 k_\sigma(t)=\frac{t^2+2\sigma t-1}{t^2+1},\qquad
 k_\sigma(t)^2=1+\sigma f(t),\quad\sigma\in\{\pm1\}.
\]

A global rational contact map is a rational function `h in Q(t)` such that, for fixed `epsilon,eta in {±1}`, some `y in Q(t)` satisfies

\[
 y^2=1+\varepsilon f(t)+\eta f(h(t))
     =k_\varepsilon(t)^2+k_\eta(h(t))^2-1.             \tag{1}
\]

**Theorem 1.** If `deg h=1`, then

\[
 f(h(t))=-\varepsilon\eta f(t),\qquad y^2=1.           \tag{2}
\]

Thus no degree-one map produces a distinct positive-spread petal contact on an interval: if both spreads are positive, (2) forces opposite ports and identical spreads, hence the same unordered primitive Pythagorean petal. Constant maps give only the zero-spread boundary.

**Theorem 2.** If `deg h=2`, write `h=P/Q` with coprime rational polynomials of maximum degree two. A necessary condition for (1) is that

\[
 P^2+Q^2\text{ is reducible over }\mathbb Q.             \tag{3}
\]

More precisely it is a product of two quadratic factors, with multiplicity allowed, and all its roots belong to `Q(i)`. Equivalently, `P+iQ` splits into linear factors over `Q(i)`.

This narrows the degree-two case to a Gaussian split-pole classification. The complete nonexistence proof is assembled in the linked degree-minimality theorem.

## 2. A residue lemma

**Lemma 3.** Let `x,z,y in Q(t)` satisfy `y²=x²+z²−1`. Suppose the finite poles of `x` and `z` are disjoint, all simple, with each pole set one irreducible Galois orbit over `Q`; suppose also all three functions are finite at infinity. Then

\[
 y=\alpha x+\beta z+c,\quad \alpha,\beta\in\{\pm1\},
 \quad c\in\mathbb Q.                                \tag{4}
\]

**Proof.** At a pole of `x`, the leading coefficient of `y²` is the square of the residue of `x`; hence the residue of `y` is `+` or `-` that of `x`. Galois conjugation preserves the quotient of these residues, which is one of the rational numbers `±1`. Irreducibility makes the sign the same at every pole of `x`. The same argument applies to `z`. Subtracting the corresponding signed functions removes every finite pole. The difference is bounded at infinity and therefore constant, necessarily rational. QED.

If (4) holds, set `X=alpha x`, `Z=beta z`. Squaring gives

\[
 2XZ+2c(X+Z)+c^2+1=0,
 \qquad
 Z=T_c(X):=-\frac{2cX+c^2+1}{2(X+c)}.                 \tag{5}
\]

For `c²=1`, the right side is constant; otherwise it is a degree-one fractional linear transformation.

## 3. Degree one with disjoint pole pairs

For a degree-one rational `h`, the poles of `x=k_epsilon(t)` are `±i`; the poles of `z=k_eta(h(t))` are the two conjugate points `h^{-1}(±i)`, another irreducible quadratic pair. Rational coefficients ensure the pairs are either equal or disjoint. All are simple, and infinity is not a pole because a rational fractional linear map takes infinity to a rational number or infinity, never to `±i`.

In the disjoint case Lemma 3 applies. Both `X` and `Z` have degree two, with exactly the same unordered branch-value set

\[
 \{-\sqrt2,+\sqrt2\}.                                \tag{6}
\]

This follows either by differentiating `k_sigma`, or by solving `k_sigma(t)=v`: the quadratic discriminant is `4(2−v²)`. Precomposition by `h`, or multiplication by a sign, does not change the unordered branch values.

The constant case `c²=1` in (5) is impossible. Otherwise, because `Z=T_c o X`, the transformation `T_c` must preserve (6). For either `u²=2`, this requires

\[
 (2cu+c^2+1)^2-8(u+c)^2=0.
\]

Reducing modulo `u²−2` gives

\[
 (c^2-3)\bigl(4cu+c^2+5\bigr).                        \tag{7}
\]

Vanishing at both values `u=±sqrt2` forces `c²=3`, which has no rational solution. Therefore the disjoint-pole case cannot occur.

## 4. Degree one with coincident pole pairs

Now `h` preserves the unordered pair `±i`. Put

\[
 x_0=\frac{t^2-1}{t^2+1},\quad z_0=\frac{2t}{t^2+1},
 \quad x_0^2+z_0^2=1,
 \quad f=2x_0z_0,
 \quad g=x_0^2-z_0^2.
\]

A rational fractional linear transformation preserving `±i` induces a rational orthogonal transformation of the circle coordinates `(x_0,z_0)`. Consequently

\[
 f(h(t))=A f(t)+B g(t),\qquad A,B\in\mathbb Q,
 \qquad A^2+B^2=1.                                   \tag{8}
\]

For an elementary verification, write `h=(at+b)/(ct+d)`. Preservation of the two roots of `t²+1` means `(at+b)²+(ct+d)²=lambda(t²+1)`, so `a²+c²=b²+d²=lambda` and `ab+cd=0`. Substituting these relations into the circle-coordinate formulas gives an orthogonal matrix, and its induced action on the two quadratic harmonics gives (8).

The square in (1) has poles of order at most two only at `±i`, so `y` has poles of order at most one there and is finite at infinity. Every such rational function has the unique form

\[
 y=w+r x_0+s z_0,\qquad w,r,s\in\mathbb Q.             \tag{9}
\]

Comparing the first and second circle harmonics in (1) gives `wr=ws=0`. If `y` is nonconstant, then `w=0`, and comparison of the remaining coefficients yields

\[
 r^2+s^2=2,\qquad rs=\varepsilon+\eta A,\qquad
 (r^2-s^2)/2=\eta B.
\]

Hence

\[
 1=(rs)^2+((r^2-s^2)/2)^2
   =(\varepsilon+\eta A)^2+B^2
   =2+2\varepsilon\eta A.
\]

This forces `A=−epsilon eta/2` and `B²=3/4`, impossible for rational `B`. Thus `y` is constant. Equation (1), or harmonic comparison again, now gives `y²=1`, `A=−epsilon eta`, and `B=0`. This is precisely (2), proving Theorem 1.

For a constant rational `h`, write `c=eta f(h)`. The equation is `y²=k_epsilon(t)²+c`. If `c!=0`, both zeros `±sqrt(-c)` would have to be branch values of the degree-two map `k_epsilon`, or else their preimages would be simple zeros of a rational square. Its branch values are `±sqrt2`, forcing `c=−2`. But `|f(h)|<=1` for every real constant parameter, so this is impossible. Thus `c=0`, giving only `h in {0,1,−1,infinity}`.

## 5. Degree-two maps with an irreducible pole quartic are impossible

Suppose `h=P/Q` has degree two and `P²+Q²` is irreducible over `Q`. It has degree exactly four because the leading coefficients are real and not both zero. Its roots are simple, disjoint from `±i`, and form a single Galois orbit. The function `z=k_eta(h)` has a simple pole at every one: at `P=±iQ`, its numerator `P²+2eta PQ−Q²` is nonzero. The degree of `z` is exactly four, by multiplication of degrees under composition.

The hypotheses of Lemma 3 hold for `x=k_epsilon(t)` and this `z`. Equation (5) would make the degree-four function `Z` a fractional linear function of the degree-two function `X`, or a constant. Both are impossible. Therefore `P²+Q²` is reducible.

It has no rational root, since `P,Q` are coprime and a sum of their real squares cannot vanish at a rational point. A reducible quartic without rational roots is a product of two irreducible quadratic factors, with multiplicity allowed. If `a` is any root, then `h(a)=±i`, so `Q(i)` is a subfield of the quadratic field `Q(a)`. Thus `Q(a)=Q(i)`. Every root belongs to that field, proving Theorem 2, including repeated-pole cases.

## 6. What remains

The degree-two case now reduces to Gaussian-rational preimages of the contact poles. The residue argument alone does not exclude it: separate quadratic pole pairs can carry different residue signs, and shared or repeated poles require a separate local analysis. These were the remaining cases at the reduction stage. Sections 7–8 below finish the structural classification, and the shared/adjacent/opposite companion proofs now eliminate every case, establishing cubic minimality.

## 7. Repeated quadratic poles are also impossible

**Proposition 4.** A degree-two contact map cannot have `P²+Q²` proportional to the square of a quadratic.

**Proof.** Let `x=k_epsilon(t)` and `z=k_eta(h(t))` as above. Suppose first that the repeated quadratic `q` differs from `t²+1`. Its two roots are disjoint from `±i`. At both roots, `z` has poles of order two. The identity `y²=z²+x²−1` forces the complete polar part of `y` to equal `beta z`, with the same sign `beta` at both roots by Galois conjugation. More precisely,

\[
 y-\beta z=\frac{x^2-1}{y+\beta z}
\]

has a zero of order two at each root of `q`, since `x²−1=epsilon f(t)` is finite and nonzero there. But `y−beta z` can have poles only at `±i`, each at most simple, and is finite at infinity. A nonzero rational function with at most two poles counted with multiplicity cannot have four zeros. It cannot be identically zero, since `x²−1` is not identically zero. This is a contradiction.

If the repeated quadratic is `t²+1`, then `z` has double poles at `±i`. The same quotient identity now makes `y−beta z` finite at both points, so it is a rational constant `c`. It is nonzero, and

\[
 z=\frac{\varepsilon}{2\beta c}f(t)-\frac{c}{2\beta}.   \tag{10}
\]

The degree-two map `h` is totally ramified at its two preimages of `±i`, which are exactly `±i` in some order. These are all its ramification points. Therefore the finite branch values of `k_eta o h` are exactly `±sqrt2`; its remaining branch value is infinity. On the other hand, `f` has branch values `-1,+1,infinity`: differentiate to see four simple critical points over `±1`, together with its two double poles. An affine transform `a f+b` has finite branch values `b−a,b+a`. Comparing with (10) forces `b=0`, hence `c=0`, a contradiction. QED.

Thus the two quadratic pole pairs of any remaining degree-two candidate are distinct and unramified.

## 8. Exact divisor patterns for the remaining two cases

Write `r=y−z`, `s=y+z`. Then

\[
 rs=\varepsilon f(t),\qquad z=(s-r)/2.                \tag{11}
\]

The divisor of `f` consists of simple zeros at `0,1,-1,infinity` and double poles at `±i`. The distinct Gaussian quadratic pole pairs of `z` are denoted by `q_A,q_B`.

### Both pole pairs disjoint from the input pair

At either new pair, exactly one of `r,s` has simple poles and the other has simple zeros, because their product is finite and nonzero there. The choices must be opposite at the two pairs; otherwise the residue lemma gives `y=alpha x+beta z+c`, contradicting degrees as in Section 5. Both `r,s` have simple poles at `±i`. Hence each has total pole degree four, already has two zeros at one Gaussian pair, and must receive exactly two of the four real zeros of `f`.

Consequently, after labeling the two pairs, every candidate has the exact divisor form

\[
 r=c\,\frac{R_S(t)q_A(t)}{(t^2+1)q_B(t)},\qquad |S|=2,
                                                               \tag{12}
\]

where `S` is one of the six two-element subsets of `{0,1,-1,infinity}`. Here `R_S` denotes the product of the corresponding homogeneous linear factors, so a factor at infinity is retained in degree counting even when its affine representative is `1`. The other factor is then forced by (11). Equation (12) is a necessary reduction, not a sufficient contact criterion: `z=(s-r)/2` must still have the form `k_eta(h)`.

### Exactly one pole pair shared with the input

Now let the other pair be `q`, distinct from `t²+1`. Both `r,s` have simple poles at `±i`: a cancellation of the leading residue in either factor would contradict the nonzero leading square of `x` in `y²-z²=x²−1`. At the new pair, one factor has simple zeros and the other simple poles. The factor having the new zeros has total pole degree two, so it receives none of the four real zeros. Up to swapping `r,s`, necessarily

\[
 r=c\,\frac{q(t)}{t^2+1},\qquad
 s=\frac{4\varepsilon}{c}\frac{t(t^2-1)}{q(t)(t^2+1)}. \tag{13}
\]

Thus all four real zeros belong to the second factor. Further coefficient conditions in (12) and (13), coming from `z=k_eta(h)`, form the precise reduced algebraic problem. The shared-pole, adjacent-pole, and opposite-pole companion notes now exclude all of them; see the complete degree-minimality theorem for the exhaustive assembly.
