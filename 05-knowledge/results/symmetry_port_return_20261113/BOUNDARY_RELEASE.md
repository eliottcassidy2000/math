# Preserve the useful port while releasing the old module boundary

**PROVED:** the odd-wrapper context divisor and the explicit release
below. **FINITE-EXACT:** the three displayed profile controls. This note
explains why the new positive constructions are allowed to escape older
uniform-context obstructions. It does not infer universal target coverage.

## 1. An asymmetric odd wrapper can have an unavoidable context divisor

For a nonempty tournament T of odd order, write

    W(T)=C3[T,x,y],     T→x→y→T.

The inherited slot polynomial has integer monomial coefficients and
Q_T(−t)=−Q_T(t); this follows directly from its odd-cycle expansion.
Consequently P_T(2)≡−H(T) modulo3. Also P_T(3)≡0 modulo3, either by
cycling its three labelled slots or by Q_T(0)=0. The wrapper law gives

    H(W(T))=P_T(2)+H(T)≡0 (mod3),
    P_W(T)(2)=2P_T(3)+4P_T(2)−2H(T)≡0 (mod3).

Let c_r count covers in r nonempty ordered slots. Then c_1=H and
c_2=P2−2H are divisible by3. For r≥3, the free permutation of the r
nonempty path slots proves r! divides c_r, hence3 divides c_r as well.
Thus **every** complete nonempty-run coefficient of W(T) is divisible
by3. Every Hamilton-count observation in a uniform exterior is an integer
linear combination of these coefficients and is consequently divisible
by3. More generally the same argument applies to a fixed number of global
path slots with the module kept uniform.

This need not come from a global automorphism. The rigid A of the current
clock has(H,P2,P3)=(95,1276,8709), and W(A) is rigid by the canonical
cyclic-module argument. Its first two slot counts are

    H(W(A))=1371, P_W(A)(2)=22332, gcd=3.

Together with the all-r argument this proves its full ordered-cover
content is exactly3. Global rigidity does not remove the uniform-context
obstruction.

## 2. The new refill releases exactly the hypothesis causing that divisor

After one zero-wrap cycle, the right module is W(A→A). Its values are

    H=1637201, P2=158188016, gcd=1.

It therefore has full ordered-cover content1. The old W(A) remains
induced, with every old arc untouched. It is **no longer a uniform module**
in the enlarged graph: a vertex in the fresh A loses to the old A and to
y, but beats x. This selective incidence crosses the old wrapper boundary.

The complete outer counts at the two endpoints are

    H(G(TT1,A))=4021131,
    H(G(TT2,A→A))=85423165841,

which are coprime. The first is divisible by3 and the second is not.
The construction does not contradict a uniform-context divisor theorem:
its nonneutral installation step deliberately changes that context type.
Within the fresh B replacement, the smaller B leaf is still uniform and
unmarked, so the neutral identity and protected-port preservation remain
valid. The useful port and the obstruction-producing module are different
physical subsets.

The determinant extensions have the same concrete distinction. Every
successful new mask is nonconstant on its old eight-vertex module, so
those old vertices remain induced but cease to be a uniform module.
The old balanced root and its old star survive. This is exactly the kind
of boundary change left open by the
[all-frame retained-module obstruction](../prime_terminal_20261113_1610/absorption/TWO_RUN_CONE.md),
not a claim that an arbitrary enlarged exterior repairs a retained module.

## 3. A two-coordinate divisor can still overstate the full obstruction

Higher observations matter even in this explicit family. At T=A^191
(ordinal product), exact arithmetic gives

    gcd(H(W(T)),P_W(T)(2))=1371=3·457,
    gcd(H(W(T)),P_W(T)(2),P_W(T)(3))=3,
    P_W(T)(3)≡312 (mod457).

The all-r parity argument still forces3. The third observation excludes
the additional457, so the full ordered-cover content is exactly3. Thus
the gcd of the first two slot counts is not always the full context
divisor. The index191 is an explicit witness, not an asserted universal
threshold or an independently significant constant.

The formula used here is just

    P_W(A^n)(3)=3Q_A(4)^n+9Q_A(3)^n−3Q_A(2)^n,
    Q_A(4)=40760.

`boundary_check.py` reproduces the three exact controls. The all-order
divisor claim is proved in§1, not extrapolated from them.
