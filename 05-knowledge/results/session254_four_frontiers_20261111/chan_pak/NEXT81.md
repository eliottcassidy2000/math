# Mod81 continuation: aliased derivative gauge and one nine-valued sidecar

**PROVED HERE.** This supplements PROOF.md without changing its mod9/mod27
claims. There is no large census and no integer-coverage conclusion.

## 1. Euler-log interpretation, with its precision boundary

Write theta=t d/dt and H=Q(1), D=Q'(1), E=Q^[2](1), F=Q^[3](1).
Since H is a3-unit, log derivatives at1 are rational with denominators
invertible3-adically. Precisely,

    (theta² log Q)(1)=2E/H+D/H-D²/H²=(2HE-D²+HD)/H².

Reducing this identity modulo3 gives the previous K/H²; the numerator
above is the integer lift, rather than its residue-class notation.

For every j>=2, (theta^j log(t^r Q))(1)=(theta^j log Q)(1).
Indeed log(t^r Q)=r log t+log Q locally at1, theta log t=1, and the
next Euler derivative vanishes. Thus padding is exactly a logarithmic
first-derivative gauge, with higher Euler-log derivatives invariant.
This is a genuine differential-operation identity. It is not by itself an
identification with any JC operator or a guarantee that these rational
invariants are independent coordinates of a finite modular function.

That last distinction matters already at81.

## 2. Cubic alias and the complete stable interface

Assume at least4 universal sinks have been padded, so Q(3u)=0 mod81.
At a unit chart,

    Q(1+3u)=H+3uD+9u²E+27u³F (mod81).

Since u³=u mod3, the last term aliases into the linear term. Define

    A = D+9F (mod27).

Then

    Q(1+3u)=H+3uA+9u²E (mod81).                (1)

Together with parity, this handles the -1 chart. Thus the complete stable
mod81 evaluation function, and hence its substitution jet, is determined
by H mod81, A mod27, E mod9, and order parity. Conversely the evaluations
at1,4,-2 recover these: the difference at4,-2 is6A mod81, and their sum
minus2H is18E mod81. Dividing by3 and9 respectively leaves units2 modulo27
and9. For a3-unit H, evaluation at-1 recovers parity.

The factorial jet now ends at8: v3(8!)=2 and v3(9!)=4, so all F_S(k) with
k>=9 vanish modulo81. This uses the same function/finite-difference bridge
as PROOF.md; no assertion that formal Taylor tuples are all realizable.

## 3. A lawful normalization with bounded vertex cost

Under r sink paddings, the Hasse coordinates transform as

    D_r=D+rH,
    E_r=E+rD+binom(r,2)H,
    F_r=F+rE+binom(r,2)D+binom(r,3)H.

Put A_r=D_r+9F_r modulo27. These identities give

    A_(r+3)=A_r+12H (mod27),
    A_(r+27)=A_r (mod27),
    A_r=D+rH (mod3).

For the first equality, F_(r+3)-F_r=H mod3, because
binom(r+3,3)-binom(r,3)=3binom(r,2)+3r+1. For each fixed r mod3,
12H generates the nine possibilities with that residue modulo3, since H
is a unit. The three r residues mod3 have three different A residues.
Consequently r -> A_r is a BIJECTION of Z/27, although it is not generally
a linear map modulo27.

There is a unique r in {4,...,57} with A_r=0 mod27 and any prescribed
final parity. This follows from the preceding bijection and coprimality
of27 and2. The operation retains every old vertex and arc and exact H.

**Mod81 normal-form theorem.** Every3-unit source can be normalized with
at most57 sink vertices so that its complete mod81 substitution jet is
represented by

    Q_r(3u)=0,
    Q_r(s(1+3u))=s^epsilon [H+9u² J] (mod81),
    J=E_r mod9, s in {1,-1}.

The normalized jet is determined exactly by H mod81, chosen parity, and
the single nine-valued coordinate J. Given fixed H and parity, different J
are distinguished by the evaluation at4. We do NOT claim all nine values
occur among tournaments of each fixed H.

J is invariant under further padding FOLLOWED BY RENORMALIZATION: replacing
the source by t^sQ just shifts the unique solution r modulo27. E_(r+27)=E_r
mod9, so the resulting J is unchanged. This is a complete stable padding
orbit coordinate at81, not merely a conjectured higher-derivative feature.
Modulo3 it recovers the previous curvature: D_r=0 mod3, hence
K=2H J mod3.

The minimal normalized r>=4 with chosen parity lies in4..57 and minimizes
the exact two-ear output H+2^r Q(2) within this stable padding progression.
The output bound is H+2^57 Q(2), final order at most |S|+59. Large integer
cost remains explicit even though H is preserved during normalization.

## 4. A sharp pure-parity derivative ambiguity

Consider the signed integer polynomial

    P(t)=t^4[(t²-1)^3-9(t²-1)].

It is even and vanishes as an integer function modulo81. If3 divides t,
t^4 supplies81; otherwise t²-1=3w and the bracket is27(w³-w), a multiple
of81. At1 its Hasse data are

    P(1)=0, P'(1)=-18, P^[2](1)=-81, P^[3](1)=-136.

Thus D changes by9 mod27 and F changes by2 mod3, while A=D+9F and E mod9
stay unchanged. One cannot recover raw D mod27 or F mod3 separately from
the mod81 evaluation function, even with pure parity and vanishing at
multiples of3. This is a hostile in the ambient integer-polynomial module;
P is signed and is NOT claimed to be a tournament polynomial or a positive
fibre operation. The preceding normal-form theorem does apply to actual
tournaments because it only uses lawful padding and their inherited parity.

This is the precise caution for a possible differential transfer: higher
logarithmic derivatives are genuine gauge invariants, but modular integer
functions can identify derivative coordinates through integer-valued
polynomial aliases. The correct observer at81 is A, not bare D.

## 5. Reproduction

`python3 probe81.py` produces output81.json; optimized mode is identical.
It checks all residue evaluations and all54 normalization choices on five
explicit positive tournament polynomials, the padding bijection, the
mod27 curvature reduction, and the signed alias on all81 residue classes.
The all-source result is proved above, not extrapolated from these controls.
