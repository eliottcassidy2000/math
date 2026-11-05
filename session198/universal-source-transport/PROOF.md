# Positive kernel witnesses forbid universal source-word transport

Session 198. **PROVED, architecture-specific obstruction.** The statements
below concern identities valid for every pair of actual sources. They do not
exclude target-dependent rearrangements, equalities on a restricted physical
fiber, or universal scalar coverage.

## 1. The physical carrier

Use the complete finite source grammar from the Session 198 bilinear-descent
packet. Every nontrivial start-one source has nonnegative coordinates

    y=(a,f,e,p)=(D-E-d,U-D+E-d,E,d-E),
    D=a+2e+p, U=a+f+2e+2p, d=e+p>=1.

The seed is t=(0,0,0,1). The actual one-vertex operations are

    O(a,f,e,p)=(a+e,a+2f+e+2p,0,e+p),
    Z(a,f,e,p)=(2a+f+e+p,f,e+p,0),
    J(a,f,e,p)=(f+p,0,a+e+p,e).

O adds a sink. Z adds a vertex beating only the old root, retaining that
root. J performs the same birth and marks the newborn. Every finite word
maps every actual source to another actual source. Letters are applied in
chronological order, so a word whose first letter is T has matrix W=...T.

The complete C2 forms are N,M. The six C3 forms A,B,C,T,J1,J2, together with
these two C2 forms, will be denoted by the set F. Their 4-by-4 matrices P_F
are those in the positive-source packet. The separate J0,J3 forms have all
entries divisible by three and are deliberately excluded below.

## 2. Three lost directions are positively visible

The source moves have the respective exact kernel vectors

    k_O=(1,-1,-1,1),
    k_Z=(0,0,1,-1),
    k_J=(1,1,0,-1).                                (1)

These are algebraic directions, not asserted physical sources. Direct
multiplication of the reviewed matrices gives:

    F    P_F k_O       P_F k_Z       P_F k_J
    N    (1,0,2,2)     (1,1,1,1)     (1,1,3,2)
    M    (1,1,3,3)     (2,1,3,2)     (1,1,3,2)
    A    (3,2,3,3)     (1,1,2,2)     (4,2,6,6)
    B    (1,1,2,2)     (3,2,3,3)     (1,1,3,2)
    C    (1,1,3,3)     (4,3,7,6)     (1,1,5,4)
    T    (2,1,5,4)     (1,1,1,1)     (6,4,12,10)
    J1   (0,0,1,1)     (0,0,2,2)     (3,6,3,4)
    J2   (0,0,2,2)     (0,0,1,1)     (6,9,8,10).

Every displayed column is entrywise at least (0,0,1,1). Consequently, for
every actual nontrivial source y and every F in F,

    y^t P_F k_T >= e+p = d >= 1.                   (2)

This is a strict physical positivity statement. It does not treat the
nonnegative cone as a source realization criterion.

## 3. No universal transfer, at arbitrary word length

**Theorem 1.** Let F belong to F. For every source word X, including the
empty word, every nonempty source word Y, and every real 4-by-4 matrix Q,
there is no matrix identity

    X^t P_F = Q Y.                                (3)

In particular no universal scalar identity

    F(Xv,w)=v^t Q(Yw)                              (4)

holds for all actual nontrivial source pairs v,w.

Proof of (3). Let T be the first letter of Y. By (1), Yk_T=0. For any
actual source v, Xv is still an actual source, so (2) gives

    v^t X^t P_F k_T >= 1.

Equation (3) would make this zero, a contradiction. There is no bound on the
lengths of X or Y in this argument.

To justify the scalar consequence without assuming that arbitrary cone
points are physical, actual source vectors span R^4. Four explicit source
recipes are

    t=(0,0,0,1),
    Zt=(1,0,1,0),
    Ot=(0,2,0,1),
    J(Zt)=(0,0,2,1).

Their column determinant has absolute value four. If (4) held on all
physical source pairs, linear extension in each variable would force (3).
Thus the matrix obstruction is genuinely a physical universal-identity
obstruction. The transpose versions follow because F is closed under
transpose: B and T exchange, as do J1 and J2; the other four are symmetric.

This theorem allows an arbitrary replacement bilinear matrix Q, not merely
one of the ten inherited scalar families. It rules out transporting a
whole nonempty source word from one independent source to the other by a
fixed universal formula. It does not forbid swapping the two entire
sources, and does not forbid equality for selected pairs.

## 4. No fixed universal affine radix rule

**Theorem 2.** Let F belong to F, let P be any real 4-by-4 matrix, and let
X,Y be source words with at least one nonempty. There are no real constants
lambda!=0 and c such that

    (Xv)^t P(Yw)=lambda F(v,w)+c                    (5)

for every actual nontrivial source pair v,w.

First c=0. The triangle vector q=(1,0,1,0) is physical, and so is 3q:
apply the ordinary L=2 gadget to the triangle, giving old coordinates
(D,E,U,d)=(9,3,9,3). Fix any physical v and compare (5) at w=q and w=3q.
Both bilinear terms triple, forcing c=3c and hence c=0.

By the physical spanning argument, (5) then forces

    X^t P Y=lambda P_F.                            (6)

If Y is nonempty, apply both sides to its first-letter kernel k_T. The
left vanishes and P_F k_T is nonzero by the displayed table. This
contradicts lambda!=0. If only X is nonempty, transpose (6) and use the
transpose closure of F. This proves the theorem.

Thus no fixed pair of source extensions can implement an identity such as
2F+1, 2F-1, 3F+2 or 3F+4 simultaneously on every source pair in these
families. This leaves open target-dependent choices and identities whose
hypotheses restrict the sources to smaller physical fibers. It is not a
proof that a finite collection of restricted cases cannot cover all targets.

## 5. The exceptional families and the useful surviving question

The omitted forms expose the boundary sharply:

    P_J0 k_O=0,       P_J3 k_Z=0.

They really have universal word-transfer identities, for example

    Z^t P_J0=P_J0 O,
    J^t P_J0=P_J3 J,
    O^t P_J3=P_J3 Z.

These follow by direct integer multiplication and concern actual operations.
Their outputs are always divisible by three. The failure of Theorem 1's
strict positivity is therefore meaningful, rather than a missing test.

For the remaining forms, a constructive exact-value proof must retain
additional source information in a transfer, or choose the transfer only
on a specified physical fiber. The separately derived rank-two J-swap
defect in the reassociation packet is such a fiber-level statement; this
all-word obstruction does not replace its exact equality classification.

The accompanying checker verifies every kernel, every displayed vector,
the physical spanning recipes by literal small tournament counts, the
physical q and3q witnesses, and all three surviving identities. Its finite
matrix certificate supports the arbitrary-word proof, not extrapolation
from a word-length cutoff. No literature-priority claim is made.
