# An exact eight-step reversal in the even-order reader plane

Use the notation and triangle-wrap operator from README.md. For an
even-order core the mod-7 state is q=(Q(1),Q(2),Q(3)) and the wrap matrix is

    A = [[1,1,0],[-2,4,2],[0,-3,5]].

Its characteristic polynomial over F_7 is (x+1)(x^2+3x-1). The quadratic
has nonsquare discriminant 6, so the state splits as a real algebraic
direct sum over F_7: a one-dimensional sign line and a two-dimensional
plane, identifiable with F_49. On the line A acts by -1. In the plane
it multiplies by a root lambda of that quadratic. The other root is
lambda^7 (Frobenius); the product is -1. Therefore lambda^8=-1 and
lambda has exact order 16. In particular A^8 preserves the line,
negates the plane, and A^16 is identity. Any state with a nonzero
plane projection has exact state period 16, and its observed counts
have the same period by the first-three-observations triangularity.

The line projector is P=2(A^2+3A-I) modulo 7. Reading its first
coordinate gives L=2Q(1)+2Q(2)+4Q(3) mod 7. Take h_m=H(U^m(T)).
Then

    h_{m+8} = 2(-1)^m L - h_m (mod 7).

Thus eight steps reflect the output about an alternating centre, and
another eight return it. For the even-order 343-path source of README,
q=(0,0,5), hence L=6, h_0=0, and h_8=5 modulo 7. The sign line and
plane are genuine finite-field analogues of an alternating baseline
and paired phase coordinates. The identity is exact but modular;
it does not supply inexpensive absolute placement. Its field is
F_49, not F_8, and the resemblance to an eight-step Bott shift is an
analogy, not a link to topological K-theory. The F_8 bridge in README
concerns the separate 1,2,4 indexing of nonzero signs.

A short direct matrix check (standard Python, no dependencies):

    A=((1,1,0),(-2,4,2),(0,-3,5))
    I=((1,0,0),(0,1,0),(0,0,1))
    P=((2,2,4),(3,3,6),(5,5,3))
    def dot(M,N):
        return tuple(tuple(sum(M[i][k]*N[k][j] for k in range(3))%7
                           for j in range(3)) for i in range(3))
    B=I
    for _ in range(8): B=dot(B,A)
    assert dot(P,P)==P
    assert dot(B,B)==I and B!=I
    assert B==tuple(tuple((2*P[i][j]-I[i][j])%7 for j in range(3))
                    for i in range(3))
