#!/usr/bin/env python3
"""
signed_pencil_spectral_monad.py -- THE UNIFICATION (monad-S7):

For a tournament, A^T = J - I - A, so EVERY "signed adjacency" matrix
  M_omega = omega*A + conj(omega)*A^T = (omega - conj(omega))*A + conj(omega)*(J - I)
is AFFINE in A modulo the rank-1 all-ones matrix J.  More generally the whole
3-parameter pencil
  P(alpha,beta,gamma) = alpha*A + beta*(J - I) + gamma*I
contains: A itself (1,0,0); the skew S=A-A^T=2A-(J-I) (2,-1,0); every M_omega.

CLAIM (this script): the characteristic polynomial of P(alpha,beta,gamma) is
SPECTRAL (a function of charA) for tournaments, for ALL scalars alpha!=0,beta,gamma.

PROOF (clean, unconditional reduction): with y=x+beta-gamma,
  det(xI - P) = det(yI - alpha*A - beta*J)
             = det(yI - alpha*A) * (1 - beta * 1^T (yI - alpha A)^{-1} 1)   [MDL]
det(yI-alpha A)=alpha^n charA(y/alpha) is spectral; the bracket is the WALK-
GENERATING FUNCTION 1^T(yI-alpha A)^{-1}1 = sum_k alpha^k w_k y^{-k-1}, spectral
because the walk counts w_k=1^T A^k 1 are spectral (monad-S7 / Moon / known
complement=converse equivalence).  So det(xI-P) is spectral.  QED (modulo w_k).

Tests:
  (1) verify M_omega = (omega-omega_bar)A + omega_bar(J-I) exactly (omega=i case:
      M_i = i*S; check S = 2A-(J-I)).
  (2) verify the pencil char poly is spectral for many integer (alpha,beta,gamma):
      group all tournaments by charA, check char(P) constant on each class.
Author: monad-explorer-2026-06-15-S7.  Pure-python exact integer arithmetic.
"""
def charpoly_int(M):
    n=len(M)
    def trace(X):return sum(X[i][i] for i in range(n))
    def mm(X,Y):return [[sum(X[i][k]*Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    def ad(X,s):return [[X[i][j]+(s if i==j else 0) for j in range(n)] for i in range(n)]
    Mp=[row[:] for row in M];cc=[1,-trace(Mp)];cp=cc[1]
    for k in range(2,n+1):
        Mk=mm(M,ad(Mp,cp));num=-trace(Mk);assert num%k==0
        ck=num//k;cc.append(ck);Mp=Mk;cp=ck
    poly=[0]*(n+1)
    for k in range(n+1):poly[n-k]=cc[k]
    return tuple(poly)
def tour_from_bits(bits,n):
    A=[[0]*n for _ in range(n)];pos=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>pos)&1:A[i][j]=1
            else:A[j][i]=1
            pos+=1
    return A
def pencil(A, alpha, beta, gamma):
    n=len(A)
    return [[ alpha*A[i][j] + (beta if i!=j else 0) + (gamma if i==j else 0)
              for j in range(n)] for i in range(n)]

# (1) M_omega = (om-omb)A + omb(J-I) ; omega=i: M_i = i*S, S=2A-(J-I)
def check_affine_identity(n=5):
    import random;random.seed(1)
    fails=0
    for _ in range(200):
        A=tour_from_bits(random.randrange(1<<(n*(n-1)//2)),n)
        S=[[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]
        # S should equal 2A-(J-I): 2A[i][j] - (1 if i!=j else 0)
        S2=pencil(A,2,-1,0)
        if S!=S2: fails+=1
    print(f"  (1) S = A-A^T == 2A-(J-I) [pencil(2,-1,0)] : {'HOLDS' if fails==0 else f'{fails} FAIL'} (n={n})")

def run(n, params, sample_limit=None):
    total=1<<(n*(n-1)//2)
    if sample_limit is None or total<=sample_limit:
        it=range(total);exh=True
    else:
        import random;random.seed(9);it=[random.randrange(total) for _ in range(sample_limit)];exh=False
    # acp -> {(alpha,beta,gamma): set of char(P)}
    acp_pencil={p:{} for p in params}
    for bits in it:
        A=tour_from_bits(bits,n)
        acp=charpoly_int(A)
        for p in params:
            cp=charpoly_int(pencil(A,*p))
            acp_pencil[p].setdefault(acp,set()).add(cp)
    print(f"  n={n} (exhaustive={exh}, {min(total,sample_limit or total)} tournaments):")
    for p in params:
        split=sum(1 for a,s in acp_pencil[p].items() if len(s)>1)
        print(f"     P{p}: {'SPECTRAL' if split==0 else f'NON-SPECTRAL (splits {split})'}")

if __name__=="__main__":
    check_affine_identity(5)
    # a spread of pencils: A(1,0,0), skew(2,-1,0), (3,1,0), (1,2,-1), (5,-3,2), (2,7,0)
    params=[(1,0,0),(2,-1,0),(3,1,0),(1,2,-1),(5,-3,2),(2,7,0)]
    print("  (2) char poly of pencil P(alpha,beta,gamma)=alpha*A+beta*(J-I)+gamma*I spectral?")
    for n in [4,5,6]:
        run(n, params)
    run(7, params, sample_limit=20000)
