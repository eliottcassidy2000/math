#!/usr/bin/env python3
"""
complex_signing_filter_monad.py -- A COMPLEX (Hermitian) signing of a tournament
that FILTERS cycle lengths by residue mod r, and the question: does its
(real, Hermitian) spectrum carry NON-spectral information beyond charA?

Setup.  For a primitive r-th root of unity omega, define the Hermitian matrix
  M_omega[i][j] = omega   if i->j
                = conj(omega) if j->i
                = 0  if i==j.
M_omega is Hermitian (M*[i][j]=conj(M[j][i])=M[i][j]); its eigenvalues are REAL.
A directed cycle of length L (following arcs) gets Coates weight omega^L; its
reversal gets conj(omega)^L=omega^{-L}; together 2cos(2*pi*L/r).  So det(xI-M)
weights cycle structure by LENGTH MOD r:
  r=4 (omega=i):  cos(pi*L/2) = 0 for odd L -> KILLS ODD = the skew matrix S
                  (M_i = i*S). [SPECTRAL, shown in monad-S7.]
  r=3:  cos(2pi L/3): L=0 mod3 -> +1, else -1/2.
  r=6:  cos(pi L/3):  L=6->+1, L=3->-1, L=+-1->1/2, L=+-2->-1/2.  <- targets the
        LENGTH-6 non-spectrality onset.
  r=2 (omega=-1):  every cycle weight (-1)^L; M=-(J-I) on arcs both ways -> trivial.

For r in {3,4,6} we have omega+conj(omega) in {-1,0,1}, so the (real) char poly
of M_omega has INTEGER coefficients -> exact arithmetic in Z[omega].

THE QUESTION: is det(xI - M_omega) SPECTRAL (a function of det(xI - A))?  If NOT,
M_omega is a Hermitian-spectral invariant carrying non-spectral info -- a new tool,
possibly one that sees H / the length-6 structure.

Author: monad-explorer-2026-06-15-S7.  Pure-python exact Z[omega] arithmetic.
"""
from itertools import combinations

# ---------- Z[omega] arithmetic: element = (a,b) meaning a + b*omega ----------
# r=4: omega=i, omega^2=-1.    mult: (a+bi)(c+di)=(ac-bd)+(ad+bc)i
# r=3: omega^2=-1-omega.       (a+bw)(c+dw)=(ac-bd)+(ad+bc-bd)w
# r=6: omega^2=omega-1.        (a+bw)(c+dw)=(ac-bd)+(ad+bc+bd)w
def make_ring(r):
    if r==4:
        def mul(x,y):
            a,b=x;c,d=y;return (a*c-b*d, a*d+b*c)
        conj=lambda x:(x[0],-x[1])      # conj(a+bi)=a-bi
        omega=(0,1)
    elif r==3:
        def mul(x,y):
            a,b=x;c,d=y;return (a*c-b*d, a*d+b*c-b*d)
        conj=lambda x:(x[0]-x[1], -x[1])  # conj(a+bw)=a-b(w_bar); w_bar=-1-w => a+b(-1-w)=(a-b)+(-b)w
        omega=(0,1)
    elif r==6:
        def mul(x,y):
            a,b=x;c,d=y;return (a*c-b*d, a*d+b*c+b*d)
        conj=lambda x:(x[0]+x[1], -x[1])  # w_bar=1-w => conj(a+bw)=a+b(1-w)=(a+b)+(-b)w
        omega=(0,1)
    else:
        raise ValueError("r must be 3,4,6")
    add=lambda x,y:(x[0]+y[0],x[1]+y[1])
    sub=lambda x,y:(x[0]-y[0],x[1]-y[1])
    zero=(0,0); one=(1,0)
    return dict(mul=mul,add=add,sub=sub,conj=conj,omega=omega,zero=zero,one=one)

def charpoly_ring(M, R):
    """char poly det(xI-M) for M with Z[omega] entries, via Faddeev-LeVerrier.
    Returns list of REAL integer coeffs [c_0..c_n] (asserts omega-part vanishes &
    division exact). poly index = power of x."""
    n=len(M); mul=R['mul'];add=R['add'];sub=R['sub'];zero=R['zero'];one=R['one']
    def trace(X):
        t=zero
        for i in range(n): t=add(t,X[i][i])
        return t
    def matmul(X,Y):
        Z=[[zero]*n for _ in range(n)]
        for i in range(n):
            Xi=X[i]
            for k in range(n):
                xik=Xi[k]
                if xik==zero: continue
                Yk=Y[k]
                Zi=Z[i]
                for j in range(n):
                    if Yk[j]!=zero:
                        Zi[j]=add(Zi[j],mul(xik,Yk[j]))
        return Z
    def adddiag(X,s):
        return [[add(X[i][j],s) if i==j else X[i][j] for j in range(n)] for i in range(n)]
    def negscalar(s): return (-s[0],-s[1])
    Mp=[[M[i][j] for j in range(n)] for i in range(n)]
    c1=negscalar(trace(Mp)); cc=[one,c1]; cp=c1
    for k in range(2,n+1):
        Mk=matmul(M, adddiag(Mp, cp))
        num=negscalar(trace(Mk))
        # divide by k (exact, real integer): num should be (k*val, 0)
        assert num[1]==0, f"omega-part nonzero at k={k}: {num}"
        assert num[0]%k==0, f"non-integer FL at k={k}"
        ck=(num[0]//k,0); cc.append(ck); Mp=Mk; cp=ck
    # to real-int poly
    poly=[0]*(n+1)
    for k in range(n+1):
        assert cc[k][1]==0, f"non-real coeff at {k}: {cc[k]}"
        poly[n-k]=cc[k][0]
    return tuple(poly)

# ---- plain integer char poly of A ----
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

def M_omega(A, R):
    n=len(A);om=R['omega'];omb=R['conj'](om);z=R['zero']
    return [[ (om if A[i][j]==1 else (omb if A[j][i]==1 else z)) for j in range(n)] for i in range(n)]

def ham_path_count(A):
    n=len(A)
    if n==1:return 1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    full=(1<<n)-1
    for mask in range(1<<n):
        for v in range(n):
            d=dp[mask][v]
            if d==0:continue
            for w in range(n):
                if mask&(1<<w):continue
                if A[v][w]==1:dp[mask|(1<<w)][w]+=d
    return sum(dp[full][v] for v in range(n))

def run(n, rs=(3,4,6), sample_limit=None):
    print(f"\n========== n={n} ==========")
    total=1<<(n*(n-1)//2)
    if sample_limit is None or total<=sample_limit:
        it=range(total);exh=True
    else:
        import random;random.seed(5);it=[random.randrange(total) for _ in range(sample_limit)];exh=False
    rings={r:make_ring(r) for r in rs}
    # per r: charA -> set of charM ; and (charA,charM)->set(H)
    data={r:{'acp_to_mcp':{}, 'acpmcp_to_H':{}, 'acp_to_H':{}} for r in rs}
    cnt=0
    for bits in it:
        A=tour_from_bits(bits,n)
        acp=charpoly_int(A)
        H=ham_path_count(A)
        for r in rs:
            mcp=charpoly_ring(M_omega(A,rings[r]), rings[r])
            d=data[r]
            d['acp_to_mcp'].setdefault(acp,set()).add(mcp)
            d['acpmcp_to_H'].setdefault((acp,mcp),set()).add(H)
            d['acp_to_H'].setdefault(acp,set()).add(H)
        cnt+=1
    print(f"  processed {cnt} (exhaustive={exh})")
    for r in rs:
        d=data[r]
        n_acp=len(d['acp_to_mcp'])
        split=sum(1 for a,s in d['acp_to_mcp'].items() if len(s)>1)
        n_pairs=len(d['acpmcp_to_H'])
        H_amb_a=sum(1 for k,v in d['acp_to_H'].items() if len(v)>1)
        H_amb_am=sum(1 for k,v in d['acpmcp_to_H'].items() if len(v)>1)
        verdict = "SPECTRAL (no new info)" if split==0 else f"NON-SPECTRAL: splits {split} cospectral classes"
        print(f"  r={r}: charM_omega -> {verdict}")
        print(f"        #distinct charA={n_acp}, #distinct (charA,charM)={n_pairs};  "
              f"H ambiguous: charA={H_amb_a} -> (charA,charM)={H_amb_am}")

if __name__=="__main__":
    for n in [3,4,5,6]:
        run(n)
    run(7, sample_limit=30000)
