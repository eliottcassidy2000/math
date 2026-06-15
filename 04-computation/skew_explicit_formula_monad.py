#!/usr/bin/env python3
"""
skew_explicit_formula_monad.py -- Explicit spectral formulas for the walk
counts w_k = 1^T A^k 1 of a tournament, and supporting tests for the theorem
"det(xI - S) is spectral" (S = A - A^T).

Tests:
  (A) Is the SCORE SEQUENCE (sorted multiset) spectral?  (constant on cospectral
      classes?)  -- a sharper companion fact.
  (B) Are sum_i s_i^p (score power moments) spectral, for p=1,2,3,4?
  (C) Explicit fit: w_2 = (n-1)C(n,2) - sum s_i^2; and sum s_i^2 = 2C(n,3) -
      2*c3 + C(n,2) with c3 = tr(A^3)/3.  Verify EXACTLY.  Then fit w_3, w_4 as
      rational-coefficient polynomials in (tr A^3, tr A^4, tr A^5,...) and report.
  (D) w_k spectral for k beyond n (up to 2n)?  -> supports the "all k" conjecture.
  (E) Verify the matrix-determinant-lemma identity
        det(xI - S) = det((x-1)I - 2A) - bordered((x-1)I-2A)
      exactly for several tournaments (x evaluated at integer points).

Author: monad-explorer-2026-06-15-S7.  Pure-python exact arithmetic.
"""
from itertools import combinations
from fractions import Fraction

def charpoly_int(M):
    n=len(M)
    def trace(X): return sum(X[i][i] for i in range(n))
    def mm(X,Y): return [[sum(X[i][k]*Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    def ad(X,s): return [[X[i][j]+(s if i==j else 0) for j in range(n)] for i in range(n)]
    Mp=[row[:] for row in M]; cc=[1,-trace(Mp)]; cp=cc[1]
    for k in range(2,n+1):
        Mk=mm(M,ad(Mp,cp)); num=-trace(Mk); assert num%k==0
        ck=num//k; cc.append(ck); Mp=Mk; cp=ck
    poly=[0]*(n+1)
    for k in range(n+1): poly[n-k]=cc[k]
    return poly

def det_int(M):
    n=len(M)
    if n==0: return 1
    A=[row[:] for row in M];sign=1;prev=1
    for k in range(n-1):
        if A[k][k]==0:
            sw=None
            for r in range(k+1,n):
                if A[r][k]!=0: sw=r;break
            if sw is None: return 0
            A[k],A[sw]=A[sw],A[k];sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                A[i][j]=(A[i][j]*A[k][k]-A[i][k]*A[k][j])//prev
        prev=A[k][k]
    return sign*A[n-1][n-1]

def tour_from_bits(bits,n):
    A=[[0]*n for _ in range(n)];pos=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>pos)&1:A[i][j]=1
            else:A[j][i]=1
            pos+=1
    return A
def skew(A):
    n=len(A);return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]
def matvec(M,v): return [sum(M[i][j]*v[j] for j in range(len(v))) for i in range(len(M))]
def powtraces(A,kmax):
    n=len(A); P=[[1 if i==j else 0 for j in range(n)] for i in range(n)]; tr=[n]
    cur=P
    for k in range(1,kmax+1):
        cur=[[sum(cur[i][t]*A[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        tr.append(sum(cur[i][i] for i in range(n)))
    return tr  # tr[k]=tr(A^k)
def walk_counts(A,kmax):
    n=len(A); v=[1]*n; w=[n]
    for k in range(1,kmax+1):
        v=matvec(A,v); w.append(sum(v))
    return w
def scores(A):
    return sorted(sum(row) for row in A)

def run(n, sample_limit=None, kmax_mult=2):
    print(f"\n===== n={n} =====")
    total=1<<(n*(n-1)//2)
    if sample_limit is None or total<=sample_limit:
        it=range(total); exh=True
    else:
        import random;random.seed(11);it=[random.randrange(total) for _ in range(sample_limit)];exh=False
    kmax=kmax_mult*n
    acp_score={}; acp_smoment={}; acp_wk={}
    # for explicit fit collect (tr3, w3) etc per distinct charA
    rows=[]
    for bits in it:
        A=tour_from_bits(bits,n)
        acp=tuple(charpoly_int(A))
        sc=tuple(scores(A))
        acp_score.setdefault(acp,set()).add(sc)
        sm=tuple(sum(s**p for s in sc) for p in range(1,5))
        acp_smoment.setdefault(acp,set()).add(sm)
        w=tuple(walk_counts(A,kmax)[1:])
        acp_wk.setdefault(acp,set()).add(w)
        tr=powtraces(A,kmax)
        rows.append((tr,w,sc))
    sscore_split=sum(1 for a,s in acp_score.items() if len(s)>1)
    print(f"  exhaustive={exh}; #cospectral classes={len(acp_score)}")
    print(f"  (A) score sequence SPECTRAL? split in {sscore_split} classes -> {'SPECTRAL' if not sscore_split else 'NOT spectral'}")
    for p in range(1,5):
        bad=sum(1 for a,s in acp_smoment.items() if len({t[p-1] for t in s})>1)
        print(f"  (B) sum s_i^{p} : {'spectral' if bad==0 else f'NOT spectral ({bad} classes)'}")
    for k in range(1,kmax+1):
        bad=sum(1 for a,s in acp_wk.items() if len({t[k-1] for t in s})>1)
        tag = '' if bad==0 else f'  <-- NON-spectral ({bad})'
        if k in (1,2,3,n,n+1,2*n) or bad>0:
            print(f"  (D) w_{k}: {'spectral' if bad==0 else 'NON-SPECTRAL'}{tag}")
    # (C) explicit formulas via exact rational least-data fit using the cospectral structure:
    #   w_2 closed form check
    import math
    C2=n*(n-1)//2; C3=n*(n-1)*(n-2)//6
    ok2=all( r[1][1] == (n-1)*C2 - (2*C3 - 2*(r[0][3]//3) + C2) for r in rows if r[0][3]%3==0)
    print(f"  (C) w_2 == (n-1)C(n,2) - [2C(n,3) - 2*tr(A^3)/3 + C(n,2)] : {'EXACT all' if ok2 else 'FAILS'}")
    # fit w_3 = a*tr3 + b   (tr1=tr2=0 for tournaments; only tr3 nontrivial up to deg3)
    # solve using two rows with different tr3
    fit3=fit_linear(rows, ycol=2, xcols=[3], n=n)
    print(f"  (C) w_3 linear in tr(A^3): {fit3}")
    fit4=fit_linear(rows, ycol=3, xcols=[3,4], n=n)
    print(f"  (C) w_4 linear in (tr A^3, tr A^4): {fit4}")
    fit5=fit_linear(rows, ycol=4, xcols=[3,4,5], n=n) if kmax>=5 else "n/a"
    print(f"  (C) w_5 linear in (tr A^3,A^4,A^5): {fit5}")

def fit_linear(rows, ycol, xcols, n):
    """Try to fit w_{ycol+1} = sum coeff_j * tr(A^{xcols[j]}) + const exactly over
    rationals, using exact solve on a subset; then VERIFY on all rows. Returns
    the formula string or 'no exact affine fit'."""
    # Build design: features = [1] + [tr^{c} for c in xcols]; target = w_{ycol+1}=r[1][ycol]
    feats=[]; targs=[]
    seen=set()
    for r in rows:
        f=tuple([1]+[r[0][c] for c in xcols])
        if f in seen:
            continue
        seen.add(f)
        feats.append([Fraction(x) for x in f]); targs.append(Fraction(r[1][ycol]))
        if len(feats)>len(xcols)+1+3: break
    # solve least-squares-exact: pick first (len+1) independent rows
    m=len(xcols)+1
    # Gaussian elimination to find a solution to feats @ coeff = targs (use first m indep rows)
    A=[]; b=[]
    for f,t in zip(feats,targs):
        A.append(f[:]); b.append(t)
        if len(A)>m+2: break
    coeff=solve_exact(A,b,m)
    if coeff is None:
        return "no exact affine fit (rank-deficient/inconsistent)"
    # verify on ALL rows
    for r_idx,r in enumerate(__rows_cache(rows)):
        pass
    # verify over the deduped feats fully
    good=True
    for r in rows:
        pred=coeff[0]+sum(coeff[j+1]*r[0][xcols[j]] for j in range(len(xcols)))
        if pred!=r[1][ycol]:
            good=False;break
    if not good:
        return "no exact affine fit"
    terms=[f"{coeff[0]}"]+[f"{coeff[j+1]}*tr(A^{xcols[j]})" for j in range(len(xcols))]
    return " + ".join(terms)+"  [EXACT on all rows]"

def __rows_cache(rows): return []

def solve_exact(A,b,m):
    """Solve for m unknowns using up to len(A) equations; return solution if a
    consistent unique solution exists on the first m independent rows, else None."""
    rows=[A[i][:]+[b[i]] for i in range(len(A))]
    ncol=m
    piv=0; pivcols=[]
    r=0
    for c in range(ncol):
        # find pivot
        sel=None
        for i in range(r,len(rows)):
            if rows[i][c]!=0: sel=i;break
        if sel is None: continue
        rows[r],rows[sel]=rows[sel],rows[r]
        pv=rows[r][c]
        rows[r]=[x/pv for x in rows[r]]
        for i in range(len(rows)):
            if i!=r and rows[i][c]!=0:
                fct=rows[i][c]
                rows[i]=[rows[i][j]-fct*rows[r][j] for j in range(ncol+1)]
        pivcols.append(c); r+=1
        if r==ncol: break
    if r<ncol: return None
    sol=[Fraction(0)]*ncol
    for idx,c in enumerate(pivcols):
        sol[c]=rows[idx][ncol]
    # check consistency of remaining rows
    for i in range(r,len(rows)):
        if any(rows[i][j]!=0 for j in range(ncol)) :
            continue
        if rows[i][ncol]!=0: return None
    return sol

def verify_mdl_identity(n, ntest=20):
    """det(xI - S) = det((x-1)I - 2A) - bordered((x-1)I - 2A), at integer x."""
    import random; random.seed(3)
    total=1<<(n*(n-1)//2)
    fails=0
    for _ in range(ntest):
        A=tour_from_bits(random.randrange(total),n)
        S=skew(A)
        for x in range(-3,5):
            lhs=det_int([[ (x if i==j else 0)-S[i][j] for j in range(n)] for i in range(n)])
            M=[[ ((x-1) if i==j else 0)-2*A[i][j] for j in range(n)] for i in range(n)]
            detM=det_int(M)
            # bordered det([[M,1],[1^T,0]])
            B=[row[:]+[1] for row in M]+[[1]*n+[0]]
            detB=det_int(B)
            rhs=detM - detB
            if lhs!=rhs: fails+=1
    print(f"  (E) n={n}: MDL identity det(xI-S)=det((x-1)I-2A)-bordered : {'HOLDS' if fails==0 else f'{fails} FAILS'}")

if __name__=="__main__":
    for n in [3,4,5,6]:
        run(n)
    run(7, sample_limit=50000)
    print("\n--- (E) matrix-determinant-lemma identity check ---")
    for n in [4,5,6,7]:
        verify_mdl_identity(n)
