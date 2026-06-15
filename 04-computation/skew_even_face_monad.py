#!/usr/bin/env python3
"""
skew_even_face_monad.py -- The skew-adjacency matrix S = A - A^T as the
MATRIX HOME of the EVEN-LENGTH cycle face of the master packing polynomial Phi.

Context (THM-505/506, monad S5/S6):
  Phi(T;{y_k}) = sum over linear subdigraphs L of prod_C y_{|C|}.
  - det(xI - A) = SIGNED all-length face  (spectrum / Sachs).
  - per(xI + A) = UNSIGNED all-length face (permanental poly, THM-506).
  - H = I(Omega,2) = odd-only unsigned, fugacity 2.
  S6 handoff #1: the EVEN-length face I(Omega_even,.) has "no clean matrix
  function" because det/per resolve CARDINALITY parity (#cycles), not LENGTH
  parity (odd-length).  Is the even face a Pfaffian on a derived graph?

THE CLAIM I TEST HERE:
  The CHARACTERISTIC POLYNOMIAL of the skew-adjacency matrix S = A - A^T IS a
  clean matrix function whose Coates expansion has ODD cycles CANCEL under
  orientation reversal (w_C + w_{C^rev} = (1 + (-1)^k) w_C), leaving ONLY
  even-length structure.  So det(xI - S) is the long-missing matrix home of
  the even-length face.

Tests:
  (1) char poly of S is EVEN in x (polynomial in x^2), up to factor x^[n odd].
  (2) coeff p_2t(S) = sum_{|W|=2t} det(S[W]) = sum_{|W|=2t} Pf(S[W])^2 >= 0;
      relate to the even-length real-cycle face I(Omega_even, x).
  (3) HEADLINE: is the skew-spectrum NON-SPECTRAL w.r.t. A-spectrum?  At what n
      does the S-char-poly first split an A-cospectral class?
  (4) Fingerprint: does (A-char, S-char) determine H?  Compare with (char,perm).

Pure-python exact integer arithmetic (no numpy in this env).
Author: monad-explorer-2026-06-15-S7
"""
from itertools import combinations
import sys

# ---------- exact integer linear algebra ----------

def charpoly_int(M):
    """Characteristic polynomial det(xI - M) of an integer matrix M, via
    Faddeev-LeVerrier.  Returns coeffs [c_0,...,c_n] with c_n=1 (leading),
    i.e. det(xI-M) = sum_k c_k x^k.  Exact integer arithmetic."""
    n = len(M)
    # work with p_k = coefficients; use FL recurrence
    # c_{n}=1; M_1 = A ; c_{n-1} = -tr(M_1); ...
    I = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    Mk = [row[:] for row in I]   # M_0 = I (we'll build)
    coeffs = [0]*(n+1)
    coeffs[n] = 1
    Acur = [[0]*n for _ in range(n)]
    # FL: N_1 = A, c1 = tr(N_1); N_{k} = A (N_{k-1} - c_{k-1} I); c_k = tr(N_k)/k
    # char poly det(xI-A) = x^n - c1 x^{n-1} - c2 x^{n-2} - ... (one convention)
    # We'll use the standard FL producing p(x)=det(xI-A).
    N = [[0]*n for _ in range(n)]   # N_0
    cs = [1]  # leading
    Nprev = [[1 if i==j else 0 for j in range(n)] for i in range(n)]  # N_0 = I
    cprev = 1
    # iterate k=1..n
    Acopy = M
    c_list = [1]  # c_0 corresponds to x^n coeff = 1
    Mprod = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    c = 1
    res = [1]  # res[k] = coefficient, will be c_k below
    # Use the clean form:
    # M_1 = A;            c_1 = -tr(M_1)
    # M_k = A(M_{k-1}) + c_{k-1} I ; c_k = -tr(A M_k)/k    (gives det(xI-A)=x^n + c_1 x^{n-1}+...+c_n)
    Mk = [[M[i][j] for j in range(n)] for i in range(n)]  # M_1 = A
    cc = [1]  # cc[0]=1 (coeff of x^n)
    def trace(X):
        return sum(X[i][i] for i in range(n))
    def matmul(X,Y):
        return [[sum(X[i][k]*Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    def adddiag(X, s):
        return [[X[i][j] + (s if i==j else 0) for j in range(n)] for i in range(n)]
    c1 = -trace(Mk)
    cc.append(c1)
    cprev = c1
    Mprev = Mk
    for k in range(2, n+1):
        # M_k = A(M_{k-1} + c_{k-1} I)
        tmp = adddiag(Mprev, cprev)
        Mk = matmul(M, tmp)
        ck_num = -trace(Mk)
        assert ck_num % k == 0, f"FL non-integer at k={k}"
        ck = ck_num // k
        cc.append(ck)
        Mprev = Mk
        cprev = ck
    # cc = [1, c1, c2, ..., cn] with det(xI-A) = x^n + c1 x^{n-1} + ... + cn
    # return as list of length n+1, index = power of x (low to high)
    poly = [0]*(n+1)
    for k in range(n+1):
        poly[n-k] = cc[k]
    return poly   # poly[j] = coeff of x^j

def det_int(M):
    """Exact integer determinant via fraction-free Bareiss."""
    n=len(M); A=[row[:] for row in M]; sign=1; prev=1
    for k in range(n-1):
        if A[k][k]==0:
            sw=None
            for r in range(k+1,n):
                if A[r][k]!=0: sw=r;break
            if sw is None: return 0
            A[k],A[sw]=A[sw],A[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                A[i][j]=(A[i][j]*A[k][k]-A[i][k]*A[k][j])//prev
        prev=A[k][k]
    return sign*A[n-1][n-1]

# ---------- tournaments ----------

def tour_from_bits(bits, n):
    A=[[0]*n for _ in range(n)]
    pos=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>pos)&1: A[i][j]=1
            else: A[j][i]=1
            pos+=1
    return A

def skew(A):
    n=len(A)
    return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]

# ---------- H via OCF / Held-Karp Hamiltonian path count ----------

def ham_path_count(A):
    """Number of Hamiltonian paths (directed) in tournament A = H(T). Held-Karp."""
    n=len(A)
    if n==1: return 1
    # dp[mask][v] = number of ham paths over vertices in mask ending at v
    from collections import defaultdict
    dp=[ [0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v]=1
    full=(1<<n)-1
    for mask in range(1<<n):
        for v in range(n):
            d=dp[mask][v]
            if d==0: continue
            for w in range(n):
                if mask&(1<<w): continue
                if A[v][w]==1:
                    dp[mask|(1<<w)][w]+=d
    return sum(dp[full][v] for v in range(n))

# ---------- even-length real-cycle face: I(Omega_even, x) ----------

def directed_cycles_by_parity(A):
    """Enumerate all simple directed cycles of T (real, following arcs),
    return dict length->list of frozenset(vertices) ... but a vertex-set may
    support multiple directed cycles, so return list of (frozenset, ) per
    distinct directed cycle.  We track vertex-sets for packing (disjointness).
    Returns: list of (length, frozenset_vertices) for EVEN cycles only, AND
    full list for all cycles (for H cross-check)."""
    n=len(A)
    cycles=[]  # (length, frozenset)
    # enumerate simple cycles via DFS fixing the smallest vertex as start
    for start in range(n):
        # path stack
        stack=[(start,[start],1<<start)]
        while stack:
            v,path,mask=stack.pop()
            for w in range(n):
                if A[v][w]!=1: continue
                if w==start and len(path)>=3:
                    cycles.append((len(path), frozenset(path)))
                elif w>start and not (mask&(1<<w)):
                    stack.append((w,path+[w],mask|(1<<w)))
    # each directed cycle counted once (start = min vertex, direction fixed by arcs)
    return cycles

def even_face_poly(A):
    """I(Omega_even, x) = sum over packings of vertex-disjoint EVEN directed
    cycles of x^{#cycles}.  Returns list coeff[k] = # packings with k cycles."""
    cyc=[c for c in directed_cycles_by_parity(A) if c[0]%2==0]
    # build conflict structure: pack disjoint vertex-sets
    sets=[c[1] for c in cyc]
    m=len(sets)
    # count independent sets by size in the intersection graph (disjoint = compatible)
    # DP over cycles
    coeff=[1]  # empty packing
    # recursive enumeration
    res=[0]*(m+1)
    res[0]=1
    def rec(idx, used, cnt):
        for j in range(idx, m):
            if sets[j] & used: continue
            res[cnt+1]+=1
            rec(j+1, used|sets[j], cnt+1)
    rec(0, frozenset(), 0)
    # trim
    while len(res)>1 and res[-1]==0: res.pop()
    return res

def H_from_ocf(A):
    """H = I(Omega,2) over ALL odd cycles, as cross-check."""
    cyc=[c for c in directed_cycles_by_parity(A) if c[0]%2==1]
    sets=[c[1] for c in cyc]
    m=len(sets)
    total=[0]  # will accumulate sum 2^cnt
    acc=[0]
    def rec(idx, used, cnt):
        acc[0]+= (1<<cnt)  # but we want to count each independent set once with weight 2^cnt; the empty set counted below
        for j in range(idx, m):
            if sets[j]&used: continue
            rec(j+1, used|sets[j], cnt+1)
    # the standard: I(G,2)=sum over independent sets 2^{|S|}. include empty.
    acc[0]=0
    def rec2(idx, used, cnt):
        # count THIS independent set (the one formed so far) once
        for j in range(idx, m):
            if sets[j]&used: continue
            acc[0]+=(1<<(cnt+1))  # adding cycle j -> independent set of size cnt+1, weight 2^{cnt+1}? careful
            rec2(j+1, used|sets[j], cnt+1)
    # simpler: enumerate all independent sets, sum 2^size
    acc[0]=1  # empty set weight 2^0=1
    def rec3(idx, used, cnt):
        for j in range(idx, m):
            if sets[j]&used: continue
            acc[0]+=(1<<(cnt+1))
            rec3(j+1, used|sets[j], cnt+1)
    rec3(0, frozenset(), 0)
    return acc[0]

# ---------- main ----------

def poly_is_even_in_x(poly):
    """poly[j]=coeff x^j. Check all ODD-power coeffs (other than possibly a
    global factor x for odd n) are zero. Return (is_even_up_to_x_factor, low)."""
    n=len(poly)-1
    low=0
    while low<=n and poly[low]==0: low+=1
    # after factoring x^low, remaining should be even
    ok=True
    for j in range(low, n+1):
        if (j-low)%2==1 and poly[j]!=0:
            ok=False;break
    return ok, low

def run_n(n, sample_limit=None):
    print(f"\n========== n={n} ==========")
    total = 1<<(n*(n-1)//2)
    bits_iter = range(total) if (sample_limit is None or total<=sample_limit) else None
    if bits_iter is None:
        # sample
        import random
        random.seed(12345)
        bits_list=[random.randrange(total) for _ in range(sample_limit)]
    else:
        bits_list=list(bits_iter)

    even_x_fail=0
    # group by A-charpoly
    data=[]       # (acp, scp, H)
    acp_to_scp={}
    acp_to_H={}
    acpscp_to_H={}
    for bits in bits_list:
        A=tour_from_bits(bits,n)
        S=skew(A)
        acp=tuple(charpoly_int(A))
        scp=tuple(charpoly_int(S))
        ok,low=poly_is_even_in_x(scp)
        if not ok: even_x_fail+=1
        H=ham_path_count(A)
        data.append((acp,scp,H))
        acp_to_scp.setdefault(acp,set()).add(scp)
        acp_to_H.setdefault(acp,set()).add(H)
        acpscp_to_H.setdefault((acp,scp),set()).add(H)

    print(f"  tournaments processed: {len(bits_list)} (exhaustive={sample_limit is None or total<=sample_limit})")
    print(f"  S-char-poly EVEN-in-x (odd cycles cancel): {'ALL PASS' if even_x_fail==0 else f'{even_x_fail} FAIL'}")

    # cospectral analysis
    n_acp=len(acp_to_scp)
    cospectral_classes=[acp for acp,scps in acp_to_scp.items() if len(scps)>1]
    print(f"  #distinct A-char-polys: {n_acp}")
    print(f"  #A-char-polys split by S-char-poly (skew is non-spectral here): {len(cospectral_classes)}")
    # how many distinct fingerprints
    n_acpscp=len(set((d[0],d[1]) for d in data))
    print(f"  #distinct (A-char):              {n_acp}")
    print(f"  #distinct (A-char, S-char):      {n_acpscp}")
    # does (A-char,S-char) determine H?
    undet=[k for k,v in acpscp_to_H.items() if len(v)>1]
    undet_a=[k for k,v in acp_to_H.items() if len(v)>1]
    print(f"  (A-char) alone leaves H ambiguous in {len(undet_a)} classes")
    print(f"  (A-char,S-char) leaves H ambiguous in {len(undet)} classes")
    return

if __name__=="__main__":
    for n in [3,4,5,6]:
        run_n(n)
    run_n(7, sample_limit=40000)
