# Verify: near-join of T1 = TT_n + x_w (port x) and T2 = converse(TT_m + x'_u) (port x'), back arc x' -> x,
# has H = (1+u)(1+v) + g(u) g(v) - e(u) e(v), where v = val(w), u = val(w2), e(v) = [v==0 or v even].
import itertools, sys
sys.setrecursionlimit(10000)
def hcount(n,out):
    full=(1<<n)-1
    end=[[0]*n for _ in range(1<<n)]
    for v in range(n): end[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            c=end[S][v]
            if not c: continue
            nb=out[v]&~S
            while nb:
                w=(nb&-nb).bit_length()-1; nb&=nb-1
                end[S|(1<<w)][w]+=c
    return sum(end[full])
def g(v):
    if v==0: return 1
    if v&1: return g(v>>1)
    u=v>>1; return g(u)+u+(u&1)
def e(v): return 1 if (v==0 or v%2==0) else 0
bad=0; cnt=0
for n in range(0,5):
  for m in range(0,5):
    for w in itertools.product([0,1],repeat=n):
      for w2 in itertools.product([0,1],repeat=m):
        # T1: chain c_0..c_{n-1} (i->j for i<j), x; x->c_i iff w[i]
        # T2 = converse of (chain d_0..d_{m-1}, x'; x'->d_i iff w2[i]): in T2, d_j -> d_i for i<j; d_i -> x' iff w2[i]
        N=n+1+m+1
        X=n; Xp=n+1+m
        idx1=list(range(n))+[X]
        idx2=list(range(n+1,n+1+m))+[Xp]
        out=[0]*N
        def arc(a,b): out[a]|=1<<b
        for i in range(n):
            for j in range(i+1,n): arc(i,j)
            if w[i]: arc(X,i)
            else: arc(i,X)
        for i in range(m):
            for j in range(i+1,m): arc(n+1+j,n+1+i)   # converse chain
            if w2[i]: arc(n+1+i,Xp)                    # converse of x'->d_i
            else: arc(Xp,n+1+i)
        # T1 => T2 except back arc x' -> x
        for a in idx1:
            for b in idx2:
                if a==X and b==Xp: arc(Xp,X)
                else: arc(a,b)
        H=hcount(N,out)
        v=int(''.join(map(str,w)),2) if n else 0
        u=int(''.join(map(str,w2)),2) if m else 0
        pred=(1+u)*(1+v)+g(u)*g(v)-e(u)*e(v)
        cnt+=1
        if H!=pred:
            bad+=1
            if bad<5: print('BAD',w,w2,H,pred)
print('near-join two one-apex checks',cnt,'failures',bad)
