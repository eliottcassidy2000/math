"""Standalone exact shared-core cut-exchange certificate; no repo imports."""
import itertools,json,sys

def need(ok,label):
    if not ok: raise RuntimeError(label)

def reverse(a):
    return [((1<<len(a))-1)^row^(1<<i) for i,row in enumerate(a)]

def tables(a):
    n=len(a); full=1<<n
    end=[[0]*full for _ in a]; h=[0]*full; h[0]=1
    for mask in range(1,full):
        for j in range(n):
            if not(mask>>j&1):continue
            old=mask^(1<<j)
            end[j][mask]=1 if old==0 else sum(end[i][old] for i in range(n) if old>>i&1 and a[i]>>j&1)
            h[mask]+=end[j][mask]
    return h,end

def coordinates(a):
    h,end=tables(a);_,start=tables(reverse(a));n=len(a);full=(1<<n)-1
    av=[sum(h[m]*start[i][full^m] for m in range(full+1)) for i in range(n)]
    W=[[0]*n for _ in a]
    for i in range(n):
        for j in range(i):
            W[i][j]=W[j][i]=sum(end[i][m]*start[j][full^m]+end[j][m]*start[i][full^m] for m in range(full+1))
    return h[-1],av,W

def closure(a,cut):
    n=len(a);g=[row|(1<<n) for row in a]+[0]
    return [row|((1<<(n+1)) if not(cut>>i&1) else 0) for i,row in enumerate(g)]+[cut]

def strong(a):
    for i in range(len(a)):
        seen=1<<i
        while True:
            new=seen
            for j in range(len(a)):
                if seen>>j&1:new|=a[j]
            if new==seen:break
            seen=new
        if seen!=(1<<len(a))-1:return False
    return True

def scalar(h,av,W,cut):
    ids=[i for i in range(len(av)) if cut>>i&1]
    return h+2*sum(av[i] for i in ids)-sum(W[i][j] for i,j in itertools.combinations(ids,2))


# New small certificate only; not a historical-result audit.
def literal_count(g):
    return sum(all(g[p[i]]>>p[i+1]&1 for i in range(len(p)-1))
               for p in itertools.permutations(range(len(g))))

def safe_bounds(av,W,p,q,k):
    ds=sorted(W[q][j]-W[p][j] for j in range(len(av)) if j not in (p,q))
    m=k-1
    return (2*(av[q]-av[p])-sum(ds[len(ds)-m:]),
            2*(av[q]-av[p])-sum(ds[:m]))

core=[24,17,3,6,12]
h,av,W=coordinates(core)
need(h==15 and av==[21]*5,'physical coordinates')
rows=[]
for S in itertools.combinations(range(5),2):
    cut=sum(1<<i for i in S);g=closure(core,cut)
    actual=literal_count(g);formula=scalar(h,av,W,cut)
    need(actual==formula,'literal permutation check')
    need(strong(g),'strong closure')
    rows.append({'cut':list(S),'H':actual})
need(sorted(set(r['H'] for r in rows))==[89,91],'physical slice')
for p,q in itertools.combinations(range(5),2):
    lo,hi=safe_bounds(av,W,p,q,2)
    need(lo>=-2 and hi<=2,'safe edge')
# Gauge W-8, a-2 for this k=2 slice.
for S in itertools.combinations(range(5),2):
    need(h+sum(2*(av[i]-2) for i in S)-sum(W[i][j]-8 for i,j in itertools.combinations(S,2))==99-W[S[0]][S[1]],'gauge')
# Hostile: transitive core has an unsafe swap.
tt=[sum(1<<j for j in range(i+1,5)) for i in range(5)]
ht,at,wt=coordinates(tt)
bounds=safe_bounds(at,wt,0,4,2)
need(bounds[0]<-2 or bounds[1]>2,'hostile unsafe')
# Dense model. It is not asserted physical.
n,k,L=8,4,100
values=sorted({1+sum(L*(k-1)+2*i for i in S)-2*L*(k*(k-1)//2) for S in itertools.combinations(range(n),k)})
need(values==list(range(1+k*(k-1),2+k*(2*n-k-1),2)),'dense model interval')
print(json.dumps({'core':core,'h':h,'a':av,'W':W,'closures':rows,'safe_graph':'K5','hostile_transitive_bounds':bounds,'formal_dense_interval':[values[0],values[-1]],'scope':'10 literal seven-vertex closures; one formal dense model; no global coverage'},indent=2))
