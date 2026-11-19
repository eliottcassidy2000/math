def need(ok, detail):
    if not ok:
        raise AssertionError(detail)

def subpaths(rows, tags=None):
    n=len(rows)
    D=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        if tags is None or tags>>i&1:
            D[1<<i][i]=1
    for s,row in enumerate(D):
        for i,c in enumerate(row):
            if c:
                w=rows[i]&~s
                while w:
                    b=w&-w;w-=b
                    D[s|b][b.bit_length()-1]+=c
    return [sum(x) for x in D]

def H(rows):
    n=len(rows)
    need(all(0<=a<1<<n and a>>i&1==0 and
        all((a>>j&1)+(rows[j]>>i&1)==1 for j in range(i))
        for i,a in enumerate(rows)),rows)
    return subpaths(rows)[-1] if rows else 1

def attach(rows,tags):
    n=len(rows)
    return [q+((1<<n) if tags>>i&1==0 else 0)
            for i,q in enumerate(rows)]+[tags]

from functools import lru_cache

@lru_cache(None)
def cached_profile(rows,tags):
    return profile(list(rows),tags)

def profile(rows,tags):
    hp=subpaths(rows);hp[0]=1
    D=subpaths(rows,tags);F=(1<<len(rows))-1
    n=H(rows);k=H(attach(rows,tags))
    return (n,k,*[sum(D[s]*hp[F^s] for s in range(1,F+1)
                       if s.bit_count()%2==t) for t in (0,1)])

def star(A,a,B,b,m):
    nb,nA=len(B),len(A);u=nb+m+nA
    ans=[0]*(u+1)
    for i,r in enumerate(B):
        ans[i]=r+((1<<u)-(1<<nb))+((1<<u) if b>>i&1==0 else 0)
    for j in range(m):
        ans[nb+j]=(1<<(u+1))-(1<<(nb+j+1))
    for i,r in enumerate(A):
        ans[nb+m+i]=(r^((1<<nA)-1)^(1<<i))<<(nb+m)
        if a>>i&1:ans[nb+m+i]|=1<<u
    ans[u]=b+sum(1<<(nb+m+i) for i in range(nA) if a>>i&1==0)
    return ans

def F(u):
    h=1
    while u:
        q=u//2
        if not u&1:h+=q+(q&1)
        u=q
    return h

def apex(v):
    n=v.bit_length()
    return [(1<<n)-(1<<(i+1))+(1<<n if v>>i&1 else 0)
            for i in range(n)]+[((1<<n)-1)^v]

def f1(b,e):
    L=apex(b-1);C=apex(e-1)
    R=[a^((1<<len(C))-1)^(1<<i) for i,a in enumerate(C)]
    nb=len(L);top=nb+len(R)-1;out=L+[q<<nb for q in R]
    out[nb-1]|=1<<top
    for i in range(nb,len(out)):
        out[i]|=(1<<(nb-1))-1 if i==top else (1<<nb)-1
    return out

for j,(b,e) in enumerate([(2,15),(0,0),(2,68),(2,109),(3,164),
                         (2,526),(2,907),(2,2196),(2,3675)]):
    N=(38<<j)+1
    if j==1:
        rows=[74,32,3,6,15,29,62]
    else:
        rows=f1(b,e)
        need(b*e+F(b-1)*F(e-1)-(b&1)*(e&1)==N,(j,b,e))
    need(H(rows)==N,(j,rows))
A=[6,12,8,1];a=12
B=[4,5,24,19,67,31,175,63];b=76
need(profile(A,a)==(5,9,4,7),A)
need(profile(B,b)==(79,525,294,314),B)
for m in range(1,4):
    T=star(A,a,B,b,m)
    need(H(T)==1+9728*2**m,(m,T))
print('PASS: defect-14 bridge, literal seed profiles, and all nine smaller heights')

A=[0,21,25,3,9];a=2
B=[2,4,9,19,7,31];b=56
need(profile(A,a)==(5,19,7,7),A)
need(profile(B,b)==(9,101,63,65),B)
for m in range(1,4):
    need(H(star(A,a,B,b,m))==1+2432*2**m,m)

from itertools import product
from pathlib import Path
bank = {}
for row in Path(__file__).with_name('star_seam_bank.tsv').read_text().splitlines():
    if row.startswith('#') or not row.strip():
        continue
    q,c,a,A,b,B=row.split()
    q,c,a,b=map(int,(q,c,a,b))
    A=tuple(map(int,A.split(','))) if A!='-' else ()
    B=tuple(map(int,B.split(','))) if B!='-' else ()
    ha,ka,ea,oa=cached_profile(A,a)
    hb,kb,eb,ob=cached_profile(B,b)
    need(ha*(kb-2*eb)+hb*(ka-ha)==1,(q,A,B))
    need(c==(ha+ea+oa)*(eb+ob) and c//(c&-c)==q,(q,c))
    need(q not in bank,q)
    bank[q]=c

bound=2055
total=0
max_small=0
for n in range(1,(2*bound+1).bit_length()):
    for p in product('01',repeat=n-1):
        h,e,t=1,0,2
        for b in (*p,'0','1'):
            h,e,t=(h+t-e,1,2*t-1) if b=='1' else (h,0,2*t)
        defect=t-h
        if defect>bound:
            continue
        total+=1
        c=t-1;q=c//(c&-c)
        need(q in bank,(defect,p,q))
        # All ray values after width 1 are exact; if a smaller ray term
        # precedes the star construction, it lies below the prior cap.
        if c<2*bank[q]:
            max_small=max(max_small,bank[q]+1)
        need(bank[q]+1<2**32,(q,bank[q]))
print('PASS: checked',len(bank),'literal module pairs on',total,
      'canonical ray prefixes; earlier terms <=',max_small)
