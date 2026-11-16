from itertools import product
from fractions import Fraction

def figs(w,k=1,ab=True):
  a=[int(c)//2 for c in w];b=[int(c)%2 for c in w];m=len(w)
  acc=[0]*3
  for z in product(range(3),repeat=m):
    I=tuple(i for i,j in enumerate(z) if j==1);J=tuple(i for i,j in enumerate(z) if j==2)
    ni=len(I);nj=len(J);empty=m-ni-nj
    A=ni>0 and not ni%2 and a[I[0]] and not a[I[-1]]
    B=nj>0 and not nj%2 and b[J[0]] and not b[J[-1]]
    if (A and not J) or (B and not I):acc[0]+=2*k**(empty+2)
    if A and B:acc[1]+=4*k**(empty+2)
    if (ni+nj)%2 and (a[I[0]] and not b[I[-1]] if I else ab) and (b[J[0]] and not a[J[-1]] if J else not ab):
      acc[2]+=2*k**(empty+1)
  return acc
def fit(vals):
  A=[[Fraction(x) for x in (3**t,2**t,1,(-1)**t,vals[t-1])] for t in range(1,5)]
  for j in range(4):
    i=next(i for i in range(j,4) if A[i][j]);A[j],A[i]=A[i],A[j]
    a=A[j][j];A[j]=[z/a for z in A[j]]
    for i in range(4):
      if i!=j:
        a=A[i][j]; A[i]=[x-a*y for x,y in zip(A[i],A[j])]
  return tuple(r[-1] for r in A)
# Assignment-to-run automaton. Each run stores (unused, open odd/even,
# closed odd/even) encoded 0..4. k=1 throughout the linear certificate.
EPS=(((0,0),),((1,1),),((0,0),(1,1)),((0,1),(1,0)))
ST=[list(product(range(5),repeat=len(ep))) for ep in EPS]
OFF=[0,5,10,35]
D=60
def mat(symbol):
    bits=(symbol//2,symbol%2)
    rows=[[] for _ in range(D)]
    for ep,S,off in zip(EPS,ST,OFF):
        idx={v:i+off for i,v in enumerate(S)}
        for v,i in idx.items():
            rows[i].append(i)
            for j,(source,sink) in enumerate(ep):
                p=v[j]
                if p>=3 or p==0 and not bits[source]:
                    continue
                q=2 if p==1 else 1
                repl=lambda z:v[:j]+(z,)+v[j+1:]
                rows[i].append(idx[repl(q)])
                if not bits[sink]:
                    rows[i].append(idx[repl(q+2)])
    return rows
M=list(map(mat,range(4)))
def push(v,m):
    out=[0]*D
    for i,x in enumerate(v):
        if x:
            for j in m[i]:out[j]+=x
    return tuple(out)
def back(m,v):
    return tuple(sum(v[j] for j in t) for t in m)
def prefix(w):
    v=tuple(int(i in OFF) for i in range(D))
    for x in w:v=push(v,M[int(x)])
    return v
def observ(arc):
    r=[[0]*D for _ in range(3)]
    r[0][4]=r[0][9]=2
    r[1][34]=4
    for i,(a,b) in enumerate(ST[-1],35):
        if ((a>=3 or a==0 and arc) and (b>=3 or b==0 and not arc)
                and ((a==3)+(b==3))==1):
            r[2][i]=2
    return tuple(map(tuple,r))
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def advance(v,w):
    # Actual paired transfer plus virtual coordinates (rho,p,beta,1).
    return (*push(v[:D],M[w]),*push(v[D:2*D],M[w]),
            v[-4]+w*v[-3],2*v[-3],2*v[-2]+w*v[-1],v[-1])
def guard_controls():
    # Verify four table rows via direct ternary assignment enumeration.
    want={('1230',True):((0,16,-4,0),(18,-8,-6,0),(18,-2,-2,0)),
          ('1230',False):((0,16,-4,0),(18,-8,-6,0),(18,-6,0,0)),
          ('1320',True):((0,18,-4,0),(18,-8,-6,0),(18,2,-2,0)),
          ('1320',False):((0,18,-4,0),(18,-8,-6,0),(18,-4,0,0))}
    for (w,arc),cs in want.items():
        vals=[figs(w+'0'*t,ab=arc) for t in range(1,5)]
        assert tuple(fit([a[i] for a in vals]) for i in range(3))==cs
    print("Exact four-state square table checked by literal run assignments.")

def span_certificate():
    ci=[Fraction(x) for x in (*prefix('1320'),*prefix('1230'),0,1,0,1)]
    rows={}
    queue=[ci]
    while queue:
        v=queue.pop()
        for i in sorted(rows):
            a=v[i]
            if a:v=[x-a*y for x,y in zip(v,rows[i])]
        try:i=next(j for j,z in enumerate(v) if z)
        except StopIteration:continue
        v=[x/v[i] for x in v]
        rows[i]=v
        queue.extend(advance(v,w) for w in (0,1))
    ao,bo=observ(True),observ(False)
    for t in range(1,5):
        ao=tuple(back(M[0],v) for v in ao)
        bo=tuple(back(M[0],v) for v in bo)
        for v in rows.values():
            r,p,beta,o=v[-4:]
            actual=[dot(v[:D],b)-dot(v[D:2*D],a) for a,b in zip(ao,bo)]
            expect=[2**(t+1)*p-2*r,0,2*o-6*r-2**(t+1)*p]
            assert actual==expect,(t,actual,expect)
            same=[sum(dot(v[:D],r0)-dot(v[D:2*D],r0) for r0 in obs0)
                for obs0 in (ao,bo)]
            assert same==[(6*p+beta)*2**t-2*r,(4*p+beta)*2**t-2*r]

    print(f"Invariant Q-span dimension {len(rows)}: all binary w and t>=1 pass.")

def subset_h(w,arc):
    n=len(w)+2;a=n-2;b=n-1
    out=[sum(1<<j for j in range(i+1,a)) for i in range(a)]+[0,0]
    for i,z in enumerate(map(int,w)):
        for v,flag in ((a,z//2),(b,z%2)):
            if flag:out[v]|=1<<i
            else:out[i]|=1<<v
    if arc:out[a]|=1<<b
    else:out[b]|=1<<a
    Dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):Dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,mass in enumerate(Dp[mask]):
            if not mass:continue
            row=out[i]&~mask
            while row:
                bit=row&-row;row-=bit
                Dp[mask|bit][bit.bit_length()-1]+=mass
    return sum(Dp[-1])

def graph_controls():
    for w in ('','0','1','01','10','11'):
        rho=sum(int(x)*(1<<i) for i,x in enumerate(w))
        for t in (1,2):
            for pre in ('1230','1320'):
                for arc in (True,False):
                    z=pre+w+'0'*t
                    assert subset_h(z,arc)==1+sum(figs(z,ab=arc))
            assert subset_h('1320'+w+'0'*t,False)-subset_h('1230'+w+'0'*t,True)==2-8*rho
    print("Literal adjacency subset DP checked on 48 physical tournaments.")

if __name__=='__main__':
    guard_controls()
    span_certificate()
    graph_controls()
