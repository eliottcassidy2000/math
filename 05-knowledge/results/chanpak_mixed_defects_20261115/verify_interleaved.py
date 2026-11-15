from itertools import product

def c(d): return 1 if d==1 else 1<<(d-2)

def formula(x):
    w=1
    while x:
        x,b=divmod(x,2)
        if not b: w += x+(x&1)
    return w

def brute(x,m):
    n=m+1
    arcs=[0]*n
    for i in range(m):
        arcs[i]=sum(1<<j for j in range(i+1,m))
        if (x>>(m-i-1))&1: arcs[m]|=1<<i
        else: arcs[i]|=1<<m
    tab=[[0]*n for _ in range(1<<n)]
    for i in range(n):tab[1<<i][i]=1
    for s in range(1,1<<n):
        for i,a in enumerate(tab[s]):
            if not a:continue
            dst=arcs[i]&~s
            while dst:
                j=dst&-dst
                tab[s|j][j.bit_length()-1]+=a
                dst-=j
    return sum(tab[-1])

def need(b):
    if not b:raise RuntimeError('exact kernel mismatch')

def rev(x,m): return int(f'{x:0{m}b}'[::-1],2)
def rc(x,m): return rev(x,m)^((1<<m)-1)

for x,v in zip((1158,1190,1114,1146),(721,723,725,727)):
    need(brute(x,11)==formula(x)==v)

for m in range(1,11):
    N=1<<m
    hs=list(map(formula,range(N)))
    old={}
    for x,v in enumerate(hs):old.setdefault(v,[]).append(x)
    for x,v in enumerate(hs):
        for y,z in product(old.get(v+2,()),old.get(v+4,())):
            need(((x^y)&(x^z)) or hs[x^y^z]!=v+6)

for a,b in product(range(20),repeat=2):
    ws=[]
    for x in (1158,1190,1114,1146):
        full=((a<<11)|x)<<max(b.bit_length(),1) | rc(b,max(b.bit_length(),1))
        ws.append(formula(full))
    need(ws[1]-ws[0]==ws[3]-ws[2]==2+32*(b-a))
    need(ws[2]-ws[0]==ws[3]-ws[1]==4+56*a-44*b)
    if a==b:need(ws[0]==2048*a*a + 2428*a+719+2*formula(a))

for m in range(1,10):
    fibers={}
    for x in range(1<<m):
        key=x-rev(x,m)
        sides=tuple(((x>>i)&1)-((x>>(m-i-1))&1) for i in range(m//2))
        need(fibers.setdefault(key,sides)==sides)
        free=[i for i in range((m+1)//2) if ((x>>i)&1)==((x>>(m-i-1))&1)]
        aligned=x
        for i in free:aligned &= ~((1<<i)|(1<<(m-i-1)))
        wt=0
        for j,i in enumerate(free):
            for k in free[:j]:
                if ((x>>i)^(x>>k))&1:
                    wt+=2*c(i-k)
                    if i!=m-i-1:wt+=2*c(m-1-i-k)
        need(formula(x)-formula(aligned)==wt)
print('literal rectangle DP, smaller masks, exterior responses, and exact cut fibers verified')
