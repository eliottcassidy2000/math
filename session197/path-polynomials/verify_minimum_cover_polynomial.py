from verify_matching_amplitude import *
from itertools import permutations,combinations

def induced(a,vs):return tuple(sum((a[u]>>v&1)<<j for j,v in enumerate(vs)) for u in vs)
def sidepoly(a,side,C,prefix,mandatory):
    n=len(side);allpaths=[]
    for l in range(1,n+1):
        for p in permutations(range(n),l):
            if all(a[side[u]]>>side[v]&1 for u,v in zip(p,p[1:])):
                allpaths.append((sum(1<<i for i in p),side[p[-1] if prefix else p[0]]))
    states={0:1}
    for c in C:
        opts={}
        if c not in mandatory:opts[0]=1
        for mask,v in allpaths:
            ok=(a[v]>>c&1) if prefix else (a[c]>>v&1)
            if ok:opts[mask]=opts.get(mask,0)+1
        nxt={}
        for used,w in states.items():
            for mask,count in opts.items():
                if not used&mask:nxt[used|mask]=nxt.get(used|mask,0)+w*count
        states=nxt
    out=[0]
    for used,w in states.items():out=add(out,scale(bin_polynomial(induced(a,[v for i,v in enumerate(side) if not used>>i&1])),w))
    return out

def check(a,r):
    O=[i for i in range(len(a)) if a[r]>>i&1];I=[i for i in range(len(a)) if a[i]>>r&1]
    edges=[(u,v) for u in O for v in I if a[u]>>v&1]
    nu,mu=cross_matching(a,r);op=skeleton_operator(a,r);top=[0]*len(a)
    for (power,shift),c in op.items():
        if shift==nu:top[power]+=c
    top=trim(top);n=0
    for C in combinations(O+I,nu):
        C=set(C)
        if not all(u in C or v in C for u,v in edges):continue
        L=sidepoly(a,[u for u in O if u not in C],sorted(C),True,C&set(I))
        R=sidepoly(a,[v for v in I if v not in C],sorted(C),False,C&set(O))
        rise=[1]
        for j in range(nu):rise=mul(rise,[j,1])
        got=mul(rise,mul(L,R))
        require(got==top,('factor mismatch',a,r,C,top,L,R,got))
        require(all(v>=0 for v in L+R),'ordinary coefficient positivity')
        require(len(L)-1==len(O)-len(C&set(O))-len(C&set(I)),'left degree')
        require(len(R)-1==len(I)-len(C&set(I))-len(C&set(O)),'right degree')
        require(L[-1]*R[-1]==mu,'maximum matching multiplicity factors')
        require(value(top,1)==factorial(nu)*value(L,1)*value(R,1),'scalar amplitude factorization');n+=1
    return n
cases=[]
for n in range(1,5):
    for mask in range(1<<(n*(n-1)//2)):
        a=tournament(n,mask)
        for r in range(n):cases.append((a,r))
for bits in range(16):cases.append((balanced(2,bits,0),0))
cases.append((balanced(3,0b100010001,0),0))
count=0
for a,r in cases:count+=check(a,r)
print(json.dumps({'contexts':len(cases),'all_minimum_cover_factorizations':count,'status':'PASS'},sort_keys=True,indent=2))
