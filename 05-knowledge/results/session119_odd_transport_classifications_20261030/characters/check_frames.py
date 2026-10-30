"""Literal digit maps, node automorphisms and Hamiltonian path controls."""
from itertools import product
import json

def require(ok,why):
    if not ok:raise RuntimeError(why)
def chi(x,p):return 1 if pow(x%p,(p-1)//2,p)==1 else -1

def sign(d,p,a,frame):
    d%=p**a
    if not d:return 0
    v=0
    while d%p==0:d//=p;v+=1
    return frame[v]*chi(d,p)
def digits(x,p,a):
    out=[]
    for _ in range(a):out.append(x%p);x//=p
    return out
def number(ds,p):return sum(x*p**i for i,x in enumerate(ds))
def check_map(mp,p,a,src,dst):
    n=p**a
    require(len(set(mp))==n,'bijection')
    for x in range(n):
        for y in range(n):
            require(sign(mp[y]-mp[x],p,a,dst)==sign(y-x,p,a,src),('arc map',p,a,x,y))

def hamiltonian(p,a,frame):
    n=p**a;D=[[0]*n for _ in range(1<<n)]
    for x in range(n):D[1<<x][x]=1
    for mask in range(1,1<<n):
        for x in range(n):
            value=D[mask][x]
            if not value:continue
            for y in range(n):
                if not(mask>>y&1) and sign(y-x,p,a,frame)==1:
                    D[mask|1<<y][y]+=value
    return sum(D[-1])

rows=[];maps=automorphisms=0
for p,a in [(3,1),(3,2),(3,3),(3,4),(7,1),(7,2),(11,1),(11,2),(19,1)]:
    n=p**a; frames=list(product((-1,1),repeat=a));nonsquare=next(u for u in range(1,p) if chi(u,p)==-1)
    for src in frames:
        for dst in frames:
            slopes=[1 if s==t else nonsquare for s,t in zip(src,dst)]
            mp=[number([(d*slopes[i])%p for i,d in enumerate(digits(x,p,a))],p) for x in range(n)]
            check_map(mp,p,a,src,dst);maps+=1
    # Verify literal independent actions at every rooted low-digit prefix.
    # Translation1 and all square slopes generate each displayed local affine group.
    nodes=0
    for depth in range(a):
        for prefix in product(range(p),repeat=depth):
            nodes+=1
            local=[(1,1)]+[(u,0) for u in range(1,p) if chi(u,p)==1]
            for slope,b in local:
                mp=[]
                for x in range(n):
                    ds=digits(x,p,a)
                    if tuple(ds[:depth])==prefix:ds[depth]=(slope*ds[depth]+b)%p
                    mp.append(number(ds,p))
                check_map(mp,p,a,(1,)*a,(1,)*a);automorphisms+=1
    require(nodes==(n-1)//(p-1),'tree node count')
    base_order=p*(p-1)//2
    group_order=base_order**nodes
    row={'p':p,'a':a,'vertices':n,'labelled_frames':2**a,'isomorphism_classes':1,'certified_subgroup_order':group_order}
    if n<=11:
        H=hamiltonian(p,a,(1,)*a)
        require(H%group_order==0,'free Hamiltonian action divisibility')
        require(H%2==1,'odd Hamiltonian control')
        row['Hamiltonian_paths']=H;row['quotient_by_subgroup_order']=H//group_order
    rows.append(row)
print(json.dumps({'status':'PASS','literal_frame_isomorphisms':maps,'literal_node_automorphisms':automorphisms,'rows':rows,'scope':'explicit subgroup only; no full automorphism or new Ramsey-number claim'},indent=2,sort_keys=True))
