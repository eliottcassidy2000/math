#!/usr/bin/env python3
"""Standalone standard-library exact checks; all guards active under -O."""
from itertools import combinations, permutations, product
from collections import Counter
import json

def require(x, msg):
    if not x: raise RuntimeError(msg)

def tournament(n, code):
    a=[[0]*n for _ in range(n)]
    for b,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=(code>>b)&1; a[j][i]=1-a[i][j]
    return a

def triangles(a):
    return tuple((i,j,k) for i,j,k in combinations(range(len(a)),3)
                 if a[i][j]==a[j][k]==a[k][i])

def paths(a):
    n=len(a); dp=[[0]*n for _ in range(1<<n)]; H=[1]+[0]*((1<<n)-1)
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if mask>>v&1:
                prev=mask^(1<<v)
                if prev: dp[mask][v]=sum(dp[prev][u] for u in range(n) if prev>>u&1 and a[u][v])
        H[mask]=sum(dp[mask])
    return H,dp[-1]

def isos(a,b):
    n=len(a); da=list(map(sum,a)); db=list(map(sum,b)); p=[-1]*n; used=set(); count=0
    order=sorted(range(n),key=lambda i:db.count(da[i]))
    def visit(k):
        nonlocal count
        if k==n: count+=1; return
        i=order[k]
        for j in range(n):
            if j in used or da[i]!=db[j]: continue
            if any(a[i][u]!=b[j][p[u]] for u in order[:k]): continue
            p[i]=j; used.add(j); visit(k+1); used.remove(j)
    visit(0); return count

def seed(bit):
    a=[[0,1,1,1],[0,0,1,0],[0,0,0,1],[0,1,0,0]]
    return [[a[j][i] if bit else a[i][j] for j in range(4)] for i in range(4)]

def compose(children):
    m=len(children[0]); a=[[0]*(3*m) for _ in range(3*m)]
    for i in range(3*m):
        for j in range(3*m):
            x,u=divmod(i,m); y,v=divmod(j,m)
            a[i][j]=children[x][u][v] if x==y else int((y-x)%3==1)
    return a

def build(word):
    if len(word)==1:return seed(word[0])
    m=len(word)//3; return compose([build(word[i*m:(i+1)*m]) for i in range(3)])

def kind(word):
    if len(word)==1:return str(word[0]),0
    m=len(word)//3; children=[kind(word[i*m:(i+1)*m]) for i in range(3)]
    types=tuple(x[0] for x in children)
    return '('+','.join(min(types[i:]+types[:i] for i in range(3)))+')',sum(x[1] for x in children)+int(len(set(types))==1)

def next_poly(p):
    cube=Counter()
    for i,a in p.items():
        for j,b in p.items():
            for k,c in p.items():cube[i+j+k]+=a*b*c
    for i,a in p.items():cube[3*i]-=a
    require(all(a%3==0 for a in cube.values()),'orbit integrality')
    q=Counter({i:a//3 for i,a in cube.items() if a})
    for i,a in p.items():q[3*i+1]+=a
    return q

def main():
    census=[]
    for n in range(1,7):
        bank={}
        for code in range(1<<(n*(n-1)//2)):
            a=tournament(n,code); c=triangles(a); h,_=paths(a)
            if c in bank: require(bank[c]==h,'C3/profile equivalence')
            else:bank[c]=h
        census.append({'n':n,'tournaments':1<<(n*(n-1)//2),'labelled_C3_images':len(bank)})
    # Direct run-reversal bijection for every tournament/module through order four.
    moves=0; checked_paths=0
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            a=tournament(n,code)
            for mask in range(1,1<<n):
                inside={i for i in range(n) if mask>>i&1}
                if any(len({a[v][u] for u in inside})!=1 for v in set(range(n))-inside):continue
                b=[r[:] for r in a]
                for u in inside:
                    for v in inside:b[u][v]=a[v][u]
                moves+=1
                for size in range(1,n+1):
                    for path in permutations(range(n),size):
                        if not all(a[u][v] for u,v in zip(path,path[1:])):continue
                        out=[]; i=0
                        while i<size:
                            j=i+1
                            if path[i] in inside:
                                while j<size and path[j] in inside:j+=1
                                out.extend(reversed(path[i:j]))
                            else:out.append(path[i])
                            i=j
                        require(all(b[u][v] for u,v in zip(out,out[1:])), 'run map')
                        require(set(path)==set(out),'label subset'); checked_paths+=1
    # All eight depth-one colourings: full H profile and global groups.
    base=build((0,0,0)); baseline,_=paths(base); groups=[]
    for word in product(range(2),repeat=3):
        a=build(word); _,e=kind(word)
        require(triangles(a)==triangles(base),'family triples')
        require(paths(a)[0]==baseline,'family profile')
        g=isos(a,a); require(g==3**(3+e),'global group')
        groups.append(g)
    # Exhaust every binary word through depth two, without global group enumeration at 36.
    p=Counter({0:2}); table=[]
    for h in range(4):
        if h<=2:
            classes={}
            for word in product(range(2),repeat=3**h):
                t,e=kind(word);classes[t]=e
            require(Counter(classes.values())==p,'class polynomial')
        require(set(p)==set(range((3**h-1)//2+1)), 'gapless exponent support')
        table.append({'h':h,'order':4*3**h,'classes':sum(p.values()),'upper_exponent_counts':dict(sorted(p.items()))})
        p=next_poly(p)
    # Endpoint counts fail already on the four-vertex seed.
    ha,ea=paths(seed(0)); hb,eb=paths(seed(1))
    require(ha==hb and ea!=eb,'endpoint hostile')
    # Nonmodule reversal changes even C3/H: reverse one edge in a C3.
    a=tournament(3,5); b=[r[:] for r in a]; b[0][1],b[1][0]=b[1][0],b[0][1]
    require(triangles(a)!=triangles(b) and paths(a)[0]!=paths(b)[0],'nonmodule hostile')
    print(json.dumps({'census':census,'module_moves':moves,'direct_paths_checked':checked_paths,
      'depth_one_global_group_orders':groups,'family_tables':table,'seed_endpoint_counts':[ea,eb],
      'seed_H':ha[-1],'depth_one_H':baseline[-1]},sort_keys=True,indent=2))
if __name__=='__main__':main()
