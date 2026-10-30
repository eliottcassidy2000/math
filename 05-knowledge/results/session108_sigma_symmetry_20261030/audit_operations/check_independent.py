"""Independent relation-object and fiber-action audit; standard library only."""
from itertools import permutations, product, combinations
from fractions import Fraction
from math import factorial, gcd, prod
import json

def check(v,msg):
    if not v: raise RuntimeError(msg)
def span(vs):
    s={0}
    for v in vs: s |= {x^v for x in list(s)}
    return frozenset(s)
def spaces(n):
    seen={frozenset({0})}; pending=list(seen)
    while pending:
        s=pending.pop()
        for v in range(1<<n):
            if v not in s:
                t=frozenset(set(s)|{x^v for x in s})
                if t not in seen: seen.add(t);pending.append(t)
    return seen
def coord(v,p):
    return sum(((v>>i)&1)<<p[i] for i in range(len(p)))
def act(o,p):
    if isinstance(o,tuple):return (frozenset(coord(x,p) for x in o[0]),coord(o[1],p))
    return frozenset(coord(x,p) for x in o)
def orbit_count(objects,aut):
    unseen=set(objects); count=0
    while unseen:
        x=unseen.pop();unseen.difference_update(act(x,p) for p in aut);count+=1
    return count
def tours(n):
    pairs=list(combinations(range(n),2)); ps=list(permutations(range(n)))
    def edges(mask):return frozenset((b,a) if mask>>i&1 else (a,b) for i,(a,b) in enumerate(pairs))
    unseen={edges(m) for m in range(1<<len(pairs))};out=[]
    while unseen:
        e=unseen.pop();aut=[]
        for p in ps:
            f=frozenset((p[a],p[b]) for a,b in e)
            unseen.discard(f)
            if f==e:aut.append(p)
        source=any(all((i,j) in e for j in range(n) if j!=i) for i in range(n))
        out.append((e,aut,source))
    return out
GL=[]
for bs in product(range(1,8),repeat=3):
    if len(span(bs))==8:
        GL.append(tuple((bs[0] if x&1 else 0)^(bs[1] if x&2 else 0)^(bs[2] if x&4 else 0) for x in range(8)))
AGL=[tuple(x^t for x in g) for g in GL for t in range(8)]
marked=[g for g in GL if all((g[x]&1)==(x&1) for x in range(8))]
groups={'GL':GL,'AGL':AGL,'marked':marked}
def parts(n,lo=1):
    if not n:yield ()
    for i in range(lo,n+1,2):
        for p in parts(n-i,i):yield (i,)+p
def burn(n,group):
    ans=Fraction()
    for g in group:
        for p in parts(n):
            e=sum((i-1)//2 for i in p)+sum(gcd(a,b) for a,b in combinations(p,2))
            z=prod(i**p.count(i)*factorial(p.count(i)) for i in set(p))
            fixed=[]
            for l in p:
                count=0
                for x in range(8):
                    y=x
                    for _ in range(l):y=g[y]
                    count+=y==x
                fixed.append(count)
            ans+=Fraction((1<<e)*prod(fixed),z)
    ans/=len(group);check(ans.denominator==1,'Burnside integrality');return int(ans)
counts={k:[] for k in groups};sources={k:[] for k in groups};split={}
for n in range(5):
    ss=spaces(n);one=(1<<n)-1
    objects={'GL':[w for w in ss if len(w)<=8], 'AGL':[w for w in ss if one in w and len(w)<=16], 'marked':[(w,x) for w in ss if len(w)<=8 for x in w if x or len(w)<=4]}
    for name,gg in groups.items():
        entries=[(len(aut),orbit_count(objects[name],aut),has_source) for _,aut,has_source in tours(n)]
        total=sum(e[1] for e in entries); src=sum(e[1] for e in entries if e[2])
        check(total==burn(n,gg),f'relation vs Burnside {name} {n}')
        counts[name].append(total);sources[name].append(src)
        if n==3:split[name]=sorted(entries)
        if n:
            remaining=set(range(8));pred=0
            while remaining:
                c=remaining.pop();remaining.difference_update(g[c] for g in gg)
                pred+=burn(n-1,[g for g in gg if g[c]==c])
            check(src==pred,f'source stabilizer {name} {n}')
check(counts=={'GL':[1,2,5,24,188],'AGL':[1,1,2,8,48],'marked':[1,3,11,70,756]},'published counts')
# Independently generate each exact-observer permutation from four independent fiber flips.
obs=lambda x:x>>1
fiber_maps=[tuple(x^(bits[x>>1]) for x in range(8)) for bits in product(range(2),repeat=4)]
def affine(f):return all((f[x^y]^f[0])==(f[x]^f[y]) for x in range(8) for y in range(8))
def commutes(f):return all(f[x^(1<<i)]==(f[x]^(1<<i)) for x in range(8) for i in range(3))
ops=[len(fiber_maps),sum(affine(f) for f in fiber_maps),sum(commutes(f) for f in fiber_maps)]
check(ops==[16,8,2],'operation groups')
grid=list(product(range(3),range(2)))
grid_ok=[]
for p in permutations(range(6)):
    if all((grid[i][d]==grid[j][d])==(grid[p[i]][d]==grid[p[j]][d]) for i in range(6) for j in range(6) for d in range(2)):grid_ok.append(p)
check(len(grid_ok)==12,'grid automorphisms')
# Divisor sums by divisor incidence sieve, not multiplicative prime-factor code.
N=6000;sigma=[0]*(N+1)
for d in range(1,N+1):
    for n in range(d,N+1,d):sigma[n]+=d
probes=0
for a in range(2,31):
    for n in range(1,201):
        if gcd(a,n)==1:
            check(sigma[a*n]-a*n==(sigma[a]-a)*n+sigma[a]*(sigma[n]-n),'coprime probe');probes+=1
check(sigma[15]-15==9 and sigma[57]-57==23,'examples')
check((sigma[3],sigma[81],sigma[343])==(4,121,400),'neutral prime-power values')
check(pow(sigma[3]%11,5,11)==1 and sigma[81]%11==0,'character zero guard')
check(sigma[9]-9 != 3+4*(sigma[3]-3),'shared prime hostile')
check(sigma[14]==sigma[23]==24 and sigma[42]==sigma[69]==96,'sigma noninjective')
# Direct two-vertex refresh: tuple entries are winner, loser, avoiding vertex
# canonicalization and using neither U/D nor the predicted spectrum to build K.
def signature(pair,name):
    rows=[sum(((c>>b)&1)<<i for i,c in enumerate(pair)) for b in range(3)]
    w=span(rows)
    if name=='AGL':return span(rows+[3])
    if name=='marked':return (w,rows[0])
    return w
def rank(mat):
    a=[list(map(Fraction,row)) for row in mat];r=0
    for c in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];v=a[r][c];a[r]=[x/v for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][c]:
                v=a[i][c];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        r+=1
    return r
# Focused final integration check: derive spaces from actual color lists.
projection_dimensions={}
for name,pairs in [('GL',[(1,0),(1,1)]),('AGL',[(0,1,1),(0,0,1)])]:
    dims=[]
    for colors in pairs:
        rows=[sum(((color>>b)&1)<<i for i,color in enumerate(colors)) for b in range(3)]
        if name=='AGL':rows.append((1<<len(colors))-1)
        w=span(rows);projected={x>>1 for x in w}
        dims.append((len(w).bit_length()-1,len(projected).bit_length()-1))
    check(dims==([(1,0),(1,1)] if name=='GL' else [(2,1),(2,2)]),'rank quotient source-deletion obstruction')
    projection_dimensions[name]=dims
spectra={}
for name in groups:
    reps={signature(pair,name):pair for pair in product(range(8),repeat=2)}
    keys=list(reps);index={key:i for i,key in enumerate(keys)};size=len(keys)
    K=[[Fraction() for _ in keys] for _ in keys]
    for key,pair in reps.items():
        row=K[index[key]]
        for removed in range(2):
            survivor=pair[1-removed]
            for color in range(8):
                for out in [(survivor,color),(color,survivor)]:
                    row[index[signature(out,name)]]+=Fraction(1,32)
    mult={str(ev):size-rank([[K[i][j]-(ev if i==j else 0) for j in range(size)] for i in range(size)]) for ev in [Fraction(0),Fraction(1,2),Fraction(1)]}
    expected={'GL':{'0':3,'1/2':1,'1':1},'AGL':{'0':1,'1/2':0,'1':1},'marked':{'0':8,'1/2':2,'1':1}}[name]
    check(mult==expected,'direct refresh eigenspaces');spectra[name]=mult
print(json.dumps({'status':'PASS','method':'relation-subspace decorations over tournament automorphism groups; independent fiber flips; divisor incidence sieve; direct refresh matrices','counts':counts,'source_counts':sources,'degree3_splits_autorder_count_source':split,'operational_orders':ops,'neutral_grid_order':len(grid_ok),'coprime_probes':probes,'rank_projection_dimensions':projection_dimensions,'degree2_refresh_eigenspace_dimensions':spectra,'h_GL':[counts['GL'][i]-counts['GL'][i-1] for i in range(1,5)]},indent=2,sort_keys=True))
