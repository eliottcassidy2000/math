import itertools, math, pickle, os
def strong(bits,n):
    a=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if bits[k]=='1': a[i]|=1<<j
            else: a[j]|=1<<i
            k+=1
    full=(1<<n)-1
    def reach(adj):
        seen=1; fr=1
        while fr:
            nf=0; x=fr
            while x:
                b=x&-x; v=b.bit_length()-1; x^=b; nf|=adj[v]
            fr=nf&~seen; seen|=nf
        return seen==full
    ra=[0]*n
    for i in range(n):
        for j in range(n):
            if a[i]>>j&1: ra[j]|=1<<i
    return reach(a) and reach(ra)
if os.path.exists('sp.pkl'): Sp=pickle.load(open('sp.pkl','rb'))
else:
    Sp={1:{1},3:{3},4:{5}}
    for n in range(5,10):
        Sp[n]={int(l.split()[1]) for l in open(f'hl{n}.txt') if strong(l.split()[0],n)}
    pickle.dump(Sp,open('sp.pkl','wb'))
def S35(r):  # is rational r = p/q a {3,5}-unit?
    p,q=r
    for x in (p,q):
        for b in (3,5):
            while x%b==0: x//=b
        if x!=1: return False
    return True
def test(total,comps):
    prods={c:{} for c in comps}
    for c in comps:
        a,b=c
        for x in Sp[a]:
            for y in Sp[b]: prods[c].setdefault(x*y,[]).append((x,y))
    out=[]
    for c1,c2 in itertools.combinations(comps,2):
        common=set(prods[c1])&set(prods[c2])
        expl=0; unexpl=[]
        for v in sorted(common):
            ok=False
            for (x,y) in prods[c1][v]:
                for (x2,y2) in prods[c2][v]:
                    g=math.gcd(x,x2)
                    if S35((x//g,x2//g)): ok=True; break
                if ok: break
            if ok: expl+=1
            else: unexpl.append((v,prods[c1][v][0],prods[c2][v][0]))
        out.append((c1,c2,len(common),expl,unexpl))
    return out
for total,comps in ((11,[(4,7),(5,6),(3,8)]),(12,[(4,8),(5,7),(6,6),(3,9)])):
    print(f"== {total} vertices, two strong components")
    for c1,c2,n,e,u in test(total,comps):
        print(f"  {c1} vs {c2}: shared counts {n}, brick-move explained {e}, unexplained {len(u)}  e.g. {u[:3]}")
