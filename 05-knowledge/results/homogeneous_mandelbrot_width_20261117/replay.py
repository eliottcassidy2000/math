from itertools import combinations,permutations

def need(condition):
    if not condition:raise ArithmeticError("failed")

def literal(T):
    n=len(T); out=[0]*(1<<n);out[0]=1
    for mask in range(1,1<<n):
        vs=[v for v in range(n) if mask>>v&1]
        out[mask]=sum(all(T[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(vs))
    return out

def profile(T,r):
    H=literal(T);full=len(H)-1
    scores=[0,0]
    for s in range(1,full+1):
        if not s>>r&1:continue
        vs=[v for v in range(len(T)) if v!=r and s>>v&1]
        for p in permutations(vs):
            w=(r,)+p
            if all(T[a]>>b&1 for a,b in zip(w,w[1:])):
                scores[s.bit_count()%2]+=H[full^s]
    return H[-1],H[full^(1<<r)],*scores

def run(T,a,b):
    old=[i for i in range(len(T)) if i not in (a,b)]
    # Enumerate trace words independently, without a gap convolution.
    d={}
    def go(word,U,i,j):
        if not U and i and j:d[i,j]=d.get((i,j),0)+1
        for v in U+[a,b]:
            if word and not (T[word[-1]]>>v&1):continue
            go(word+[v],[x for x in U if x!=v],i+(v==a),j+(v==b))
    go([],old,0,0)
    return d

core=(24,1,3,6,14)
for n in range(1,6):
    edges=list(combinations(range(n),2))
    found=[]
    for m in range(1<<len(edges)):
        t=[0]*n
        for i,(a,b) in enumerate(edges):
            if m>>i&1:t[a]|=1<<b
            else:t[b]|=1<<a
        if literal(t)[-1]==11:found.append(tuple(t))
    if n<5:need(not found)
    else:
        orbit={tuple(sum((core[a]>>b&1)<<i for i,b in enumerate(p)) for a in p)
               for p in permutations(range(n))}
        need(len(orbit)==120 and orbit==set(found))
for a,b in combinations(range(5),2):
    tab=run(core,a,b)
    if 0 in (a,b):need(any(max(ij)>2 for ij in tab))
    else:need(tab[2,2]==(3 if (a,b) in ((1,2),(3,4)) else 4)
             and tab[2,1]==tab[1,2]==6 and max(map(max,tab))==2)
print("H11 census and traces passed")
for p,T,r,z in [(19,(2,4,1),0,95011519192),(23,(2,4,25,19,35,15),2,2330279290103)]:
    H,g,E,O=profile(T,r)
    need(O==E+g)
    while z:
        if z&1:H,E,O=H+2*E,E+O-g,E+O
        else:E,O=E+O,E+O
        z>>=1
    need(H==(2**p+1)**2+2)
    print("period",3*p,"Q",H,"H",H*(2**(p-1)-1))
