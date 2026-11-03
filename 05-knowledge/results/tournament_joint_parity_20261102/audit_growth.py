from collections import Counter
from itertools import permutations, combinations
from math import comb


def ck(ok, text):
    if not ok: raise RuntimeError(text)


def join(blocks):
    # Blocks consist of named vertices, with a directed triangle specified by order.
    n=sum(len(b) for b in blocks)
    arcs=set()
    for i,b in enumerate(blocks):
        if len(b)==3:
            arcs.update(zip(b,b[1:]+b[:1]))
        else: ck(len(b)==1, 'unsupported block')
        for c in blocks[i+1:]: arcs.update((u,v) for u in b for v in c)
    ck(len(arcs)==comb(n,2),'tournament size')
    return n,arcs


def paths(vertices, arcs):
    return [p for p in permutations(vertices) if all(e in arcs for e in zip(p,p[1:]))]


def stats(n,arcs):
    verts=tuple(range(n))
    hp=paths(verts,arcs)
    start=Counter(p[0] for p in hp)
    end=Counter(p[-1] for p in hp)
    P=[]
    C=Counter()
    for k in range(1,n+1):
        count=0
        for p in permutations(verts):
            for cuts in combinations(range(1,n),k-1):
                if all((p[i-1],p[i]) in arcs for i in range(1,n) if i not in cuts):
                    count+=1
                    if k==2: C[(p[cuts[0]-1],p[cuts[0]])]+=1
        P.append(count)
    hist=Counter()
    vector=[]
    for b in range(1<<n):
        extended=arcs|{(i,n) if b>>i&1 else (n,i) for i in verts}
        h=len(paths(tuple(range(n+1)),extended))
        formula=sum(start[i]*(1-(b>>i&1))+end[i]*(b>>i&1) for i in verts)
        formula+=sum(c*(b>>i&1)*(1-(b>>j&1)) for (i,j),c in C.items())
        ck(h==formula,('endpoint expansion',b))
        hist[h]+=1
        vector.append(h)
    aut=sum(all(((p[i],p[j]) in arcs)==((i,j) in arcs) for i in verts for j in verts if i!=j) for p in permutations(verts))
    return {'H':len(hp),'aut':aut,'covers':P,'hist':dict(sorted(hist.items())),
            'moment1':sum(vector),'moment2':sum(x*x for x in vector),
            'start':[start[i] for i in verts],'end':[end[i] for i in verts],
            'C':[[C[i,j] for j in verts] for i in verts], 'vector':vector}

U=stats(*join(((0,),(1,),(2,3,4))))
V=stats(*join(((0,),(1,2,3),(4,))))
for name,d in [('U',U),('V',V)]: print(name,d)
ck(U['covers']==V['covers']==[3,42,162,240,120],'cover profile')
ck(U['aut']==V['aut']==3,'aut')
ck(U['moment1']==V['moment1']==432,'mean')
ck([U['moment2'],V['moment2']]==[8064,8640],'moments')
# Exhaust all labelled n<=4, group by H. No canonicalizer is involved.
for n in range(5):
    hist_by_H={}
    es=list(combinations(range(n),2))
    for mask in range(1<<len(es)):
        arcs={(i,j) if mask>>k&1 else (j,i) for k,(i,j) in enumerate(es)}
        h=len(paths(tuple(range(n)),arcs))
        histogram=Counter()
        for cut in range(1<<n):
            child=arcs|{(i,n) if cut>>i&1 else (n,i) for i in range(n)}
            histogram[len(paths(tuple(range(n+1)),child))]+=1
        signature=tuple(sorted(histogram.items()))
        if h in hist_by_H: ck(hist_by_H[h]==signature,('smaller hostile',n,h,mask))
        hist_by_H[h]=signature
    print('sharpness',n,hist_by_H)
# Independently transport the literal joins into the parent producer's mask gauge.
for blocks,transport,code in [(((0,),(1,),(2,3,4)),(4,3,1,0,2),2),
                              (((0,),(1,2,3),(4,)),(4,1,3,2,0),32)]:
    n,arcs=join(blocks)
    transported={(transport[u],transport[v]) for u,v in arcs}
    es=list(combinations(range(n),2))
    encoded=sum((1<<k) for k,(i,j) in enumerate(es) if (i,j) in transported)
    ck(encoded==code,('named transport',code))
print('PASS: literal joins, exact named transports, profile, histograms, moments, sharp order-five boundary')
