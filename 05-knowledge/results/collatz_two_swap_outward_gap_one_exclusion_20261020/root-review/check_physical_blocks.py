import math,json
shapes=[(2,3),(5,8),(7,11),(12,19),(41,65),(53,84),(306,485),(665,1054)]
expected=[(0,None),(4,1),(6,1),(10,1),(34,551),(44,11),(254,1),(552,37)]
rows=[]
for (p,q),(count,mx) in zip(shapes,expected):
    w=[((i+1)*p)//q-(i*p)//q for i in range(q)]
    gs=[]
    for i in range(q):
        idx=[(i+j)%q for j in range(4)]
        if len(set(idx))<4 or [w[t] for t in idx]!=[0,1,1,0]:continue
        for order in [(0,1),(1,0)]:
            v=w[:]
            edges=[idx[:2],idx[2:]]
            for e in order:
                a,b=edges[e]
                if v[a]==v[b]:raise RuntimeError('illegal edge')
                v[a],v[b]=v[b],v[a]
            B=0
            for j,bit in enumerate(v):
                if bit:B=3*B+(1<<j)
            D=abs((1<<q)-3**p)
            g=math.gcd(B,D)
            if g==D:raise RuntimeError('integer output')
            gs.append(g)
    actual=(len(gs),max(gs,default=None))
    if actual!=(count,mx):raise RuntimeError((p,q,actual,(count,mx)))
    rows.append(dict(p=p,q=q,ordered_paths=len(gs),maximum_g=max(gs,default=None)))
print(json.dumps(dict(status='PASS',total=sum(r['ordered_paths'] for r in rows),rows=rows),sort_keys=True,indent=2))
