from functools import lru_cache

def check(ok, msg):
    if not ok: raise RuntimeError(msg)

def hp(rows):
    n=len(rows)
    @lru_cache(None)
    def walks(end,remaining):
        if not remaining:return 1
        choices=rows[end]&remaining; total=0
        while choices:
            bit=choices&-choices; choices-=bit
            total += walks(bit.bit_length()-1,remaining^bit)
        return total
    full=(1<<n)-1
    return sum(walks(i,full^(1<<i)) for i in range(n))

def orient(rows,mask):
    n=len(rows)
    return tuple(row | ((1<<n) if not (mask>>i&1) else 0) for i,row in enumerate(rows))+(mask,)

def marked(rows,port):
    p=tuple(i for i in range(len(rows)) if i!=port)+(port,)
    r=tuple(sum(1<<j for j,v in enumerate(p) if rows[u]>>v&1) for u in p)
    h=hp(r);q=orient(r,1<<(len(r)-1));hq=hp(q);qq=orient(q,1<<(len(q)-1));hqq=hp(qq)
    check((hq-h)%2==0 and (hqq-hq)%2==0, 'parity')
    return h,(hq-h)//2,(hqq-hq)//2

def count_formula(rows):
    # Endpoint/start subset path counts, calculated independently of the target code's
    # singleton/pair interpolation. Partition each extended Hamilton path at new z.
    n=len(rows);size=1<<n;end=[[0]*n for _ in range(size)]
    start=[[0]*n for _ in range(size)]
    for i in range(n):end[1<<i][i]=start[1<<i][i]=1
    for mask in range(1,size):
        for v in range(n):
            if mask>>v&1 and mask!=(1<<v):
                previous=mask^(1<<v)
                end[mask][v]=sum(end[previous][u] for u in range(n) if previous>>u&1 and rows[u]>>v&1)
                start[mask][v]=sum(start[previous][u] for u in range(n) if previous>>u&1 and rows[v]>>u&1)
    full=size-1
    covers=[[0]*n for _ in range(n)]
    for mask in range(1,full):
        other=full^mask
        for u in range(n):
            if end[mask][u]:
                for v in range(n):
                    covers[u][v]+=end[mask][u]*start[other][v]
    def evaluate(d):
        return (sum(start[full][v] for v in range(n) if d>>v&1)
                +sum(end[full][u] for u in range(n) if not(d>>u&1))
                +sum(covers[u][v] for u in range(n) if not(d>>u&1) for v in range(n) if d>>v&1))
    return evaluate

sources={
 'U':((112,97,115,103,74,16,32),(31,14,25)),
 'U2':((252,249,194,228,204,148,32,64),(95,41,75)),
 'U4':((496,497,499,487,456,400,288,64,128),(97,40,75)),
 'R2':((1008,1017,1011,997,776,912,816,592,128,256),(285,123,225)),
 'R4':((996,993,994,999,975,912,800,576,128,256),(291,120,225))}
expected={'U':[(1,97),(33,95),(117,95),(125,97)],'U2':[(4,289)],'U4':[],'R2':[],'R4':[]}
for name,(rows,source) in sources.items():
    n=len(rows);check(marked(rows,n-1)==source,('source state',name))
    targets={(3*source[0]+2*j,3*source[1]-j,3*source[2]) for j in (1,2)}
    formula=count_formula(rows)
    hits=[]
    for m in range(1<<n):
        h=formula(m)
        if h in {v[0] for v in targets}:
            g=orient(rows,m)
            check(hp(g)==h,('literal target hit',name,m))
            hits.append((m,h))
            for p in range(n+1):check(marked(g,p) not in targets,('target state found',name,m,p))
    check(hits==expected[name],('hits',name,hits))
    print(name,'PASS all',1<<n,'masks, independently derived cut formula; hits',hits,'all ports checked')
