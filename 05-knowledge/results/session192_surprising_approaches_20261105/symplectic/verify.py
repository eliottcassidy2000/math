#!/usr/bin/env python3
"""Integer symplectic compiler, sharp rank-two hostile, and smooth witnesses.
No third-party packages. Every verification gate remains active with python -O.
"""
import hashlib,itertools,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def require(test, message):
    if not test:raise ValueError(message)

def omega(u,v):
    # Coordinates are (q1,p1,q2,p2), NOT (q1,q2,p1,p2).
    return u[0]*v[1]-u[1]*v[0]+u[2]*v[3]-u[3]*v[2]

def compile_word(word):
    w=list(map(int,word));K=3*len(w)+1
    points=[(1,0,K,0),(1,0,0,K)]
    for i,c in enumerate(w,1):
        s=-1 if c&2 else 1;t=1 if c&1 else -1
        points.append((1,3*i,s,t))
    return points

def word_graph(word):
    w=list(map(int,word));n=len(w)+2;r=[0]*n;r[0]=2
    for i,c in enumerate(w,2):
        for a in (0,1):
            u,v=(a,i) if c>>a&1 else(i,a);r[u]|=1<<v
        for j in range(2,i):r[j]|=1<<i
    return r

def point_graph(points):
    n=len(points);r=[0]*n;least=None
    for i in range(n):
        for j in range(i+1,n):
            p=omega(points[i],points[j]);require(p!=0,'nonzero pairings')
            least=abs(p) if least is None else min(least,abs(p))
            u,v=(i,j) if p>0 else(j,i);r[u]|=1<<v
    return r,least

def path_count(rows):
    n=len(rows)
    if n==0:return 1
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,1<<n):
        for u,count in enumerate(d[mask]):
            if not count:continue
            targets=rows[u]&~mask
            while targets:
                bit=targets&-targets;targets-=bit
                d[mask|bit][bit.bit_length()-1]+=count
    return sum(d[-1])

def permutation_count(rows):
    return sum(all(rows[u]>>v&1 for u,v in zip(p,p[1:]))
               for p in itertools.permutations(range(len(rows))))

def transitive_on(rows,mask):
    names=[i for i in range(len(rows)) if mask>>i&1]
    return not any((rows[a]>>b&1) and (rows[b]>>c&1) and (rows[c]>>a&1)
                   for a,b,c in itertools.permutations(names,3))

def local_order(rows):
    full=(1<<len(rows))-1
    return all(transitive_on(rows,rows[v]) and
               transitive_on(rows,full^rows[v]^(1<<v)) for v in range(len(rows)))

def full_graph(n,mask):
    r=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            u,v=(i,j) if mask>>k&1 else(j,i);r[u]|=1<<v;k+=1
    return r

def primes_to(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,math.isqrt(p)+1))]

BANK={r['H']:r['word'] for r in json.loads((ROOT/'inherited_prefix_witnesses.json').read_text())}

def factor_smooth(N):
    require(N>0 and N%2==1,'positive odd input')
    out={};r=N
    for p in primes_to(139):
        while r%p==0:out[p]=out.get(p,0)+1;r//=p
    require(r==1,'input must be 139-smooth')
    require(N not in (7,21),'excluded targets 7 and 21')
    e=out.pop(7,0);f=[]
    if e==1:
        p=next((p for p,k in sorted(out.items()) if p>=5 and k),None)
        if p is not None:
            f.append(7*p);out[p]-=1
        else:
            require(out.get(3,0)>=2,'the only uncovered cases are 7 and 21')
            f.append(63);out[3]-=2
    elif e>=2:
        if e%2:f.append(343);e-=3
        f.extend([49]*(e//2))
    for p,k in sorted(out.items()):f.extend([p]*k)
    require(math.prod(f)==N,'factorization reconstructs target')
    require(all(h in BANK for h in f),'all generators have certified words')
    return f

def ordered_join(points_blocks):
    if not points_blocks:return [(1,0,0,0)]
    require(all(p[0]==1 for b in points_blocks for p in b),'common affine chart')
    B=max(abs(p[1]) for b in points_blocks for p in b)
    C=max(max(abs(p[2]),abs(p[3])) for b in points_blocks for p in b)
    L=2*B+2*C*C+1
    points=[];block_ids=[]
    for j,block in enumerate(points_blocks):
        points.extend((p[0],p[1]+j*L,p[2],p[3]) for p in block)
        block_ids.extend([j]*len(block))
    rows,margin=point_graph(points)
    offset=0
    for j,block in enumerate(points_blocks):
        old,_=point_graph(block)
        for u in range(len(block)):
            for v in range(len(block)):
                require(((rows[offset+u]>>(offset+v))&1)==((old[u]>>v)&1),'join preserves block')
        offset+=len(block)
    require(all(rows[u]>>v&1 for u in range(len(points)) for v in range(len(points))
                if block_ids[u]<block_ids[v]),'all cross arcs point forward')
    return points

def main():
    receipt={'coordinate_order':['q1','p1','q2','p2'],
       'inherited_prefix_sha256':hashlib.sha256((ROOT/'inherited_prefix_witnesses.json').read_bytes()).hexdigest(),
       'status':'All arithmetic, sign, and finite count checks passed. All-order assertions use the written proofs.'}
    require(set(BANK)==set(range(1,1002,2))-{7,21},'exact inherited target prefix')
    base=[]
    for h,w in sorted(BANK.items()):
        p=compile_word(w);r,margin=point_graph(p)
        require(r==word_graph(w),'compiled word incidences')
        require(path_count(r)==h,'independent subset DP on compiled geometry')
        base.append({'H':h,'word':w,'order':len(p),'minimum_margin':margin})
    receipt['base_witnesses_checked']=len(base)
    receipt['maximum_base_order']=max(r['order'] for r in base)
    receipt['global_minimum_base_margin']=min(r['minimum_margin'] for r in base)
    compiler_checks=0
    for m in range(6):
        for w in itertools.product(range(4),repeat=m):
            w=''.join(map(str,w));r,_=point_graph(compile_word(w))
            require(r==word_graph(w),'complete small alphabet sign control')
            compiler_checks+=1
    receipt['complete_small_word_sign_checks']=compiler_checks
    # An independent full labelled universe, not only cut-word presentations.
    rounds={};round_checks=0
    for n in range(1,6):
        values=set();count=0
        for mask in range(1<<(n*(n-1)//2)):
            r=full_graph(n,mask)
            if local_order(r):
                h=path_count(r);require(h==permutation_count(r),'literal permutation control')
                values.add(h);count+=1
            round_checks+=1
        rounds[str(n)]={'labelled_local_orders':count,'H_values':sorted(values)}
    require(all(13 not in x['H_values'] for x in rounds.values()),'H13 rank-two finite gate')
    require(5**5>13**3,'strong-floor order cutoff for H13')
    receipt['all_labelled_small_tournaments']=round_checks
    receipt['rank_two_finite_gate']=rounds
    p=compile_word('120');r,margin=point_graph(p)
    require(path_count(r)==permutation_count(r)==13,'explicit rank-four H13 witness')
    require(not local_order(r),'H13 is not locally transitive')
    # {a,b,x1,x2} principal Pfaffian certifies actual rank at least four.
    pf=omega(p[0],p[1])*omega(p[2],p[3])-omega(p[0],p[2])*omega(p[1],p[3])+omega(p[0],p[3])*omega(p[1],p[2])
    require(pf==120,'rank-four nonzero Pfaffian')
    receipt['H13']={'points':p,'rows':r,'margin':margin,'pfaffian_0_1_2_3':pf}
    generators=sorted(set(primes_to(139))-{2,7}|{49,63,343}|{7*p for p in primes_to(139) if p>=5 and p!=7})
    require(all(g in BANK for g in generators),'finite smooth generator bank')
    receipt['smooth_generators']=[{'H':g,'word':BANK[g]} for g in generators]
    # Exhausts input integers, not just random exponent vectors.
    integers=0
    for N in range(1,100002,2):
        t=N
        for p0 in primes_to(139):
            while t%p0==0:t//=p0
        if t!=1 or N in (7,21):continue
        factor_smooth(N);integers+=1
    receipt['smooth_factorization_checks_through_100001']=integers
    joins=[]
    for N in (1,3,13,39,49,63,147,343,973,4459,7**11*3**9*139**3*137**4):
        factors=factor_smooth(N);blocks=[compile_word(BANK[h]) for h in factors]
        points=ordered_join(blocks);rows,margin=point_graph(points)
        actual=path_count(rows) if len(points)<=14 else None
        if actual is not None:require(actual==N,'direct joined-geometry path count')
        joins.append({'H':N,'factors':factors,'order':len(points),'margin':margin,
                      'direct_subset_count':actual,'coordinate_height':max(map(abs,itertools.chain.from_iterable(points)))})
    receipt['join_controls']=joins
    receipt['constructor_output']={'139_smooth_cutoff':139,'excluded':[7,21],
       'positive_dimension_bound':4,'first_globally_attainable_target_requiring_dimension_4':13}
    (ROOT/'base_receipt.json').write_text(json.dumps(base,indent=2)+'\n')
    out=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    (ROOT/'receipt.json').write_text(out);print(out,end='')

if __name__=='__main__':main()
