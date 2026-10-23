#!/usr/bin/env python3
"""Independent all-size dyadic/switch checks and sparse blue-kernel projection.
Only stdlib is needed for projection; SymPy is used for small exact determinants.
No inherited research producer is imported.
"""
import argparse, gzip, hashlib, itertools as it, json, time
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(p,s):
    if not p: raise RuntimeError(s)

def components(n,edges):
    adj=[[] for _ in range(n)]
    for u,v,w in edges:
        if w and u!=v: adj[u].append(v);adj[v].append(u)
    result=[];seen=set()
    for u in range(n):
        if u in seen:continue
        seen.add(u);todo=[u];block=[]
        while todo:
            v=todo.pop();block.append(v)
            for w in adj[v]:
                if w not in seen:seen.add(w);todo.append(w)
        result.append(sorted(block))
    return result

def rank2(rows):
    basis={}
    for z in rows:
        while z:
            p=z.bit_length()-1
            if p in basis:z^=basis[p]
            else:basis[p]=z;break
    return len(basis)

def nullbasis2(rows,n):
    piv={}
    for z in rows:
        while z:
            p=(z&-z).bit_length()-1
            if p in piv:z^=piv[p]
            else:piv[p]=z;break
    for p in sorted(piv,reverse=True):
        for q in list(piv):
            if q<p and (piv[q]>>p)&1:piv[q]^=piv[p]
    out=[]
    for j in range(n):
        if j not in piv:
            z=1<<j
            for p,row in piv.items():
                if (row>>j)&1:z|=1<<p
            out.append(z)
    return out

def bits(x):
    while x:
        b=x&-x;yield b.bit_length()-1;x^=b

def lap(n,edges):
    A=[[0]*n for _ in range(n)]
    for u,v,w in edges:
        if u!=v:A[u][u]+=w;A[v][v]+=w;A[u][v]-=w;A[v][u]-=w
    return A

def projection(n,blue,black):
    """black means paired-black, and loops remain in input but vanish in A."""
    gc=components(n,blue+black);bc=components(n,blue)
    roots={c[-1] for c in gc};keep=[v for v in range(n) if v not in roots];ix={v:i for i,v in enumerate(keep)}
    rows=[0]*len(keep)
    for u,v,w in blue:
        if u==v or w%2==0:continue
        delta=(1<<ix[u] if u in ix else 0)^(1<<ix[v] if v in ix else 0)
        if u in ix:rows[ix[u]]^=delta
        if v in ix:rows[ix[v]]^=delta
    basis=nullbasis2(rows,len(keep));dim=len(basis)
    membership=[0]*n
    for j,z in enumerate(basis):
        for v in bits(z):membership[keep[v]]|=1<<j
    q=[0]*dim
    # Black contribution is a binary incidence Gram matrix.
    for u,v,w in black:
        if w%2:
            d=membership[u]^membership[v]
            for i in bits(d):q[i]^=d
    # Blue half-energy requires integer sums before division by two.
    blue_raw={}
    for u,v,w in blue:
        if u==v:continue
        a=membership[u];b=membership[v];d=a^b
        ds=[(i,1 if a>>i&1 else -1) for i in bits(d)]
        for i,di in ds:
            for j,dj in ds:blue_raw[i,j]=blue_raw.get((i,j),0)+w*di*dj
    qb=[0]*dim
    for (i,j),v in blue_raw.items():
        need(v%2==0,'blue harmonic half pairing is integral')
        if v//2%2:qb[i]^=1<<j;q[i]^=1<<j
    need(all(((q[i]>>j)&1)==((q[j]>>i)&1) for i in range(dim) for j in bits(q[i])),'q symmetry')
    rank=rank2(q)
    return {'vertices':n,'full_components':len(gc),'blue_components':len(bc),
      'blue_bicycle_dimension':dim-len(bc)+len(gc),'critical_two_rank':dim,
      'half_form_rank':rank,'factors_divisible_by_four':dim-rank,
      'blue_half_form_rank':rank2(qb),'blue_half_form_alternating':not any(z>>i&1 for i,z in enumerate(qb)),
      'full_half_form_alternating':not any(z>>i&1 for i,z in enumerate(q)),
      '_q':q,'_members':membership,'_basis':basis,'_keep':keep}

def public(p):return {k:v for k,v in p.items() if not k.startswith('_')}

def fixed_path_census(n):
    free=[(i,j) for i in range(n) for j in range(i+2,n)];pairs=list(it.combinations(range(n),2));top=(1<<len(free))-1
    reflection=[free.index((n-1-j,n-1-i)) for i,j in free];keys=[];col=[];cache={}
    def canon(rows):
        if rows in cache:return cache[rows]
        deg=[r.bit_count() for r in rows];blocks=[[v for v in range(n) if deg[v]==d] for d in sorted(set(deg))]
        best=min(sum(((rows[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(pairs)) for pp in it.product(*(it.permutations(b) for b in blocks)) for p in [sum(pp,())])
        key=(tuple(sorted(deg)),best);cache[rows]=key;return key
    for x in range(top+1):
        rows=[0]*n
        for i in range(n-1):rows[i]|=1<<(i+1)
        for k,(i,j) in enumerate(free):
            u,v=(i,j) if x>>k&1 else (j,i);rows[u]|=1<<v
        rev=tuple(((1<<n)-1)^r^(1<<i) for i,r in enumerate(rows));keys.append(min(canon(tuple(rows)),canon(rev)))
        col.append(all((x>>i&1)==(x>>j&1) for i,j in enumerate(reflection)))
    nodes=sorted(set(keys));ix={key:i for i,key in enumerate(nodes)};blue=Counter();black=Counter()
    for x in range(top+1):
        if x<top^x:(blue if col[x] else black)[tuple(sorted((ix[keys[x]],ix[keys[top^x]])))]+=1
    need(all(v%2==0 for v in black.values()),'literal twins')
    return len(nodes),[[u,v,w] for (u,v),w in sorted(blue.items())],[[u,v,w//2] for (u,v),w in sorted(black.items())]

def exhaustive_tower(A,maxk=3):
    """Every residue lift is enumerated, independently of Smith or elimination."""
    n=len(A);layers=[];bases=[]
    for k in range(1,maxk+2):
        mod=1<<k;lifts={}
        for x in it.product(range(mod),repeat=n):
            if all(sum(a*b for a,b in zip(row,x))%mod==0 for row in A):lifts.setdefault(tuple(a%2 for a in x),x)
        bases.append(lifts)
        if k>maxk:break
        vals=list(lifts);q={}
        for u in vals:
            for v in vals:
                raw=sum(lifts[u][i]*A[i][j]*lifts[v][j] for i in range(n) for j in range(n));need(raw%mod==0,'tower integral')
                q[u,v]=(raw//mod)%2
        radical={u for u in vals if all(q[u,v]==0 for v in vals)}
        layers.append((radical,{'k':k,'V_size':len(vals),'radical_size':len(radical)}))
    for k,(rad,_) in enumerate(layers):need(rad==set(bases[k+1]),'higher radical equals next liftable space')
    return [p for _,p in layers]

def controls():
    import sympy as sp
    from sympy.matrices.normalforms import smith_normal_form
    reports=[];switches=[]
    for n in range(3,7):
        count,blue,black=fixed_path_census(n);p=projection(count,blue,black)
        A=sp.Matrix(lap(count,blue+[[u,v,2*w] for u,v,w in black]));snf=smith_normal_form(A[:-1,:-1],domain=sp.ZZ)
        ds=[abs(int(snf[i,i])) for i in range(count-1)];need(all(ds),'finite exact connectedness')
        need(sum(d%2==0 for d in ds)==p['critical_two_rank'],'2 rank vs Smith')
        need(sum(d%4==0 for d in ds)==p['factors_divisible_by_four'],'4 rank vs Smith')
        reports.append(dict(n=n,**public(p)))
        if n!=5:continue
        Q=p['_q'];mem=p['_members'];base_tree=int(A[:-1,:-1].det());tested=detected=0;invisible=None
        for e,f in it.combinations(range(len(black)),2):
            a,b,w=black[e];c,d,z=black[f]
            if len({a,b,c,d})<4:continue
            for c,d in [(c,d),(d,c)]:
                new=[(u,v,z-(i==e)-(i==f)) for i,(u,v,z) in enumerate(black) if z-(i==e)-(i==f)>0]+[(a,c,1),(b,d,1)]
                if len(components(count,blue+new))!=1:continue
                p2=projection(count,blue,new);need(p2['_members']==mem,'fixed kernel basis')
                x=mem[a]^mem[d];y=mem[b]^mem[c];pred=[0]*len(Q)
                for i in bits(x):pred[i]^=y
                for i in bits(y):pred[i]^=x
                actual=[u^v for u,v in zip(Q,p2['_q'])];need(actual==pred,'exact switch form')
                test=bool(x and y and x!=y);need(bool(any(actual))==test,'independence iff detection')
                tested+=1;detected+=test
                if not test and invisible is None:
                    A2=sp.Matrix(lap(count,blue+[(u,v,2*z) for u,v,z in new]));tree=int(A2[:-1,:-1].det())
                    if tree!=base_tree:
                        invisible={'remove':[[a,b,1],[c,d,1]],'add':[[a,c,1],[b,d,1]],'same_full_half_form':True,'same_blue_graph':True,'same_color_degrees':True,'original_tree':base_tree,'changed_tree':tree,'blue_edges':blue,'original_paired_black_edges':black,'changed_paired_black_edges':new}
        switches.append({'n':5,'connected_unit_switches_checked':tested,'half_form_detected':detected,'half_form_invisible':tested-detected,'invisible_control':invisible})
    mats=[[[2]],[[4]],[[8]],[[0,2],[2,0]],[[2,2],[2,6]],[[3,-1],[-1,3]],[[2,0,0],[0,4,0],[0,0,8]],[[6,-2,-2],[-2,4,0],[-2,0,4]]]
    towers=[{'matrix':A,'layers':exhaustive_tower(A)} for A in mats]
    return {'status':'FINITE-EXACT independent controls','census':reports,'switches':switches,'tower_controls':towers}

def n8_projection(base):
    with gzip.open(base/'nodes.tsv.gz','rt') as f:
        header=f.readline().strip().split('\t');lines=[line.rstrip().split('\t') for line in f if line.strip()]
    with gzip.open(base/'edges.tsv.gz','rt') as f:
        head=f.readline().strip().split('\t');need(head==['source','target','blue','black'],'edge schema');rows=[list(map(int,line.split())) for line in f if line.strip()]
    n=len(lines);blue=[];black=[]
    for u,v,b,k in rows:
        need(k%2==0,'black twin multiplicity')
        if b:blue.append((u,v,b))
        if k:black.append((u,v,k//2))
    p=projection(n,blue,black)
    result={'status':'FINITE-EXACT independent projection from full n8 ledger','n':8,'source_sha256':{name:hashlib.sha256((base/name).read_bytes()).hexdigest() for name in ['nodes.tsv.gz','edges.tsv.gz']},**public(p)}
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n8',type=Path);ap.add_argument('--out',type=Path);args=ap.parse_args();start=time.time()
    result=n8_projection(args.n8) if args.n8 else controls();result['source_sha256']=result.get('source_sha256',{});result['source_sha256']['check_electrical.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=args.out or HERE/('n8_electrical_receipt.json' if args.n8 else 'receipt.json');out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));print('ALL CHECKS PASSED',round(time.time()-start,3),'seconds')
