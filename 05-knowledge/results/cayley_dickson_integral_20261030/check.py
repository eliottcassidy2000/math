"""Exact checks for the scoped Cayley--Dickson integral packet. No assert gates."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from functools import lru_cache
from math import gcd, prod
from pathlib import Path
import argparse
import importlib.util
import json
import random
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
GATES = 0

def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)

def signs(d):
    if d == 0:
        return [[1]]
    old = signs(d-1)
    h = len(old)
    t = [[0]*(2*h) for _ in range(2*h)]
    for i in range(2*h):
        for j in range(2*h):
            u, v = i % h, j % h
            conjugate = 1 if u == 0 else -1
            if i < h:
                t[i][j] = old[u][v] * (1 if j < h else conjugate)
            else:
                t[i][j] = old[v][u] * (1 if j < h else -conjugate)
    return t

def add(x,y): return [a+b for a,b in zip(x,y)]
def neg(x): return [-a for a in x]
def conj(x): return [x[0]]+[-a for a in x[1:]]

def recursive_product(x,y):
    """Independent vector-level doubling, not an XOR-table multiplication."""
    if len(x) == 1:
        return [x[0]*y[0]]
    h = len(x)//2
    a,b,c,d = x[:h],x[h:],y[:h],y[h:]
    return (add(recursive_product(a,c), neg(recursive_product(d,conj(b))))
            + add(recursive_product(conj(a),d), recursive_product(c,b)))

def mul(t,x,y):
    out = [0]*len(t)
    for i,a in enumerate(x):
        for j,b in enumerate(y):
            out[i^j] += a*b*t[i][j]
    return out

def basis(n,i): return [int(i==j) for j in range(n)]
def left(t,weights):
    n = len(t)
    m = [[0]*n for _ in range(n)]
    for i,a in weights.items():
        for j in range(n): m[i^j][j] += a*t[i][j]
    return Matrix(m)

def pair_count(d,i,j):
    if i == j: raise ValueError('distinct pair required')
    if not i or not j: return 0
    h = 1 << (d-1)
    u,v = i % h,j % h
    if (i<h) == (j<h): return 2*pair_count(d-1,u,v)
    if not u or not v or u == v: return 0
    return h-4-2*pair_count(d-1,u,v)

def blocks(t,i,j):
    step=i^j
    out=[]
    for b in range(len(t)):
        c=b^step
        if b<c:
            s=t[j][c]*t[i][b]
            u=t[j][b]*t[i][c]
            out.append((b,c,s,u))
    return out

def histogram(d):
    if d==1: return Counter({0:1})
    old=histogram(d-1);h=1<<(d-1);out=Counter({0:3*h-2})
    for k,v in old.items():
        out[2*k]+=2*v
        cross=v-((h-1) if k==0 else 0)
        if cross: out[h-4-2*k]+=2*cross
    return +out

def valuation(a,p):
    if not a: return None
    k=0
    while a%p==0: a//=p;k+=1
    return k

def canonical_diagonal(ds):
    ds=[abs(int(d)) for d in ds]; nonzero=[d for d in ds if d]
    # Canonical invariant factors by elementary gcd/lcm operations.
    for i in range(len(nonzero)):
        for j in range(i+1,len(nonzero)):
            g=gcd(nonzero[i],nonzero[j])
            nonzero[j]=nonzero[i]//g*nonzero[j];nonzero[i]=g
    return nonzero+[0]*(len(ds)-len(nonzero))

def predicted_diagonal(n,k,a,b):
    g=gcd(a,b)
    if not g: raise ValueError('zero multiplier excluded')
    return canonical_diagonal([g, (a*a-b*b)//g]*k
                              +[g,(a*a+b*b)//g]*(n//2-k))

def rank_mod(m,p):
    a=[[int(v)%p for v in row] for row in m.tolist()];r=0
    for c in range(m.cols):
        pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None: continue
        a[r],a[pivot]=a[pivot],a[r];u=pow(a[r][c],-1,p)
        a[r]=[(v*u)%p for v in a[r]]
        for i in range(r+1,len(a)):
            if a[i][c]:
                u=a[i][c];a[i]=[(v-u*w)%p for v,w in zip(a[i],a[r])]
        r+=1
    return r

def sigma(n):
    out=1;p=2
    while p*p<=n:
        power=1;total=1
        while n%p==0: n//=p;power*=p;total+=power
        out*=total;p+=1
    if n>1: out*=n+1
    return out

def affine_blocks(t,weights):
    i0=min(weights);span={0}
    for i in weights: span |= {u^i^i0 for u in list(span)}
    seen=set();out=[]
    for b in range(len(t)):
        if b in seen: continue
        columns=sorted(b^u for u in span);seen.update(columns)
        rows=sorted(c^i0 for c in columns)
        m=left(t,weights)
        out.append((rows,columns,m.extract(rows,columns)))
    return out

@lru_cache(None)
def sign_step(code,a,b):
    f,g,za,zb=code&1,(code>>1)&1,(code>>2)&1,(code>>3)&1
    forward=(f+b*(1-za) if a==0 else g+b*za)&1
    reverse=(g+a*(1-zb) if b==0 else f+a*zb)&1
    return forward+2*reverse+4*(za and not a)+8*(zb and not b)

def span_add(mask,g):
    return mask | sum(1<<(v^g) for v in range(mask.bit_length()) if (mask>>v)&1)

def frame_histogram(d,r,shape):
    queries=list(product(shape,range(1<<r)))
    states={(tuple([12]*len(queries)),1):1}
    sizes=[]
    for _ in range(d):
        nxt=Counter()
        for (codes,span),count in states.items():
            for origin,base,g in product(range(2),range(2),range(1<<r)):
                updated=tuple(sign_step(code,origin^((s&g).bit_count()%2),base^((u&g).bit_count()%2))
                              for code,(s,u) in zip(codes,queries))
                nxt[updated,span_add(span,g)]+=count
        if len(nxt)==len(states):
            need(set(nxt)==set(states),'full state-set closure, not only equal cardinality')
        states=nxt;sizes.append(len(states))
    out=Counter();full=(1<<(1<<r))-1
    for (codes,span),count in states.items():
        if span==full:
            mask=sum((code&1)<<i for i,code in enumerate(codes));out[mask]+=count
    return out,sizes

def direct_frames(d,r,shape):
    t=signs(d);n=len(t);out=Counter()
    for generators in product(range(n),repeat=r):
        image=[0]*(1<<r)
        for s in range(1<<r):
            for i,g in enumerate(generators):
                if (s>>i)&1: image[s]^=g
        if len(set(image))!=1<<r: continue
        for origin,base in product(range(n),repeat=2):
            mask=0
            for index,(s,u) in enumerate(product(shape,range(1<<r))):
                if t[origin^image[s]][base^image[u]]<0: mask |= 1<<index
            out[mask]+=1
    return out

def minimal_frame_machine(r,shape):
    queries=list(product(shape,range(1<<r)))
    letters=list(product(range(2),range(2),range(1<<r)))
    start=(tuple([12]*len(queries)),1)
    states=[start];index={start:0};edges=[];outputs=[]
    for codes,span in states:
        outputs.append(sum((code&1)<<i for i,code in enumerate(codes))
                       if span==(1<<(1<<r))-1 else -1)
        row=[]
        for origin,base,g in letters:
            updated=tuple(sign_step(code,origin^((s&g).bit_count()%2),base^((u&g).bit_count()%2))
                          for code,(s,u) in zip(codes,queries))
            target=(updated,span_add(span,g))
            if target not in index: index[target]=len(states);states.append(target)
            row.append(index[target])
        edges.append(row)
    need(len(states)==(38 if r==1 else 4534),'BFS independent reachable closure')
    labels={v:i for i,v in enumerate(sorted(set(outputs)))}
    part=[labels[v] for v in outputs];sizes=[len(labels)]
    while True:
        ids={};new=[]
        for i in range(len(states)):
            sig=(part[i],tuple(part[j] for j in edges[i]))
            if sig not in ids: ids[sig]=len(ids)
            new.append(ids[sig])
        sizes.append(len(ids))
        if len(ids)==len(set(part)):
            part=new;break
        part=new
    representatives={}
    for i,c in enumerate(part):
        sig=(outputs[i],tuple(part[j] for j in edges[i]))
        if c in representatives: need(representatives[c]==sig,'minimal quotient transition consistency')
        else: representatives[c]=sig
    # Explicit current-reading collision separated by one further legal bit append.
    groups=defaultdict(list)
    for i,o in enumerate(outputs):
        if o>=0: groups[o].append(i)
    witness=None
    for o,group in groups.items():
        if witness: break
        for i,j in combinations(group,2):
            for a,letter in enumerate(letters):
                oi,oj=outputs[edges[i][a]],outputs[edges[j][a]]
                if oi!=oj:
                    witness={'state_a':[list(states[i][0]),states[i][1]],
                             'state_b':[list(states[j][0]),states[j][1]],
                             'common_reading':o,'next_letter':list(letter),'next_readings':[oi,oj]}
                    break
            if witness: break
    need(witness is not None,'phase/zero guard erasure hostile')
    return {'r':r,'reachable':len(states),'minimal_states':len(set(part)),
            'refinement_sizes':sizes,'reading_collision':witness}

def run():
    rng=random.Random(20261030);census={};pair_total=0;snf_cases=0
    for d in range(1,9):
        t=signs(d);n=len(t);observed=Counter()
        for i,j in combinations(range(n),2):
            bs=blocks(t,i,j);k=sum(s*u==1 for _,_,s,u in bs)
            need(k==pair_count(d,i,j),f'pair recursion {d,i,j}')
            observed[k]+=1;pair_total+=1
        need(observed==histogram(d),'literal complete census')
        if d>=3:
            need(observed[0]==(6*d-5)*2**(d-2)+2,'closed zero count')
            h=n//2
            need(sum(k*v for k,v in observed.items())==(h-4)*(h-1)*(h-2),'first moment')
            need(set(observed)==set(range(0,max(0,h-4)+1,4)),'annihilator dimensions')
        census[str(d)]={str(k):v for k,v in sorted(observed.items())}
        if d<=5:
            for i,j in product(range(n),repeat=2):
                need(recursive_product(basis(n,i),basis(n,j))==
                     [t[i][j]*v for v in basis(n,i^j)],'independent product convention')
            for _ in range(12):
                x=[rng.randrange(-2,3) for _ in range(n)]
                y=[rng.randrange(-2,3) for _ in range(n)]
                need(recursive_product(x,y)==mul(t,x,y),'dense product convention')
        if d<=5:
            pairs=list(combinations(range(n),2))
            chosen=pairs if d<=3 else [next((i,j) for i,j in pairs if pair_count(d,i,j)==k)
                                       for k in observed]+rng.sample(pairs,8)
            for i,j in chosen:
                for a,b in [(1,1),(1,-1),(1,2),(2,3),(2,2),(2,4),(0,3)]:
                    m=left(t,{i:a,j:b});k=pair_count(d,i,j)
                    ds=predicted_diagonal(n,k,a,b)
                    raw=smith_normal_form(m,domain=ZZ)
                    need([abs(int(raw[z,z])) for z in range(n)]==ds,'independent integer Smith form')
                    snf_cases+=1
                    for p in (2,3,5,7):
                        need(rank_mod(m,p)==sum(z%p!=0 for z in ds),'independent modular rank')
                    # Verify the complete transported block, not just its determinant.
                    transported=left(t,{i:1}).T*m
                    for r,c,s,u in blocks(t,i,j):
                        need(transported.extract([r,c],[r,c])==Matrix([[a,b*s],[b*u,a]]),'named phase transport')
                    if d<=2:
                        for modulus in (2,3,4):
                            rows=m.tolist()
                            count=sum(all(sum(int(v)*w for v,w in zip(row,x))%modulus==0 for row in rows)
                                      for x in product(range(modulus),repeat=n))
                            need(count==prod(gcd(modulus,z) for z in ds),'literal composite-modulus kernel')
    # All-size formulas checked well beyond the direct census.
    for d in range(9,15):
        h=1<<(d-1);hist=histogram(d)
        need(sum(hist.values())==(2*h)*(2*h-1)//2,'large recursive total')
        need(hist[0]==(6*d-5)*2**(d-2)+2,'large closed count')
        need(sum(k*v for k,v in hist.items())==(h-4)*(h-1)*(h-2),'large first moment')
    # Independent coefficient aggregation for the dense norm-defect law.
    energies={}
    for d in range(2,7):
        t=signs(d);n=len(t);coeff=defaultdict(int)
        for target in range(n):
            for a,c in combinations(range(n),2):
                b,e=a^target,c^target
                coeff[a,c,min(b,e),max(b,e)]+=2*t[a][b]*t[c][e]
        h=n//2;energy=sum(v*v for v in coeff.values())
        expected=0 if d<3 else 16*(h-4)*(h-1)*(h-2)
        need(energy==expected,'independent norm polynomial energy')
        need(set(coeff.values())<={-4,0,4},'all mixed coefficients')
        energies[str(n)]=str(Fraction(energy,n**4))
    # Sparse compiler: exact rows/columns, including four or more support points.
    compiler_cases=0
    for d in range(3,7):
        t=signs(d);n=len(t)
        for size in (1,2,3,4,5):
            support=rng.sample(range(n),size)
            weights={i:rng.choice([-2,-1,1,2]) for i in support}
            m=left(t,weights);reconstructed=Matrix.zeros(n)
            predicted=[]
            for rows,cols,block in affine_blocks(t,weights):
                for u,r in enumerate(rows):
                    for v,c in enumerate(cols): reconstructed[r,c]=block[u,v]
                if d<=5:
                    s=smith_normal_form(block,domain=ZZ)
                    predicted += [int(s[z,z]) for z in range(s.rows)]
            need(reconstructed==m,'arbitrary-support compiler exact transport')
            if d<=5:
                s=smith_normal_form(m,domain=ZZ)
                need(canonical_diagonal(predicted)==[abs(int(s[z,z])) for z in range(n)],'sparse compiler Smith')
            compiler_cases+=1
    # Explicit hostile: inverse identity is weaker than cancellation.
    t=signs(4);x=add(basis(16,1),basis(16,10))
    r,c,s,u=next(v for v in blocks(t,1,10) if v[2]*v[3]==1)
    y=add(basis(16,r),[-u*v for v in basis(16,c)])
    need(mul(t,x,y)==[0]*16,'explicit zero divisor')
    need(mul(t,conj(x),x)==[2]+[0]*15,'conjugate inverse still exists')
    m=left(t,{1:1,10:2})
    need(rank_mod(m,3)==12 and rank_mod(m,5)==12,'new prime obstruction 3 versus 5')
    m2=left(t,{1:1,10:1})
    need(rank_mod(m2,2)==8 and rank_mod(m2.T*m2,2)==0,'Gram loses binary multiplication')
    # Characteristic-two ring law and Frobenius; named coordinates remain observable.
    for d in range(1,8):
        t=signs(d);n=len(t)
        for _ in range(10):
            x=[rng.randrange(2) for _ in range(n)]
            square=[v%2 for v in mul(t,x,x)]
            need(square==[sum(x)%2]+[0]*(n-1),'binary Frobenius augmentation')
            for i in range(n):
                need(mul(t,basis(n,i),x)[0]%2==x[i],'all named coordinates observable')
    # Original tournament transporter hostile, inherited rather than reconstructed here.
    frame_results=[]
    for r,shape,maxd in [(1,(0,1),8),(2,(0,1,2),7)]:
        for d in range(r,maxd+1):
            hist,sizes=frame_histogram(d,r,shape);n=1<<d
            need(sum(hist.values())==n*n*prod(n-2**j for j in range(r)),'framed transfer universe')
            if d<=4:
                need(hist==direct_frames(d,r,shape),'independent literal joint sign histogram')
            if r==2 and d>=6:
                need(sizes[5]==4534 and len(hist)==3840,'closed rank-two universe')
                for mask in range(4096):
                    forbidden=all((mask&h).bit_count()%2==v for h,v in [(15,1),(240,1),(3840,1),(854,0)])
                    need((mask not in hist)==forbidden,'complete three-support realization syndrome')
            frame_results.append({'d':d,'r':r,'states':sizes[-1],
                                  'sign_patterns':len(hist),'framed_mass':sum(hist.values())})
    minimal_machines=[minimal_frame_machine(1,(0,1)),minimal_frame_machine(2,(0,1,2))]
    arithmetic_bridge=[];primes=(7,11,13);t=signs(3)
    for source,expected_rank in [(1001,8),(343,8),(57,1)]:
        response=[]
        for e in range(8):
            row=[]
            for h in range(8):
                n=source*prod(p**(((e>>i)&1)+((h>>i)&1)) for i,p in enumerate(primes))
                row.append((sigma(n)-n)%2)
            response.append(row)
        H=Matrix(response);need(rank_mod(H,2)==expected_rank,'literal source response rank')
        for i in range(8):
            L=left(t,{i:1}).applyfunc(lambda x:int(x)%2)
            need((H*L-L*H).applyfunc(lambda x:int(x)%2)==Matrix.zeros(8),'binary arithmetic/algebra intertwiner')
        arithmetic_bridge.append({'source':source,'rank':expected_rank})
    spec=importlib.util.spec_from_file_location('hostile_runtime',ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py')
    hostile=importlib.util.module_from_spec(spec);spec.loader.exec_module(hostile)
    control=hostile.hostile_checks()
    return {'scope':'all-pair census through dimension 256; all-size proofs in PROOFS.md',
            'pair_cases':pair_total,'smith_cases':snf_cases,'compiler_cases':compiler_cases,
            'frame_transfers':frame_results,
            'minimal_machines':minimal_machines,
            'arithmetic_bridge':arithmetic_bridge,
            'census':census,'normalized_norm_variances':energies,
            'zero_divisor':{'x':{'1':1,'10':1},'y':{str(i):v for i,v in enumerate(y) if v}},
            'hostile_81_27':control,'inherited_gates':hostile.GATES,'local_gates':GATES}

def main():
    p=argparse.ArgumentParser();p.add_argument('--record',action='store_true');args=p.parse_args()
    result=run();path=HERE/'receipt.json'
    if args.record: path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    else: need(result==json.loads(path.read_text(encoding='utf-8')),'recorded receipt')
    print(json.dumps({k:result[k] for k in ('pair_cases','smith_cases','compiler_cases','local_gates','inherited_gates')}))

if __name__=='__main__': main()
