#!/usr/bin/env python3
"""Standalone exact checks; no removable assertions or external dependencies."""
from itertools import product
from collections import Counter
import json

def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)

def graph(bits, rules=(), fixed=0):
    m=len(bits); k=len(rules); n=m+2+k
    r=[0]*n; r[0]=2
    for i,b in enumerate(bits,2):
        r[i]|=1
        if b: r[1]|=1<<i
        else: r[i]|=2
        for j in range(2,i): r[j]|=1<<i
    for j,rule in enumerate(rules):
        z=m+2+j
        for v in [0,1]+list(range(m+2,z)):
            if fixed&1: r[z]|=1<<v
            else: r[v]|=1<<z
            fixed>>=1
        for i,b in enumerate(bits,2):
            if rule>>b&1: r[z]|=1<<i
            else: r[i]|=1<<z
    return r

def induced(rows, vs):
    return [sum(1<<j for j,v in enumerate(vs) if rows[u]>>v&1) for u in vs]

def hp(rows):
    n=len(rows)
    if not n:return 1
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            a=d[mask][v]
            if not a:continue
            can=rows[v]&~mask
            while can:
                bit=can&-can;can-=bit; w=bit.bit_length()-1
                d[mask|bit][w]+=a
    return sum(d[-1])

def max_module_runs(rows, module):
    n=len(rows); M=set(module)
    d=[[-1000]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=int(i in M)
    for mask in range(1,1<<n):
        for v in range(n):
            a=d[mask][v]
            if a<0:continue
            can=rows[v]&~mask
            while can:
                bit=can&-can;can-=bit;w=bit.bit_length()-1
                val=a+int(w in M and v not in M)
                d[mask|bit][w]=max(d[mask|bit][w],val)
    return max(d[-1])

def g(u):
    if u==1:return 1
    if u%2==0:return g(u//2)
    v=(u+1)//2
    return g(v)+v-v%2

def g_bits(bits):
    m=len(bits)
    return bits[0]+1-bits[-1]+sum((1-bits[i])*bits[i+1] for i in range(m-1))+sum(bits[j]*(1-bits[i])*2**(i-j-1) for j in range(m) for i in range(j+1,m))

def bitword(u):return [int(c) for c in bin(u-1)[2:]]

def insert_path(rows,path):
    for v in range(len(rows)):
        if v in path:continue
        if rows[v]>>path[0]&1:path.insert(0,v)
        elif rows[path[-1]]>>v&1:path.append(v)
        else:
            for i in range(len(path)-1):
                if rows[path[i]]>>v&1 and rows[v]>>path[i+1]&1:
                    path.insert(i+1,v);break
            else:raise RuntimeError('insertion failed')
    require(len(set(path))==len(rows),'Hamilton permutation')
    require(all(rows[a]>>b&1 for a,b in zip(path,path[1:])),'Hamilton edges')
    return path

def check_witness(rows,module,edge1,edge2):
    a,b=edge1;c,d=edge2
    require(len({a,b,c,d})==4,'disjoint edges')
    require(rows[a]>>b&1 and rows[c]>>d&1,'cross arcs')
    for x in module:
        require(rows[x]>>a&1 and rows[x]>>c&1 and rows[b]>>x&1 and rows[d]>>x&1,'cross partition')
    path=insert_path(rows,[module[0],a,b,module[1],c,d,module[2]])
    positions=[path.index(x) for x in module]
    require(all(y-x>1 for x,y in zip(positions,positions[1:])),'separated module visits')
    require(max_module_runs(rows,module)>=3,'independent max-run DP')
    require(hp(rows)>=3**3-3*2**3+3,'onto count bound')

def classify(rules,fixed):
    # Construct the exact context in the all-order proof, with pump length 3.
    for j,rule in enumerate(rules):
        if rule==1: # eta(0)=1, eta(1)=0
            bits=[1,0,1,1,1]; rows=graph(bits,rules,fixed); z=7+j
            check_witness(rows,[4,5,6],(0,1),(z,3))
            return 'nonuniform opposite'
        if rule==2:
            bits=[1,1,1,0];rows=graph(bits,rules,fixed);z=6+j
            check_witness(rows,[2,3,4],(0,1),(5,z))
            return 'nonuniform copied'
    # Every attachment is now uniform. Use context 1,0,0,0.
    bits=[1,0,0,0];rows=graph(bits,rules,fixed)
    A=[6+j for j,r in enumerate(rules) if r==0]
    B=[6+j for j,r in enumerate(rules) if r==3]
    for a in [0]+A:
        for b in B:
            if rows[a]>>b&1:
                check_witness(rows,[3,4,5],(1,2),(a,b))
                return 'uniform forbidden cross'
    # Surviving structure and independent quadrilateral checks.
    base=graph([1],rules,fixed)
    Av=[3+j for j,r in enumerate(rules) if r==0]
    Bv=[3+j for j,r in enumerate(rules) if r==3]
    kap=hp(induced(base,Bv))*hp(induced(base,[0]+Av))
    for k in range(2,6):
        q=2**k
        hs=[hp(graph(bitword(q+j),rules,fixed)) for j in range(1,5)]
        require(hs[3]-hs[2]-hs[1]+hs[0]==-2*kap,('quadrilateral',rules,fixed,k,hs,kap))
    return 'uniform ordered survivor'

def main():
    counts={};total=0
    for k in range(3):
        c=Counter()
        for rules in product(range(4),repeat=k):
            for fixed in range(1<<(2*k+k*(k-1)//2)):
                c[classify(rules,fixed)]+=1;total+=1
        counts[str(k)]=dict(sorted(c.items()))
    source_checks=0
    for m in range(1,8):
        for tail in product(range(2),repeat=m-1):
            bits=(1,)+tail;u=1+sum(b*2**(m-1-i) for i,b in enumerate(bits))
            rows=graph(bits);W=induced(rows,list(range(1,m+2)))
            require(hp(W)==g(u)==g_bits(bits),('source g',bits,u))
            require(hp(rows)==u+g(u)-u%2,('source H',bits,u))
            # Independently enumerate labelled two-bin covers ending at c.
            ell=0;n=len(W)
            for mask in range(1<<n):
                if not mask&1:continue
                left=[i for i in range(n) if mask>>i&1]
                right=[i for i in range(n) if not mask>>i&1]
                # Reserved path ends at c: private chain followed by c.
                pp=[i for i in left if i!=0]
                ok=not pp or bool(W[max(pp)]&1)
                if ok:ell+=hp(induced(W,right))
            require(ell==1+sum((1-b)*2**i for i,b in enumerate(bits)),('L source',bits,ell))
            source_checks+=1
    print(json.dumps({'status':'PASS','labelled_grafts_checked':total,'grafts_by_added_order':counts,'source_bitwords_checked':source_checks,'quadrilaterals_per_survivor':4,'normalization':'r=0,c=1,private vertices in binary order','theorem_scope':'fixed finite letter-local extensions retaining E(u), arbitrary eventual affine slope'},sort_keys=True,indent=2))
if __name__=='__main__':main()
