#!/usr/bin/env python3
"""Independent finite controls for signed primitive classification.
No imported producer code. Tests labelled tournaments through order five,
then canonical degree-seven and imprimitive degree-nine examples.
"""
import argparse, hashlib, itertools, json
from pathlib import Path

def require(test, message):
    if not test:
        raise RuntimeError(message)

def adj_from_mask(n, mask):
    a=[[False]*n for _ in range(n)]
    for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
        a[i][j]=bool(mask>>k&1); a[j][i]=not a[i][j]
    return a

def cyclic(n, steps):
    return [[i!=j and (j-i)%n in steps for j in range(n)] for i in range(n)]

def signed_groups(a):
    n=len(a); degree=[sum(row) for row in a]; pos=[]; neg=[]
    pairs=list(itertools.combinations(range(n),2))
    for p in itertools.permutations(range(n)):
        if all(degree[p[i]]==degree[i] for i in range(n)) and all(a[p[i]][p[j]]==a[i][j] for i,j in pairs):
            pos.append(p)
        if n>=2 and all(degree[p[i]]==n-1-degree[i] for i in range(n)) and all(a[p[i]][p[j]]!=a[i][j] for i,j in pairs):
            neg.append(p)
    return pos,neg

def image(s,p):
    return sum(1<<p[i] for i in range(len(p)) if s>>i&1)

def regular_subsets(group,n):
    identity=tuple(range(n))
    return [s for s in range(1<<n) if all(p==identity or image(s,p)!=s for p in group)]

def transitive(group,n):
    return n>0 and len({p[0] for p in group})==n

def primitive(group,n):
    if not transitive(group,n): return False
    for s in range(1,1<<n):
        if not s&1 or s.bit_count() in (1,n): continue
        if all((t:=image(s,p))==s or not t&s for p in group): return False
    return True

def modules(a):
    n=len(a); found=[]
    for s in range(1,1<<n):
        if s.bit_count() in (1,n): continue
        if all(len({a[x][y] for y in range(n) if s>>y&1})<=1 for x in range(n) if not s>>x&1):
            found.append(s)
    return found

def strongly_connected(a):
    n=len(a)
    for start in range(n):
        seen={start}; todo=[start]
        while todo:
            x=todo.pop()
            for y in range(n):
                if a[x][y] and y not in seen: seen.add(y); todo.append(y)
        if len(seen)!=n: return False
    return True

def hamilton_data(a, g, neg):
    n=len(a)
    paths=[p for p in itertools.permutations(range(n)) if all(a[p[i]][p[i+1]] for i in range(n-1))]
    remaining=set(paths); representatives=[]
    while remaining:
        p=min(remaining); representatives.append(p)
        remaining.difference_update(tuple(t[v] for v in p) for t in g)
    representative_set=set(representatives)
    fixed_orbits=None
    if neg:
        tau=neg[0]
        fixed_orbits=0
        for p in representatives:
            q=tuple(tau[v] for v in reversed(p))
            qrep=min(tuple(t[v] for v in q) for t in g)
            require(qrep in representative_set,'Converse failed to map Hamilton orbits')
            fixed_orbits+=qrep==p
    # Independent subset dynamic programming count.
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1<<n):
        for i in range(n):
            if dp[s][i]:
                for j in range(n):
                    if a[i][j] and not s>>j&1: dp[s|1<<j][j]+=dp[s][i]
    require(sum(dp[-1])==len(paths),'Hamilton count implementations disagree')
    require(len(paths)==len(g)*len(representatives),'Ordinary Hamilton action not free')
    return {'ordered_hamilton_paths':len(paths), 'ordinary_hamilton_orbits':len(representatives),
            'converse_fixed_hamilton_orbits':fixed_orbits,
            'signed_hamilton_orbits':(len(representatives)+fixed_orbits)//2 if neg else len(representatives)}

def data(a):
    n=len(a); g,neg=signed_groups(a); e=g+neg
    rg=regular_subsets(g,n); re=regular_subsets(e,n)
    return {'n':n, 'ordinary_order':len(g), 'signed_order':len(e),
            'ordinary_regular_subsets':len(rg), 'signed_regular_subsets':len(re),
            'ordinary_regular_orbits':len(rg)//len(g),
            'signed_regular_orbits':len(re)//len(e),
            'ordinary_primitive':primitive(g,n), 'signed_primitive':primitive(e,n),
            'module_masks':modules(a), 'strong':strongly_connected(a),
            'outneighbors':[[j for j,v in enumerate(row) if v] for row in a],
            'hamilton_data':hamilton_data(a,g,neg)}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',default='receipt.json'); args=parser.parse_args()
    census=[]; indecomp_nonprimitive=None
    for n in range(2,6):
        row={'n':n,'labelled_tournaments':0,'signed_primitive':0,'exceptional':0,'exceptional_primitive':0,'non_strong':0}
        for mask in range(1<<(n*(n-1)//2)):
            a=adj_from_mask(n,mask); g,neg=signed_groups(a); e=g+neg
            prim=primitive(e,n); pg=primitive(g,n); rg=regular_subsets(g,n); re=regular_subsets(e,n)
            require(bool(rg),'Gluck ordinary binary failure')
            require(n==2 or pg==prim,'Primitive action equivalence failure')
            require(not prim or neg,'Primitive tournament not self-converse')
            require(not prim or re or (n in (3,5) and transitive(g,n)),'Unexpected primitive exception')
            strong=strongly_connected(a)
            require(strong or re,'Non-strong signed exception')
            require(transitive(g,n) or re,'Intransitive ordinary group yielded signed exception')
            if n==5 and not modules(a) and not pg and indecomp_nonprimitive is None:
                indecomp_nonprimitive={'mask':mask,**data(a)}
            row['labelled_tournaments']+=1; row['signed_primitive']+=prim
            row['exceptional']+=not re; row['exceptional_primitive']+=prim and not re
            row['non_strong']+=not strong
        census.append(row)
    named={
        'C3':data(cyclic(3,{1})),
        'C5_12':data(cyclic(5,{1,2})),
        'P7':data(cyclic(7,{1,2,4})),
        'C7_123':data(cyclic(7,{1,2,3})),
    }
    for key in ('C3','C5_12','P7'):
        require(named[key]['signed_regular_subsets']==0,key+' should be exceptional')
    require(named['C7_123']['signed_regular_subsets']>0,'Cyclic seven control should not be exceptional')
    # C3[C3] uses vertex (outer,inner) encoded 3*outer+inner.
    a=[[False]*9 for _ in range(9)]
    for i in range(9):
        for j in range(9):
            if i!=j:
                a[i][j]=((j//3-i//3)%3==1) if i//3!=j//3 else ((j%3-i%3)%3==1)
    named['C3_of_C3']=data(a)
    require(named['C3_of_C3']['ordinary_order']==81,'Incorrect wreath group order')
    require(named['C3_of_C3']['signed_order']==162,'Incorrect signed wreath group order')
    require(named['C3_of_C3']['signed_regular_subsets']==0,'Triangle tower should be exceptional')
    require(not named['C3_of_C3']['signed_primitive'],'Triangle tower should be imprimitive')
    receipt={'status':'PASS','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'census':census,'named_examples':named,
             'indecomposable_but_not_permutation_primitive':indecomp_nonprimitive,
             'scope':'All labelled tournaments 2<=n<=5, plus four named odd-order controls and C3[C3]. No all-order completeness claim.'}
    Path(args.output).write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt,indent=2,sort_keys=True))

if __name__=='__main__': main()
