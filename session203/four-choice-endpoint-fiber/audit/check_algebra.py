"""Independent direct two-path core partitions, operators, and labeled-bin eigenvalues."""
from itertools import permutations,product
from collections import Counter
from pathlib import Path
from functools import lru_cache
import json
Z=(0,0,0,0);ONE={Z:1};X=[{tuple(int(i==j) for i in range(4)):1} for j in range(4)]
A=[102,0,10,83,135,30,182,47];V=5;SEEDS=(10,14,17,21)

def need(b,s):
    if not b:raise RuntimeError(s)

def add(*terms):
    c=Counter()
    for k,p in terms:
        for m,v in p.items():c[m]+=k*v
    return {m:v for m,v in c.items() if v}

def mul(p,q):
    c=Counter()
    for m,v in p.items():
        for n,w in q.items():c[tuple(a+b for a,b in zip(m,n))]+=v*w
    return {m:v for m,v in c.items() if v}

def core_to_body(v,color):return bool(color) if v==3 else bool(A[v]>>V&1)

def path_poly(word,root_first=False,root_last=False):
    if not word:return add((1,ONE),*((1,x) for x in X))
    if root_first and word[0]!=0:return {}
    if root_last and word[-1]!=0:return {}
    p=ONE
    for gap in range(len(word)+1):
        l=word[gap-1] if gap else None;r=word[gap] if gap<len(word) else None
        options=ONE if l is None or r is None or A[l]>>r&1 else {}
        if not ((root_first and gap==0) or (root_last and gap==len(word))):
            for f,t in product((0,1),repeat=2):
                if (l is None or core_to_body(l,f)) and (r is None or not core_to_body(r,t)):
                    options=add((1,options),(1,X[2*f+t]))
        p=mul(p,options)
        if not p:break
    return p

def tensors():
    names=('h','s','e','L','R','d','beta');out={n:{} for n in names};orders=cuts=0
    for word in permutations([i for i in range(8) if i!=V]):
        orders+=1;r=word.index(0)
        for name,rf,rl in [('h',False,False),('s',True,False),('e',False,True)]:out[name]=add((1,out[name]),(1,path_poly(word,rf,rl)))
        out['R']=add((1,out['R']),(1,mul(path_poly(word[:r]),path_poly(word[r:],True))))
        out['L']=add((1,out['L']),(1,mul(path_poly(word[:r+1],False,True),path_poly(word[r+1:]))));cuts+=2
    for word in permutations([i for i in range(8) if i not in (0,V)]):
        orders+=1;out['d']=add((1,out['d']),(1,path_poly(word)))
        for k in range(len(word)+1):out['beta']=add((1,out['beta']),(1,mul(path_poly(word[:k]),path_poly(word[k:]))));cuts+=1
    out['B']=add((1,out['h']),(1,out['L']),(-1,out['s']));out['V']=add((1,out['R']),(1,out['beta']))
    return out,orders,cuts

def derivative(p,i):
    c={}
    for m,v in p.items():
        if m[i]:
            n=list(m);n[i]-=1;c[tuple(n)]=v*m[i]
    return c

def D(p):return add((1,mul(add((1,X[0]),(1,X[2])),derivative(p,2))),(1,mul(add((1,ONE),(1,X[1]),(1,X[3])),derivative(p,3))))
def E(p):return add((1,mul(add((1,ONE),(1,X[0]),(1,X[1])),derivative(p,0))),(1,mul(add((1,X[2]),(1,X[3])),derivative(p,2))))
SUBS=[add((1,X[1]),(-1,ONE),(-1,X[3])),X[3],add((1,X[0]),(-1,X[1]),(-1,X[2]),(1,ONE),(1,X[3])),add((1,X[2]),(-1,ONE),(-1,X[3]))]
@lru_cache(None)
def transform_monomial(m):
    p=ONE
    for i,n in enumerate(m):
        for _ in range(n):p=mul(p,SUBS[i])
    return p

def transform(p):return add(*((v,transform_monomial(m)) for m,v in p.items()))
def diag_D(p):return {m:v*(m[0]+m[2]) for m,v in p.items() if m[0]+m[2]}
def diag_E(p):return {m:v*(m[0]+m[1]) for m,v in p.items() if m[0]+m[1]}

@lru_cache(None)
def value(m,word,n,diagonal):
    slots=[i for i,power in enumerate(m) for _ in range(power)]
    if not slots:return int(n==0)
    total=0
    for assignment in product(range(len(slots)),repeat=n):
        masks=[0]*len(slots)
        for v,slot in enumerate(assignment):masks[slot]|=1<<v
        ok=True
        for kind,mask in zip(slots,masks):
            if mask:
                f=(word>>((mask&-mask).bit_length()-1))&1;t=(word>>(mask.bit_length()-1))&1
                yes=(kind==0 or (kind==1 and f==0) or (kind==2 and t==1) or (kind==3 and f==0 and t==1)) if diagonal else (kind==2*f+t)
            else:yes=(diagonal and kind!=3)
            if not yes:ok=False;break
        total+=ok
    return total

def eval_poly(p,word,n,diagonal=False):return sum(v*value(m,word,n,diagonal) for m,v in p.items())
def encode(p):return [[list(m),v] for m,v in sorted(p.items())]

def body_functions(rows,colors):
    n=len(rows);f=[[0]*(1<<n) for _ in range(4)]
    def dfs(first,last,mask):
        f[2*colors[first]+colors[last]][mask]+=1
        todo=rows[last]&~mask
        while todo:
            bit=todo&-todo;todo-=bit;dfs(first,bit.bit_length()-1,mask|bit)
    for v in range(n):dfs(v,v,1<<v)
    return f

def body_values(rows,colors,basis):
    n=len(rows);f=body_functions(rows,colors);values={}
    for m in basis:
        slots=[i for i,c in enumerate(m) for _ in range(c)]
        if not slots:values[m]=int(n==0);continue
        total=0
        for assignment in product(range(len(slots)),repeat=n):
            masks=[0]*len(slots)
            for v,j in enumerate(assignment):masks[j]|=1<<v
            count=1
            for kind,mask in zip(slots,masks):count*=f[kind][mask]
            total+=count
        values[m]=total
    return values

def arbitrary_operator_controls(basis):
    from itertools import combinations
    cases=checks=0
    for n in range(4):
        pairs=list(combinations(range(n),2))
        for code in range(1<<len(pairs)):
            a=[0]*n
            for bit,(i,j) in enumerate(pairs):
                if code>>bit&1:a[i]|=1<<j
                else:a[j]|=1<<i
            for colors in product((0,1),repeat=n):
                old=body_values(a,colors,basis)
                source=a+[(1<<n)-1];sink=[row|(1<<n) for row in a]+[0]
                src=body_values(source,colors+(1,),basis);snk=body_values(sink,colors+(0,),basis)
                for m in basis:
                    p={m:1}
                    need(sum(c*old[k] for k,c in D(p).items())==src[m],'arbitrary-body source operator')
                    need(sum(c*old[k] for k,c in E(p).items())==snk[m],'arbitrary-body sink operator')
                    checks+=2
                cases+=1
    return cases,checks


def main():
    producer=json.loads(Path(__file__).with_name('PRODUCER_TENSORS.json').read_text());raw,orders,cuts=tensors();report={}
    basis=[m for m in product(range(5),repeat=4) if sum(m)<=4];operator_checks=0
    for m in basis:
        p={m:1};q=transform(p)
        need(transform(D(p))==diag_D(q) and transform(E(p))==diag_E(q),'diagonal generator action')
        need(D(E(p))==E(D(p)),'commuting physical additions')
        for n in range(4):
            for word in range(1<<n):
                need(eval_poly(D(p),word,n)==value(m,(word<<1)|1,n+1,False),'literal source addition')
                need(eval_poly(E(p),word,n)==value(m,word,n+1,False),'literal sink addition')
                operator_checks+=2
    for name,p in raw.items():
        archived=producer[name]
        need(encode(p)==archived['tensor'],('original tensor',name))
        diag=transform(p);need(encode(diag)==archived['diagonal_tensor'],('full diagonal tensor',name))
        values={}
        for seed in SEEDS:
            grouped=Counter()
            for m,c in diag.items():grouped[(m[0]+m[2])*(m[0]+m[1])]+=c*value(m,seed,5,True)
            grouped={k:v for k,v in grouped.items() if v};values[seed]=grouped
            need(grouped=={int(k):v for k,v in archived['eigenvalues'][str(seed)].items()},('all eigencomponents',name,seed))
            need(sum(grouped.values())==eval_poly(p,seed,5),'zero-wrap sum, including lambda zero')
        for seed,offset in zip(SEEDS,(0,8,14,22)):
            delta={k:values[seed].get(k,0)-values[10].get(k,0) for k in set(values[seed])|set(values[10])};delta={k:v for k,v in delta.items() if v}
            need(delta==({2:offset} if offset and name in ('h','s') else {}),('exact four offsets',name,seed))
        report[name]={'tensor':encode(p),'diagonal_tensor':encode(diag),'eigenvalues':values}
    table={}
    for line in Path(__file__).with_name('REVIEWED_PROOF.md').read_text().splitlines():
        if line.startswith('| '):
            cells=[cell.strip() for cell in line.strip('|').split('|')]
            if cells[0] in raw:table[cells[0]]=[int(v) for v in cells[1:]]
    need(set(table)==set(raw),'all displayed coefficient rows')
    for name,coefficients in table.items():
        observed={lam:c for lam,c in zip((0,1,2,4,6,9,12,16),coefficients) if c}
        need(observed==report[name]['eigenvalues'][10],('displayed coefficient table',name))
    need(raw['e']=={},'retained sink forces zero endings')
    arbitrary_cases,arbitrary_checks=arbitrary_operator_controls(basis)
    print(json.dumps({'status':'PASS','core_orders':orders,'core_partition_cuts':cuts,'monomial_basis_through_degree4':len(basis),'literal_source_sink_operator_checks':operator_checks,'arbitrary_colored_operator_inputs_through3':arbitrary_cases,'arbitrary_body_source_sink_checks':arbitrary_checks,'displayed_coefficient_table_rows':len(table),'profiles':report,'zero_eigenvalue_convention':'lambda^0=1, including lambda=0; zero for positive wrap count'},indent=2,sort_keys=True))
if __name__=='__main__':main()
