"""Independent audit: subset support DP, set reachability, no producer imports."""
from pathlib import Path
from itertools import combinations,permutations
import json,random,ast

def demand(p,msg):
    if not p:raise RuntimeError(msg)

def encode(n,mask,fixed=False):
    rows=[0]*n;pos=0
    for a,b in combinations(range(n),2):
        if fixed and b==a+1:forward=True
        else:forward=bool(mask&1<<pos);pos+=1
        if forward:rows[a]+=2**b
        else:rows[b]+=2**a
    return tuple(rows)

def trace_dp(t):
    n=len(t);full=2**n-1;ends=[0]*(full+1);starts=ends.copy()
    incoming=[sum(1<<i for i in range(n) if t[i]>>j&1) for j in range(n)]
    for mask in range(1,full+1):
        bits=mask
        while bits:
            bit=bits&-bits;bits-=bit;v=bit.bit_length()-1;rest=mask^bit
            if not rest or ends[rest]&incoming[v]:ends[mask]|=bit
            if not rest or starts[rest]&t[v]:starts[mask]|=bit
    d=[0]*n
    for mask in range(1,full):
        left=ends[mask];right=starts[full^mask]
        while left:
            bit=left&-left;left-=bit;v=bit.bit_length()-1
            d[v]|=t[v]&right
    return tuple(d)

def closure(t,vertices=None):
    vs=set(range(len(t))) if vertices is None else set(vertices)
    reach={u:{u} for u in vs}
    for u in vs:
        todo=[u]
        while todo:
            v=todo.pop()
            for w in vs-reach[u]:
                if t[v]>>w&1:reach[u].add(w);todo.append(w)
    return reach

def decode(d):
    n=len(d);reach=closure(d);answer=[0]*n
    for u,v in combinations(range(n),2):
        if d[u]>>v&1:a,b=u,v
        elif d[v]>>u&1:a,b=v,u
        elif v not in reach[u] or u not in reach[v]:
            demand((v in reach[u])!=(u in reach[v]),'component order')
            a,b=(u,v) if v in reach[u] else (v,u)
        else:
            forward=[w for w in range(n) if d[u]>>w&1 and d[w]>>v&1]
            backward=[w for w in range(n) if d[v]>>w&1 and d[w]>>u&1]
            demand((len(forward)>1)!=(len(backward)>1),'two-route uniqueness')
            a,b=(u,v) if len(forward)>1 else (v,u)
        answer[a]|=1<<b
    return tuple(answer)

def check(t):
    d=trace_dp(t);demand(decode(d)==t,'independent reconstruction')
    return sum(r.bit_count() for r in t)-sum(r.bit_count() for r in d)

def upset(n):
    # Transitive i->j for i<j, reverse the path n-1,n-3,...,0 into increasing arrows.
    # Equivalently start with all larger->smaller and reverse consecutive even slots.
    t=[(1<<i)-1 for i in range(n)]
    for i in range(0,n-2,2):t[i]|=1<<(i+2);t[i+2]^=1<<i
    return tuple(t)

def cycle_substitution(blocks):
    starts=[0]
    for b in blocks:starts.append(starts[-1]+len(b))
    rows=[]
    for i,b in enumerate(blocks):
        j=(i+1)%3;ext=((1<<len(blocks[j]))-1)<<starts[j]
        rows.extend((r<<starts[i])|ext for r in b)
    return tuple(rows)

def hcount(t):
    n=len(t);full=(1<<n)-1;dp=[[0]*n for _ in range(full+1)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,full+1):
        for v,c in enumerate(dp[mask]):
            if not c:continue
            bits=t[v]&(full^mask)
            while bits:
                b=bits&-bits;bits-=b;dp[mask|b][b.bit_length()-1]+=c
    return sum(dp[full])

def main():
    labelled=[]
    for n in range(1,6):
        total=1<<(n*(n-1)//2);missing=0
        for mask in range(total):missing+=check(encode(n,mask))
        labelled.append({'n':n,'states':total,'missing':missing})
    fixed_missing=sum(check(encode(6,m,True)) for m in range(1024))
    root=Path(__file__).resolve().parent.parent
    atlas=json.loads((root/'strong_atom_continuation_vitali_20261026/strong-atlas.json').read_text(encoding='utf-8'))
    sevens=[a for a in atlas if a['n']==7]
    demand(len(sevens)==353,'atlas size')
    seven_missing=sum(check(tuple(a['source'])) for a in sevens)
    large=[]
    for n in range(5,18,2):
        t=upset(n);missing=check(t)
        demand(missing==(n*n-4*n+3)//8,'upset invisible count')
        large.append({'n':n,'missing':missing})
    rng=random.Random(20261026729);samples=[]
    for n in range(8,13):
        for _ in range(20):
            t=encode(n,rng.getrandbits(n*(n-1)//2));check(t)
        samples.append({'n':n,'states':20})
    plus=(10,12,9,0);minus=(4,1,2,7)
    twins=[cycle_substitution([plus]*3),cycle_substitution([minus,plus,plus])]
    triples=[];counts=[]
    for t in twins:
        check(t);counts.append(hcount(t))
        triples.append({s for s in combinations(range(12),3) if all(sum(t[v]>>w&1 for w in s)==1 for v in s)})
    demand(triples[0]==triples[1] and len(triples[0])==67,'81/27 labelled triple control')
    demand(counts==[462915,462915],'81/27 Hamilton counts')
    # Opposite guard outcomes are witnessed by a path using an arc different in the twins.
    demand(trace_dp(twins[0])!=trace_dp(twins[1]),'81/27 trace supports distinct')
    source=root.parent.parent/'04-computation/vitali_holonomy_composition.py'
    node=next(n for n in ast.parse(source.read_text(encoding='utf-8')).body if isinstance(n,ast.FunctionDef) and n.name=='count_directed_k_cycles')
    env={'combinations':combinations,'permutations':permutations}
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),env)
    old=[]
    for mask in [401896,440806,828520,801894]:
        t=encode(7,mask);A=[[t[i]>>j&1 for j in range(7)] for i in range(7)]
        by_k=[]
        for k in range(3,8):
            directed=sum(all(t[p[i]]>>p[(i+1)%k]&1 for i in range(k)) for p in permutations(range(7),k))//k
            demand(env['count_directed_k_cycles'](A,7,k)==directed,'legacy routine correction')
            by_k.append(directed)
        old.append({'mask':mask,'H':hcount(t),'simple_cycles_3_to_7':by_k})
    pe={(0,1),(0,4),(0,5),(1,2),(1,6),(2,3),(2,7),(3,4),(3,8),(4,9),(5,7),(5,8),(6,8),(6,9),(7,9)}
    pm=[frozenset(s) for s in combinations(sorted(pe),5) if len({v for e in s for v in e})==10]
    demand(len(pm)==6 and all(len(a&b)==1 for a,b in combinations(pm,2)),'independent Petersen matching census')
    print(json.dumps({'all_labelled':labelled,'fixed_path_six':{'states':1024,'missing':fixed_missing},'strong_seven':{'classes':353,'missing':seven_missing},'upset_family':large,'random_seed':20261026729,'random_samples':samples,'81_27_control':{'H':counts,'same_labelled_cyclic_triples':67,'trace_supports_differ':True},'legacy_cycle_correction':old,'Petersen_matchings':len(pm)},sort_keys=True,indent=2))

if __name__=='__main__':main()
