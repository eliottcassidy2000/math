"""Independent degree-permutation canonicalization and literal adjacency checks.
Uses no assertions (the same tests remain active under python -O).
"""
from collections import Counter, defaultdict
from itertools import permutations, product
from pathlib import Path
import csv,gzip,json,random,sys,time
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'session67_metagraph_enumeration_information_20261026/session67-root/output'
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def canon(rows):
    buckets=defaultdict(list)
    for i,r in enumerate(rows):buckets[r.bit_count()].append(i)
    best=None
    for choices in product(*(permutations(vs) for d,vs in sorted(buckets.items()))):
        p=sum(choices,());code=0;bit=1
        for i,u in enumerate(p):
            for v in p[i+1:]:
                if rows[u]>>v&1:code|=bit
                bit<<=1
        if best is None or code<best:best=code
    return best
def decode(n,x):
    rows=[0]*n;b=0
    for i in range(n):
        for j in range(i+1,n):
            q=True if j==i+1 else bool(x>>b&1)
            if j>i+1:b+=1
            if q:rows[i]|=1<<j
            else:rows[j]|=1<<i
    return rows
def converse(rows):return [((1<<len(rows))-1)^r^(1<<i) for i,r in enumerate(rows)]
def H(rows):
    n=len(rows)
    if not n:return 1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for S in range(1,1<<n):
        for v,t in enumerate(dp[S]):
            if not t:continue
            q=rows[v]&~S
            while q:
                b=q&-q;q-=b;dp[S|b][b.bit_length()-1]+=t
    return sum(dp[-1])
def child(rows,S):
    n=len(rows);return [r|((0 if S>>i&1 else 1)<<n) for i,r in enumerate(rows)]+[S]
def tail(rows,S,t):
    n=len(rows);v=n+t;out=list(rows)+[0]*(t+1)
    for i in range(n):
        for z in range(n,v):out[i]|=1<<z
        if S>>i&1:out[v]|=1<<i
        else:out[i]|=1<<v
    for z in range(n,v):
        for w in range(z+1,v+1):out[z]|=1<<w
    return out
rs=[json.loads(s) for s in gzip.open(HERE/'records.jsonl.gz','rt')]
receipt={'canonical_frames':0,'literal_children':0,'literal_tails':0,'relabelings':0,'line_rows':0,'orders':[]}
maxn=int(sys.argv[1]) if len(sys.argv)>1 else 7
for n in range(1,maxn+1):
    group=[r for r in rs if r['n']==n];M=(n-1)*(n-2)//2
    canonical={canon(r['rows']):r for r in group};need(len(canonical)==len(group),'inherited class distinctness')
    words=[];count=Counter()
    for x in range(1<<M):
        c=canon(decode(n,x));need(c in canonical,'full frame coverage');words.append(c);count[c]+=1
    need(all(count[c]==r['tc'] for c,r in canonical.items()),'presentation fibre exact')
    for c,r in canonical.items():need((canon(converse(r['rows']))==c)==bool(r['sc']),'converse flag')
    actual=Counter();edges=[(i,j) for i in range(n) for j in range(i+2,n)];index={p:i for i,p in enumerate(edges)}
    for x in range((1<<M)//2 if M else 0):
        y=x^((1<<M)-1);sx=0
        for b,(i,j) in enumerate(edges):sx|=((x>>b)&1)<<index[n-1-j,n-1-i]
        a=canonical[words[x]]['merged'];b=canonical[words[y]]['merged']
        actual[min(a,b),max(a,b),int(sx==x)]+=1
    expected=Counter()
    for row in csv.DictReader((SOURCE/f'lines_{n}.tsv').open(),delimiter='\t'):
        u,v,b,m=map(int,row.values());expected[u,v,b]=m
    need(actual==expected,'independent complement line table')
    receipt['canonical_frames']+=1<<M;receipt['line_rows']+=len(actual)
    receipt['orders'].append({'n':n,'ordinary':len(canonical),'frames':1<<M,'line_rows':len(actual)})
    print('Canonical/line audit order',n,'PASS',flush=True)
# Independent direct child DP for all attachments through base order five;
# all attachments of the first/last source at each later order; seeded masks elsewhere.
rng=random.Random(33117)
for r in rs:
    n=r['n']
    if n<=5:choices=range(1<<n)
    elif n<=7:choices=sorted({0,(1<<n)-1,rng.randrange(1<<n),rng.randrange(1<<n)})
    else:choices=[]
    for s in choices:need(H(child(r['rows'],s))==r['extension'][s],'literal child');receipt['literal_children']+=1
for n in (6,7,8):
    group=[r for r in rs if r['n']==n]
    for r in (group[0],group[-1]):
        for s in range(1<<n):need(H(child(r['rows'],s))==r['extension'][s],'whole response control');receipt['literal_children']+=1
for n in range(1,9):
    candidates=[(r,s,c,e) for r in rs if r['n']==n for s,c,e in r['seeds']]
    choices=candidates if n<=6 else candidates[::max(1,len(candidates)//40)]
    for r,s,c,e in choices:
        for t in (1,2):need(H(tail(r['rows'],s,t))==1+c*(1<<(e+t-1)),'literal seed tail');receipt['literal_tails']+=1
    group=[r for r in rs if r['n']==n]
    for r in group[::max(1,len(group)//16)]:
        p=list(range(n));rng.shuffle(p);renamed=[sum(((r['rows'][p[i]]>>p[j])&1)<<j for j in range(n)) for i in range(n)]
        need(H(renamed)==r['H'],'relabel H');need(H(converse(r['rows']))==r['H'],'converse H');receipt['relabelings']+=1
# Direct action proof audit for the entire saved corpus: a positive unit selection
# count must be divisible by Aut, and rooted attachment histograms sum correctly.
for r in rs:
    need(len(r['seeds'])%r['aut']==0,'unit orbit count');need(r['H']==r['tc']*r['aut'],'fibre')
    for s,o in enumerate(r['attachment_orbit']):need(r['extension'][s]==r['extension'][o],'attachment orbit invariant')
# Exact recovered small-H grammars, not an interpolation claim.
for n in range(1,9):
    g=[r for r in rs if r['n']==n]
    need(sum(r['H']==11 for r in g)==max(0,n-4),'H11 grammar')
    need(sum(r['H']==33 for r in g)==(n*n-5*n-1 if n>=6 else 0),'H33 grammar')
    need(not any(r['H']==7 or r['tc']==7 for r in g),'inherited gap controls')
receipt.update(status='PASS',record_count=len(rs),random_seed=33117,
    scope='complete independent frame/line audit through '+str(maxn)+'; all statistics use inherited complete n<=8 universe; specified literal child/tail controls')
(HERE/'verification.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps(receipt,sort_keys=True))
