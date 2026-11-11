"""Cheap structured probes beyond the exhaustive nine-vertex boundary."""
from pathlib import Path
from collections import Counter
import json
HERE=Path(__file__).resolve().parent
def response(rows):
    n=len(rows);N=1<<n;full=N-1
    end=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)]
    h=[0]*N;h[0]=1
    for S in range(1,N):
        for v in range(n):
            if not S>>v&1:continue
            R=S^(1<<v)
            if not R:end[S][v]=start[S][v]=1
            else:
                end[S][v]=sum(end[R][u] for u in range(n) if rows[u]>>v&1)
                start[S][v]=sum(start[R][u] for u in range(n) if rows[v]>>u&1)
        h[S]=sum(end[S])
    uv=[sum(start[S][j]*h[full^S] for S in range(1,N)) for j in range(n)]
    C=[[0]*n for _ in range(n)]
    for S in range(1,full):
        for i,a in enumerate(end[S]):
            if a:
                for j,b in enumerate(start[full^S]):C[i][j]+=a*b
    U=[0]*N;W=[0]*N;E=[0]*N;F=[h[full]]*N;seeds=[]
    for S in range(1,N):
        bit=S&-S;j=bit.bit_length()-1;R=S^bit
        U[S]=U[R]+uv[j];E[S]=E[R]+end[full][j]
        W[S]=W[R]+sum(C[i][j]+C[j][i] for i in range(n) if R>>i&1)
        F[S]+=U[S]-E[S]-W[S]
        if W[S]==h[full]-1 and U[S]>0:
            c=U[S];e=1
            while c%2==0:c//=2;e+=1
            seeds.append([S,c,e])
    return h[full],F,seeds
def dual(rows):return [((1<<len(rows))-1)^r^(1<<i) for i,r in enumerate(rows)]
def cycle_substitution(blocks):
    out=[];offset=0;total=sum(map(len,blocks))
    for i,B in enumerate(blocks):
        target=(i+1)%3;begin=sum(len(A) for A in blocks[:target]);targetmask=((1<<len(blocks[target]))-1)<<begin
        out += [(r<<offset)|targetmask for r in B];offset+=len(B)
    return out
def cayley(p,S):return [sum(1<<((i+s)%p) for s in S) for i in range(p)]
B=[14,4,8,2] # singleton -> directed triangle
cases=[('decoration_81_27',cycle_substitution([B,B,B]),cycle_substitution([dual(B),B,B])),
       ('uniform_incidence_13',cayley(13,{2,7,8,9,10,12}),cayley(13,{1,2,4,7,8,10}))]
out=[]
for name,A,B in cases:
    a,fa,sa=response(A);b,fb,sb=response(B);ea=set(fa);eb=set(fb)
    out.append({'name':name,'order':len(A),'rows':[A,B],'H':[a,b],
        'support_sizes':[len(ea),len(eb)],'equal_support':ea==eb,
        'equal_distribution':Counter(fa)==Counter(fb),
        'symmetric_difference_size':len(ea^eb),'first_separating_value':min(ea^eb) if ea^eb else None})
    print(name,'H',a,b,'support equal',ea==eb,flush=True)
# A concrete converse hostile to seed-coefficient folding, from the full n6 corpus.
import gzip
rs=[json.loads(s) for s in gzip.open(HERE/'records.jsonl.gz','rt')]
r=next(r for r in rs if r['n']==6 and r['rep']==1)
a,fa,sa=response(r['rows']);b,fb,sb=response(dual(r['rows']))
out.append({'name':'seed_direction_is_not_converse_invariant','order':6,'rows':[r['rows'],dual(r['rows'])],
            'H':[a,b],'seeds':[sa,sb],'equal_extension_support':set(fa)==set(fb)})
# Padding theorem controls, checking U exponents as well as odd coefficients.
for r in rs:
    if r['n']>5 or not r['seeds']:continue
    n=r['n'];padded=[((1<<(n+1))-1)^1]+[row<<1 for row in r['rows']]
    hp,fp,sp=response(padded);expected={(S<<1,c,e) for S,c,e in r['seeds']}
    if not expected.issubset(set(map(tuple,sp))):raise RuntimeError('source padding failed')
(HERE/'hostiles.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print('source-padding controls PASS')
