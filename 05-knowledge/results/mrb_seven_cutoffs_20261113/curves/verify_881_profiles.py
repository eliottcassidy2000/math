"""Independent literal construction/count audit of both inherited 881 sources."""
from pathlib import Path
from hashlib import sha256
import json,subprocess
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
checks=0

def require(ok,msg):
    global checks
    checks+=1
    if not ok:raise RuntimeError(msg)

def count(t,start=None):
    n=len(t);dp=[{} for _ in range(1<<n)]
    for i in range(n) if start is None else [start]:dp[1<<i][i]=1
    for mask,d in enumerate(dp):
        for u,value in d.items():
            rem=t[u]&~mask
            while rem:
                b=rem&-rem;v=b.bit_length()-1;target=dp[mask|b]
                target[v]=target.get(v,0)+value
                rem-=b
    return sum(dp[-1].values())

def restrict(t,vs):return [sum(((t[u]>>v)&1)<<j for j,v in enumerate(vs)) for u in vs]

def funnel(t,r,word):
    n=len(t);l=len(word);q=n+l
    out=[row|sum(1<<j for j in range(n,n+l))|(0 if i==r else 1<<q) for i,row in enumerate(t)]
    for j,bit in enumerate(word):
        out.append(sum(1<<v for v in range(n+j+1,n+l))|((1<<q) if bit=='1' else 0))
    out.append((1<<r)|sum(1<<(n+j) for j,bit in enumerate(word) if bit=='0'))
    return out,q

def ready(words):
    t,r=[0],0
    for word in words:t,r=funnel(t,r,word)
    n=len(t)
    return [row|((1<<n) if i==r else 0) for i,row in enumerate(t)]+[((1<<n)-1)^(1<<r)],n,r

def inflate(t,v,k):
    owner=[u for u in range(len(t)) for _ in range(k if u==v else 1)]
    return [sum((int(i<j) if u==w else ((t[u]>>w)&1))<<j for j,w in enumerate(owner) if i!=j) for i,u in enumerate(owner)]

inputs=[(['11','1','11'],[881,1529,1038,678,678,657,462,462]),(['1','','00101'],[770,913,824,617,430,430,501,339,673])]
rows=[]
for words,expected in inputs:
    t,z,q=ready(words)
    require(all(((t[i]>>j)&1)+((t[j]>>i)&1)==1 for i in range(len(t)) for j in range(i+1,len(t))),'literal source tournament type')
    require([i for i in range(len(t)) if (t[i]>>z)&1]==[q],'sole predecessor')
    h=count(t);parent=restrict(t,[i for i in range(len(t)) if i!=z])
    require(h==881,'full source H881')
    require(count(parent)-count(t,z)==1,'ready defect one')
    require(count(parent,q)==1,'unique-start inverse')
    profiles=[]
    for v in range(len(t)):
        if v in (z,q):continue
        f=[count(inflate(t,v,k)) for k in [1,2,3]]
        require(f[0]==881,'collapsed count')
        E=f[2]-3*f[1]+2*f[0]
        require((f[1]-f[0])%2==0 and E>=0 and E%12==0,'coefficient type/sign/divisibility')
        A2=(f[1]-f[0])//2
        require(A2!=405,'next-ray coefficient absent')
        profiles.append({'vertex':v,'F1_F2_F3':f,'A1':881,'A2':A2,'A3':E//6})
    require(sorted(p['A2'] for p in profiles)==sorted(expected),'all 17 profile values independently match')
    rows.append({'words':words,'masks':t,'root':z,'predecessor':q,'order':len(t),'H':h,'parent_H':count(parent),'profiles':profiles})

commit='a03df4746338daf016df8819bcce373ec35c95e8'
paths=['05-knowledge/results/session296_multicore_metagraph_20261113/cores/PROOF.md','05-knowledge/results/session297_parity_cutoff_20261113/BRIDGES.md']
pins=[]
for path in paths:
    data=subprocess.check_output(['git','show',commit+':'+path],cwd=ROOT).decode('utf-8').replace('\r\n','\n')
    pins.append({'commit':commit,'path':path,'text_lf_sha256':sha256(data.encode('utf-8')).hexdigest()})
receipt={'status':'FINITE-EXACT independent literal 17-vertex-profile audit; completeness inherited, not inferred','checks':checks,'universe':'both classified rooted 881 source representatives, every non-root/non-predecessor vertex, TT_k inflation at k=1,2,3','source_completeness_dependency':'Session296 cores/PROOF.md Section4 in pinned incoming commit','controls':'graph type, root-star, full H, direct defect and unique-start inverse, coefficient integrality/divisibility','sources':rows,'incoming_source_pins':pins}
(HERE/'881_profiles.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'checks':checks,'profiles':sum(len(row['profiles']) for row in rows),'status':'PASS'}))
