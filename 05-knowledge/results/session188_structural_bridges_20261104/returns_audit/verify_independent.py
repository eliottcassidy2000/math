"""Independent literal-path audit; imports no producer module."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import hashlib
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


SOURCE = Path(__file__).resolve().parent.parent / 'returns'
BASE = (2,20,9,3,45,79,31)
OUTER = (1,2,4,5,6)


def inflate(m):
    labels = [('clone',i) for i in range(m)] + [('old',u) for u in OUTER]
    rows = [0]*len(labels)
    for i,j in combinations(range(len(labels)),2):
        a,b = labels[i],labels[j]
        if a[0] == b[0] == 'clone':
            direction = a[1] < b[1]
        else:
            aa = 3 if a[0] == 'clone' else a[1]
            bb = 3 if b[0] == 'clone' else b[1]
            direction = bool(BASE[aa] >> bb & 1)
        if direction: rows[i] |= 1 << j
        else: rows[j] |= 1 << i
    W={m+OUTER.index(u) for u in (1,5)}
    return rows,W


def path_response(rows):
    # Enumerate every literal directed simple path, without subset-DP aggregation.
    n=len(rows);N=(1<<n)-1
    H=Counter({0:1});ends={};starts={};visits=0
    def visit(mask,start,end):
        nonlocal visits
        visits+=1;H[mask]+=1
        ends.setdefault(mask,Counter())[end]+=1
        starts.setdefault(mask,Counter())[start]+=1
        allowed=rows[end]&~mask
        while allowed:
            bit=allowed&-allowed;allowed-=bit
            visit(mask|bit,start,bit.bit_length()-1)
    for v in range(n):visit(1<<v,v,v)
    F=[0]*n;G=[0]*n;J=[[0]*n for _ in range(n)]
    p=0
    for A in range(N+1):
        B=N^A;p+=H[A]*H[B]
        for u,count in ends.get(A,{}).items():F[u]+=count*H[B]
        for v,count in starts.get(B,{}).items():G[v]+=H[A]*count
        for u,countu in ends.get(A,{}).items():
            for v,countv in starts.get(B,{}).items():J[u][v]+=countu*countv
    return (H[N],p,F,G,J,[ends[N][v] for v in range(n)],
            [starts[N][v] for v in range(n)]),visits


def signature(R,W):
    h,p,F,G,J,e,s=R
    return (h,p,sum(F[v] for v in W),sum(G[v] for v in W),
            sum(J[u][v] for u in W for v in W),sum(e[v] for v in W),sum(s[v] for v in W))


def ordered(a,b):
    if a is None:return b
    if b is None:return a
    h,p,f,g,j,e,s=a;H,P,F,G,J,E,S=b
    return (h*H,p*P,f*H+p*F,g*P+h*G,j*H+h*J+e*S+g*F,h*E,s*H)


def contexts():
    out=[None]
    for n in (1,2,3):
        edges=list(combinations(range(n),2))
        for mask in range(1<<len(edges)):
            rows=[0]*n
            for k,(u,v) in enumerate(edges):
                if mask>>k&1:rows[u]|=1<<v
                else:rows[v]|=1<<u
            R,_=path_response(rows)
            for ports in range(1<<n):out.append(signature(R,{i for i in range(n) if ports>>i&1}))
    return out


def run():
    prior=json.loads((SOURCE/'receipt.json').read_text())
    for f,h in prior['sha256'].items():
        require(hashlib.sha256((SOURCE/f).read_bytes()).hexdigest()==h,'producer source hash: '+f)
    require((SOURCE/'normal.stdout').read_bytes()==(SOURCE/'optimized.stdout').read_bytes(),
            'producer normal/optimized receipts differ')
    samples=[];survivors=[];masks=0
    for m in (1,2,3,4):
        rows,W=inflate(m);R,paths=path_response(rows)
        h,p,F,G,J,e,s=R
        require(h==9*2**m-1,'Hamilton family')
        require(signature(R,W)[4]==6*2**m-4,'external-port joint')
        for mask in range(1,1<<m):
            Z=W|{i for i in range(m) if mask>>i&1}
            sig=signature(R,Z);c=sig[0]-sig[4];masks+=1
            if m>=2:
                require(sig[4]>=35*2**(m-1),'all-port lower bound')
                require((sig[4]==35*2**(m-1))==(mask==1),'sharp equality mask')
                if c>=-103:survivors.append(sig)
            if m in (2,3):
                old=next(x for x in prior['small_candidates'] if x['m']==m and x['module_port_mask']==mask)
                require(old['rows']==rows and old['signature']==list(sig),'literal small receipt')
        for i in range(m):
            j=m-1-i
            pred=6*2**m-4+(18*3**i+5)*2**j+8*2**i-int(j==0)-4*int(i==0)
            require(signature(R,W|{i})[4]==pred,'singleton position formula')
        if m>=2:
            expected=(9*2**m-1,66*3**m+28*2**m-30,66*2**m-12,
                      24*3**m+15*2**m-9,35*2**(m-1),13,7*2**(m-1))
            require(signature(R,W|{0})==expected,'seven-coordinate first-clone formula')
        samples.append({'m':m,'vertices':m+5,'literal_directed_paths':paths,
                        'first_clone_signature':signature(R,W|{0})})
    ctx=contexts();tests=0
    require(len(survivors)==3,'analytic small-case survivor count')
    for B in survivors:
        for X in ctx:
            for Y in ctx:
                Z=ordered(ordered(X,B),Y);c=Z[0]-Z[4]
                require(c not in (-103,-101),'finite ordered context hostile')
                tests+=1
    result={'status':'PASS','method':'literal enumeration of every directed simple path; no producer imports',
            'named_masks_checked':masks,'literal_families':samples,
            'marked_context_presentations_including_empty':len(ctx),
            'ordered_context_pairs_checked':tests,
            'producer_source_hashes':prior['sha256'],
            'limitations':'All-order reasoning is audited in AUDIT.md; these finite tests are corroboration.'}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':run()
