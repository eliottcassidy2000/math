"""Independent exact full-comb replay. Imports no producer or source modules."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from math import ceil
from pathlib import Path
import json, sys, hashlib
R=F(1,14)
def need(ok, label):
    if not ok: raise RuntimeError(label)
@lru_cache(None)
def comb(v,r=R):
    return tuple(((j+r)/v,(j+1-r)/v) for j in range(v))
def intersect(a,b):
    i=j=0; out=[]
    while i<len(a) and j<len(b):
        l=max(a[i][0],b[j][0]); h=min(a[i][1],b[j][1])
        if l<=h:
            if out and l<=out[-1][1]: out[-1]=(out[-1][0],max(h,out[-1][1]))
            else: out.append((l,h))
        if a[i][1]<b[j][1]: i+=1
        elif a[i][1]>b[j][1]: j+=1
        else: i+=1; j+=1
    return tuple(out)
def safe(speeds,r=R):
    out=((F(0),F(1)),)
    for v in sorted(set(speeds)): out=intersect(out,comb(v,r))
    return out
def width(a): return max((h-l for l,h in a),default=F(0))
def run(path):
    raw=Path(path).read_bytes(); d=json.loads(raw); families=d['families']
    need([tuple(x['core']) for x in families]==list(combinations(range(1,14),10)),'complete unique core census')
    need(d['radius']==str(R),'radius')
    first=second=leaves=0; mb=mv=mt=0
    for index, rec in enumerate(families):
        C=tuple(rec['core']); I=safe(C); W=width(I)
        need(I==tuple(tuple(F(x) for x in row) for row in rec['safe_intervals']),'base cover')
        need(W>0 and W==F(rec['width']),'base positive width'); B=ceil(1/W); mb=max(mb,B)
        need(B==rec['triple_cutoff'],'base cutoff')
        need([b['u'] for b in rec['branches']]==list(range(1,B)),'first coverage')
        for b in rec['branches']:
            u=b['u']; J=intersect(I,comb(u)); W1=width(J); first+=1
            need(W1>0 and W1==F(b['width']),'first width'); V=ceil(1/W1); mv=max(mv,V)
            need(V==b['pair_cutoff'],'pair cutoff')
            need([row[0] for row in b['finite_branches']]==list(range(u,V)),'second coverage')
            for v,w2,t,rejects in b['finite_branches']:
                K=intersect(J,comb(v)); W2=width(K); second+=1
                need(W2>0 and W2==F(w2),'second width'); T=ceil(2*R/W2); mt=max(mt,T)
                need(T==t,'single cutoff'); actual=[]
                for w in range(v,T):
                    L=intersect(K,comb(w)); leaves+=1
                    if not L: actual.append(w)
                    else:
                        witness=L[0][0]
                        need(all(R<= (s*witness)%1 <=1-R for s in (*C,u,v,w)),'original-clock witness')
                need(actual==rejects==[],'finite rejection')
        if (index+1)%25==0: print('audited',index+1,'cores',second,'pairs',leaves,'leaves',file=sys.stderr,flush=True)
    need(d['core_count']==286 and d['all_triples_accepted'] is True,'receipt summary')
    need((first,leaves,mb,mv,mt)==(d['independent_first_joint_checks'],d['finite_queries'],d['maximum_triple_cutoff'],d['maximum_pair_cutoff'],d['maximum_one_tail_cutoff']),'summary counts')
    # Hostile to marginal-only acceptance and loss of closed singleton endpoints.
    Q=F(1,4); I=safe((1,2),Q)
    need(bool(intersect(I,comb(3,Q))) and bool(intersect(I,comb(12,Q))),'individual acceptance')
    need(not safe((1,2,3,12),Q),'joint hostile')
    need(intersect(((F(1,3),F(1,3)),),comb(1,F(1,3)))==((F(1,3),F(1,3)),),'closed endpoint')
    need(not intersect(((F(1,3),F(1,3)),),comb(3,F(1,3))),'unsafe endpoint')
    # Sufficient single-tail threshold is inclusive; strict-width test would lose this point.
    need(intersect(((F(0),F(1,7)),),comb(1,R))==((R,F(1,7)),),'one-tail equality')
    # Tight interval around bad arc at 1/2, v=2, W=2r/v.
    need(intersect(((F(1,2)-R/2,F(1,2)+R/2),),comb(2,R))==((F(1,2)-R/2,F(1,2)-R/2),(F(1,2)+R/2,F(1,2)+R/2)),'equality only endpoints')
    boundary=0
    for k in (2,3,4):
        r=F(1,4*k-2)
        for W in (F(1,2),F(1,3),F(2,5)):
            B=ceil(1/W)
            for left in (F(0),F(1,17),F(2,7)):
                for tails in combinations_with_replacement(range(B,B+6),k):
                    I=((left,left+W),)
                    for v in tails: I=intersect(I,comb(v,r))
                    need(sum(h-l for l,h in I)>0,'critical conductor boundary')
                    boundary+=1
    return {'status':'PASS','core_count':286,'first_branches':first,'second_branches':second,'finite_leaves':leaves,'max_cutoffs':[mb,mv,mt],'conductor_boundary_controls':boundary,'hostile_endpoint_controls':6,'receipt_sha256':hashlib.sha256(raw).hexdigest()}
if __name__=='__main__':
    path=sys.argv[1] if len(sys.argv)>1 else str(Path(__file__).parent.parent/'lrc/triple_receipt.json')
    print(json.dumps(run(path),sort_keys=True))
