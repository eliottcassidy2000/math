"""Independent complete replay and compact binding, using full-comb geometry.

No producer code, copied source geometry, or producer receipt fields determine
which finite branches are searched. Every per-core tree is reconstructed.
"""
import sys,json,hashlib
from pathlib import Path
from itertools import combinations
from fractions import Fraction as F
from math import ceil
from geometry import R,need,comb,intersection,safe,width,check_witness
COEFFICIENT={1:F(1,7),2:F(1,3),3:F(3,5),4:F(1)}
def encode(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def validate_header(d):
    need(d['status']=='PASS' and d['radius']==str(R),'status/radius')
    need(d['bounded_pool']==list(range(1,14)),'bounded pool')
    need(d['core_size']==9 and d['tail_count']==4 and d['core_count']==715,'family scope')
    need([tuple(rec['core']) for rec in d['families']]==list(combinations(range(1,14),9)),'complete unique core universe')
    need(d['coefficient']=={str(k):str(v) for k,v in COEFFICIENT.items()},'analytic coefficients')
    need(d['all_four_tail_families_accepted'] is True,'successful summary')
    need(d['inherited_theorem']=='session121_near_frontier_closures_20261031/lrc/LRC.md T3','inherited scope')
def reconstruct(C,with_tree=False):
    I=safe(C); counts=[0]*5; cuts=[0]*5; inherited=0
    def visit(J,k,lo,tails):
        nonlocal inherited
        counts[k]+=1; W=width(J); need(W>0,'positive branch width')
        B=ceil(COEFFICIENT[k]/W);cuts[k]=max(cuts[k],B);children=[]
        for v in range(lo,B):
            if 1<=v<=13 and v not in C:
                need(len(set(C)|{v})==10,'inherited ten-core condition')
                # Original row has nine bounded speeds and four arbitrary
                # tails. Promoting this v leaves exactly three arbitrary
                # attachments, including the earlier fixed tail prefix.
                children.append([v,'inherit_ten_core']);inherited+=1;continue
            K=intersection(J,comb(v))
            if k==1:
                counts[0]+=1;need(bool(K),'literal finite leaf')
                t=K[0][0];check_witness(t,(*C,*tails,v))
                children.append([v,'witness',str(t)])
            else:children.append([v,visit(K,k-1,v,(*tails,v))])
        return {'width':str(W),'cutoff':B,'children':children}
    tree=visit(I,4,1,());raw=encode(tree)
    rec={'core':list(C),'safe_intervals':[[str(a),str(b)] for a,b in I],
            'width':tree['width'],'cutoff':tree['cutoff'],'node_counts':counts,
            'max_cutoffs':cuts,'inherited_branches':inherited,
            'tree_bytes':len(raw),'tree_sha256':hashlib.sha256(raw).hexdigest()}
    return (rec,tree) if with_tree else rec
def verify_family(rec,expected=None):
    if expected is None:expected=reconstruct(tuple(rec['core']))
    need(rec==expected,'complete family tree/geometry commitment')
    return expected

def run(path):
    raw=Path(path).read_bytes();d=json.loads(raw);validate_header(d)
    totals=[0]*5;cuts=[0]*5;inherited=0;tree_bytes=0
    for ix,rec in enumerate(d['families']):
        expected=verify_family(rec)
        totals=[a+b for a,b in zip(totals,expected['node_counts'])]
        cuts=[max(a,b) for a,b in zip(cuts,expected['max_cutoffs'])]
        inherited+=expected['inherited_branches'];tree_bytes+=expected['tree_bytes']
        if (ix+1)%25==0:print('audited',ix+1,'cores',totals,file=sys.stderr,flush=True)
    need(totals==d['node_counts'] and cuts==d['max_cutoffs'] and inherited==d['inherited_branches'],'aggregate counters')
    return {'status':'PASS','core_count':715,'node_counts':totals,'max_cutoffs':cuts,
            'inherited_branches':inherited,'tree_commitments_checked':715,
            'full_tree_bytes_reconstructed':tree_bytes,'receipt_sha256':hashlib.sha256(raw).hexdigest()}
if __name__=='__main__':
    print(json.dumps(run(sys.argv[1]),sort_keys=True,indent=2))
