"""Complete exact eight-core/five-tail certificate with inherited nine-core ports.

All arithmetic is rational. The compact receipt commits every deterministic
finite tree; --full-directory optionally writes the corresponding full trees.
"""
from fractions import Fraction as F
from itertools import combinations
from math import ceil, floor
from pathlib import Path
import argparse, hashlib, json
import source_lrc_banks as source
R=F(1,14)
COEFFICIENT={1:F(1,7),2:F(2,7),3:F(3,7),4:F(4,7),5:F(10,7)}
def encode(x): return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def need(ok,message):
    if not ok: raise RuntimeError(message)
def attach(intervals,v):
    out=[]
    for l,h in intervals:
        for j in range(ceil(v*l+R-1),floor(v*h-R)+1):
            a=max(l,(j+R)/v); b=min(h,(j+1-R)/v)
            need(a<=b,'clipped tooth inverted')
            out.append((a,b))
    return out
def width(intervals): return max((h-l for l,h in intervals),default=F(0))
def first_arrival(intervals,v):
    for l,h in intervals:
        p=(v*l)%1
        d=R-p if p<R else F(0) if p<=1-R else 1+R-p
        t=l+d/v
        if t<=h: return t
    return None

def produce(full_directory=None):
    records=[]; totals=[0]*6; maximum=[0]*6; inheritance=0
    for ix,C in enumerate(combinations(range(1,14),8)):
        I=source.safe_set(C,R); counts=[0]*6; cuts=[0]*6; inherited=0
        def visit(J,k,lo,tails):
            nonlocal inherited
            counts[k]+=1; W=width(J)
            need(W>0,('no positive interval',C,tails))
            B=ceil(COEFFICIENT[k]/W); cuts[k]=max(cuts[k],B)
            children=[]
            for v in range(lo,B):
                if 1<=v<=13:
                    children.append([v,'inherit_nine_core']); inherited+=1
                    continue
                if k==1:
                    counts[0]+=1; t=first_arrival(J,v)
                    literal=source.intersection(J,source.safe_set((v,),R))
                    need(t is not None and bool(literal),('finite leaf fails',C,tails,v))
                    need(all(R<=s*t%1<=1-R for s in (*C,*tails,v)),'original-clock witness')
                    children.append([v,'witness',str(t)])
                else:
                    K=attach(J,v)
                    children.append([v,visit(K,k-1,v,(*tails,v))])
            return {'width':str(W),'cutoff':B,'children':children}
        tree=visit(I,5,1,()); raw=encode(tree)
        if full_directory:
            full_directory.mkdir(parents=True,exist_ok=True)
            (full_directory/('%03d.json'%ix)).write_bytes(raw)
        rec={'core':list(C),'safe_intervals':[[str(a),str(b)] for a,b in I],
             'width':tree['width'],'cutoff':tree['cutoff'],'node_counts':counts,
             'max_cutoffs':cuts,'inherited_branches':inherited,
             'tree_bytes':len(raw),'tree_sha256':hashlib.sha256(raw).hexdigest()}
        records.append(rec)
        totals=[a+b for a,b in zip(totals,counts)]
        maximum=[max(a,b) for a,b in zip(maximum,cuts)]; inheritance+=inherited
        if (ix+1)%25==0: print('certified',ix+1,'cores',totals,flush=True)
    return {'status':'PASS','radius':str(R),'bounded_pool':list(range(1,14)),
            'core_size':8,'tail_count':5,'core_count':len(records),
            'coefficient':{str(k):str(v) for k,v in COEFFICIENT.items()},
            'node_counts':totals,'max_cutoffs':maximum,'inherited_branches':inheritance,
            'inherited_theorem':'session122/lrc/LRC.md P3 nine-core/four-tail',
            'inheritance_rule':'all bounded tails; pad by one missing bounded speed for a repeated tail',
            'all_five_tail_families_accepted':True,'families':records}
if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--output',type=Path,required=True)
    p.add_argument('--full-directory',type=Path); a=p.parse_args()
    result=produce(a.full_directory); a.output.write_bytes(encode(result))
    print(json.dumps({k:v for k,v in result.items() if k!='families'},sort_keys=True))
