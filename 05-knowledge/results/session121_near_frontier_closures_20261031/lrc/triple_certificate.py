"""Exact complete finite stopping certificate for three arbitrary attached speeds."""
from fractions import Fraction as F
from math import ceil,floor
from itertools import combinations
from pathlib import Path
import json,argparse
import source_lrc_banks as src
R=F(1,14)
def attach(intervals,v,radius=R):
    r=radius
    out=[]
    for l,h in intervals:
        for k in range(ceil(v*l+r-1),floor(v*h-r)+1):
            out.append((max(l,(k+r)/v),min(h,(k+1-r)/v)))
    return out

def width(intervals):
    return max((h-l for l,h in intervals),default=F(0))
def one_tail(intervals,v):
    # Separate arrival calculation: first safe point after each actual left endpoint.
    for l,h in intervals:
        phase=(v*l)%1
        gap= R-phase if phase<R else F(0) if phase<=1-R else 1+R-phase
        t=l+gap/v
        if t<=h: return t
    return None

def run(limit=None):
    records=[]; queries=0; joint_checks=0; maxB=maxV=maxW=0
    for ix,C in enumerate(combinations(range(1,14),10)):
        if limit and ix>=limit: break
        I=src.safe_set(C,R); W0=width(I)
        if W0<=0: raise RuntimeError(('Core has no positive interval',C))
        B=ceil(1/W0); maxB=max(maxB,B); branches=[]
        for u in range(1,B):
            J=attach(I,u)
            # Independent full comb intersection, not the clipped-tooth implementation.
            if J!=src.intersection(I,src.safe_set((u,),R)): raise RuntimeError('First joint mismatch')
            joint_checks+=1
            W1=width(J)
            if W1<=0: raise RuntimeError(('First branch has no positive interval',C,u))
            V=ceil(1/W1); maxV=max(maxV,V); sub=[]
            for v in range(u,V):
                K=attach(J,v); W2=width(K)
                if W2<=0: raise RuntimeError(('Second branch has no positive interval',C,u,v))
                T=ceil(2*R/W2); maxW=max(maxW,T); rejects=[]
                for w in range(v,T):
                    t=one_tail(K,w); truth=src.intersection(K,src.safe_set((w,),R))
                    if (t is not None)!=bool(truth): raise RuntimeError('Independent finite-head disagreement')
                    if t is None: rejects.append(w)
                    elif not all(R<= (s*t)%1 <=1-R for s in (*C,u,v,w)):
                        raise RuntimeError('Invalid original-clock witness')
                    queries+=1
                sub.append([v,str(W2),T,rejects])
            branches.append({'u':u,'width':str(W1),'pair_cutoff':V,'finite_branches':sub})
        rec={'core':C,'safe_intervals':[[str(a),str(b)] for a,b in I],'width':str(W0),'triple_cutoff':B,'branches':branches}
        records.append(rec)
        print('core',ix+1,C,'queries',queries,flush=True)
    success=all(not row[3] for rec in records for b in rec['branches'] for row in b['finite_branches'])
    return {'radius':str(R),'core_count':len(records),'finite_queries':queries,'independent_first_joint_checks':joint_checks,'maximum_triple_cutoff':maxB,'maximum_pair_cutoff':maxV,'maximum_one_tail_cutoff':maxW,'all_triples_accepted':success,'families':records}
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--limit',type=int); ap.add_argument('--compact',action='store_true',help='Recompute the complete tree in memory and publish only its compact committed receipt'); args=ap.parse_args()
    data=run(args.limit)
    if args.compact:
        from compact_receipt import compact_data
        data=compact_data(data)
    args.output.write_text(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in data.items() if k!='families'}))
