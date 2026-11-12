"""All order10 extreme-gated quotients via their complete order8 interiors.

Use three corners to solve candidate power-of-two crossing thresholds;
only candidates require the doubly inflated fourth corner. This includes
arbitrary frozen interior structure, not just extensions of the old J2 core.
"""
from pathlib import Path
from collections import Counter
import gzip,json,time
import numpy as np
from numba import njit
from census_nine import need,save,tables,strong,defect
from joint_nine import cones,key,kernel,quotient
from stable_census import certify
HERE=Path(__file__).resolve().parent

@njit(cache=True)
def candidates(h0,w0,h1,w1,h2,w2,g,h):
    out=[];d0=np.zeros(256,dtype=np.int64);d1=np.zeros(256,dtype=np.int64);d2=np.zeros(256,dtype=np.int64)
    for m in range(1,256):
        b=m&-m;v=0
        while 1<<v!=b:v+=1
        rest=m^b;d0[m]=d0[rest];d1[m]=d1[rest];d2[m]=d2[rest]
        for u in range(v+1,8):
            if rest>>u&1:d0[m]+=w0[v,u];d1[m]+=w1[v,u];d2[m]+=w2[v,u]
        if not m>>g&1 or not m>>h&1:continue
        a=h0-d0[m]-1;b1=h1-d1[m]-1;b2=h2-d2[m]-1
        if not a or a==b1 or a==b2:continue
        if 2*a%(a-b1) or 2*a%(a-b2):continue
        X=2+2*a//(a-b1);Y=2+2*a//(a-b2)
        if X<4 or Y<4 or X&(X-1) or Y&(Y-1):continue
        A=B=0;xx=X;yy=Y
        while xx>1:A+=1;xx//=2
        while yy>1:B+=1;yy//=2
        if a%((2-X)*(2-Y)):continue
        out.append((m,A,B,a//((2-X)*(2-Y))))
    return out

def run():
    # The known crossing, with an unselected ordinal source and the eight
    # exterior vertices placed first, controls the complete corner solver.
    q=[2,100,49,7,75,409,173,31,223]+[511];order=[9,1,2,3,4,5,6,7,0,8]
    control=[sum(1<<j for j,v in enumerate(order) if q[u]>>v&1) for u in order];ct=[]
    for sizes in [[1,1],[2,1],[1,2]]:
        rr,_,_=cones.expand(control,[],[8,9],sizes);h,u,w,_=tables(np.array(rr,dtype=np.int64));ct.extend([h,w])
    need((50,2,2,1) in candidates(*ct,1,5),'known physical crossing positive control')
    interiors=[]
    with gzip.open(HERE.parent/'metagraph_aspect_sequences_20261111/records.jsonl.gz','rt') as f:
        for line in f:
            r=json.loads(line)
            if r['n']==8:interiors.append(r['rows'])
    need(len(interiors)==6880,'complete ordinary order8 interiors');stats=Counter();out=[];false=[];seen=set();start=time.time()
    for idx,R in enumerate(interiors):
        for g in range(8):
            for h in range(8):
                rows=[r|((1<<8) if v!=g else 0)|((1<<9) if v==h else 0) for v,r in enumerate(R)]+[1<<g,511^(1<<h)]
                stats['graph_presentations']+=1
                if not strong(rows):continue
                stats['strong_presentations']+=1;H0,U0,W0,_=tables(np.array(rows,dtype=np.int64))
                rr,_,_=cones.expand(rows,[],[8,9],[2,1]);H1,U1,W1,_=tables(np.array(rr,dtype=np.int64))
                rr,_,_=cones.expand(rows,[],[8,9],[1,2]);H2,U2,W2,_=tables(np.array(rr,dtype=np.int64))
                cand=candidates(H0,W0,H1,W1,H2,W2,g,h)
                if not cand:continue
                stats['three_corner_candidates']+=len(cand);rr,_,_=cones.expand(rows,[],[8,9],[2,2]);H3,U3,W3,_=tables(np.array(rr,dtype=np.int64))
                for O,A,B,c in cand:
                    d=defect(H3,W3,[v for v in range(8) if O>>v&1])
                    if d!=1+c*(4-2**A)*(4-2**B):
                        false.append({'rows':rows,'selected':O,'thresholds':[A,B],'predicted_mixed_coefficient':c,'predicted_double_defect':1+c*(4-2**A)*(4-2**B),'actual_double_defect':d});continue
                    stats['crossing_presentations']+=1;ckey=key(rows,O)
                    if ckey in seen:continue
                    seen.add(ckey);_,_,comps=quotient(rows,O);K=kernel(rows,O,[8,9]);need(K.get((2,2),0)==c,'independent full mixed coefficient')
                    roots=[]
                    for lengths in [[1,B],[A,1]]:
                        ss=[v for v in range(8) if O>>v&1];root,ss,_=cones.expand(rows,ss,[8,9],lengths);m=sum(1<<v for v in ss)
                        if len(root)>18:roots.append({'order':len(root),'axis_lengths':lengths,'not_literally_audited':True});continue
                        pp,cert,_,_=certify(root,[m]);_,_,rootblocks=quotient(root,m);eligible=[bl for bl in rootblocks if len(bl)>1 and set(bl)&set(pp[m])]
                        roots.append({'rows':root,'selected':m,'stable_ports':pp[m],'eligible_blocks':eligible,'axis_lengths':lengths})
                    terminal=all('eligible_blocks' in r and not r['eligible_blocks'] for r in roots);nonconverse=None
                    if terminal:
                        r0,r1=roots;full=(1<<len(r1['rows']))-1;rev=[full^r^(1<<i) for i,r in enumerate(r1['rows'])];nonconverse=key(r0['rows'],r0['selected']) not in [key(r1['rows'],r1['selected']),key(rev,r1['selected'])]
                    rec={'interior_index':idx,'gates':[g,h],'quotient_rows':rows,'selected':O,'ports':[8,9],'thresholds':[A,B],'coefficient':c,'crossing_order':8+A+B,'twin_blocks':[bl for bl in comps if len(bl)>1],'terminal_roots':terminal,'roots_nonconverse':nonconverse,'roots':roots,'Stirling_terms':[[*a,v] for a,v in sorted(K.items())]}
                    out.append(rec);print('TEN_EXTREME_CROSSING '+json.dumps(rec),flush=True)
        if idx%100==0:print(json.dumps({'interiors':idx,'total':6880,'statistics':dict(stats),'distinct_marked_quotients':len(out),'seconds':round(time.time()-start,1)}),flush=True)
    need(stats['graph_presentations']==440320,'complete declared gate choices')
    result={'scope':'Complete ordinary order8 interior bank, all64 choices of unique selected out/in gates, strong order10 quotients, every selection containing both gates. All two-run crossing thresholds>=2 are solved exactly from three corners and verified at the fourth; complete kernels on hits. This covers every strong order10 quotient with the two specified extreme-port types, up to presentations.','statistics':dict(stats),'known_crossing_positive_control':True,'false_three_corner_candidates':false,'distinct_marked_crossing_quotients':len(out),'records':out}
    save('extreme_ten.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
