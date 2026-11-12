"""All crossing thresholds for extreme-gated two-port quotients through9.

No assumption that the crossing is at lengths(2,2): the exact two-run
polynomial determines both roots, which are tested for positive powers2.
"""
from pathlib import Path
from collections import Counter
import gzip,csv,json,time
import numpy as np
from census_nine import need,save,tables,defect,strong
from joint_nine import cones,key,kernel,quotient
from two_block_probe import decode
HERE=Path(__file__).resolve().parent

def graphs():
    with gzip.open(HERE.parent/'metagraph_aspect_sequences_20261111/records.jsonl.gz','rt') as f:
        for line in f:yield json.loads(line)['rows']
    with gzip.open(HERE.parent/'blue_quotient_20261020b/n9_classes.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):yield decode(int(r['code']),9)

def pow2(x):return x.bit_length()-1 if x>=2 and not x&(x-1) else None

def run():
    stats=Counter();out=[];seen=set();start=time.time()
    for rows in graphs():
        n=len(rows);stats[f'classes_{n}']+=1
        if not strong(rows):continue
        plus=[v for v in range(n) if rows[v].bit_count()==1];minus=[v for v in range(n) if rows[v].bit_count()==n-2]
        for p in plus:
            for q in minus:
                if p==q or rows[p]>>q&1:continue
                gp=rows[p].bit_length()-1;gq=next(v for v in range(n) if v!=q and rows[v]>>q&1)
                need(gp not in [p,q] and gq not in [p,q],'fixed singleton selected gates');stats['extreme_port_pairs']+=1
                corners={}
                for a,b in [(1,1),(2,1),(1,2),(2,2)]:
                    rr,_,labs=cones.expand(rows,[],[p,q],[a,b]);h,u,w,_=tables(np.array(rr,dtype=np.int64));corners[a,b]=(int(h),w,[labs.index(v) for v in range(n)])
                for O in range(1<<n):
                    if O>>p&1 or O>>q&1 or not O>>gp&1 or not O>>gq&1:continue
                    stats['mark_presentations']+=1;chosen=[v for v in range(n) if O>>v&1];ds={}
                    for ab,(h,w,positions) in corners.items():ds[ab]=defect(h,w,[positions[v] for v in chosen])
                    a,b,c,d=ds[1,1],ds[2,1],ds[1,2],ds[2,2];k=(d-b-c+a)//4
                    need((d-b-c+a)%4==0,'integral bilinear coefficient')
                    if not k:continue
                    ell=(b-a)//2-2*k;emm=(c-a)//2-2*k;const=a-4*k-2*ell-2*emm
                    need((b-a)%2==0 and (c-a)%2==0,'integral linear coefficients')
                    if const*k-ell*emm!=k or (-emm)%k or (-ell)%k:continue
                    A=pow2((-emm)//k);B=pow2((-ell)//k)
                    if A is None or B is None:continue
                    ckey=(key(rows,O),p,q) # exact port labels retained in this presentation count
                    if ckey in seen:continue
                    seen.add(ckey);K=kernel(rows,O,[p,q]);_,_,comps=quotient(rows,O)
                    rec={'rows':rows,'selected':O,'ports':[p,q],'thresholds':[A,B],'coefficient':k,'crossing_order':n-2+A+B,'marked_quotient_key':key(rows,O),'twin_blocks':[c for c in comps if len(c)>1],'Stirling_terms':[[*r,v] for r,v in sorted(K.items())]}
                    out.append(rec);print('EXTREME_CROSSING '+json.dumps(rec),flush=True)
        if sum(v for k,v in stats.items() if k.startswith('classes_'))%50000==0:print(json.dumps({'sources':sum(v for k,v in stats.items() if k.startswith('classes_')),'pairs':stats['extreme_port_pairs'],'crossings':len(out),'seconds':round(time.time()-start,1)}),flush=True)
    result={'scope':'Every ordinary tournament through9, restricted to strong quotients with one unique-out gate and one unique-in gate at distinct unselected ports. Every selection containing the two singleton gates. The uniform two-run identity classifies all positive crossing thresholds, with no length cutoff. Port presentations are retained, not asserted to be isomorphism classes of port-marked quotients.','statistics':dict(stats),'crossing_presentations':len(out),'records':out}
    save('extreme_crossings.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
