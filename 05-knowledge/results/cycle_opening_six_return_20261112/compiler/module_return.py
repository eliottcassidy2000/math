"""Literal full-count and named-arc audit of the uniform module return."""
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent.parent/'berggren_target_return_20261111'/'reset'))
import funnel_six_state as f
P=(78,12,24,16,35,79,30)
BASE_MAP={1:3,8:8,9:9,0:0,2:5,3:6,4:7,5:2,6:4,7:1}
FIXED_FLIPS={(0,4),(0,7),(4,5),(4,7),(4,8),(5,8)}
MODULE_FLIPS={4,6,8}

def require(ok,message):
    if not ok: raise RuntimeError(message)

def birth(T,r):
    n=len(T)
    return tuple(row|((1<<n) if i==r else 0) for i,row in enumerate(T))+(((1<<n)-1)^(1<<r),)

def inflate(T,v,M):
    owners=[(i,j) for i in range(len(T)) for j in range(len(M) if i==v else 1)]
    rows=[0]*len(owners)
    for i,(a,u) in enumerate(owners):
        for j,(b,w) in enumerate(owners):
            if i==j: continue
            edge=((M[u]>>w)&1) if a==b else ((T[a]>>b)&1)
            if edge: rows[i]|=1<<j
    return tuple(rows),owners

def raw(M):
    T,_=inflate(P,1,M);m=len(M)
    return birth(birth(birth(T,m+3),m+6),m+7)

def target(M):
    parent,r=f.a.funnel((2,4,1),[0],(1,0,0,0,1))
    return inflate(birth(parent,r),3,M)[0]

def counts(T,r):
    h=f.a.tables(T)[0][-1]
    starts=f.a.tables(f.a.opposite(T))[1][-1][r]
    deletion=f.a.tables(f.a.induced(T,[i for i in range(len(T)) if i!=r]))[0][-1]
    return {'h':h,'start':starts,'deletion':deletion,'defect':deletion-starts,'x':h-deletion}

def main():
    cases=[('TT'+str(k),f.a.transitive(k)) for k in range(1,8)]
    cases += [('C3',(2,4,1)),('C4',(2,12,9,1)),('C3_then_single',(10,12,9,0))]
    receipts=[]
    for name,M in cases:
        pm,_,_=f.a.bin_values(M,0,2);h,p=pm[1],pm[2]
        R=raw(M);Q=target(M);m=len(M);r=m+8
        cr,cq=counts(R,r),counts(Q,r)
        require(cr['h']==cq['h']==135*p+23*h,'universal full target '+name)
        require(cr['defect']==24*p+5*h and cq['defect']==h,'defect repair '+name)
        require(cr['deletion']==69*p+13*h and cq['deletion']==27*p+27*h,'deletion state '+name)
        require(f.a.strong(R) and f.a.strong(Q),'strong endpoints '+name)
        _,oldnames=inflate(tuple([0]*10),1,M)
        _,newnames=inflate(tuple([0]*10),3,M)
        newindex={x:i for i,x in enumerate(newnames)}
        mapping={i:newindex[(BASE_MAP[base],j)] for i,(base,j) in enumerate(oldnames)}
        require(mapping[r]==r and mapping[r-1]==r-1,'same physical root and sole predecessor')
        differences=[]
        for i in range(len(R)):
            for j in range(i+1,len(R)):
                if ((R[i]>>j)&1)!=((Q[mapping[i]]>>mapping[j])&1):
                    differences.append((i,j))
                    bi,bj=oldnames[i][0],oldnames[j][0]
                    require((1 in (bi,bj) and (bj if bi==1 else bi) in MODULE_FLIPS) or tuple(sorted((bi,bj))) in FIXED_FLIPS,'exact named arc support')
        require(len(differences)==6+3*m,'exact arc count')
        require([i for i in range(len(R)) if R[i]>>r&1]==[r-1],'raw incoming port')
        require([i for i in range(len(Q)) if Q[i]>>r&1]==[r-1],'target incoming port')
        receipts.append({'module':name,'module_rows':list(M),'module_h':h,'module_p2':p,
                         'order':len(R),'raw':cr,'returned':cq,'raw_rows':list(R),
                         'returned_rows':list(Q),'old_to_new':mapping,'changed_arcs':differences})
        f.a.tables.cache_clear()
    result={'status':'PASS','scope':'literal TT1..7 plus three nontransitive module controls',
            'base_mapping':BASE_MAP,'fixed_flips':sorted(FIXED_FLIPS),'module_boundary_flips':sorted(MODULE_FLIPS),'cases':receipts}
    (HERE/'module_return.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','cases':len(receipts),'largest_order':max(x['order'] for x in receipts)}))
if __name__=='__main__':main()
