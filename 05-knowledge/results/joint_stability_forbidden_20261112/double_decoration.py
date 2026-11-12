"""Two fixed exterior vertices around the physical crossing's extreme ports.

All16384 incidences with the old seven exterior vertices, with new9->10.
Relabeling the two new vertices covers the opposite orientation. Both named
ports retain their unique selected gate, so the uniform two-run theorem
turns the four-corner check into a complete all-length decision.
"""
from collections import Counter
import json,time
import numpy as np
from census_nine import need,save,tables,defect,strong
from joint_nine import cones,key,kernel,quotient
from stable_census import certify
Q=[2,100,49,7,75,409,173,31,223];P=[0,8]

def unit_marks(rows,lengths):
    rr,_,labels=cones.expand(rows,[],P,lengths);h,u,w,marks=tables(np.array(rr,dtype=np.int64));out={}
    for mask in marks:
        chosen=[labels[i] for i in range(len(rr)) if mask>>i&1]
        if any(v in P for v in chosen) or not {1,5}<=set(chosen):continue
        old=sum(1<<v for v in chosen);out[old]=(int(mask),rr,int(h),sum(int(u[i]) for i in range(len(rr)) if mask>>i&1))
    return out

def run():
    stats=Counter();out=[];seen=set();start=time.time()
    for bits in range(1<<14):
        rows=list(Q)+[1,1];rows[8]|=(1<<9)|(1<<10);rows[9]|=1<<10
        for i,v in enumerate(range(1,8)):
            for j in range(2):
                w=9+j
                if bits>>(2*i+j)&1:rows[w]|=1<<v
                else:rows[v]|=1<<w
        if not strong(rows):continue
        stats['strong_quotients']+=1;left=unit_marks(rows,[1,2])
        if not left:continue
        stats['one_unit_line_candidates']+=len(left);right=unit_marks(rows,[2,1]);common=left.keys()&right.keys()
        if not common:continue
        both=unit_marks(rows,[2,2]);common&=both.keys()
        for mask in sorted(common):
            h,u,w,_=tables(np.array(rows,dtype=np.int64));d=defect(h,w,[i for i in range(11) if mask>>i&1]);stats['three_unit_corners']+=1
            if d==1:continue
            need((d-1)%4==0,'integral two-run crossing coefficient');ckey=key(rows,mask)
            if ckey in seen:continue
            seen.add(ckey);K=kernel(rows,mask,P);need(K.get((2,2),0)==(d-1)//4,'full tensor crossing coefficient')
            roots=[]
            for m,rr,hh,uu in [left[mask],right[mask]]:
                ps,cert,_,_=certify(rr,[m]);_,_,comps=quotient(rr,m);eligible=[c for c in comps if len(c)>1 and set(c)&set(ps[m])]
                roots.append({'rows':rr,'selected':m,'H':hh,'U':uu,'stable_ports':ps[m],'eligible_blocks':eligible})
            r0,r1=roots;full=(1<<len(r1['rows']))-1;rev=[full^r^(1<<i) for i,r in enumerate(r1['rows'])]
            terminal=all(not r['eligible_blocks'] for r in roots);different=key(r0['rows'],r0['selected']) not in [key(r1['rows'],r1['selected']),key(rev,r1['selected'])]
            rec={'extension_bits':bits,'quotient_rows':rows,'quotient_selected':mask,'crossing_coefficient':(d-1)//4,'terminal_roots':terminal,'roots_different_modulo_converse':different,'roots':roots,'Stirling_terms':[[*a,c] for a,c in sorted(K.items())]}
            out.append(rec);print('CROSSING '+json.dumps(rec),flush=True)
        if bits%1000==0:print(json.dumps({'at':bits,'stats':dict(stats),'types':len(out),'seconds':round(time.time()-start,1)}),flush=True)
    result={'scope':'Complete two-fixed-vertex decoration of the nine-vertex crossing quotient with both extreme-port incidences preserved. All16384 exterior incidence choices with new9->10, all selections including the two singleton selected gates. Exact uniform two-run classification; not all order11 quotients.','graph_presentations':16384,'statistics':dict(stats),'distinct_marked_crossings':len(out),'records':out}
    save('double_decoration.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
