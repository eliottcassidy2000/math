"""Compile every single birth clock through10 into its exact slope modes.

For unselected p, a taxed U-cover has at most one run beginning a path
and at most two ending paths. Thus R<=min(dminus+1,dplus+2), uniformly
in its own length. R values determine every positive-length slope.
"""
from collections import Counter
import json,time
import numpy as np
from census_nine import read,save,need,tables,ib

def onto(m,r):
    from math import comb
    return sum((-1)**(r-j)*comb(r,j)*j**m for j in range(r+1))

def run():
    bank=read('birth_union_ten.json.gz')['records'];out=[];binary=[];stats=Counter();start=time.time()
    for idx,r in enumerate(bank):
        rows=r['rows'];n=len(rows);O=r['selected']
        for v in r['ports']:
            R=min(n-rows[v].bit_count(),rows[v].bit_count()+2);values=[];coeff=[]
            for m in range(1,R+1):
                labels,rr=ib.inflate(rows,v,m);ss=[i for i,q in enumerate(labels) if O>>q&1];h,u,w,_=tables(np.array(rr,dtype=np.int64));value=sum(int(u[i]) for i in ss);values.append(value)
                rem=value-sum(coeff[j-1]*onto(m,j) for j in range(1,m));den=onto(m,m);need(rem%den==0,'integer slope skeleton coefficient');coeff.append(rem//den)
            need(all(c>=0 for c in coeff),'positive U skeleton counts')
            rec={'birth_id':idx,'rows':rows,'selected':O,'port':v,'R':R,'values':values,'Stirling_coefficients':coeff}
            out.append(rec);stats['axes']+=1
            if any(coeff[2:]):stats['higher_modes']+=1;continue
            q=coeff[1] if len(coeff)>1 else 0;b=coeff[0]-2*q
            if q==0 or b==0:stats['constant_odd_part_axes']+=1;continue
            j=(abs(b)&-abs(b)).bit_length()-1;vq=(q&-q).bit_length()-1;a0=max(1,j-vq+1)
            A=(q*2**(a0-1))//2**j;B=b//2**j
            need(A>0 and B%2==1 and q*2**a0%2**(j+1)==0,'odd coefficient normalization')
            heads=[{'s':s,'offset':1+B*2**s,'coefficient':A*2**s} for s in range(1,j+1)]
            item={**rec,'slope_binary':[b,q],'valuation_after_start':j,'start_port_length':a0,'odd_coefficient_A':A,'odd_coefficient_B':B,'coefficient_law':'A*2^k+B, k>=1','source_order_constant':n+a0-2,'unpaid_heads':heads}
            binary.append(item);stats['binary_nonconstant_axes']+=1
        if idx%100==0:print(json.dumps({'births':idx,'total':len(bank),'stats':dict(stats),'seconds':round(time.time()-start,1)}),flush=True)
    result={'scope':'Every single unselected stable birth clock from the complete1150-birth bank through10. Exact full U modes certified by the one-sided selected-start run bound. Binary affine normalization only when higher Stirling coefficients vanish. No general coverage conclusion.','statistics':dict(stats),'distinct_binary_laws':len({(r['odd_coefficient_A'],r['odd_coefficient_B']) for r in binary}),'records':out,'binary_families':binary}
    save('affine_births.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k not in ['records','binary_families']}),flush=True)
if __name__=='__main__':run()
