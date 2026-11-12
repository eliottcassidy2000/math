"""Finite seed coefficient band, complete lower heads, and a new renewable ray.

The seed universe is exactly the ordinary strong marked unit seeds of order10.
Each coefficient stores an actual source with smallest dyadic valuation.
All finite heads are checked against the inherited complete order11 spectrum.
"""
from pathlib import Path
from collections import Counter
import hashlib,json
import numpy as np
from census_nine import read,save,need,tables,defect
from joint_nine import cones
from c3_pair_family import exponentials
from construct import tail
HERE=Path(__file__).resolve().parent

def run():
    source=read('unit_ten.json.gz');best={}
    for r in source['records']:
        for m in r['marks']:
            u=m['U'];need(u>0,'positive selected slope');v=(u&-u).bit_length()-1;c=u>>v
            candidate={'c':c,'v2U':v,'U':u,'rows':r['rows'],'selected':m['selected'],'H':r['H'],'code':r['code']}
            if c not in best or (v,r['code'],m['selected'])<(best[c]['v2U'],best[c]['code'],best[c]['selected']):best[c]=candidate
    first=next(c for c in range(1,max(best)+3,2) if c not in best);need(first==53559,'complete coefficient prefix boundary')
    spectrum=HERE.parent/'path_spectrum_two_angles_metagraph_20261106/data/least_order_11.txt'
    least={}
    for line in spectrum.read_text().splitlines():
        if line.startswith('#'):continue
        fields=line.split();least[int(fields[0])]=int(fields[1]) if fields[1].isdigit() else None
    heads=[];vh=Counter();witness_checks=0
    for c in range(1,first,2):
        r=best[c];vh[r['v2U']]+=1;sel=[v for v in range(10) if r['selected']>>v&1]
        h,u,w,_=tables(np.array(r['rows'],dtype=np.int64))
        need(int(h)==r['H'] and defect(h,w,sel)==1 and sum(int(u[v]) for v in sel)==r['U'],'every stored seed endpoint certificate');witness_checks+=1
        for s in range(1,r['v2U']+1):
            n=1+c*2**s;need(n in least,'finite head within independently inherited complete spectrum')
            need(least[n] is not None or n in [7,21],'all lower heads covered except proved omissions')
            heads.append({'c':c,'s':s,'H':n,'least_order':least[n]})
    need(max(h['H'] for h in heads)==43131 and len(heads)==3954,'exact lower-head invoice')
    need([h['H'] for h in heads if h['least_order'] is None]==[7,21],'only all-order exceptions in head list')
    # A renewable coefficient1039 seed with a single stable port.
    rows=[4,17,26,19,673,15,63,111,255,495];O=[2,3,4];P=[5]
    D,U=cones.kernel(rows,P,O);need(D=={(1,):1} and exponentials(U)=={(1,):98,(2,):990},'full new 1039 family tensor')
    controls=[]
    for a in range(1,6):
        rr,ss,_=cones.expand(rows,O,P,[a]);h,u,w,_=tables(np.array(rr,dtype=np.int64));slope=98+990*2**a
        need(defect(h,w,ss)==1 and sum(int(u[v]) for v in ss)==slope,'literal renewable family')
        output=tail(rr,ss,1);oh,_,_,_=tables(np.array(output,dtype=np.int64));need(int(oh)==1+2*slope,'literal selected-tail target')
        controls.append({'port_length':a,'c':49+495*2**a,'U':slope,'tail1_H':int(oh)})
    need(least[2079]==9,'1039 below-ray head')
    # Close the whole infinite family's s=1 head by a physical product.
    # At a=1 its factor21 is forbidden, so retain an actual separate witness.
    five=[2,4,9,19,7];ninety_nine=[42,44,9,16,39,76,31]
    exception=[320,25,35,389,141,27,62,103,246]
    fh,fu,fw,_=tables(np.array(five,dtype=np.int64))
    need(defect(fh,fw,[0,1])==1 and int(fu[0]+fu[1])==20,'coefficient5 seed and exact start s=3')
    need(int(tables(np.array(ninety_nine,dtype=np.int64))[0])==99,'literal factor99')
    need(int(tables(np.array(exception,dtype=np.int64))[0])==2079,'literal exceptional head')
    head_controls=[]
    for a in range(2,6):
        second=tail(five,[0,1],a-1);n=len(ninety_nine);m=len(second)
        joined=[r|(((1<<m)-1)<<n) for r in ninety_nine]+[r<<n for r in second]
        target=99*(1+5*2**(a+1));need(target==1+2*(49+495*2**a),'exact head factorization')
        need(int(tables(np.array(joined,dtype=np.int64))[0])==target,'literal product head')
        head_controls.append({'a':a,'H':target,'order':len(joined)})
    out={'scope':'Complete odd-part supply from all ordinary strong marked order-ten unit seeds, minimized by v2(U). Literal endpoint recertification of every prefix witness. All lower dyadic heads checked against the inherited complete ordinary order-eleven spectrum. This gives infinite rays for a finite coefficient band, not Chan-Pak completeness.',
         'unit_input_sha256':hashlib.sha256((HERE/'unit_ten.json.gz').read_bytes()).hexdigest(),
         'spectrum_input':str(spectrum.relative_to(HERE.parents[2])).replace('\\','/'),'spectrum_input_lf_sha256':hashlib.sha256(spectrum.read_text().encode()).hexdigest(),
         'distinct_odd_coefficients':len(best),'first_absent_coefficient_in_declared_bank':first,'covered_odd_prefix_through':first-2,'prefix_seed_certificates':witness_checks,
         'prefix_minimum_valuation_histogram':dict(sorted(vh.items())),'lower_heads':len(heads),'maximum_lower_head':max(h['H'] for h in heads),'head_omissions':[7,21],
         'renewable_1039':{'rows':rows,'selected':28,'port':5,'c_a':'49+495*2^a, a>=1','recurrence':'c_next=2*c-49','source_order':'a+9','U':'2*c','tail_exponents':'s>=2','all_heads_closed':True,'head_factorization':'1+2*c_a=99*(1+5*2^(a+1))','head_factorization_valid_a':'a>=2; a=1 uses separate2079 witness because21 is forbidden','factor99_rows':ninety_nine,'coefficient5_seed_rows':five,'coefficient5_selected':3,'exception2079_rows':exception,'lower_head_H':2079,'lower_head_least_order':9,'literal_controls':controls,'literal_head_controls':head_controls},
         'heads':heads,'records':[best[c] for c in sorted(best)]}
    save('coefficient_band.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k not in ['records','heads']}),flush=True)
if __name__=='__main__':run()
