"""Odd-slope axes in all stable states through10, including descendants."""
from common import *
from affine_births import onto
from collections import Counter
from hamilton_clocks import trim
from fractions import Fraction
import time

def run():
    old=json.loads((HERE.parent/'five_stress_connection_atlas_20261112/cyclic_port_receipt.json').read_text())['positive_screen_certificates']
    sources=[{'rows':r['rows'],'selected':r['selected'],'stable_ports':r['individually_stable']} for r in old if r['individually_stable']]
    sources+=read('stable_nine.json.gz',True)['records']
    sources+=read('stable_ten_selected.json.gz',True)['records']
    out=[];stats=Counter();start=time.time()
    for idx,r in enumerate(sources):
        rows=r['rows'];n=len(rows);O=r['selected'];h,u,w,_=tables(np.array(rows,dtype=np.int64));U=sum(int(u[v]) for v in range(n) if O>>v&1)
        if U%2==0:continue
        _,_,comps=quotient(rows,O);ports=[next(v for v in c if v in r['stable_ports']) for c in comps if set(c)&set(r['stable_ports'])]
        stats['odd_sources']+=1
        for p in ports:
            R=min(n-rows[p].bit_count(),rows[p].bit_count()+2);values=[];coeff=[]
            for a in range(1,R+1):
                labs,rr=ib.inflate(rows,p,a);hh,uu,ww,_=tables(np.array(rr,dtype=np.int64));ss=[v for v,q in enumerate(labs) if O>>q&1]
                need(defect(hh,ww,ss)==1,'inherited stable axis replay');val=sum(int(uu[v]) for v in ss);values.append(val)
                rem=val-sum(coeff[j-1]*onto(a,j) for j in range(1,a));den=onto(a,a);need(rem%den==0,'integer full U coefficients');coeff.append(rem//den)
            need(all(c>=0 for c in coeff) and coeff[0]%2==1,'positive odd full slope kernel')
            out.append({'rows':rows,'selected':O,'port':p,'R':R,'values':values,'Stirling_coefficients':coeff,'exponential_coefficients':list(trim(expcoeff(coeff)))})
        if idx%100==0:print(json.dumps({'sources':idx,'axes':len(out),'seconds':round(time.time()-start,1)}),flush=True)
    lookup={tuple(r['exponential_coefficients']):i for i,r in enumerate(out)};matches=[]
    for r in read('even_valuations.json.gz')['records']:
        if r['type']!='bounded_constant':continue
        j=r['eventual_valuations_by_parity'][0];e=r['exponential_coefficients']
        for shift in range(-6,7):
            f=[Fraction(c,2**j)*Fraction(i+1)**shift for i,c in enumerate(e)]
            if any(c.denominator!=1 for c in f):continue
            ix=lookup.get(trim([int(x) for x in f]))
            if ix is not None:matches.append({'birth_id':r['birth_id'],'port':r['port'],'division_valuation':j,'original_parameter_shift':shift,'odd_axis_index':ix})
    result={'scope':'Every strong marked stable state through10, odd U only, one representative per stable twin block. Full exact single-axis slope coefficients from the proved selected-start run bound. Matches test shifts-6..6 for eventual constant-valuation even families; finite prefix obligations are separate.','source_states':len(sources),'odd_sources':stats['odd_sources'],'axes':len(out),'distinct_profiles':len(lookup),'matches':matches,'records':out}
    save('odd_axis_bank.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
