"""All odd-slope higher-mode birth clocks have no below-ray head debt."""
import json
from census_nine import read,save,need

def run():
    allr=read('affine_births.json.gz')['records'];higher=[r for r in allr if any(r['Stirling_coefficients'][2:])]
    odd=[r for r in higher if r['values'][0]%2];even=[r for r in higher if not r['values'][0]%2]
    need((len(higher),len(odd),len(even))==(798,320,478),'exact mode-parity classification')
    need(len({tuple(r['Stirling_coefficients']) for r in odd})==320,'distinct full coefficient sequences')
    # r!S(a,r) is even for every r>=2, so the first coefficient controls
    # parity for every positive a. This uses the complete certified tensor.
    for r in odd:need(r['Stirling_coefficients'][0]%2==1,'odd constant parity mode')
    constant=[]
    for r in allr:
        c=r['Stirling_coefficients'];q=c[1] if len(c)>1 else 0
        if not any(c[2:]) and (q==0 or c[0]==2*q):
            oddpart=c[0]//(c[0]&-c[0]);need(oddpart<=53557,'constant coefficient in complete band');constant.append({'birth_id':r['birth_id'],'port':r['port'],'c':oddpart})
    out={'scope':'Every higher-mode birth clock through10 whose base U is odd. Complete positive Stirling tensors prove U(a) odd at all positive lengths; the physical tail therefore realizes every1+U(a)*2^s for a,s>=1. No unproved assertion about dyadic valuations of the even-slope higher-mode clocks.','higher_mode_axes':len(higher),'complete_odd_slope_families':len(odd),'even_slope_higher_mode_axes_remaining':len(even),'constant_odd_part_axes':len(constant),'constant_coefficients':sorted({r['c'] for r in constant}),'records':odd,'constant_axes':constant}
    save('higher_odd_families.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k not in ['records','constant_axes']}),flush=True)
if __name__=='__main__':run()
