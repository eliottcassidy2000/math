"""All-order affine TT-chart barriers from small scalar starts and growth rank."""
from common import *
from hamilton_clocks import trim
from math import gcd

def run():
    bank=read('hamilton_clocks_9.json.gz')['records'];lookup={tuple(r['exponential_coefficients']) for r in bank};out=[]
    for bid,p in [(3,0),(4,0)]:
        src=next(r for r in read('even_valuations.json.gz')['records'] if (r['birth_id'],r['port'])==(bid,p))
        for shift in [0,1]:
            e=[c*(i+1)**shift for i,c in enumerate(src['exponential_coefficients'])];e[0]+=1;g=abs(gcd(*e));divisors=[d for d in range(1,g+1) if g%d==0]
            need(value(e,1)<125 and e[2]>0 and not any(e[3:]),'dominant3 and small strong cutoff')
            for d in divisors:need(trim([c//d for c in e]) not in lookup,'every possible fixed SCC scalar factor excluded')
            out.append({'birth_id':bid,'port':p,'original_parameter_shift':shift,'target_profile':e,'first_H':value(e,1),'fixed_SCC_factor_divisors':divisors})
    save('small_head_barriers.json.gz',{'scope':'No fixed affine TT chart, defined at every parameter a>=1 (or every a>=2 for shift1), realizes the specified head family. Uses the proved growth-rank reduction, inherited Busch minimum b(10)=125, and complete ordinary pointed bank through9. Does not exclude later-start charts or individual values.','records':out});print(json.dumps(out),flush=True)
if __name__=='__main__':run()
