"""Algebra controls for Seed-2, not a tournament seed census."""
from pathlib import Path
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
checks=0
for h in range(1,32):
    for w in range(h+6):
        for u in range(1,32):
            a1=h-w+2*u;a2=h-w+4*u
            if min(a1,a2)<0:continue
            for c in range(1,32,2):
                for j in range(5):
                    m=c*2**j;cap=2+4*m
                    positive=(min(cap,a1)==1+2*m and min(cap,a2)==1+4*m)
                    seed=(h-w==1 and u==m)
                    need(positive==seed,'joint positive seed equivalence');checks+=1
h,w,u=1000,999,7
cap=30
need(min(cap,h)==min(cap,w),'signed overflow hostile control')
need((h-w+2*u,h-w+4*u)==(15,29),'positive controls retain cancelled offset')
data={'status':'PASS','algebra_cases':checks,'scope':'integer affine-law controls only, not a physical seed census',
      'cancellation_control':{'H':h,'W':w,'U':u,'positive_counts':[15,29],'cap':cap}}
Path(__file__).with_name('seed_probe_receipt.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(data))
