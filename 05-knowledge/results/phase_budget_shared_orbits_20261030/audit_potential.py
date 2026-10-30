#!/usr/bin/env python3
"""Independent vectorized audit of exported little-endian integer potentials."""
import hashlib,json,subprocess,tempfile
from pathlib import Path
import numpy as np

def need(b,msg):
    if not b:raise RuntimeError(msg)

root=Path(__file__).parent
reports=[]
with tempfile.TemporaryDirectory(prefix='phase-budget-') as tmp:
    tmp=Path(tmp);exe=tmp/'producer'
    subprocess.run(['c++','-O2','-std=c++17',str(root/'epoch_potential.cpp'),'-o',str(exe)],check=True)
    for E in range(1,9):
        out=tmp/f'V{E}.bin'
        receipt=json.loads(subprocess.check_output([str(exe),str(E),str(out)]))
        q=5**E;scale=8*q;cost=5*q;cap=receipt['cap']
        raw=out.read_bytes();V=np.frombuffer(raw,dtype='<i8');s=np.arange(2*q,dtype=np.int64)
        z=s//2;active=s%2;h=np.zeros(q,dtype=np.int64)
        for e in range(1,E+1):h[np.arange(q)%(5**e)==0]+=1
        need(len(V)==2*q and int(V.max())<=0,'finite nonpositive vector')
        powers=np.array([5**i for i in range(E+1)],dtype=np.int64)
        zero=0
        for a in range(1,cap+1):
            two=pow(2,a,q);inv=pow(two,-1,q)
            t=2*((3*z+two-2)*inv%q)+(active^(a%2==0))
            if a==1:delta=np.where(active!=0,-cost//powers[h[z]],-cost)
            else:
                k=sum((2**a-2)%(5**e)==0 for e in range(1,E+1))
                delta=scale*(t%2)*np.maximum(h[t//2],k)-cost*a
                seed=2*((two-2)*inv%q)+(1^(a%2==0))
                need(int(V[seed])>=int(delta[1]),'first reset seed')
            slack=V[t]-V-delta
            need(int(slack.min())>=0,'every exact transition inequality')
            zero+=int(np.count_nonzero(slack==0))
        need(cost*(cap+1)>=scale*E-int(V.min()),'all larger letters and seeds')
        receipt.update(sha256=hashlib.sha256(raw).hexdigest(),audited_edges=2*q*cap,
                       zero_slack_edges=zero,unbounded_letter_tail=True)
        reports.append(receipt)
        print('audited',E,2*q*cap,flush=True)
(root/'potential_audit.json').write_text(json.dumps(reports,indent=2)+'\n')
