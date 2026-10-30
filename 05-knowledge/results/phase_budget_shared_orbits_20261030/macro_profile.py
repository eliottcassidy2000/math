#!/usr/bin/env python3
"""All-period address formula versus independent C++ modular orbit enumeration."""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import hashlib,json,subprocess,tempfile
import numpy as np
from collatz_audit import affine,address,v5,need,analyze
root=Path(__file__).parent
U=[3,1,3]+[1]*6+[3,1,3]
R=[5]+[1]*8+[5];W=R+[1]*233+R
reports=[]
with tempfile.TemporaryDirectory(prefix='macro-profile-') as tmp:
    tmp=Path(tmp);exe=tmp/'profile'
    subprocess.run(['c++','-std=c++17','-O2',str(root/'macro_profile.cpp'),'-o',str(exe)],check=True)
    for name,u in (('U3',U),('W5',W)):
        p,q,b=affine(u);Z=b+q-p;h=v5(Z);cap=[];signals=Counter()
        pv=qv=1;bv=0
        for k in range(len(u)+1):
            if k:pv,qv,bv=3*pv,2**u[k-1]*qv,3*bv+qv
            C=p*(bv+qv-pv)+qv*Z
            pp,qq,bb=affine(u[:k]+u)
            need(C==bb+qq-pp and C>0,'entry-preserving divisor independently composed')
            c=v5(C);cap.append(c)
            if k and p*pv*pow(q*qv,-1,5)%5 in (1,4):
                for e,r in enumerate(address(p*pv,q*qv,c),1):signals[(2*5**(e-1),r)]+=1
        K=max(h,*cap);T=2*5**(K-1);L=np.arange(T,dtype=np.int64)
        prefix=analyze(u)['charge'];profile=np.full(T,prefix,dtype=np.int64)
        active=p*pow(q,-1,5)%5 in (1,4)
        gain=F(5*T*(5**h-1),8*5**h) if active else F(0)
        need(gain.denominator==1,'integer period charge')
        if active:
            for e,r in enumerate(address(p,q,h),1):
                t=2*5**(e-1);profile+=np.maximum(0,1+(L-r)//t)-int(r==0)
        for (t,r),multiplicity in signals.items():profile[r::t]+=multiplicity
        out=tmp/f'{name}.bin'
        direct=json.loads(subprocess.check_output([str(exe),name,str(out)]))
        raw=out.read_bytes();direct_profile=np.frombuffer(raw,dtype='<i4')
        need(np.array_equal(profile,direct_profile),'all residues: addresses vs modular trajectory')
        score=T*profile-int(gain)*L;best=int(score.max());where=int(score.argmax())
        need(best==direct['max_numerator'] and where==direct['first_maximizer'],'optimal uniform constant')
        # Actual unbounded integer coefficients at chosen short lengths, including both known hostiles.
        for ell in (0,1,8,231):
            need(int(profile[ell])==analyze(u+[1]*ell+u)['charge'],'literal integer gcd control')
        reports.append(dict(name=name,entry_length=len(u),entry_clock=sum(u),entry_charge=prefix,
                            entry_height=h,K=K,phase_cap_histogram=dict(sorted(Counter(cap).items())),
                            period=T,gain=int(gain),slope=str(gain/T),
                            optimal_additive_constant=str(F(best,T)),
                            optimal_clock_saving=str(F(2*sum(u))*gain/T-F(best,T)),
                            first_maximizer=where,maximizers=int(np.count_nonzero(score==best)),
                            profile_sha256=hashlib.sha256(raw).hexdigest(),checked_residues=T))
        print(json.dumps(reports[-1]),flush=True)
(root/'macro_profile.json').write_text(json.dumps(reports,indent=2)+'\n')
