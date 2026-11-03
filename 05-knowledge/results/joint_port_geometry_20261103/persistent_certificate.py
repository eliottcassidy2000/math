"""Actual-source Hankel lower bound and modular algebra for the fixed spare."""
import json
import numpy as np
import sympy as s
from return_probe import HERE,junction,observer,need
import effects as e
from algebra import mod_rank
from algebra_probe import algebra

def main():
    p=1009;basis={};witnesses=[];examined=0
    for n in (3,4,5):
        for number,T in enumerate(junction.all_tournaments(n)):
            if not junction.strong(T):continue
            for r in range(n):
                U=junction.ear(junction.ear(T,1<<r),1<<n)
                for b in range(n):
                    if b==r:continue
                    V=observer(U,(n,b));v=[int(V.get(k,0)) for k in e.STATES];row=[x%p for x in v];examined+=1
                    for pivot,base in sorted(basis.items()):
                        c=row[pivot]
                        if c:row=[(x-c*y)%p for x,y in zip(row,base)]
                    pivot=next((i for i,x in enumerate(row) if x),None)
                    if pivot is not None:
                        inv=pow(row[pivot],-1,p);basis[pivot]=[x*inv%p for x in row]
                        witnesses.append({'rows':U,'ports':(n,b),'cap_vector':v,'provenance':[n,number,r,b]})
            if len(basis)==27:break
        if len(basis)==27:break
    d=json.loads((HERE/'persistent_port.json').read_text());W=s.Matrix([x['cap_vector'] for x in witnesses])
    need(len(basis)==27==mod_rank(W*s.Matrix(d['query_basis']),p)==d['observer_rank'],'actual fixed-spare observer is minimal27')
    ops={name:np.array([[int(s.Rational(x).p)*pow(int(s.Rational(x).q),-1,p)%p for x in row] for row in rows],dtype=np.int64) for name,rows in d['matrices'].items()}
    alg=algebra(ops,p);need(alg['algebra_rank']==95,'fixed-spare modular word algebra')
    out={'status':'EXACT_MINIMAL_OBSERVER_WITH_MODULAR_ALGEBRA','observer_rank':27,'examined':examined,
         'witnesses':witnesses,'modular_algebra':alg,'rational_algebra_classification':'OPEN'}
    (HERE/'persistent_certificate.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: actual minimal persistent-spare observer27; modular algebra95',flush=True)
if __name__=='__main__':main()
