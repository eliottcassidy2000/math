"""Sharp cyclic drift controls for the inherited binary address recurrence."""
from pathlib import Path
from fractions import Fraction
import sys,json
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'changing_core_returns_20261113'/'returns'))
import inverse as old

def bits(r,L):return tuple((r>>i)&1 for i in range(L))

def data(r,L):
    b=bits(r,L);Q=(1<<L)-1
    rev=sum(b[i]<<(L-1-i) for i in range(L))
    a=Q-rev;beta=r&1;firstzero=1-b[-1];c=old.g(r)-1
    D=Fraction(firstzero*beta+c)-Fraction(a*r,Q)
    Z=tuple(sum((1-b[i])*b[(i+d)%L] for i in range(L)) for d in range(1,L))
    energy=Fraction(Z[0] if Z else 0)-Fraction(sum((1<<i)*z for i,z in enumerate(Z)),Q)
    return dict(D=D,energy=energy,Z=Z,a=a,b=firstzero,c=c,beta=beta)

def main():
    levels=[];count=0;minimum_controls=0
    for L in range(1,15):
        Q=(1<<L)-1;lower=Fraction((1<<(L//2))+(1<<((L+1)//2))-2,Q)
        minimum=None;minimizers=[]
        for r in range(1<<L):
            z=data(r,L);D=z['D'];count+=1
            old.need(D==z['energy'],'inherited drift versus cyclic correlations')
            old.need((D==0)==(r in (0,Q)),'zero drift iff constant period')
            if r not in (0,Q):
                old.need(D>=lower,'sharp mixed-period bound')
                if minimum is None or D<minimum:minimum=D;minimizers=[r]
                elif D==minimum:minimizers.append(r)
                Z1=z['Z'][0]
                old.need(all(v<=min(d,L-d)*Z1 for d,v in enumerate(z['Z'],1)),'cyclic transition bound')
                old.need(D<=Z1,'drift run-count upper bound')
            # Rotation, reversal, and complement preserve the actual drift.
            rot=(r>>1)|((r&1)<<(L-1))
            rev=sum(v<<(L-1-i) for i,v in enumerate(bits(r,L)))
            old.need(data(rot,L)['D']==D and data(rev,L)['D']==D and data(Q-r,L)['D']==D,'dihedral and complement symmetry')
            if L<=8:
                for prefix in (0,1,2,5):
                    old.need(old.g((prefix<<L)+r)==old.g(prefix)+z['a']*prefix+z['b']*(prefix&1)+z['c'],'literal suffix law')
        if L>=2:
            old.need(minimum==lower,'attained sharp lower bound')
            for r in minimizers:
                b=bits(r,L);runs=sum((1-b[i])*b[(i+1)%L] for i in range(L))
                old.need(runs==1 and sum(b) in (L//2,(L+1)//2),'all equality cases')
            minimum_controls+=len(minimizers)
        levels.append(dict(length=L,mixed_minimum=None if minimum is None else str(minimum),minimizers=minimizers))
    # Actual ready tournaments on the mixed 01-period address family.
    core=old.cores(3)[1];n=1;physical=[]
    for j in range(1,5):
        n=4*n+1;T,y,V,x,w=old.physical(core,n)
        expected=(64*4**j+6*j+23)//3
        old.need((64*4**j+6*j+23)%3==0,'integral physical mixed ray')
        old.need(old.direct(T)==expected and old.direct(V,x)==1,'literal mixed-period H and unit start')
        old.need(old.direct(V)-old.direct(T,y)==1,'literal mixed-period ready defect')
        physical.append(dict(j=j,address=n,word=''.join(map(str,w)),H=expected,rows=T,root=y,predecessor=x))
    dataout=dict(status='PASS',all_word_controls=count,max_length=14,equality_controls=minimum_controls,levels=levels,physical_mixed_ray=physical)
    suffix='' if __debug__ else '.optimized'
    (HERE/('cyclic_drift_receipt'+suffix+'.json')).write_text(json.dumps(dataout,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS',count,'periods and four literal unit-ready mixed-ray graphs')

if __name__=='__main__':main()
