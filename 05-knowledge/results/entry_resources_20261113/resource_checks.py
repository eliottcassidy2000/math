"""Expected-defect identity and a literal non-delta-matroid witness."""
from common import *
from math import factorial,comb
from fractions import Fraction
from independent_checks import signed_cut

def inflate(rr,O,ports,sizes):
    labels=[]
    for v in range(len(rr)):labels.extend([v]*(sizes[ports.index(v)] if v in ports else 1))
    rows=[0]*len(labels)
    for i,a in enumerate(labels):
        for j,b in enumerate(labels):
            if (a==b and i<j) or (a!=b and rr[a]>>b&1):rows[i]|=1<<j
    return rows,sum(1<<i for i,v in enumerate(labels) if O>>v&1)

def run():
    import content_checks as cc
    records=load(HERE/'budget_complexes.json.gz')['records'];checks=[]
    for n in range(2,9):
        sums=[0]*(n+1);count=0
        for r in records:
            if r['n']!=n:continue
            rr=cc.decode(n,r['code']);h,u,w,_=tables(np.array(rr,dtype=np.int64));total=sum(int(w[i,j]) for i in range(n) for j in range(i+1,n))
            mass=factorial(n)//r['automorphisms'];count+=mass
            for k in range(n+1):sums[k]+=mass*(comb(n,k)*int(h)-(comb(n-2,k-2)*total if k>=2 else 0))
        EH=Fraction(factorial(n),2**(n-1));labelled=1<<(n*(n-1)//2)
        need(count==labelled,'full labelled universe')
        for k in range(n+1):
            expected=labelled*comb(n,k)*EH*(1-Fraction(2*k*(k-1),n))
            need(sums[k]==expected,'exact average defect by marked-set size')
        checks.append({'n':n,'sum_defect_over_all_labelled_marked_pairs':sums})
    rr=[16,1,67,7,142,31,187,303,127];O=10;ports=[5,8];physical=[]
    for X in range(8):
        a=1+(X&3).bit_count();b=1+int(bool(X&4));rows,mark=inflate(rr,O,ports,[a,b]);d=int(signed_cut(np.array(rows,dtype=np.int64),mark))
        need(d==1-(2**a-2)*(2**b-2),'literal J1 pool inflation')
        physical.append({'optional_mask':X,'lengths':[a,b],'defect':d})
    feasible={r['optional_mask'] for r in physical if r['defect']==1}
    X,Y,e=3,4,4
    need(X in feasible and Y in feasible and all((X^(e|f)) not in feasible for f in [1,2,4]),'symmetric exchange fails for every proposed second toggle')
    save('resource_checks',{'scope':'Exact expected-defect identity checked over all labelled tournaments through8 by complete unmarked orbit masses, and direct signed-cut counts on all8 optional clone subsets of the physical J1 example. All-order statements follow from the separate proofs.','expected_defect_controls':checks,'physical_delta_matroid_hostile':physical,'exchange_witness':{'X':X,'Y':Y,'e':e,'feasible':sorted(feasible)}})
    print('expected-defect levels',len(checks),'J1 optional subsets',len(physical),flush=True)
if __name__=='__main__':run()
