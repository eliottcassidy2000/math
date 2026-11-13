"""Mixed terminal-parent forests and a literal two-ancestry target family."""
from collections import Counter
from functools import cache
from pathlib import Path
import json
import sys
import check as base
from targets import one_apex,induced

@cache
def mixed(n):
    if n in base.SEEDS:return 0,0,((n,1),),('seed',)
    choices=[]
    if base.atoms(n)==[n]:
        digits=(2,) if n%3==2 else (4,-2)
        for d in digits:
            h=(n-d)//3
            if h<3 or h in (7,21):continue
            g,depth,leaves,_=mixed(h)
            choices.append((g+1,depth+1,leaves,('gate',h,d)))
    else:
        for u in base.divisors(n):
            v=n//u
            if not(1<u<=v) or u in (7,21) or v in (7,21):continue
            g,d,l,_=mixed(u);G,D,L,_=mixed(v)
            leaves=Counter(dict(l));leaves.update(dict(L))
            choices.append((g+G,max(d,D),tuple(sorted(leaves.items())),('product',u,v)))
    base.need(choices,'guarded mixed predecessor')
    return min(choices,key=lambda x:(x[0],x[1]))

def main():
    improvements=0;maxsave=(0,None);checked=0
    for n in range(3,base.LIMIT+1,2):
        if n in (7,21):continue
        g,depth,leaves,_=mixed(n);weight=3**g;count=0
        for b,c in leaves:
            base.need(b>1,'unit leaf excluded in shifted potential')
            weight*=(b-1)**c;count+=c
        base.need(weight<=n-1 and 3**g*2**count<=n-1,'shifted total forest budget')
        old=base.optimal(n)[0]
        base.need(g<=old,'additional parent cannot worsen optimum')
        if g<old:
            improvements+=1
            if old-g>maxsave[0]:maxsave=(old-g,n)
        checked+=1
    base.need(mixed(199)[0]==1 and mixed(199)[3]==('gate',67,-2),'formal alternate-parent improvement')
    physical=[]
    for m in range(1,9):
        B=1<<m;V=base.apex('101'+'0'*m+'10');T=base.apex('1001'+'0'*m+'110')
        hv,_,_=one_apex(V,len(V)-1);ht,_,_=one_apex(T,len(T)-1)
        base.need(hv==15*B+1 and ht==45*B+1,'inherited actual negative family')
        mapping=[0,1,3]+list(range(4,m+4))+[m+4,m+6,m+7]
        base.need(len(mapping)==len(V) and all(((V[i]>>j)&1)==((T[mapping[i]]>>mapping[j])&1)
                  for i in range(len(V)) for j in range(len(V))),'literal induced negative embedding')
        expected=[3*B+3,9*B+1,6*B+3]+[15 if m==1 else 15*(B//2)+1]*m+[10*B-1,5*B-1,1]
        actual=[]
        for v in range(len(V)):
            keep=[i for i in range(len(V)) if i!=v];sub=induced(V,keep)
            if v==len(V)-1:count=base.hp(sub)
            else:count=one_apex(sub,keep.index(len(V)-1))[0]
            actual.append(count)
        base.need(actual==expected and max(actual)==10*B-1<15*B-1,'uniform proper-deletion hostile')
        physical.append({'m':m,'canonical_source':15*B-1,'negative_source':hv,'target':ht,
                         'negative_source_rows':V,'target_rows':T,'induced_embedding':mapping,
                         'deletion_counts':actual,'maximum_deletion':max(actual)})
    report={'status':'PASS','scope':'Formal mixed-parent proof forests; actual inherited negative family with new complete deletion barrier',
            'allowed_nonunit_targets':checked,'improved_over_positive_only':improvements,'largest_gate_saving':maxsave,
            'formal_199':mixed(199),'physical_controls':physical}
    suffix='.optimized' if sys.flags.optimize else ''
    (Path(__file__).parent/('alternative_receipt'+suffix+'.json')).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='physical_controls'}))

if __name__=='__main__':main()
