"""Complete one-ready base fibres for the five nearest ternary offsets."""
import sys
from support import HERE,PREV,save,need,junction,clone,matching
sys.path.insert(0,str(PREV))
import moment_search as mom
from apex_carry import B

def main():
    offsets=(-4,-2,0,2,4);targets={4275+d:d for d in offsets};rows=[];total=0
    need(B(3729)==-1999 and all(B(3729)!=-3183+d for d in offsets),'clone-mark constant barrier')
    for n in range(3,16):
        for j in range(mom.COUNTS[n]):
            total+=1;h=mom.moment(n,j,1,0);F=mom.moment(n,j,2,1);H=h+F-1
            if H not in targets:continue
            d=targets[H];P=mom.unrank(n,j);T=junction.ear(P,1);r=n
            z=junction.profile(T,r);need(z[0]==H and z[3]-z[0]==1,'literal ready base')
            queries=[]
            for q in range(1,n):
                M,C=matching(T,q)
                if len(M)!=1:continue
                U,ids=clone(T,q,2);rr=ids.index((r,0));zz=junction.profile(U,rr)
                a=(zz[0]-H)//2;b=2*H-zz[0]
                queries.append({'q':q,'H2':zz[0],'defect2':zz[3]-zz[0],'coefficients':[a,b],
                    'exact_target':a==3729 and b==-3183+d and zz[3]-zz[0]==1,'matching':M,'cover':C})
                junction.endpoints.cache_clear()
            rows.append({'offset':d,'parent_n':n,'parent_rank':j,'rows':T,'root':r,'profile':z,'queries':queries})
    save('offset_fibres.json',{'status':'FINITE_EXACT_COMPLETE_ONE_READY_FIBRES_USING_CITED_ORDER_BOUND',
        'source_coefficients':[1243,-1061],'clone_mark_B3729':B(3729),'offsets':offsets,'parent_classes_examined':total,'fibres':rows,
        'uniform_hits':sum(q['exact_target'] for z in rows for q in z['queries']),
        'scope':'strong one-ready target at m=1; all one-start parents3..15; every matching-one exterior clone except the unique in-port, which cannot remain ready on cloning; clone-root cases handled by RC1/RC4'})
    print('PASS offset fibres',total,'parents',len(rows),'targets',sum(len(z['queries']) for z in rows),'clone queries',flush=True)
    print('hits',[(z['offset'],z['parent_n'],z['parent_rank'],q['q']) for z in rows for q in z['queries'] if q['exact_target']],flush=True)
if __name__=='__main__':main()
