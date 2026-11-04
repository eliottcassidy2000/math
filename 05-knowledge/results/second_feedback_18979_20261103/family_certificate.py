"""Positive exterior-ready families, literal ternary outputs and joint data."""
import json
from support import HERE,need,save,junction,clone,phases,observer
from observer import records,root_profile

def main():
    seed=json.loads((HERE/'exterior_search.json').read_text())['base_targets'][1]
    T=tuple(seed['rows']);r=len(T)-1;out=[]
    for q in (4,7,8,9,12):
        samples=[]
        for m in (1,2,3):
            U,ids=clone(T,q,m);rr=ids.index((r,0));qq=ids.index((q,0));z=junction.profile(U,rr)
            need(junction.strong(U) and z[3]-z[0]==1,'literal exterior-root readiness')
            samples.append({'m':m,'rows':U,'addresses':ids,'live_root':rr,'growth_port':qq,'profile':z,'root_phases':phases(U,rr),'growth_port_phases':phases(U,qq)})
            junction.endpoints.cache_clear()
        a=(samples[1]['profile'][0]-4277)//2;b=4277-2*a
        need(all(z['profile'][0]==a*2**z['m']+b for z in samples),'independent third clone sample')
        out.append({'growth_vertex':q,'coefficients':[a,b],'samples':samples})
        print('positive exterior law',q,a,b,flush=True)
    U=junction.ear(junction.ear(T,1<<r),1<<len(T))
    need(junction.H(U)==12833 and junction.strong(U),'literal strong 3H+2 output')
    joint=observer(T,(r,8));need(root_profile(joint,0)==junction.profile(T,r),'full eighteen-polynomial observer vs root DP')
    save('family_certificate.json',{'status':'FINITE_EXACT_POSITIVE_FAMILIES_WITH_PROVED_TWO_SAMPLE_EXTENSION','seed_rows':T,'seed_root':r,
        'families':out,'seed_joint_ports':[r,8],'seed_full_observer':records(joint),
        'actual_12833':{'rows':U,'birth_cuts':[1<<r,1<<len(T)],'H':12833,'strong':True},
        'scope':'these are explicit new source families, not replacements for the named 18979 source; all-size readiness follows from SF2, not three-sample interpolation'})
    print('PASS positive families, two physical births, full seed observer',flush=True)
if __name__=='__main__':main()
