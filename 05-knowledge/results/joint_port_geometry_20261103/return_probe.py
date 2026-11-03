"""Exact named source comparison for the 77 -> 233 ternary child."""
from pathlib import Path
from itertools import product,permutations
import sys,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PREV=HERE.parent/'lambda_emergent_20261103';sys.path.insert(0,str(PREV))
from effects import old,junction,need,observer
from controls import values
from one_hole import family

def arcs(T):return [(i,j) for i in range(len(T)) for j in range(i+1,len(T))]
def distance(U,T,mapping):
    return [(mapping[i],mapping[j]) for i,j in arcs(T) if ((U[mapping[i]]>>mapping[j])^(T[i]>>j))&1]
def main():
    w,u=family(2,1);T=old.source(w,u);a,b,ell=5,6,4
    need(junction.profile(T,a)[3]-junction.H(T)==1 and junction.H(T)==77,'actual ready77')
    X=junction.ear(T,1<<a);U=junction.ear(X,1<<7);r=7
    need(junction.profile(U,r)==(233,77,77,233),'actual ternary child233')
    targets=[]
    for m in range(2,8):
        for middle in product((0,1),repeat=m-2):
            ww=(1,)*(m-1)+(0,);uu=(1,)+middle+(0,)
            if old.gap_h(ww,uu)==233:
                S=old.source(ww,uu);targets.append({'w':ww,'u':uu,'rows':S,'profile':junction.profile(S,m)})
    print('targets',[(x['w'],x['u'],x['profile']) for x in targets],flush=True)
    candidates=[];minimum=100;tested=0
    for target in targets:
        S=tuple(target['rows']);m=len(S)-2
        free=[i for i in range(9) if i not in (r,b,ell)]
        for assignment in permutations(free,m-1):
            mapping=assignment+(ell,r,b)
            if a not in mapping:continue
            flips=distance(U,S,mapping);tested+=1
            if len(flips)<minimum:minimum=len(flips);candidates=[]
            if len(flips)==minimum:
                candidates.append({'target':target,'map_target_to_source':mapping,'deleted':[v for v in range(9) if v not in mapping],'flips':flips})
    print('minimum',minimum,'minimizers',len(candidates),'universe',tested,flush=True)
    out={'status':'FINITE_EXACT_DISCOVERY','source77':T,'source_ports':[a,b,ell],'child233':U,'child_ports':[r,b,ell],
         'targets':targets,'tested_maps':tested,'minimum_old_arc_changes_after_deletion':minimum,'minimizers':candidates}
    (HERE/'return_probe.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':main()
