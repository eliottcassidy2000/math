"""Independent carrier-cost closure and small fixed-height source-face bases."""
import json
from pathlib import Path
from tournament_kernel import *
from audit_constructor_cone import state

def ordinal(a,b):
    n,m=len(a),len(b)
    return tuple(row+(1,)*m for row in a)+tuple((0,)*n+row for row in b)

def run():
    data=json.loads(Path(__file__).with_name('cycles15.normal.json').read_text())
    components={(r['m'],r['H']) for r in data['connected_types']}
    values=[{1}]+[set() for _ in range(15)]
    for cost in range(1,16):
        for m,h in components:
            if m<=cost:values[cost].update(h*x for x in values[cost-m])
    union=set().union(*values);require(len(union)==300,'cycle-budget closure')
    require(7 not in union and 21 not in union,'bounded-budget omission control')
    minimum={h:min(m for m,v in enumerate(values) if h in v) for h in union}
    selected={h:minimum[h] for h in (1,3,5,9,11,13,15,17,19,23,25,27,29,31,35,37,41,53)}
    require(selected=={1:0,3:1,5:2,9:2,11:5,13:6,15:3,17:6,19:9,23:9,25:4,27:3,29:10,31:9,35:13,37:14,41:14,53:14},'minimal cycle costs')
    cores={3:[],5:[]}
    for n,bank in class_bank(5):
        if n<3:continue
        for code in bank:
            a=matrix(n,code);h=hp(a)
            if h in cores and strong(a):cores[h].append(a)
    one=matrix(1,0);faces={1:[1]};raw={}
    for h,bank in cores.items():
        candidates=set()
        for a in bank:
            n=len(a);end=ordinal(a,one)
            cases=[(end,1<<i) for i in range(n)]
            cases.extend([(ordinal(one,a),1),(ordinal(end,one),1<<n)])
            for parent,cut in cases:
                C,L,g,t=state(parent,cut);require(L==g==0 and t==h,'actual face base')
                candidates.add(C)
        raw[h]=sorted(candidates)
        minimal=[]
        for b in sorted(candidates):
            if any(b%a==0 and (b//a)&((b//a)-1)==0 for a in minimal):continue
            minimal.append(b)
        faces[h]=minimal
    require(faces=={1:[1],3:[3],5:[5,6,7,9,16]},'minimal dyadic banks')
    print(json.dumps({'scope':'FINITE-EXACT carrier table cost<=15 closure; fixed-height h3/5 complete strong-core enumeration n<=h and literal marked-base construction',
                      'cycle_budget_value_count':len(union),'selected_global_minimum_cycle_costs':selected,
                      'minimum_dyadic_face_banks':faces,'raw_face_candidates':raw,
                      'values_by_exact_cycle_cost':[sorted(v) for v in values]},sort_keys=True,indent=2))
if __name__=='__main__':run()
