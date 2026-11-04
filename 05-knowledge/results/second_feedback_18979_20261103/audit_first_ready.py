"""Independent matrix replay of the additional11735 return."""
import json,sys
from support import HERE,need,save
sys.path.insert(0,str(HERE.parent/'return6325_completion_20261103'))
from audit import birth,native,strong,hp,matrix

def sub(T,V):return tuple(sum(((T[a]>>b)&1)<<j for j,b in enumerate(V)) for a in V)

def main():
    d=json.loads((HERE/'first_ready.json').read_text());A=tuple(d['repaired_seed_rows']);ids=d['clone_addresses']
    T=tuple(sum(((i==j and a<b) or (i!=j and A[i]>>j&1))<<v for v,(j,b) in enumerate(ids)) for i,a in ids)
    need(T==tuple(d['cloned_source']) and hp(T)==3911,'independent literal repaired clone')
    need(sub(birth(A,d['clone_birth_cut']),d['address_to_append'])==T,'one literal birth realizes clone addresses')
    U=birth(birth(T,1<<d['cloned_root']),1<<len(T));need(U==tuple(d['source']) and hp(U)==11735,'independent source reset')
    names=d['source_names'];ns=d['retained_names'];V=[names.index(x) for x in ns];U=sub(U,V);r=d['target_root']
    need(set(names)-set(ns)==set(d['deleted_names'])=={'2','5'},'exact old losses')
    need({'6','21','31','31.1','reset-internal','reset-tip'}<=set(ns),'all named persistent pins')
    need(U==tuple(d['after_deletion']) and strong(U),'strong deletion endpoint')
    for z in d['trace']:
        before=U;X=U
        for cut in z['cuts']:X=birth(X,cut)
        labels=ns+z['helpers'];P,Q=[[labels.index(x) for x in path] for path in z['paths']]
        frame=lambda path:(path[path.index(r)-1] if path.index(r) else None,path[path.index(r)+1] if path.index(r)+1<len(path) else None)
        need(frame(P)==frame(Q),'common actual root frame')
        U=sub(native(native(X,P),Q),list(range(len(ns))))
        diff={(i,j) for i in range(len(ns)) for j in range(i+1,len(ns)) if ((U[i]>>j)^(before[i]>>j))&1}
        need(diff=={tuple(sorted(ns.index(x) for x in z['arc']))} and r not in next(iter(diff)),'independent single exterior arc')
        need(U==tuple(z['rows']) and strong(U),'strong completed repair')
    need(U==tuple(d['target']) and hp(U)==11735,'same-H actual endpoint')
    X=birth(birth(U,1<<r),1<<len(U));need(X==tuple(d['next_rows']) and hp(X)==35207 and strong(X),'independent literal next child')
    save('audit_first_ready.json',{'status':'PASS_INDEPENDENT_MATRIX_AND_HAMILTON_REPLAY','source_H':11735,'target_H':11735,'next_H':35207,
        'old_names_deleted':d['deleted_names'],'native_events':2*len(d['trace']),'completed_repairs_strong':True,'persistent_pins':['6','21','31','31.1','reset-internal','reset-tip'],
        'scope':'actual first stable clone size only; exact old-loss minimum via complete ready fibre, no native-cost minimum or uniform return claim'})
    print('PASS independent11735 native return and next35207',flush=True)
if __name__=='__main__':main()
