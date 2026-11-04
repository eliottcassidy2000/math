"""Replay the named 18979 returns and retain its full inherited joint observer."""
import sys,json,importlib.util
from support import HERE,PREV,ROOT,INCOMING,need,save,source,junction,phases
sys.path.insert(0,str(HERE.parent/'return6325_completion_20261103'))
from audit import matrix,birth,native,strong
sys.path.insert(0,str(HERE.parent/'digit_cycle_ternary_20261103'))
from hostile_embedding import tournament
from observer import root_profile
from collections import Counter

def induced(T,V):return tuple(sum(((T[a]>>b)&1)<<j for j,b in enumerate(V)) for a in V)
def triples(T):
    from itertools import combinations
    return [c for c in combinations(range(len(T)),3) if sum(bool(T[a]>>b&1) for a,b in ((c[0],c[1]),(c[1],c[2]),(c[2],c[0]))) in (0,3)]

def main():
    T,r,b,names=source();inherited=json.loads((PREV/'hostiles.json').read_text())['actual18979']
    need(tuple(inherited['rows'])==T and inherited['names']==names,'exact source and names, not same-H substitution')
    need(junction.profile(T,r)==(18979,6325,6325,18979),'fresh source profile')
    need(json.loads(json.dumps(phases(T,r)))==inherited['root_phases'] and json.loads(json.dumps(phases(T,b)))==inherited['spare_phases'],'fresh five Q4 and F49/charge phases')
    V=Counter({(tuple(tuple(p) for p in z['skeleton']),z['degree']):z['count'] for z in inherited['observer']})
    need(root_profile(V,0)==junction.profile(T,r) and root_profile(V,1)==junction.profile(T,b),'inherited full observer retained consistently')
    returns=[]
    for z in json.loads((INCOMING/'returns.json').read_text()):
        if z['source_H']!=18979:continue
        need(tuple(z['source'])==T and z['source_names']==names,'incoming exact named 18979 source')
        ns=z['names'];rr=z['root'];U=induced(T,[names.index(x) for x in ns])
        need(set(names)-set(ns)==set(z['deleted_old_names']) and U==tuple(z['after_deletion']),'explicit old-name loss')
        need(ns[rr]==43 and ns[z['spare']]==6 and 31 in ns,'retained designated pins')
        for step in z['trace']:
            before=U;local=ns+[step['helper']];P,Q=[[local.index(a) for a in path] for path in step['paths']]
            frame=lambda path:(path[path.index(rr)-1] if path.index(rr) else None,path[path.index(rr)+1] if path.index(rr)+1<len(path) else None)
            need(frame(P)==frame(Q),'native common root frame')
            U=induced(native(native(birth(U,step['cut']),P),Q),list(range(len(ns))))
            diff={(i,j) for i in range(len(ns)) for j in range(i+1,len(ns)) if ((U[i]>>j)^(before[i]>>j))&1}
            need(diff=={tuple(sorted(ns.index(a) for a in step['arc']))} and U==tuple(step['rows']) and strong(U),'literal completed native repair')
        zz=junction.profile(U,rr);need(zz[0]==18979 and zz[3]-zz[0]==z['budget'],'returned source readiness')
        X=birth(birth(U,1<<rr),1<<len(U));need(X==tuple(z['next_rows']) and junction.H(X)==3*18979+2*z['budget'],'literal ternary continuation')
        returns.append({'budget':z['budget'],'deleted_names':z['deleted_old_names'],'native_events':2*len(z['trace']),
            'rows':U,'names':ns,'root':rr,'spare':z['spare'],'profile':zz,'root_phases':phases(U,rr),'spare_phases':phases(U,z['spare']),
            'next_rows':X,'next_names':z['next_names'],'next_root':z['next_root'],'next_H':junction.H(X)})
        junction.endpoints.cache_clear()
    path=ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py';spec=importlib.util.spec_from_file_location('fresh81',path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    A,B=tournament(False),tournament(True);aut=[len(mod.iso_maps(X,X)) for X in (A,B)]
    need(aut==[81,27] and triples(A)==triples(B),'full 81/27 symmetry control')
    states=json.loads((HERE.parent/'return6325_observer_20261103/completion_hostile.json').read_text())['states'];responses=[]
    for rows in states:
        X=tuple(rows);need(junction.H(X)==462915,'same-H hostile sources');responses.append(junction.H(junction.ear(X,291)));junction.endpoints.cache_clear()
    need(responses==[2450493,2421909,2665485,2343867,2753631,2633253,2329119,2772027],'actual port distinctions')
    save('source_pins.json',{'status':'FINITE_EXACT_FRESH_NATIVE_REPLAY_AND_PINNED_FULL_OBSERVER','actual18979':inherited,'returns':returns,
        'hostile':{'automorphism_orders':aut,'common_labelled_triples':len(triples(A)),'named_cut':291,'responses':responses},
        'clocks':'native event counts and all five Q4/charge phases; no invented interval ancestry or free re-marking',
        'observer_scope':'full source observer inherited unchanged and hash pinned; returned targets retain literal relations and two root phases, not a claimed new minimal joint quotient'})
    print('PASS actual 18979 returns, sources/ports/clocks and 81/27',flush=True)
if __name__=='__main__':main()
