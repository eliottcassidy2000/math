"""Fresh source, symmetry, joint-cut and clock controls at the new checkpoint."""
import json,sys,importlib.util
from support import HERE,PREV,ROOT,need,save,source,junction,phases,observer
from observer import records
sys.path.insert(0,str(HERE.parent/'digit_cycle_ternary_20261103'))
from hostile_embedding import tournament

def triples(T):
    from itertools import combinations
    return [c for c in combinations(range(len(T)),3) if sum(bool(T[a]>>b&1) for a,b in ((c[0],c[1]),(c[1],c[2]),(c[2],c[0]))) in (0,3)]

def main():
    d=json.loads((PREV/'completion_hostile.json').read_text());states=[tuple(T) for T in d['states']]
    A,B=tournament(False),tournament(True)
    path=ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py';spec=importlib.util.spec_from_file_location('fresh81',path);v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
    aut=[len(v.iso_maps(T,T)) for T in (A,B)];need(aut==[81,27] and triples(A)==triples(B),'actual full hostile symmetry and labelled triple set')
    responses=[]
    for T in states:
        need(junction.H(T)==462915,'same scalar baseline')
        U=junction.ear(T,291);need(junction.strong(U),'actual guarded read-head source')
        responses.append(junction.H(U));junction.endpoints.cache_clear()
    need(responses==[2450493,2421909,2665485,2343867,2753631,2633253,2329119,2772027] and len(set(responses))==8,'source-separating port context retained')
    T,r,b,names=source();z=junction.profile(T,r)
    need(z==(18979,6325,6325,18979),'no accidental change to requested source')
    save('hostiles.json',{'status':'FINITE_EXACT_FRESH_REPLAY','automorphism_orders':aut,'common_triples':len(triples(A)),'eight_source_readhead_cut':291,'readhead_responses':responses,
      'actual18979':{'rows':T,'names':names,'root':r,'spare':b,'profile':z,'root_phases':phases(T,r),'spare_phases':phases(T,b),'observer':records(observer(T,(r,b)))},
      'failed55_control':'inherited pinned one_birth55_orientations; not rerun or promoted to a broader impossibility',
      'scope':'no symmetry or observer quotient of the actual source is used as a replacement tournament'})
    print('PASS hostile81/27, eight-source port selector and actual18979 clocks',flush=True)
if __name__=='__main__':main()
