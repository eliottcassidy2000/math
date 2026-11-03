"""Independent forward/backward, literal graph and source-clock controls."""
from pathlib import Path
import argparse,importlib.util,json,sys
from bank import *
from inverse import decode_state,physical_cone,norm,solve_h,source_certificate
HERE=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='verification.json');args=parser.parse_args()
    for l,m in MATS.items():
        for v in ((1,0,0),(1,16,0),(1,0,16),(1,16,16)):
            w=step(v,m)
            need(w[0]>=v[0] and 0<=w[1]<=16*w[0] and 0<=w[2]<=16*w[0],('bounding cone',l,v))
        need(physical_cone(step((1,0,0),m)),('first letter physical',l))
    states,parents,seen=closure(151)
    for v in states:
        word=decode_state(v);need(word is not None and vector(word)==v,('inverse positive',v))
        n=source_certificate(word)['vertices'];need(n==sum(COSTS[l] for l in word),'vertex clock')
        need(5**len(word)*norm(v)>=6*9**len(word),'exponential certificate bound')
    negative=0;small_states=set(closure(15)[0])
    for h in range(1,16,2):
        for x in range(2*h+1):
            for y in range(2*h+1):
                v=(h,x,y)
                if not physical_cone(v):continue
                w=decode_state(v)
                need((w is not None)==(v in small_states),('complete small inverse universe',v))
                negative+=w is None
    false=(37,18,19);need(physical_cone(false) and decode_state(false) is None,'realizable value / unrealizable state')
    scalar={}
    for h in (1,7,21,37,41,55,111,133,491):
        w=solve_h(h)
        need((w is None)==(h in (7,21)),('bounded scalar inverse',h))
        scalar[str(h)]=None if w is None else source_certificate(w)
    incoming=json.loads((HERE.parent/'golden_seam_recurrence_20261103/coverage.json').read_text())
    prime_sources={}
    for entry in incoming['selected_six']['values']:
        if entry['H'] not in (499,1009,19861,21001):continue
        v=tuple(entry['vector']);w=decode_state(v)
        need(w is not None and vector(w)==v,('large actual inverse',v))
        prime_sources[str(v[0])]=source_certificate(w)
    old=HERE.parent/'golden_hamilton_transfer_20261103';sys.path.insert(0,str(old))
    spec=importlib.util.spec_from_file_location('independent_subset_h',old/'transfer.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    literal=[]
    for word in (('r29','b'),('r43','r21'),('r29','r21','r21','b'),('r43','r21','r21','b')):
        rows=actual_target(word);value=mod.subset_h(rows)
        need(value==vector(word)[0],('literal quadratic-family target',word))
        literal.append({'word':word,'vertices':len(rows),'H':value,'rows':rows})
    for entry in prime_sources.values():
        need(mod.subset_h(entry['target_rows'])==entry['H'],('literal prime inverse',entry['H']))
    stats=decode_state.cache_info()
    result={'status':'PASS','complete_forward_states_through151':len(states),
            'rejected_small_physical_cone_states':negative,'false_state_with_allowed_H':false,
            'scalar_inverse_examples':scalar,'prime_source_examples':prime_sources,
            'literal_quadratic_examples':literal,'inverse_cache_states':stats.currsize,
            'scope':'Decision and witness construction for the stated bank, not proof of its scalar surjectivity.'}
    (HERE/args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','complete_forward_states_through151','rejected_small_physical_cone_states','inverse_cache_states')}))
if __name__=='__main__':main()
