"""Actual81/27 replay and independent endpoint counts for the701 return."""
import importlib.util,json
from common import HERE,ROOT,need,junction,values,save

def main():
    path=ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
    spec=importlib.util.spec_from_file_location('regular_parent_hostile',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    control=module.hostile_checks();checks=[];d=json.loads((HERE/'native_return.json').read_text())
    for stage,T,r in [('source',d['source_rows'],d['input_root']),('returned',d['final_rows'],d['root']),('next',d['next_rows'],d['next_root'])]:
        v,p=values(tuple(T),r);need(p==junction.profile(tuple(T),r),'independent endpoint-DP source profile')
        checks.append({'stage':stage,'profile':p})
    save('hostiles.json',{'status':'FINITE_EXACT_INDEPENDENT_CONTROLS','actual81_27':control,'return_profiles':checks})
    print('PASS: actual81/27 groups,67 equal cyclic triples, regular-mask separation; independent return profiles',flush=True)
if __name__=='__main__':main()
