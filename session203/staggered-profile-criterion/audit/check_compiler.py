#!/usr/bin/env python3
"""Exercise the frozen actual compiler, with independent graph checks."""
from pathlib import Path
import importlib.util,json
from check_certificates import need,valid,literal

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('reviewed_compiler',HERE/'REVIEWED_COMPILER.py')
C=importlib.util.module_from_spec(spec);spec.loader.exec_module(C)

def expect_failure(target,seeds,profiles,needle):
    try:C.compile_target(target,seeds,profiles)
    except ValueError as exc:
        need(needle in str(exc),'wrong rejection: '+str(exc));return str(exc)
    raise RuntimeError('invalid case compiled')

def main():
    seeds,profiles=C.read_data(HERE/'data');outputs=[];calls=0
    for target in (23,409,411,413,415,1239,2995):
        a,trace=C.compile_target(target,seeds,profiles);valid(a)
        h,_,_,work=literal(a);calls+=work
        need(h==target,'compiled graph literal count')
        outputs.append({'target':target,'order':len(a),'trace':trace})
    failures={}
    failures['missing_3001']=expect_failure(3001,seeds,profiles,'MISSING_OPEN_PROFILE_CERTIFICATE: B=999, residue=2')
    failures['missing_2999']=expect_failure(2999,seeds,profiles,'MISSING_OPEN_PROFILE_CERTIFICATE: B=999, residue=1')
    failures['even_target']=expect_failure(24,seeds,profiles,'target must be odd')
    failures['below_floor']=expect_failure(21,seeds,profiles,'target must be odd')
    bad=dict(seeds);bad[23]=[0]
    failures['wrong_seed_count']=expect_failure(23,bad,profiles,'seed path count failed')
    j,root,u=profiles[137,1]
    bad=dict(profiles);bad[137,1]=(9,root,u)
    failures['j_out_of_scope']=expect_failure(413,seeds,bad,'profile parameter outside scope')
    bad=dict(profiles);bad[137,1]=(j,len(u),u)
    failures['root_out_of_scope']=expect_failure(413,seeds,bad,'profile parameter outside scope')
    bad=dict(profiles);bad[137,1]=(0,0,[0])
    failures['wrong_profile']=expect_failure(413,seeds,bad,'profile equations failed')
    bad=dict(profiles);bad[137,1]=(0,0,[0,0])
    failures['invalid_tournament']=expect_failure(413,seeds,bad,'invalid pair')
    # A >20-vertex output is built by certified products, not passed to literal DP.
    recorded=[];original_count=C.count
    def guarded_count(a):
        recorded.append(len(a));need(len(a)<=20,'large graph falsely passed to DP')
        return original_count(a)
    C.count=guarded_count
    target=413*3**8;a,trace=C.compile_target(target,seeds,profiles);valid(a)
    base_order=next(item['order'] for item in outputs if item['target']==413)
    need(len(a)==base_order+24 and len(trace)==9 and recorded==[base_order],'large join scope')
    remainder=a
    for _ in range(8):
        n=len(remainder)-3;suffix=((1<<n)-1)<<3
        need(remainder[:3]==[2|suffix,4|suffix,1|suffix],'triangle product edges')
        need(all(not row&7 for row in remainder[3:]),'wrong cross edges')
        remainder=[row>>3 for row in remainder[3:]]
    h,_,_,work=literal(remainder);calls+=work;need(h==413,'large join base')
    print(json.dumps({'status':'PASS','literal_compiled_examples':outputs,'rejection_controls':failures,
                      'large_join':{'target':target,'order':len(a),'triangle_factors':8,'count_DP_orders_called':recorded,'base_literal_H':h,'count_method':'product bijection after independently checked factor'},
                      'literal_DFS_visits':calls,'all_order_profile_hypothesis':'OPEN'},indent=2,sort_keys=True))

if __name__=='__main__':main()
