"""A precise completed-library miss whose first count is nevertheless realized."""
from pathlib import Path
import importlib.util,json
import numpy as np
from census_nine import read,save,need,tables,defect
from affine_births import onto
from stable_census import certify
HERE=Path(__file__).resolve().parent

def run():
    c=53559;need(c not in {r['c'] for r in read('coefficient_band.json.gz')['records']},'absent from complete order10 strong coefficient bank')
    binary=read('completed_bulk_heads.json.gz')['closed_laws'];higher=read('higher_odd_families.json.gz')['records'];checks=[]
    for A,B in binary:
        numerator=c-B
        if numerator%A==0:
            q=numerator//A;need(not(q>=2 and not q&(q-1)),'no positive-parameter binary representation')
    for r in higher:
        a=1
        while True:
            value=sum(v*onto(a,j) for j,v in enumerate(r['Stirling_coefficients'],1));need(value!=c,'no higher odd-mode representation')
            if value>c:break
            a+=1
        # Positive Stirling coefficients make U nondecreasing, so the
        # first value above c ends an exact all-length membership test.
        checks.append({'birth_id':r['birth_id'],'port':r['port'],'first_length_above_target':a,'U':value})
    path=HERE.parent/'session251_creative_completion_20261111/ternary/compiler.py';sp=importlib.util.spec_from_file_location('finite_head_compiler',path);cc=importlib.util.module_from_spec(sp);sp.loader.exec_module(cc)
    core=[2,4,9,19,7];O=3;A=2728;B=0;sig=cc.signature(core,O);rr=cc.construct(core,O,A,B);H=1+2*c
    need(cc.count(sig,A,B)==H and int(tables(np.array(rr,dtype=np.int64))[0])==H,'actual missing-address count realized by another grammar')
    old=[22,36,8,19,38,77,159,63];h,u,w,_=tables(np.array(old,dtype=np.int64));need(int(h)==137 and defect(h,w,[4,7])==1 and int(u[4]+u[7])==1039,'inherited coefficient1039 seed provenance')
    need(not certify(old,[144])[0][144],'inherited1039 seed has no stable port')
    out={'scope':'53559 is the first odd coefficient outside the completed sublibrary consisting of all strong order10 unit-seed odd parts, all540 completed binary laws, and all320 higher-mode odd-slope laws. This excludes none of the other possible seeds, the478 even higher-mode axes, or arbitrary tournaments. Its first desired count107119 is explicitly realized by an inherited different grammar.',
         'coefficient':c,'binary_laws_checked':len(binary),'higher_laws_checked':len(higher),'higher_stopping_certificates':checks,
         'actual_count_witness':{'H':H,'core':core,'cut':O,'A':A,'B':B,'signature':sig,'rows':rr,'order':len(rr)},
         'inherited_1039_seed':{'rows':old,'selected':144,'H':137,'defect':1,'U':1039,'stable_ports':[]}}
    save('library_boundary.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k!='higher_stopping_certificates'}),flush=True)
if __name__=='__main__':run()
