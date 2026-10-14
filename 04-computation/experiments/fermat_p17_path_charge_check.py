#!/usr/bin/env python3
"""Exact all-class p=17 witness verification, using the published integral checker.
Standard library only. Optional argument: certificate JSON path.
"""
import importlib.util,json,math,sys
from pathlib import Path
HERE=Path(__file__).resolve()
helper=HERE.with_name('path_charge_code_obstruction_check.py')
if not helper.exists():helper=HERE.parent.parent/'fermat_path_charge_stage'/'path_charge_code_obstruction_check.py'
spec=importlib.util.spec_from_file_location('charge_integral_verifier',helper)
core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)

def main():
    name='fermat_p17_path_charge_certificates_20261014.json'
    path=Path(sys.argv[1])if len(sys.argv)>1 else HERE.with_name(name)
    if not path.exists():path=HERE.parents[2]/'05-knowledge'/'results'/name
    data=json.loads(path.read_text());rows=data['representatives'];observed=set();summaries=[]
    core.require(data['scope_primes']==[17]and len(rows)==16,'exact p17 scope')
    for row in rows:
        summary=core.verify_ring(row,data['alphabet']);p=row['p'];d=row['norm']
        core.require(p==17 and row['surjective'],'all representatives positively witnessed')
        orbit={frozenset(a*s%p for s in row['S'])for a in range(1,p)}
        core.require(len(orbit)==16 and observed.isdisjoint(orbit),'distinct full multiplier orbits');observed|=orbit
        if d>1:
            core.require(row['active_coordinates']==[0]and row['moduli']==[d]and row['zeta_coordinates'][0]==[1],'actual quotient is Z/dZ')
            z=row['zeta_coordinates'][1][0];factors=core.prime_divisors(d)
            core.require(pow(z,17,d)==1 and (1+sum(pow(z,s,d)for s in row['S']))%d==0,'integer character realizes quotient ring')
            core.require(math.prod(factors)==d,'all observed quotients are reduced')
            for ell in factors:core.require(ell%34==1 and core.order(z,ell)==17,'degree-one Fermat residue factors')
            summary.update({'zeta_image':z,'prime_factors':factors})
        else:summary.update({'zeta_image':0,'prime_factors':[]})
        summaries.append(summary)
    expected={frozenset(j+1 if mask>>j&1 else 16-j for j in range(8))for mask in range(256)}
    core.require(observed==expected,'all 256 labelled connection sets')
    core.require(sum(r['norm']for r in rows)==48160,'total residue witnesses')
    composite=[r for r in summaries if len(r['prime_factors'])>1]
    core.require(len(composite)==1 and composite[0]['norm']==24617 and composite[0]['prime_factors']==[103,239]and composite[0]['zeta_image']==1535,'exact composite product control')
    print(json.dumps({'status':'PASS','p':17,'multiplier_classes':16,'labelled_connection_sets':256,'explicit_residue_witnesses':48160,'all_classes_surjective':True,'representatives':summaries,'active_checks':core.CHECKS},indent=2,sort_keys=True))
if __name__=='__main__':main()
