#!/usr/bin/env python3
"""Mutation controls: each false repair/data claim must fail, including under -O."""
import json,pathlib,subprocess,sys,tempfile
root=pathlib.Path(__file__).resolve().parent
source=(root/'check_integral_order.py').read_text()
cases=[
 ('wrong_generator','A=[0,1,1,0,0]','A=[0,1,0,1,0]','cyclic power-basis columns'),
 ('wrong_inverse','N=[[11,-17,17,-7,-3]','N=[[11,-16,17,-7,-3]','two-sided inverse numerator'),
 ('wrong_congruence','ell=[0,1,-1,3,6]','ell=[0,1,-1,3,5]','all powers satisfy congruence'),
 ('wrong_field_inverse',"lin([-2,-3,-2,-1],field_cols[:4])==z","lin([-1,-3,-2,-1],field_cols[:4])==z",'integral cyclotomic inverse'),
 ('wrong_new_prime','field_norm([2,1,0,0])==11','field_norm([2,1,0,0])==13','new gluing norm11'),
 ('wrong_hadamard_sign','matrixA[i][j]-had[i][j]','matrixA[i][j]+had[i][j]','Hadamard repair recovers'),
 ('growing_prime_power_deficit','        w=[modulus//11*x%modulus for x in v]','        g=121\n        w=[modulus//11*x%modulus for x in v]','kernel order remains11'),
 ('false_mod5_spectral_heuristic','rank_mod(M,5)==5','rank_mod(M,5)==1','generation despite coincident mod5 eigenvalues')]
results=[]
with tempfile.TemporaryDirectory(prefix='integral-order-controls-') as temporary:
    for name,old,new,marker in cases:
        if source.count(old)!=1:raise RuntimeError('Mutation anchor not unique: '+name)
        p=pathlib.Path(temporary)/(name+'.py')
        p.write_text(source.replace(old,new))
        for flags in ([],['-O']):
            result=subprocess.run([sys.executable,*flags,str(p)],capture_output=True,text=True)
            if result.returncode==0 or marker not in result.stderr:
                raise RuntimeError('Mutation not rejected at intended mathematical check: '+name+' '+repr(flags)+' '+result.stderr)
        results.append(dict(name=name,normal_rejected=True,optimized_rejected=True))
print(json.dumps(dict(status='PASS',mutants=len(cases),runs=2*len(cases),cases=results),indent=2,sort_keys=True))
