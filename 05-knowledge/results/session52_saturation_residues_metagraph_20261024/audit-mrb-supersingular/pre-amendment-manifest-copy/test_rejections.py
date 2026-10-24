#!/usr/bin/env python3
"""Negative controls: intentionally false certificates must fail also under -O."""
from pathlib import Path
import subprocess, sys, tempfile
src=(Path(__file__).resolve().parent/'check.py').read_text()
mutations={
 'wrong_Ogg_membership':('O=[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]','O=[2,3,5,7,11,13,17,19,23,29,31,37,41,47,59,71]'),
 'wrong_CM_character':("character==(1 if p%8 in (1,3) else -1)","character==(1 if p%8 in (1,5) else -1)"),
 'forgot_parity_in_generator':("c+(coeff.get(n-15,0))==expected","c-(coeff.get(n-15,0))==expected"),
 'false_analytic_error_target':("B.b<iv.mpf('2.654e-84')","B.b<iv.mpf('1e-90')")}
with tempfile.TemporaryDirectory(prefix='session52-filter-hostiles-') as tmp:
    for name,(old,new) in mutations.items():
        if src.count(old)!=1: raise RuntimeError('mutation target not unique: '+name)
        p=Path(tmp)/(name+'.py'); p.write_text(src.replace(old,new))
        for optimized in (False,True):
            cmd=[sys.executable]+(['-O'] if optimized else [])+[str(p)]
            result=subprocess.run(cmd,capture_output=True,text=True)
            if result.returncode==0 or 'RuntimeError:' not in result.stderr:
                raise RuntimeError('mutation not explicitly rejected: '+name+' '+str(optimized))
            print('REJECTED',name,'optimized='+str(optimized),result.stderr.strip().splitlines()[-1])
print('PASS: all 8 negative-control executions rejected explicit false claims')
