#!/usr/bin/env python3
"""Optional reproduction of this audit's mutations of a supplied producer file.
Usage: python3 audit_producer_mutations.py PATH_TO_check_integral_order.py
This never edits that file; all variants run from temporary copies.
"""
import hashlib,json,pathlib,subprocess,sys,tempfile
if len(sys.argv)!=2:raise SystemExit('Provide the producer checker path')
source=pathlib.Path(sys.argv[1]).read_text();records=[]
mutations=[
 ('reciprocal-step','A=[0,1,1,0,0]','A=[0,1,0,0,1]'),
 ('bad-cyclotomic-inverse','lin([-2,-3,-2,-1],field_cols[:4])','lin([-1,-3,-2,-1],field_cols[:4])'),
 ('growing-cokernel','g=gcd(modulus,11)','g=gcd(modulus,121)'),
 ('wrong-Hadamard-sign','matrixA[i][j]-had[i][j]','matrixA[i][j]+had[i][j]'),
]
with tempfile.TemporaryDirectory(prefix='session143-producer-mutations-') as tmp:
 for name,old,new in mutations:
  count=source.count(old)
  if count<1:raise RuntimeError('mutation anchor absent '+name)
  p=pathlib.Path(tmp)/(name+'.py');p.write_text(source.replace(old,new))
  for flags in ([],['-O']):
   z=subprocess.run([sys.executable,*flags,str(p)],capture_output=True,text=True)
   if z.returncode!=1 or z.stdout or 'RuntimeError:' not in z.stderr:raise RuntimeError('mutation not rejected '+name+' '+str(flags))
   records.append({'mutation':name,'optimized':bool(flags),'anchor_count':count,'exit_code':z.returncode,'last_error':z.stderr.splitlines()[-1]})
print(json.dumps({'status':'PASS','producer_sha256':hashlib.sha256(source.encode()).hexdigest(),'tests':records},indent=2))
