#!/usr/bin/env python3
from pathlib import Path
import subprocess,tempfile,json
B=Path(__file__).resolve().parent
lines=(B/'producer/parents-nine.txt').read_text().splitlines();count=int(lines[0]);pos=1;records=[]
for _ in range(count):
 ng=int(lines[pos].split()[-1]);records.append(lines[pos:pos+ng+1]);pos+=ng+1
if pos!=len(lines):raise RuntimeError('parse')
rigid=[i for i,r in enumerate(records) if int(r[0].split()[-1])==1][:2]
nontrivial=next(i for i,r in enumerate(records) if int(r[0].split()[-1])>1)
identity=' '.join(map(str,range(9)));results=[]
with tempfile.TemporaryDirectory(prefix='parent-input-mutations-') as tmp:
 D=Path(tmp)
 for label,flags in [('O2',['-O2']),('O3-ndebug',['-O3','-DNDEBUG'])]:
  exe=D/('independent-'+label);subprocess.run(['g++','-std=c++17',*flags,str(B/'independent.cpp'),'-o',str(exe)],check=True,capture_output=True)
  for name,failure in [('duplicate_same_mass','nonisomorphic nine parents'),('omit_nontrivial_actions','full independent order-nine action equality')]:
   mutant=records.copy()
   if name=='duplicate_same_mass':mutant[rigid[1]]=records[rigid[0]]
   else:
    row=records[nontrivial][0].split();row[-1]='1';mutant[nontrivial]=[' '.join(row),identity]
   f=D/(name+'.txt');f.write_text(str(count)+'\n'+'\n'.join(s for rec in mutant for s in rec)+'\n')
   r=subprocess.run([str(exe),str(B/'producer/parents-eight.txt'),str(f),str(D/'result.json')],capture_output=True,text=True)
   if r.returncode==0 or failure not in r.stderr:raise RuntimeError(name+' not rejected: '+r.stderr)
   results.append({'mode':label,'mutation':name,'failure':r.stderr.strip(),'exit_code':r.returncode,'labelled_mass_unchanged':name=='duplicate_same_mass'})
  for name,payload in [('empty',''),('malformed_count','not-an-integer\n'),('negative_count','-1\n'),('truncated_record','1\n0 0\n')]:
   f=D/(name+'.txt');f.write_text(payload);r=subprocess.run([str(exe),str(f),str(B/'producer/parents-nine.txt'),str(D/'result.json')],capture_output=True,text=True)
   if r.returncode!=1:raise RuntimeError(name+' input failure must be checked exception')
   results.append({'mode':label,'mutation':name,'exit_code':r.returncode,'failure':r.stderr.strip()})
(B/'input-control-results.json').write_text(json.dumps({'status':'PASS','semantic_failures':len(results),'results':results},indent=2)+'\n')
print(json.dumps({'status':'PASS','semantic_failures':len(results)}))
