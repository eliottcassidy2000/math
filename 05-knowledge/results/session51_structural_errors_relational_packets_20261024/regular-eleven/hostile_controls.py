#!/usr/bin/env python3
"""Small source corruptions rejected by the same optimized production executable."""
import argparse,json,pathlib,subprocess,tempfile
here=pathlib.Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--binary',required=True);p.add_argument('--output',required=True);a=p.parse_args()
lines=(here/'sources/parents-5.txt').read_text().splitlines()
records=[]
def check(name,text):
 with tempfile.TemporaryDirectory() as d:
  src=pathlib.Path(d)/'input.txt';src.write_text(text)
  r=subprocess.run([str(pathlib.Path(a.binary).resolve()),'5',str(src),str(pathlib.Path(d)/'output')],capture_output=True,text=True)
  if r.returncode==0:raise RuntimeError(f'mutation accepted: {name}')
  records.append({'mutation':name,'rejected':True,'stderr':r.stderr.strip()})
check('truncated final action','\n'.join(lines[:-1])+'\n')
check('extra trailing token','\n'.join(lines)+'\n99\n')
bad=lines.copy();first=list(map(int,bad[1].split()));first[0]|=1;bad[1]=' '.join(map(str,first));check('self loop','\n'.join(bad)+'\n')
bad=lines.copy();bad[2]='0 0 2 3 4';check('nonpermutation action','\n'.join(bad)+'\n')
bad=lines.copy();bad[2]='1 0 2 3 4';check('orientation reversing transposition','\n'.join(bad)+'\n')
# Duplicate the first record while retaining all records; literal distinctness rejects.
first_action_count=int(lines[1].split()[-1]);bad=[str(int(lines[0])+1)]+lines[1:]+lines[1:2+first_action_count];check('duplicated source representative','\n'.join(bad)+'\n')
pathlib.Path(a.output).write_text(json.dumps({'status':'PASS','controls':records},indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','rejected_mutations':len(records)}))
