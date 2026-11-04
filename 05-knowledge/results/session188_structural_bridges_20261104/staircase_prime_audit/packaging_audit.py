"""Independent lossless receipt packaging audit with isolated negative controls."""
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile,shutil
HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'staircase_prime_quotient'
FILES=('REASSEMBLY.md','REASSEMBLY.json','reassemble.py','full_records.json.gz.b64.part000')
def need(x,msg):
 if not x:raise RuntimeError(msg)
def digest(b):return hashlib.sha256(b).hexdigest()
def main():
 originals={n:(SRC/n).read_bytes() for n in ('normal.json','optimized.json')};need(originals['normal.json']==originals['optimized.json'] and len(originals['normal.json'])==103600,'complete mathematical originals')
 records=[]
 for optimized in (False,True):
  for case in ('clean','existing_equal','corrupt_chunk','late_output_conflict','bad_output_hash','duplicate_output'):
   with tempfile.TemporaryDirectory(prefix='staircase-package-audit-') as d:
    root=Path(d)
    for name in FILES:shutil.copyfile(SRC/name,root/name)
    if case=='existing_equal':
     for n,b in originals.items():(root/n).write_bytes(b)
    if case=='corrupt_chunk':
     p=root/FILES[-1];b=p.read_bytes();p.write_bytes(bytes([b[0]^1])+b[1:])
    if case=='late_output_conflict':(root/'optimized.json').write_bytes(b'conflicting existing record\n')
    if case in ('bad_output_hash','duplicate_output'):
     p=root/'REASSEMBLY.json';meta=json.loads(p.read_text())
     if case=='bad_output_hash':meta['groups'][0]['outputs'][-1]['sha256']='0'*64
     else:meta['groups'][0]['outputs'].append(dict(meta['groups'][0]['outputs'][0]))
     p.write_text(json.dumps(meta))
    before={n:(root/n).read_bytes() if (root/n).exists() else None for n in originals};cmd=[sys.executable]+(['-O'] if optimized else [])+[str(root/'reassemble.py')];run=subprocess.run(cmd,cwd='/',capture_output=True)
    if case in ('clean','existing_equal'):
     need(run.returncode==0,'successful isolated reconstruction');report=json.loads(run.stdout);need(report['status']=='PASS','reassembler success report')
     for n,b in originals.items():need((root/n).read_bytes()==b,'byte-for-byte full original restoration')
    else:
     need(run.returncode!=0,'negative control must fail')
     after={n:(root/n).read_bytes() if (root/n).exists() else None for n in originals};need(after==before,'all validation before any output write')
    records.append({'optimized':optimized,'case':case,'result':'PASS','output_bytes_unchanged_on_failure':case not in ('clean','existing_equal')})
 pins={name:digest((SRC/name).read_bytes()) for name in FILES};pins.update({name:digest(raw) for name,raw in originals.items()})
 print(json.dumps({'status':'PASS','complete_original_sizes':{n:len(b) for n,b in originals.items()},'source_sha256':pins,'isolated_cases':records,'audit_sha256':digest(Path(__file__).read_bytes())},indent=2,sort_keys=True))
if __name__=='__main__':main()
