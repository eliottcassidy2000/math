"""Verify all curated SHA256SUMS files and their hostile-copy rejection."""
from pathlib import Path
import hashlib,json,shutil,tempfile
ROOT=Path(__file__).resolve().parent

def check(root):
 manifests=sorted(root.rglob('SHA256SUMS'))
 counts=[]
 for manifest in manifests:
  count=0
  for line in manifest.read_text().splitlines():
   if not line.strip():continue
   digest,name=line.split(None,1);name=name.lstrip('* ')
   p=(manifest.parent/name).resolve()
   if not p.is_relative_to(root.resolve()):raise RuntimeError('path escape')
   if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:raise RuntimeError('hash mismatch '+str(p.relative_to(root)))
   count+=1
  counts.append({'manifest':str(manifest.relative_to(root)),'entries':count})
 return counts

def main():
 counts=check(ROOT);negative=[]
 for relative,action in [('local_closure/check.py','append'),('gluing/width5.normal.json','append'),('ear_geometry/positive_cut_geometry.md','missing')]:
  work=Path(tempfile.mkdtemp(prefix='s136-artifact-negative-'));copy=work/'packet';shutil.copytree(ROOT,copy)
  p=copy/relative
  if action=='append':p.write_bytes(p.read_bytes()+b'\nCORRUPTED\n')
  else:p.unlink()
  try:check(copy)
  except RuntimeError:negative.append({'path':relative,'mutation':action,'rejected':True})
  else:raise RuntimeError('mutant accepted '+relative)
 return {'status':'PASS','manifests':counts,'entries':sum(c['entries']for c in counts),'negative_controls':negative}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
