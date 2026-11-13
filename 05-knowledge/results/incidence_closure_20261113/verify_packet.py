#!/usr/bin/env python3
"""Verify closed files, reviewed promotions, paired receipts and historical pins.
Heavy mathematical reproduction commands are documented in AUDIT.md.
"""
from hashlib import sha256
from pathlib import Path
import json,os,re,subprocess
P=Path(__file__).resolve().parent
REPO=P.parents[2]
IGNORED={'_scratch','__pycache__'}
PAIRS=[('closure/normal.out','closure/optimized.out'),
 ('closure/analysis.normal.out','closure/analysis.optimized.out'),
 ('closure/black.normal.json','closure/black.optimized.json'),
 ('pairings/normal.json','pairings/optimized.json'),
 ('bridges/check.normal.json','bridges/check.optimized.json'),
 ('audit/normal.json','audit/optimized.json'),
 ('audit/anchors.normal.json','audit/anchors.optimized.json'),
 ('closure/signed_review.normal.json','closure/signed_review.optimized.json'),
 ('audit_closure/normal.json','audit_closure/optimized.json'),
 ('audit_closure/root_anchors.normal.json','audit_closure/root_anchors.optimized.json'),
 ('audit_bridges/normal.json','audit_bridges/optimized.json'),
 ('audit_bridges/root.normal.json','audit_bridges/root.optimized.json')]

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def digest(b):return sha256(b).hexdigest()
def read(name):return json.loads((P/name).read_text())
def actual_files():
 return {q.relative_to(P).as_posix() for q in P.rglob('*')
         if q.is_file() and not set(q.relative_to(P).parts)&IGNORED
         and q.name!='MANIFEST.json'}
def main():
 manifest=read('MANIFEST.json')['files']
 need(actual_files()==set(manifest),('closed file set differs',actual_files()^set(manifest)))
 for name,row in manifest.items():
  q=P/name;need(not q.is_symlink(),'symlink '+name);b=q.read_bytes()
  need(len(b)==row['bytes'] and digest(b)==row['sha256'],'manifest '+name)
 historical={}
 promoted = dict(read('PROMOTIONS.json')['files'])
 promoted.update(read('BLACK_PROMOTION.json')['files'])
 promoted.update(read('SHADOW_PROMOTION.json')['files'])
 for name,row in promoted.items():
  text=(P/name).read_text();need(digest(text.encode())==row['after_sha256'],'promoted bytes '+name)
  for edit in reversed(row['edits']):
   need(text.count(edit['new'])==edit['count'],'ambiguous inverse '+name)
   text=text.replace(edit['new'],edit['old'])
  need(digest(text.encode())==row['before_sha256'],'reviewed snapshot '+name)
  historical[name]=row['before_sha256']
 for a,b in PAIRS:
  need((P/a).read_bytes()==(P/b).read_bytes(),'normal/-O mismatch '+a)
  if a.endswith('.json'):read(a)
 def bound(name,h):
  need(digest((P/name).read_bytes())==h or historical.get(name)==h,'stale binding '+name)
 for rec,field,code in [
  ('audit/normal.json','source_sha256','audit/independent.py'),
  ('audit/anchors.normal.json','source_sha256','audit/anchors.py'),
  ('pairings/normal.json','script_sha256','pairings/check.py'),
  ('bridges/check.normal.json','source_sha256','bridges/check_bridges.py'),
  ('closure/signed_review.normal.json','source_sha256','closure/check_signed_review.py'),
  ('audit_bridges/normal.json','source_sha256','audit_bridges/check_independent.py'),
  ('audit_bridges/root.normal.json','source_sha256','audit_bridges/check_root.py')]:bound(code,read(rec)[field])
 for rec in ['audit/normal.json','audit/anchors.normal.json','closure/analysis.json']:
  bound('closure/data.json',read(rec)['data_sha256'])
 bound('closure/analysis.json',read('audit/anchors.normal.json')['claims_sha256'])
 bound('closure/analyze.py',read('closure/analysis.json')['source_sha256'])
 bound('closure/data.json',read('pairings/normal.json')['producer_sha256'])
 bound('closure/check_black.py',read('closure/black.normal.json')['source_sha256'])
 bound('closure/data.json',read('closure/black.normal.json')['data_sha256'])
 pv=read('audit_bridges/promotion_verification.json')
 need(pv['promotion_manifest_sha256']==digest((P/'PROMOTIONS.json').read_bytes()),'promotion review manifest')
 bound('pairings/normal.json',read('closure/signed_review.normal.json')['input_sha256'])
 for rec in ['audit_bridges/normal.json','audit_bridges/root.normal.json']:
  for name,h in read(rec)['reviewed_sha256'].items():bound(name,h)
 for name,h in read('audit_closure/normal.json')['source_pins'].items():
  if isinstance(h,str):bound(name,h)
 for lane in ['SHA256SUMS','PAIRINGS_REVIEW_SHA256SUMS','BLACK_SHA256SUMS']:
  for line in (P/'closure'/lane).read_text().splitlines():
   h,name=line.split('  ',1);bound('closure/'+name,h)
 receipt=json.loads((P/'closure/normal.out').read_text().splitlines()[-1])
 bound('closure/check.py',receipt['source_sha256']);bound('closure/data.json',receipt['data_sha256'])
 bound('closure/orbits.cpp',receipt['orbits_source_sha256'])
 pins=read('SOURCE_PINS.json')['files'];env=dict(os.environ,GIT_NO_LAZY_FETCH='1',GIT_TERMINAL_PROMPT='0')
 for name,row in pins.items():
  b=subprocess.check_output(['git','show',row['git_commit']+':'+name],cwd=REPO,env=env)
  need(digest(b)==row['sha256'] and len(b)==row['bytes'],'historical source '+name)
 for source,key in [('closure/source-pins.json','sources'),('pairings/source-pins.json','files')]:
  for name,h in read(source)[key].items():need(pins[name]['sha256']==h,'inherited lane pin '+name)
 for md in P.rglob('*.md'):
  if set(md.relative_to(P).parts)&IGNORED:continue
  for link in re.findall(r'\]\(([^\n)]+)\)',md.read_text()):
   if '://' in link or link.startswith('#'):continue
   target=(md.parent/link.split('#')[0]).resolve()
   if target.is_relative_to(P):need(target.exists(),'broken internal link '+str(target))
   else:
    need(target.is_relative_to(REPO),'escaping link '+link)
    need(target.relative_to(REPO).as_posix() in pins,'unpinned external reference '+link)
 print('PASS: closed manifest, reversible proof promotions, 12 receipt pairs, input/code/review bindings and historical source pins.')
if __name__=='__main__':main()
