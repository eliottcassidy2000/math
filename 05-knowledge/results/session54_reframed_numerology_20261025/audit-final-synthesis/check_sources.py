from pathlib import Path
import json,hashlib,re,posixpath
A=Path(__file__).resolve().parent;mp=json.loads((A/'inputs/pre-audit-map.json').read_text());base=mp['destination']
def need(c,m):
 if not c:raise RuntimeError(m)
def sha(b):return hashlib.sha256(b).hexdigest()
rows=mp['files'];by={x['destination']:x for x in rows};need(len(by)==len(rows)==514,'map uniqueness')
roots={}
for label,pin in mp['packet_manifests'].items():
 candidates=[x for x in rows if x['destination'] in [base+'/'+label+'/'+s for s in ['MANIFEST.json','manifest.json','AUDIT_MANIFEST.json']]];need(len(candidates)==1,'manifest '+label);row=candidates[0];p=Path(row['source']);need(sha(p.read_bytes())==pin,'pin '+label);roots[label]=p.parent
 raw=json.loads(p.read_text());fs=raw if isinstance(raw,list) else raw.get('files',raw.get('sha256',raw));fs={v['path']:v['sha256'] for v in fs} if isinstance(fs,list) else fs
 need({str(q.relative_to(p.parent)) for q in p.parent.rglob('*') if q.is_file()}==set(fs)|{p.name},'closure '+label)
 mapped={x['destination'][len(base+'/'+label+'/'):] for x in rows if x['destination'].startswith(base+'/'+label+'/')};need(mapped==set(fs)|{p.name},'mapped closure '+label)
 for f,h in fs.items():
  h=h if isinstance(h,str) else h['sha256'];q=p.parent/f;need(not q.is_symlink() and q.resolve().is_relative_to(p.parent.resolve()),'path '+f);need(sha(q.read_bytes())==h,'payload '+label+'/'+f)
for row in rows:
 b=Path(row['source']).read_bytes();b.decode('utf8');need(sha(b)==row['sha256'] and len(b)==row['bytes'],'row '+row['destination']);need(not b.startswith(b'%PDF'),'article payload')
readme=(A/'inputs/README.md').read_bytes();need(sha(readme)==by[base+'/README.md']['sha256'],'summary')
syn=(A/'inputs/SYNTHESIS.md').read_bytes();need(sha(syn)==by[base+'/SYNTHESIS.md']['sha256'],'synthesis')
links=re.findall(r'\]\(([^)]+)\)',readme.decode()+syn.decode())
for link in links:need(posixpath.normpath(base+'/'+link.split('#')[0]) in by,'link '+link)
# Preserve only our research prose/manifests as portable source evidence.
source_docs={}
for label,r in roots.items():
 if label not in ['paley-endpoint-observer','audit-paley-endpoint','larger-paley-certificates','audit-larger-paley','endpoint-blocker-response','audit-endpoint-response','diagrams']:continue
 for name in ['README.md','REPORT.md','AUDIT.md','THEOREM.md','CONVERSE.md','TWO-SHEET.md','RESULT.md','PROOF.txt','AUDIT.txt','CRITERION.md']:
  p=r/name
  if p.exists():
   rel=label+'/'+name;target=A/'inputs/source-docs'/rel;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(p.read_bytes());source_docs[rel]=sha(p.read_bytes())
(A/'inputs/source-doc-pins.json').write_text(json.dumps(source_docs,indent=2)+'\n')
print(json.dumps({'status':'PASS','mapped_files':len(rows),'packet_roots':len(roots),'summary_links':links,'source_documents':source_docs,'packet_pins':mp['packet_manifests']},indent=2))
