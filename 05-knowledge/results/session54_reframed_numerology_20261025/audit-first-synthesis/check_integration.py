from pathlib import Path
import hashlib,json,re,sys,subprocess,tempfile
R=Path(__file__).resolve().parent.parent
A=Path(__file__).resolve().parent
def need(c,m):
 if not c:raise RuntimeError(m)
def sha(b):return hashlib.sha256(b).hexdigest()
pairs={'dismissed-atlas':'session54-dismissed-atlas','cross-base-fractal':'session54-cross-base-fractal','audit-cross-base':'session54-audit-cross-base','tangent-cm':'session54-tangent-cm-publication','audit-tangent-cm':'session54-audit-tangent-cm-publication'}
mp=json.loads((A/'inputs/pre-audit-map.json').read_text());dest=mp['destination'];expected=set();pins={}
for label,dn in pairs.items():
 d=R/dn;m=next(p for p in [d/'MANIFEST.json',d/'manifest.json'] if p.exists());raw=json.loads(m.read_text());fs=raw.get('files',raw)
 need({str(p.relative_to(d)) for p in d.rglob('*') if p.is_file()}==set(fs)|{m.name},'membership '+label)
 for rel,h in fs.items():
  p=d/rel;need(not p.is_symlink() and p.resolve().is_relative_to(d.resolve()),'path '+rel);need(sha(p.read_bytes())==h,'digest '+rel);expected.add(dest+'/'+label+'/'+rel)
 expected.add(dest+'/'+label+'/'+m.name);pins[label]=sha(m.read_bytes())
need(pins==mp['packet_manifests'],'map pins');expected.add(dest+'/README.md')
need({x['destination'] for x in mp['files']}==expected and len(mp['files'])==len(expected),'map coverage')
for row in mp['files']:
 b=Path(row['source']).read_bytes();b.decode('utf8');need(len(b)==row['bytes'] and sha(b)==row['sha256'],'map source '+row['source'])
need((R/'session54-root/first/README.md').read_bytes()==(A/'inputs/README.md').read_bytes(),'summary changed')
links=re.findall(r'\]\(([^)]+)\)',(A/'inputs/README.md').read_text())
for link in links:need(dest+'/'+link in expected,'summary link '+link)
# Verify every bounded atlas slice against the complete local original.
src=json.loads((R/'session54-dismissed-atlas/SOURCES.json').read_text())
for s in src['sources']:
 b=(R/s['original_path']).read_bytes();need(sha(b)==s['source_sha256'],'atlas original '+s['key']);cut=b''.join(b.splitlines(keepends=True)[s['line_start']-1:s['line_end']]);snap=(R/'session54-dismissed-atlas'/s['snapshot']).read_bytes();need(cut==snap and sha(snap)==s['snapshot_sha256'],'atlas slice '+s['key'])
# Publication mathematical payload must match the already audited private original.
eq=json.loads((R/'session54-audit-tangent-cm-publication/source-equivalence.json').read_text())
for name,key in [('session54-tangent-cm','old_manifest_sha256'),('session54-tangent-cm-publication','publication_manifest_sha256')]:need(sha((R/name/'MANIFEST.json').read_bytes())==eq[key],key)
for row in eq['unchanged']:
 old=(R/'session54-tangent-cm'/row['file']).read_bytes();new=(R/'session54-tangent-cm-publication'/row['file']).read_bytes();need(old==new and sha(new)==row['sha256'],'source equivalence '+row['file'])
private_hashes=set()
for dn in ['session54-tangent-cm','session54-audit-tangent-cm']:
 for p in (R/dn).rglob('*'):
  if p.is_file() and p.suffix in ['.pdf','.txt'] and 'sources' in p.parts:private_hashes.add(sha(p.read_bytes()))
for label in ['tangent-cm','audit-tangent-cm']:
 for p in (R/pairs[label]).rglob('*'):
  if p.is_file():need(p.suffix!='.pdf' and sha(p.read_bytes()) not in private_hashes,'third-party payload '+str(p))
# Read-only independent verification of staged navigation edit.
nav=R/'session54-publication/navigation-option/00-navigation/CURRENT-FRONTIER.md';old=(R/'session53-publication/clean-repo/00-navigation/CURRENT-FRONTIER.md').read_bytes();new=nav.read_bytes();needle=b'[S53](S53.md).';need(old.count(needle)==1,'old anchor');need(new==old.replace(needle,needle+b' [S54](S54.md).'),'exact nav edit');need(len(new)-len(old)==15 and len(new.splitlines())==450,'nav bounds')
receipt=json.loads((R/'session54-publication/navigation-receipt.json').read_text());need(sha(old)==receipt['source_sha256'] and sha(new)==receipt['result_sha256'],'nav receipt')
need(re.findall(rb'\]\(([^)]+)\)',old)==[x for x in re.findall(rb'\]\(([^)]+)\)',new) if x!=b'S54.md'],'prior routes')
print(json.dumps({'status':'PASS','mapped_files':len(expected),'summary_links':links,'atlas_slices':len(src['sources']),'source_equivalence_files':len(eq['unchanged']),'no_article_payload':True,'navigation_delta':15,'navigation_lines':450,'packet_pins':pins},indent=2))
