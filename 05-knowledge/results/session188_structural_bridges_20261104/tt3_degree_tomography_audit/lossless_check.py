"""Independent complete-byte package audit; no shared helper imports."""
from pathlib import Path
import ast,base64,gzip,hashlib,json,shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sha(b):return hashlib.sha256(b).hexdigest()
def permitted(path):
 need('MISTAKES.md' not in path.parts,'excluded source is not read');return path

def main():
 meta=json.loads((ROOT/'REASSEMBLY.json').read_text());need(meta['groups'],'nonempty groups');files=['REASSEMBLY.md','REASSEMBLY.json','reassemble.py'];originals={};records=[]
 for gi,g in enumerate(meta['groups']):
  chunks=[]
  for p in g['parts']:
   name=p['file'];need(Path(name).name==name,'local chunk name');raw=(ROOT/name).read_bytes();need(len(raw)==p['size'] and sha(raw)==p['sha256'],'independent chunk length/hash');chunks.append(b''.join(raw.split()));files.append(name)
  decoded=gzip.decompress(base64.b64decode(b''.join(chunks),validate=True))
  if all(o['file'].endswith('.json') for o in g['outputs']):json.loads(decoded)
  elif all(o['file'].endswith('.tsv') for o in g['outputs']):
   table=[list(map(int,line.split())) for line in decoded.decode('ascii').splitlines()];need(len(table)==32768 and all(len(row)==15 and row[0]==i for i,row in enumerate(table)),'complete ordered 32768-row numeric table schema')
  else:raise RuntimeError('unreviewed raw output format')
  for o in g['outputs']:
   name=o['file'];need(Path(name).name==name and name not in originals,'unique local output');original=permitted(ROOT/name).read_bytes();need(decoded==original,'complete byte equality with untouched original');need(len(original)==o['size'] and sha(original)==o['sha256'],'independent original metadata agreement');originals[name]=original;records.append({'file':name,'bytes':len(original),'sha256':sha(original),'group':gi})
 pins={f:sha((ROOT/f).read_bytes()) for f in files};tree=ast.parse((ROOT/'reassemble.py').read_text());body=tree.body[1:] if isinstance(tree.body[0],ast.Expr) and isinstance(tree.body[0].value,ast.Constant) and isinstance(tree.body[0].value.value,str) else tree.body;ast_digest=sha(ast.dump(ast.Module(body=body,type_ignores=[]),include_attributes=False).encode());writes=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in('write_bytes','write_text')];need(len(writes)==1,'single postvalidation write site')
 manifests={};manifest_records=[]
 for name in (() if ROOT.name=='even_mod16_audit' else ('SHA256SUMS','manifest.sha256')):
  p=ROOT/name
  if not p.exists():continue
  manifests[name]=sha(p.read_bytes())
  for line in p.read_text().splitlines():
   if not line.strip():continue
   h,f=line.split(maxsplit=1);f=f.lstrip('*');path=permitted(ROOT/f);need(path.is_file() and sha(path.read_bytes())==h,'existing mathematical manifest '+f);manifest_records.append({'manifest':name,'file':f,'sha256':h})
 receipt_pins={}
 for output,raw in originals.items():
  if not output.endswith('.json'):continue
  data=json.loads(raw)
  for field in('source_sha256','source_pins','dependency_sha256','dependency_pins','auditor_pins'):
   table=data.get(field,{})
   if not isinstance(table,dict):continue
   for name,h in table.items():
    if not isinstance(h,str) or len(h)!=64:continue
    rel=Path(name)
    if rel.is_absolute():path=rel
    elif len(rel.parts)>1:path=ROOT.parent/rel
    elif ROOT.name=='even_mod16_audit' and field.startswith('source'):path=ROOT.parent/'even_cyclic_mod16'/rel
    elif ROOT.name.endswith('_audit') and field.startswith('source'):path=ROOT.parent/ROOT.name[:-6]/rel
    else:path=ROOT/rel
    path=permitted(path);need(path.is_file() and sha(path.read_bytes())==h,'final mathematical receipt pin '+name);receipt_pins[str(path.relative_to(ROOT.parent)) if path.is_relative_to(ROOT.parent) else str(path)]=h
  flat_paths={'proof_sha256':ROOT/'PROOF.md','verifier_sha256':ROOT/'verify.py'}
  if ROOT.name.endswith('_audit'):
   source_root=ROOT.parent/ROOT.name[:-6]
   flat_paths.update({'source_proof':source_root/'PROOF.md','source_verifier':source_root/'verify.py'})
  for field,path in flat_paths.items():
   h=data.get(field)
   if not isinstance(h,str) or len(h)!=64:continue
   path=permitted(path);need(path.is_file() and sha(path.read_bytes())==h,'flat mathematical receipt pin '+field);receipt_pins[str(path.relative_to(ROOT.parent))]=h
 first=next(iter(originals));last=list(originals)[-1];need(first!=last,'late-output control needs distinct outputs');trials=[]
 modes=('fresh','idempotent','first_chunk_corrupt','late_group_chunk_corrupt','late_existing_conflict','late_conflict_with_good_first','late_raw_hash','late_raw_size','duplicate_output','nonlocal_output','nonlocal_chunk')
 for optimized in(False,True):
  for mode in modes:
   with tempfile.TemporaryDirectory(prefix='lossless-package-audit-') as td:
    d=Path(td)
    for name in files:shutil.copyfile(ROOT/name,d/name)
    changed=json.loads(json.dumps(meta));late=changed['groups'][-1]['outputs'][-1];before={}
    if mode=='idempotent':
     for name,raw in originals.items():(d/name).write_bytes(raw)
    if mode in('first_chunk_corrupt','late_group_chunk_corrupt'):
     part=changed['groups'][0 if mode=='first_chunk_corrupt' else -1]['parts'][0 if mode=='first_chunk_corrupt' else -1]['file'];p=d/part;raw=bytearray(p.read_bytes());raw[0]=ord('A') if raw[0]!=ord('A') else ord('B');p.write_bytes(raw)
    if mode in('late_existing_conflict','late_conflict_with_good_first'):(d/last).write_bytes(b'preexisting late conflict')
    if mode=='late_conflict_with_good_first':(d/first).write_bytes(originals[first]);before[first]=((d/first).read_bytes(),(d/first).stat().st_mtime_ns)
    if mode=='late_raw_hash':late['sha256']='0'*64
    if mode=='late_raw_size':late['size']+=1
    if mode=='duplicate_output':late['file']=first
    if mode=='nonlocal_output':late['file']='../outside.json'
    if mode=='nonlocal_chunk':changed['groups'][-1]['parts'][-1]['file']='../outside.chunk'
    if mode in('late_raw_hash','late_raw_size','duplicate_output','nonlocal_output','nonlocal_chunk'):(d/'REASSEMBLY.json').write_text(json.dumps(changed))
    cmd=[sys.executable]+(['-O'] if optimized else[])+[str(d/'reassemble.py')];p=subprocess.run(cmd,cwd='/',capture_output=True)
    if mode in('fresh','idempotent'):
     need(p.returncode==0,'isolated success');report=json.loads(p.stdout);need(report['status']=='PASS' and {r['file'] for r in report['outputs']}==set(originals),'complete success report')
     for name,raw in originals.items():need((d/name).read_bytes()==raw,'isolated exact complete restoration')
    else:
     need(p.returncode!=0,'explicit isolated failure '+mode)
     for name,raw in originals.items():
      if name==last and mode in('late_existing_conflict','late_conflict_with_good_first'):need((d/name).read_bytes()==b'preexisting late conflict','late conflict untouched')
      elif name in before:need((d/name).read_bytes()==before[name][0] and (d/name).stat().st_mtime_ns==before[name][1],'valid early output was not rewritten')
      else:need(not(d/name).exists(),'no output before complete validation')
     need(not(d.parent/'outside.json').exists(),'no nonlocal output')
    trials.append({'mode':mode,'optimized':optimized,'passed':True})
 need(pins=={f:sha((ROOT/f).read_bytes()) for f in files},'package unchanged during audit');need(all((ROOT/n).read_bytes()==raw for n,raw in originals.items()),'original mathematical data untouched');need(manifests=={n:sha((ROOT/n).read_bytes()) for n in manifests},'existing mathematical manifests untouched')
 print(json.dumps({'status':'PASS','package':ROOT.name,'complete_outputs':records,'complete_raw_bytes':sum(len(x) for x in originals.values()),'groups':len(meta['groups']),'parts':sum(len(g['parts']) for g in meta['groups']),'isolated_controls':trials,'validated_existing_manifests':manifest_records,'validated_mathematical_receipt_pins':receipt_pins,'package_sha256':pins,'reassembler_normalized_ast_sha256':ast_digest,'checker_sha256':sha(Path(__file__).read_bytes())},indent=2,sort_keys=True))
if __name__=='__main__':main()
