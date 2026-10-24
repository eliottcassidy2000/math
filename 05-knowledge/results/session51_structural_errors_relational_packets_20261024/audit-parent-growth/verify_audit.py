#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
B=Path(__file__).resolve().parent

def need(test,message):
 if not test:raise RuntimeError(message)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads((B/path).read_text())
def envelope(directory,name):
 root=B/directory;manifest=json.loads((root/name).read_text());items=manifest['files'];entries=items.items() if isinstance(items,dict) else [(r['path'],r['sha256']) for r in items]
 for path,digest in entries:need(sha(root/path)==digest,'hash mismatch: '+str(root/path))
 return len(items)
count=envelope('','manifest.json');p=envelope('producer','MANIFEST.json');s=envelope('structural-audit','manifest.json');envelope('structural-audit/producer','MANIFEST.json')
need(sha(B/'producer/MANIFEST.json')==read('producer-pin.json')['manifest_sha256'],'production manifest pin')
need(sha(B/'structural-audit/manifest.json')=='81381e8c41aac5bdb57d0b8c8a7744c45e0684e090ec03e854249ca91efd0068','structural audit pin')
need(sha(B/'structural-audit/producer/MANIFEST.json')=='905dd04456f9b4cf302e8b71c7bdd915db209ba081b4735b0f7a733378866f43','structural producer pin')
a=read('independent-O2.json');need(a==read('independent-O3.json'),'independent O2/O3 equality');pb=read('producer/build-normal.json')
for k in ['constructed_classes','pointed_states','refinement_error','upper_trace_order_nine','automorphism_order_histogram','error_histogram']:
 need(a[k]==pb[k],'atlas comparison '+k)
need(a['constructed_mass']==pb['labelled_orbit_mass']==2**36,'nine mass')
need(a['parent_mass']==2**28,'eight mass');need(a['audit_only_order_nine_canonicalization'] is True,'audit scope')
t=read('incoming-O2.json');u=read('incoming-O3.json');pt=read('producer/trace-normal.json')
for k in t:
 if k not in ['passes','singleton_early_exits']:need(t[k]==u[k],'incoming O2/O3 equality '+k)
for left,right in [('upper_trace','upper_refinement_trace'),('pointed_states','pointed_extension_states'),('prediction','predicted_tournament_count'),('stable_root_histogram','stable_root_cell_histogram')]:need(t[left]==pt[right],'trace comparison '+left)
need(u['singleton_early_exits']==0,'full independent refinement');need(read('trace-O3-ndebug-full.json')['singleton_early_exits']==0,'full production refinement')
need(t['per_parent_records_compared']==a['constructed_classes'],'complete per-parent comparison')
need(t['upper_trace']==sum(k*v for k,v in t['stable_root_histogram']),'trace weighted histogram')
need(t['pointed_states']==sum(v for k,v in t['stable_root_histogram']),'pointed histogram')
st=read('strata.json');need(st==read('strata-optimized.json'),'strata normal optimized');need(st['total_error']==t['structural_error_correction']==5508,'structural correction');need(st['corrected_trace']==t['upper_trace']-st['total_error']==10*t['prediction'],'prediction identity')
for name in ('producer-replays.json','independent-replays.json','mutation-results.json','input-control-results.json','strata.json'):need(read(name)['status']=='PASS','pass receipt '+name)
for name,digest in read('producer-replays.json')['source_hashes'].items():need(sha(B/'producer'/name)==digest,'replayed final source '+name)
need(read('mutation-results.json')['semantic_failures']==10,'semantic mutations');need(read('input-control-results.json')['semantic_failures']==12,'input mutations')
need(read('structural-replay.json')['full_error_10']==5508,'structural physical replay')
print(json.dumps({'status':'PASS','payloads':count,'producer_payloads':p,'structural_audit_payloads':s,'producer_manifest_sha256':sha(B/'producer/MANIFEST.json'),'parent_classes':a['constructed_classes'],'pointed_states':t['pointed_states'],'upper_trace':t['upper_trace'],'structural_correction':5508,'prediction':t['prediction'],'per_parent_comparisons':t['per_parent_records_compared'],'semantic_rejections':22},sort_keys=True))
