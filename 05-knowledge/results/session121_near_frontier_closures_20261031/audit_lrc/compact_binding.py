"""Independent binding audit for the publication-sized summary; not a proof replay."""
import hashlib, json, sys
from pathlib import Path
P=Path(__file__).parent.parent/'lrc'
full_path=Path(sys.argv[1]) if len(sys.argv)>1 else P/'triple_receipt.json'
compact_path=Path(sys.argv[2]) if len(sys.argv)>2 else P/'triple_compact_receipt.json'
raw=full_path.read_bytes(); full=json.loads(raw)
compact=json.loads(compact_path.read_bytes())
def need(ok,tag):
    if not ok: raise RuntimeError(tag)
need(compact['complete_ledger_sha256']==hashlib.sha256(raw).hexdigest(),'ledger hash')
need(compact['complete_ledger_bytes']==len(raw),'ledger length')
need(len(compact['families'])==len(full['families'])==286,'census')
for key,value in full.items():
    if key!='families': need(compact[key]==value,'summary '+key)
total=0
for small,large in zip(compact['families'],full['families']):
    for key in ('core','safe_intervals','width','triple_cutoff'): need(small[key]==large[key],key)
    rows=large['branches']; total+=sum(len(b['finite_branches']) for b in rows)
    need(small['first_joint_count']==len(rows),'first count')
    need(small['second_joint_count']==sum(len(b['finite_branches']) for b in rows),'second count')
    need(small['finite_query_count']==sum(max(0,t-v) for b in rows for v,width,t,rej in b['finite_branches']),'leaf count')
    need(small['all_rejection_lists_empty']==all(not rej for b in rows for v,width,t,rej in b['finite_branches']),'acceptance')
    encoded=(json.dumps(rows,separators=(',',':'),sort_keys=True)+'\n').encode()
    need(small['complete_branch_tree_sha256']==hashlib.sha256(encoded).hexdigest(),'tree hash')
need(total==compact['complete_second_joint_branches'],'total pair count')
print(json.dumps({'status':'PASS','core_tree_bindings':286,'complete_second_branches':total,'full_receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
