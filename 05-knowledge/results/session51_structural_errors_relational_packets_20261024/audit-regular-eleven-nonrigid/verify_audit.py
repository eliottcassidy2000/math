#!/usr/bin/env python3
import hashlib,json,pathlib
ROOT=pathlib.Path(__file__).resolve().parent
checks=0
def need(x,m):
 global checks;checks+=1
 if not x:raise RuntimeError(m)
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def main():
 for f in read('MANIFEST.json')['files']:need((ROOT/f['path']).stat().st_size==f['bytes'] and sha(f['path'])==f['sha256'],'audit manifest '+f['path'])
 pins=read('SOURCE_PINS.json')
 for file,key in [('producer/MANIFEST.json','producer_manifest_sha256'),('sources/correction-producer-manifest.json','correction_producer_manifest_sha256'),('sources/correction-audit-manifest.json','correction_audit_manifest_sha256')]:need(sha(file)==pins[key],'source manifest pin')
 for f in read('producer/MANIFEST.json')['files']:need(sha('producer/'+f['path'])==f['sha256'],'producer payload')
 cm=read('sources/correction-producer-manifest.json');am=read('sources/correction-audit-manifest.json')
 for local,source in [('sources/parents-ten.txt','eleven-O2-parents.txt'),('sources/states.tsv','eleven-O2-states.tsv')]:need(sha(local)==next(f['sha256'] for f in cm['files'] if f['path']==source),'holdout payload certified by source manifest')
 need(sha('sources/correction-independent.cpp')==next(f['sha256'] for f in am['files'] if f['path']=='independent.cpp'),'exact inherited independent engine')
 for a,b in [('normal.json','optimized.json'),('representatives.txt','representatives-optimized.txt'),('atlas-control-O2.json','atlas-control-O3.json'),('crosscheck-normal.json','crosscheck-optimized.json'),('hostile-normal.json','hostile-optimized.json')]:need((ROOT/a).read_bytes()==(ROOT/b).read_bytes(),'mode equality '+a)
 x=read('normal.json');control=read('atlas-control-O2.json');prod=read('producer/normal.json')
 need(x['classes']==len(x['representatives'])==prod['nonrigid_classes'],'all independent classes counted')
 delta=sum(11-r['pointed_types'] for r in x['representatives']);need(delta==x['pointing_deficit']==prod['nonrigid_unpointing_deficit'],'independent deficit sum')
 need(sum(r['pointed_types'] for r in x['representatives'])==x['nonrigid_pointed_types']==prod['nonrigid_pointed_states'],'pointed sum')
 need(x['error_histogram']==prod['nonrigid_error_histogram'],'full independent error histogram')
 need(control['candidate_classes']==x['classes'] and control['nonrigid_pointing_deficit']==delta,'independent holdout consistency')
 need(control['regular_classes_holdout']==(control['parent_classes']+delta)//11,'predicted regular count holdout')
 need(control['regular_error_from_deficit']==10*control['parent_classes']-delta,'predicted regular error holdout')
 print(json.dumps(dict(status='PASS',checks=checks,manifest_sha256=sha('MANIFEST.json'),producer_manifest_sha256=pins['producer_manifest_sha256']),sort_keys=True))
if __name__=='__main__':main()
