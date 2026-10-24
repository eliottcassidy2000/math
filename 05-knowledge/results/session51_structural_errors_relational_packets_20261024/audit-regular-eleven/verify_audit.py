#!/usr/bin/env python3
"""Check pinned bytes and reconcile independently computed results, explicitly."""
import hashlib,json,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parent
checks=0
def need(x,m):
    global checks;checks+=1
    if not x:raise RuntimeError(m)
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def main():
    for f in read('MANIFEST.json')['files']:
        p=ROOT/f['path'];need(p.stat().st_size==f['bytes'] and sha(f['path'])==f['sha256'],'audit manifest '+f['path'])
    pins=read('SOURCE_PINS.json');need(sha('producer/MANIFEST.json')==pins['producer_manifest_sha256'],'exact producer freeze')
    for f in read('producer/MANIFEST.json')['files']:
        need(sha('producer/'+f['path'])==f['sha256'],'producer manifest '+f['path'])
    dp=read('degree-normal.json');need(dp==read('degree-optimized.json'),'DP optimization independence')
    fields={'parent_classes':'almost_regular_parent_classes','eligible_source_parents':'eligible_source_parents','eligible_masks':'eligible_labeled_masks','parent_labelled_mass':'labelled_almost_regular_mass','regular_trace':'regular_trace','regular_classes':'regular_classes','regular_error':'regular_error','self_converse_parents':'self_converse_parents','self_converse_regular':'self_converse_regular_classes','parent_aut_histogram':'parent_automorphism_histogram','regular_aut_histogram':'regular_automorphism_histogram','error_histogram':'regular_error_histogram'}
    for label in ('seven','nine','eleven'):
        a=read('independent-'+label+'-normal.json');p=read('producer/'+label+'-O2.json')
        for mode in ('O2','O3'):need(a==read('independent-'+label+'-'+mode+'.json'),'normal O2 O3 independent equality')
        for ka,kp in fields.items():need(a[ka]==p[kp],'independent production match '+ka)
        k=a['source_order']+2
        d=next(x for x in dp['rows'] if x['regular_order']==k)
        need(a['parent_labelled_mass']==d['all_labelled_almost']==d['all_labelled_regular'],'literal DP labeled mass')
        need(a['regular_trace']==k*a['regular_classes'],'independent trace identity')
        need(a['regular_error']==k*(a['parent_classes']-a['regular_classes']),'independent error identity')
        need(a['target_canonicalizations']==0,'target canonicalization boundary')
    need(read('hostile-normal.json')==read('hostile-optimized.json'),'semantic hostile optimization independence')
    regular=read('independent-eleven-normal.json');nonreg=read('producer/sources/nonregular-eleven-receipt.json');combined=read('producer/combined-correction.json')
    need(combined['total_error']==regular['regular_error']+nonreg['known_total_error'],'total correction independent addition')
    need(combined['strong_error']==regular['regular_error']+nonreg['nonregular_strong_error'],'strong correction independent addition')
    print(json.dumps(dict(status='PASS',manifest_sha256=sha('MANIFEST.json'),checks=checks,producer_manifest_sha256=pins['producer_manifest_sha256']),sort_keys=True))
if __name__=='__main__':main()
