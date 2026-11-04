#!/usr/bin/env python3
"""Verify frozen evidence by default; --rerun repeats the routine checks.

The default verifies saved normal/-O evidence, receipt/input digests and exact
payload binding. It does not claim freshly executing those mathematical checks.
--rerun executes sequentially; --include-producers also rebuilds the two original
named/native payloads, whose original producers have no saved normal/-O claim.
"""
from pathlib import Path
import argparse,hashlib,json,re,subprocess,sys,time
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent
ROOT=HERE.parents[2]
JOBS=[
 ('component/check.py','component/check.json','stdout','component/check.optimized.json'),
 ('recovery/check.py','recovery/check.json','stdout','recovery/check.optimized.json'),
 ('component/audit_reader.py','component/audit_reader.json','stdout','component/audit_reader.optimized.json'),
 ('component/audit_native.py','component/audit_native.json','stdout','component/audit_native.optimized.json'),
 ('interface/check.py','interface/check.json','stdout',None),
 ('interface/phase.py','interface/phase.json','stdout',None),
 ('interface/natural_pair.py','interface/natural_pair.json','stdout',None),
 ('continuation/build.py','continuation/named_tree.json','file',None),
 ('continuation/native_tree.py','continuation/native_tree.json','file',None),
 ('continuation/verify.py','continuation/verify.json','file',None),
 ('component/audit_continuation.py','component/audit_continuation.json','stdout','component/audit_continuation.optimized.json'),
 ('continuation/lanes.py','continuation/lanes.json','stdout',None),
]
PRODUCERS=[('inverse/named.py','inverse/named_targets.json','file',None),
           ('inverse/native.py','inverse/native.json','file',None)]
BOOTSTRAP_OPT={'interface/check.json':'/tmp/zero-ready-interface-check-O.json',
               'interface/phase.json':'/tmp/zero-ready-interface-phase-O.json',
               'interface/natural_pair.json':'/tmp/zero-ready-natural-pair-O.json',
               'continuation/lanes.json':'/tmp/zero-ready-lanes-O.json'}
SCALAR_PINS={
 ('recovery/check.json','script_sha256'):'recovery/check.py',
 ('component/audit_reader.json','reader_source_sha256'):'inverse/reader.py',
 ('component/audit_reader.json','auditor_sha256'):'component/audit_reader.py',
 ('component/audit_native.json','auditor_sha256'):'component/audit_native.py',
 ('component/audit_continuation.json','auditor_sha256'):'component/audit_continuation.py',
 ('component/audit_continuation.json','elementary_auditor_sha256'):'component/audit_native.py',
 ('inverse/native.json','source_sha256'):'05-knowledge/results/adaptive_metagraph_reformulation_20261103/metagraph_probe.py',
 ('inverse/native.json','target_graph_sha256'):'inverse/named_targets.json',
}

def need(ok,message):
    if not ok:raise RuntimeError(message)

def digest(data):return hashlib.sha256(data).hexdigest()
def sha(path):return digest(path.read_bytes())
def read(rel):return json.loads((HERE/rel).read_text())
def payload_sha(obj):return digest(json.dumps(obj,sort_keys=True,separators=(',',':')).encode())
def command(script,optimized):
    return ['python3']+(['-O'] if optimized else [])+[str((HERE/script).relative_to(ROOT))]

def resolve(key,receipt):
    choices={p.resolve() for p in (ROOT/key,RESULTS/key,HERE/key,receipt.parent/key) if p.is_file()}
    need(len(choices)==1,'unique input path for '+key+' in '+str(receipt))
    return choices.pop()

def collect_pins():
    inputs={};references=0
    receipts=[HERE/job[1] for job in JOBS+PRODUCERS]+[HERE/'continuation/replay.json']
    def add(path,expected,origin):
        nonlocal references
        need(sha(path)==expected,'stale input hash: '+str(path)+' at '+origin)
        key=str(path.relative_to(ROOT));item=inputs.setdefault(key,dict(sha256=expected,referenced_by=[]))
        need(item['sha256']==expected,'conflicting inherited source hashes: '+key)
        if origin not in item['referenced_by']:item['referenced_by'].append(origin)
        references+=1
    for receipt in receipts:
        obj=json.loads(receipt.read_text());rel=str(receipt.relative_to(HERE))
        def walk(value,trail=''):
            if isinstance(value,dict):
                for key,val in value.items():
                    at=trail+'/'+key
                    if isinstance(val,str) and re.fullmatch('[0-9a-f]{64}',val) and ('/' in key or Path(key).suffix in ('.py','.json','.md','.cpp','.txt')):
                        add(resolve(key,receipt),val,rel+':'+at)
                    else:walk(val,at)
            elif isinstance(value,list):
                for i,val in enumerate(value):walk(val,trail+'/'+str(i))
        walk(obj)
        for (which,field),key in SCALAR_PINS.items():
            if which==rel:add(resolve(key,receipt),obj[field],rel+':/'+field)
    # Optional original alignment producer reads this source, but its output
    # has no source-digest field. Current input is bound explicitly here.
    explicit=RESULTS/'triangle_return_170807_20261103/interface.json'
    add(explicit,sha(explicit),'inverse/named.py explicit source input')
    for item in inputs.values():item['referenced_by'].sort()
    external={p:v for p,v in sorted(inputs.items()) if not (ROOT/p).is_relative_to(HERE)}
    internal={p:v for p,v in sorted(inputs.items()) if (ROOT/p).is_relative_to(HERE)}
    return inputs,external,internal,references

def payload_checks():
    first=read('inverse/named_targets.json');native=read('inverse/native.json');audit=read('component/audit_native.json')
    payload={'source_names':first['source_names'],'source_rows':first['source_rows'],
      'targets':[{k:t[k] for k in ('H','names','rows','component_names','component_rows','new_names','old_arc_differences')} for t in first['targets']]}
    need(payload_sha(payload)==audit['stable_named_graph_payload_sha256'],'first graph payload binds to native audit')
    need(payload_sha(native['programs'])==audit['native_program_payload_sha256'],'first native programs bind to audit')
    tree=read('continuation/named_tree.json');programs=read('continuation/native_tree.json');a=read('component/audit_continuation.json')
    keys=('branch','parent_H','H','depth','offset','names','rows','component_names','component_rows','mark','root','spare','prefix','new_names','old_arc_differences','old_name_identity')
    targets=sorted(tree['targets'],key=lambda t:(len(t['branch']),t['branch']))
    payload={'source':{k:tree['source'][k] for k in ('branch','H','names','rows')},'targets':[{k:t[k] for k in keys} for t in targets]}
    need(payload_sha(payload)==a['stable_named_tree_payload_sha256'],'full named tree binds to independent audit')
    need(payload_sha(programs['programs'])==a['stable_native_program_payload_sha256'],'full native tree binds to independent audit')
    counted=read('continuation/verify.json')
    need(len(targets)==len(programs['programs'])==len(counted['records'])==a['edges']==14,'complete14 endpoint/edge receipts')
    need(a['guarded_native_moves']==944 and a['finite_tree_depth']==3,'finite audit scope')
    need({x['H'] for x in targets}=={x['H'] for x in counted['records']},'same fourteen counted labels')
    return dict(first_named_targets=2,continuation_endpoints=14,finite_tree_depth=3,guarded_native_moves=944,
                meaning='Digest and cardinality binding to independently replayed/counting receipts; no new graph census.')

def saved_evidence(record=False):
    old={} if record else {x['script']:x for x in read('verification.json')['normal_optimized_evidence']}
    replay=read('continuation/replay.json');out=[]
    need(replay['status']=='PASS' and replay['normal_optimized_identical'],'continuation replay evidence')
    for script,receipt,mode,optimized in JOBS:
        normal=(HERE/receipt).read_bytes();how='';reference=None
        if optimized:
            need(normal==(HERE/optimized).read_bytes(),'normal/-O paired receipt '+receipt)
            how='byte-identical saved normal and optimized receipts';reference=optimized
        elif receipt.startswith('continuation/') and Path(receipt).name in replay['receipt_sha256']:
            need(sha(HERE/receipt)==replay['receipt_sha256'][Path(receipt).name],'continuation replay digest '+receipt)
            how='recorded normal/-O replay with current receipt digest';reference='continuation/replay.json'
        elif record:
            opt=Path(BOOTSTRAP_OPT[receipt]);need(normal==opt.read_bytes(),'session normal/-O stdout receipt '+receipt)
            how='session normal/-O stdout compared byte-for-byte; digest preserved in this evidence record'
            reference='session temporary receipt '+opt.name
        else:
            before=old[script]
            need(sha(HERE/receipt)==before['normal_sha256']==before['optimized_sha256'],'frozen normal/-O evidence '+receipt)
            need(sha(HERE/script)==before['script_sha256'],'producer changed after recorded normal/-O execution '+script)
            how=before['evidence_kind'];reference=before['evidence_reference']
        out.append(dict(script=script,receipt=receipt,output_semantics=mode,normal_command=command(script,False),
                        optimized_command=command(script,True),normal_sha256=digest(normal),optimized_sha256=digest(normal),
                        script_sha256=sha(HERE/script),evidence_kind=how,evidence_reference=reference))
    return out

def rerun(jobs):
    outcomes=[]
    for script,receipt,mode,_ in jobs:
        baseline=(HERE/receipt).read_bytes();observed=[];started=time.monotonic()
        for optimized in (False,True):
            cmd=[sys.executable]+(['-O'] if optimized else [])+[str(HERE/script)]
            try:
                run=subprocess.run(cmd,cwd=ROOT,capture_output=True,check=True)
                output=run.stdout if mode=='stdout' else (HERE/receipt).read_bytes()
                json.loads(output)
                need(output==baseline,'fresh producer differs from frozen receipt: '+script+(' -O' if optimized else ''))
                observed.append(digest(output))
            finally:
                if mode=='file':(HERE/receipt).write_bytes(baseline)
        need(observed[0]==observed[1],'fresh normal/-O equality '+script)
        outcomes.append(dict(script=script,normal_sha256=observed[0],optimized_sha256=observed[1],wall_seconds=round(time.monotonic()-started,3)))
        print('PASS '+script,flush=True)
    return outcomes

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',action='store_true',help='integration-only: import already compared session evidence and write the two root receipts')
    parser.add_argument('--rerun',action='store_true',help='rerun all routine producers/auditors normally and with -O; includes the bounded C++ menu')
    parser.add_argument('--include-producers',action='store_true',help='with --rerun, also rebuild inverse/named.py and inverse/native.py')
    args=parser.parse_args();need(not args.include_producers or args.rerun,'--include-producers requires --rerun')
    started=time.monotonic();evidence=saved_evidence(args.record)
    inputs,external,internal,references=collect_pins();binding=payload_checks()
    sources=dict(scope='External source inputs checked against existing receipt pins; explicit optional alignment input identified separately. This is not a packet artifact manifest.',
                 paths_relative_to='repository root',external_file_count=len(external),external_sources=external)
    if args.record:
        (HERE/'sources.json').write_text(json.dumps(sources,indent=2,sort_keys=True)+'\n')
        report=dict(status='PASS',mode='recorded evidence and current-input verification; no expensive search repeated',
          scope=f'All {len(JOBS)} routine producers/auditors have saved normal/-O equality evidence. Complete named return tree depth3 only; arbitrary-depth induction OPEN.',
          normal_optimized_evidence=evidence,payload_binding=binding,external_source_files=len(external),internal_pinned_inputs=len(internal),
          checked_hash_references=references,sources_sha256=sha(HERE/'sources.json'),reproducer_sha256=sha(Path(__file__)),
          original_producer_scope=dict(scripts=['inverse/named.py','inverse/native.py'],claim='Saved graph/program payloads independently audited normally and with -O; no original-producer normal/-O claim.',optional_reproduction='--rerun --include-producers'),
          integration_validation_wall_seconds=round(time.monotonic()-started,3))
        (HERE/'verification.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    else:
        report=read('verification.json')
        need(read('sources.json')==sources,'external source index matches current receipt pins')
        need(sha(HERE/'sources.json')==report['sources_sha256'],'source index digest')
        need(sha(Path(__file__))==report['reproducer_sha256'],'reproducer digest')
        need(report['payload_binding']==binding and report['checked_hash_references']==references,'verification binding/counts')
    if args.rerun:
        outcomes=rerun((PRODUCERS if args.include_producers else [])+JOBS)
        collect_pins();payload_checks()
        print(json.dumps(dict(status='PASS',mode='fresh normal/-O execution',jobs=outcomes),sort_keys=True))
    else:
        print(json.dumps(dict(status='PASS',mode='saved-evidence verification',routine_checks=len(evidence),external_sources=len(external),internal_pinned_inputs=len(internal),hash_references=references,finite_tree_depth=3),sort_keys=True))
if __name__=='__main__':main()
