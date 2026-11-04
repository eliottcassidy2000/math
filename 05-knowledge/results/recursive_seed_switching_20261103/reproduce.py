#!/usr/bin/env python3
"""Verify frozen evidence by default; --rerun executes selected routine jobs.

The default checks saved normal/-O evidence, all source hashes, and exact
payload binding. It does not freshly rerun mathematical searches or phase
counts. One declared prime endpoint is checked by bounded trial division.
--record integrates already-existing evidence and refuses missing or
partial receipts. --rerun is sequential and can take several minutes for the
full named-port counts; --only SCRIPT limits a requested rerun.
"""
from pathlib import Path
import argparse,hashlib,json,re,subprocess,sys,time
from math import isqrt,prod
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parent;ROOT=HERE.parents[2]
JOBS=[
 ('arithmetic/check.py','arithmetic/check.json','stdout','arithmetic/check.optimized.json'),
 ('hilbert/check.py','hilbert/check.json','stdout','hilbert/check.optimized.json'),
 ('selector/check.py','selector/check.json','stdout','selector/check.optimized.json'),
 ('selector/named.py','selector/named_repair.json','file',None),
 ('selector/followup_named.py','selector/followup_named.json','file',None),
 ('transport/compile.py','transport/native.json','file','transport/native.optimized.json'),
 ('transport/audit_native.py','transport/audit_native.json','stdout','transport/audit_native.optimized.json'),
 ('transport/control.py','transport/control.json','stdout','transport/control.optimized.json'),
 ('transport/fast_gap.py','transport/fast_gap_check.json','stdout','transport/fast_gap_check.optimized.json'),
 ('transport/fast_compare.py','transport/fast_compare.json','stdout','transport/fast_compare.optimized.json'),
 ('transport/verify.py','transport/verify.json','file','transport/verify.optimized.json'),
 ('transport/followup.py','transport/followup.json','file','transport/followup.optimized.json'),
 ('transport/followupaudit.py','transport/followupaudit.json','stdout','transport/followupaudit.optimized.json'),
]
JOB_ARGUMENTS={'transport/fast_gap.py':['--self-test']}
ALIASES={
 ('arithmetic/check.json','reader'):RESULTS/'zero_ready_component_returns_20261103/inverse/reader.py',
 ('arithmetic/check.json','sections'):HERE/'selector/inverse.py',
 ('arithmetic/check.json','s159'):RESULTS/'session159_constructions_20261103/automatic/NOTE.md',
}
EXTRA_RECEIPTS=['selector/replay.json']
GENERATED=['selector/bank.json','selector/growth.json']


def need(ok,message):
    if not ok:raise RuntimeError(message)

def digest(data):return hashlib.sha256(data).hexdigest()
def sha(path):return digest(path.read_bytes())
def read(rel):return json.loads((HERE/rel).read_text())
def payload_sha(value):return digest(json.dumps(value,sort_keys=True,separators=(',',':')).encode())
def command(script,optimized):return ['python3','-B']+(['-O'] if optimized else [])+[str((HERE/script).relative_to(ROOT))]+JOB_ARGUMENTS.get(script,[])

def resolve(key,receipt):
    matches={p.resolve() for p in (ROOT/key,RESULTS/key,HERE/key,receipt.parent/key) if p.is_file()}
    need(len(matches)==1,'Unique pinned input path required: '+key+' at '+str(receipt.relative_to(HERE)))
    return matches.pop()

def required_files():
    needed={item for job in JOBS for item in (job[0],job[1],job[3]) if item}
    needed.update(EXTRA_RECEIPTS+GENERATED)
    missing=sorted(p for p in needed if not(HERE/p).is_file())
    need(not missing,'Evidence is not complete; missing files: '+', '.join(missing))

def saved_evidence():
    replay=read('selector/replay.json');need(replay['status']=='PASS' and replay['normal_optimized_checker_identical'] and replay['normal_optimized_named_templates_identical'],'Selector replay claims')
    need(sha(HERE/'selector/check.json')==replay['check_sha256'],'Selector checker digest bound to saved replay')
    evidence=[]
    for script,receipt,mode,optimized in JOBS:
        normal=(HERE/receipt).read_bytes();obj=json.loads(normal)
        need('PARTIAL' not in str(obj.get('status','')),'Partial receipt cannot enter final verification: '+receipt)
        if optimized:
            need(normal==(HERE/optimized).read_bytes(),'Normal/-O receipts differ: '+receipt)
            how='byte-identical saved normal and optimized receipts';reference=optimized
        else:
            need(Path(receipt).name in replay['template_sha256'],'Named producer has saved replay digest')
            need(digest(normal)==replay['template_sha256'][Path(receipt).name],'Named template matches saved normal/-O digest')
            how='saved normal/-O template replay with current receipt digest';reference='selector/replay.json'
        evidence.append(dict(script=script,receipt=receipt,output_semantics=mode,normal_command=command(script,False),optimized_command=command(script,True),
            normal_sha256=digest(normal),optimized_sha256=digest(normal),script_sha256=sha(HERE/script),evidence_kind=how,evidence_reference=reference))
    return evidence

def collect_pins():
    inputs={};references=0
    def add(path,expected,origin):
        nonlocal references
        path=path.resolve();need(sha(path)==expected,'Stale pinned source: '+str(path)+' at '+origin)
        key=str(path.relative_to(ROOT));entry=inputs.setdefault(key,dict(sha256=expected,referenced_by=[]))
        need(entry['sha256']==expected,'Conflicting hashes for '+key)
        if origin not in entry['referenced_by']:entry['referenced_by'].append(origin)
        references+=1
    for rel in [j[1] for j in JOBS]+EXTRA_RECEIPTS:
        receipt=HERE/rel;obj=json.loads(receipt.read_text())
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
        for (which,key),path in ALIASES.items():
            if which==rel:add(path,obj['source_sha256'][key],rel+':/source_sha256/'+key)
        if rel=='selector/replay.json':add(HERE/'selector/check.json',obj['check_sha256'],rel+':/check_sha256')
        if rel=='transport/followup.json':add(HERE/'selector/followup_named.json',obj['input_sha256'],rel+':/input_sha256')
    for entry in inputs.values():entry['referenced_by'].sort()
    external={p:v for p,v in sorted(inputs.items()) if not(ROOT/p).is_relative_to(HERE)}
    internal={p:v for p,v in sorted(inputs.items()) if (ROOT/p).is_relative_to(HERE)}
    return external,internal,references

def payload_checks():
    selector=read('selector/check.json');bank=read('selector/bank.json');growth=read('selector/growth.json')
    need(selector['status']=='PASS' and selector['independent_miss']==dict(H=13835447,hits=0,candidate_integers=1473275,method='Inherited complete interval inversion'),'Exact declared fixed-bank miss')
    need(len(bank['seeds'])==bank['description']['eight_field_keys']==selector['bank_description']['eight_field_keys']==219,'Exact bank size')
    need(bank['description']==selector['bank_description'],'Literal bank description matches checker')
    need(len({tuple(s['fields']) for s in bank['seeds']})==219,'Distinct eight-field bank keys')
    selected={z['H']:z for z in selector['selected_targets']}
    need(selected[13835447]['hits']==0 and selected[41506343]['hits']==selected[41506345]['hits']==9,'Declared selected acceptance trace')
    need(growth['target']==13835447 and len(growth['second_menu']['witnesses'])==3,'Local growth receipt bound to failed target')
    repair=read('selector/named_repair.json');source,target=repair['source'],repair['target']
    native=read('transport/native.json');audit=read('transport/audit_native.json');counted=read('transport/verify.json')
    need(source['H']==4611815 and target['H']==3*source['H']+2==13835447,'Actual repaired parent/child labels')
    need(target['names']==source['names'] and len(target['names'])==23 and not target['new_names'],'Same23 physical names')
    need(target['prefix']==source['prefix']==[50,51,1,2] and target['mark']==4 and target['previous_mark']==31,'Actual repaired mark and clock')
    need(native['program']['parent_H']==source['H'] and native['program']['target_H']==target['H'],'Native program target labels')
    need(not native['program']['permanent_births'] and len(native['program']['changes'])==18,'Repair program scope')
    need([z['changed_pair'] for z in native['program']['changes']]==target['old_arc_differences'],'Exact repair arc-discrepancy binding')
    need(audit['status']=='PASS' and audit['native_complements']==native['receipt']['native_moves']==36 and audit['target_H']==target['H'],'Independent native audit binding')
    need(counted['H']==target['H'] and counted['names']==target['names'] and counted['component_eight']==target['component_eight_fields'],'Independent repair count binding')
    fast=read('transport/fast_gap_check.json')
    need(fast['actual_repair']['H']==target['H'] and fast['actual_repair']['P2']==counted['P2'] and fast['literal_queries']==190,'Optional compiled engine control binds actual repair H/P2')
    fast_pair=read('transport/fast_compare.json');prefix_pair=next(z for z in counted['joint_pairs'] if z['ports']==[50,51])
    need(fast_pair['status']=='PASS' and fast_pair['H']==target['H'] and fast_pair['ports']==[50,51],'Compiled prefix-port audit scope')
    need(fast_pair['fields']==[prefix_pair['H'],prefix_pair['P2'],*prefix_pair['port_data'][0],*prefix_pair['port_data'][1],*prefix_pair['ordered_W2']] and fast_pair['observer17']==prefix_pair['observer17'],'Full prefix-port payload agrees across implementations')
    pair_order=[[39,6],[31,6],[4,6],[50,51]]
    need([z['ports'] for z in counted['joint_pairs']]==pair_order and all(len(z['observer17'])==17 and z['observer17'][0]==target['H'] for z in counted['joint_pairs']),'Repair named joint observer binding')
    follow=read('selector/followup_named.json');full=read('transport/followup.json');fa=read('transport/followupaudit.json')
    need(follow['source']==target and full['source_names']==target['names'] and full['source_rows']==target['rows'],'Actual repaired source reused verbatim')
    need(full['status']=='FINITE_EXACT_NAMED_FOLLOWUP' and fa['status']=='PASS','Complete follow-up counting and independent audit')
    targets={z['H']:z for z in follow['targets']};counts={z['H']:z for z in full['records']};audits={z['target_H']:z for z in fa['records']}
    need(set(targets)==set(counts)==set(audits)=={41506343,41506345},'Only the two declared next targets')
    need(len(full['programs'])==len(full['native_receipts'])==len(full['records'])==len(follow['targets'])==len(fa['records'])==2,'Two complete follow-up programs and endpoint receipts')
    for program in full['programs']:
        t=targets[program['target_H']];c=counts[t['H']];a=audits[t['H']]
        need(t['H']==3*target['H']+t['offset'] and program['parent_H']==target['H'],'Follow-up exact parent and digit')
        need(c['rows']==t['rows'] and c['names']==t['names'],'Exact named graph bound to count receipt')
        need(a['actual_parent_graph_sha256']==payload_sha([target['names'],target['rows']]) and a['actual_target_graph_sha256']==payload_sha([t['names'],t['rows']]),'Full named graph bound to independent audit digests')
        need([z['changed_pair'] for z in program['changes']]==t['old_arc_differences'],'Follow-up exact arc discrepancies')
        need([z['name'] for z in program['permanent_births']]==t['new_names'],'Follow-up permanent birth binding')
        need(c['component_eight']==t['component_eight_fields'] and c['prefix']==[50,51] and c['active_mark']==4,'Follow-up component/frame binding')
        need([z['ports'] for z in c['joint_pairs']]==pair_order and all(len(z['observer17'])==17 and z['observer17'][0]==t['H'] for z in c['joint_pairs']),'Four named joint observers per next target')
        need(c['named_port_charts']==a['named_port_charts'],'Independently audited named field codecs')
    need(fa['native_complements']==86 and fa['permanent_births']==7 and fa['direct_guard_controls']==410,'Declared follow-up transport totals and direct-guard controls')
    prime=41506343;bound=isqrt(prime)
    need(bound==6442 and all(prime%d for d in range(2,bound+1)),'Prime endpoint by complete bounded trial division')
    factors={13835447:[37,179,2089],41506345:[5,8301269]}
    need(all(all(f>1 for f in fs) and prod(fs)==N for N,fs in factors.items()),'Composite endpoint factor-product certificates')
    return dict(fixed_bank_keys=219,certified_bank_miss=13835447,repair_target=13835447,repair_order=23,repair_native_complements=36,
        followup_targets=sorted(targets),followup_orders=[len(targets[h]['names']) for h in sorted(targets)],followup_native_complements=86,followup_permanent_births=7,
        exact_named_edges=3,total_native_complements=122,direct_guard_controls=410,joint_pair_order=pair_order,
        prime_endpoint=dict(H=prime,method='Trial division by every integer from2 through floor(sqrt(H)), inclusive.',last_trial_divisor=bound,trial_divisors_checked=bound-1),composite_endpoint_factor_products={str(N):fs for N,fs in factors.items()},
        graph_and_program_payload_sha256=dict(repair_graph=payload_sha([target['names'],target['rows']]),repair_program=payload_sha(native['program']),followup_programs=payload_sha(full['programs'])),
        meaning='Saved count/audit input hashes plus explicit graph, program, frame and cardinality binding. This is not a fresh search, count or unbounded induction.')

def checker_controls():
    """Reject changed payloads in memory, without altering any saved artifact."""
    global read,sha
    baseline_read=read;baseline_sha=sha;checks=[]
    cases=[
        ('transport/native.json','repair arc payload',lambda d:d['program']['changes'][0].update(changed_pair=[1,2]),'Exact repair arc-discrepancy binding'),
        ('transport/followup.json','counted graph payload',lambda d:d['records'][0]['rows'].__setitem__(0,d['records'][0]['rows'][0]^1),'Exact named graph bound to count receipt'),
        ('transport/followupaudit.json','audited graph digest',lambda d:d['records'][0].update(actual_target_graph_sha256='0'*64),'Full named graph bound to independent audit digests'),
        ('transport/followup.json','partial checkpoint',lambda d:d.update(status='PARTIAL_NAMED_FOLLOWUP'),'Complete follow-up counting and independent audit'),
    ]
    try:
        for receipt,label,mutate,expected in cases:
            def changed(rel):
                obj=baseline_read(rel)
                if rel==receipt:mutate(obj)
                return obj
            read=changed;rejected=False
            try:payload_checks()
            except RuntimeError as error:
                need(expected in str(error),'Payload control failed for an unrelated reason: '+label);rejected=True
            need(rejected,'Modified evidence was accepted: '+label);checks.append(label)
        read=baseline_read;bad=(HERE/'selector/inverse.py').resolve()
        sha=lambda path:'0'*64 if path.resolve()==bad else baseline_sha(path)
        rejected=False
        try:collect_pins()
        except RuntimeError as error:
            need('Stale pinned source' in str(error),'Hash control failed for an unrelated reason');rejected=True
        need(rejected,'Changed source digest was accepted');checks.append('stale source pin')
    finally:read=baseline_read;sha=baseline_sha
    return dict(rejected_in_memory_hostiles=len(checks),controls=checks,scope='Checker integrity tests only; no saved file is modified.')

def rerun(jobs):
    outcomes=[]
    for script,receipt,mode,_ in jobs:
        baseline=(HERE/receipt).read_bytes();observed=[];started=time.monotonic()
        for optimized in (False,True):
            cmd=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(HERE/script)]+JOB_ARGUMENTS.get(script,[])
            try:
                run=subprocess.run(cmd,cwd=ROOT,capture_output=True,check=True)
                output=run.stdout if mode=='stdout' else (HERE/receipt).read_bytes()
                json.loads(output);need(output==baseline,'Fresh producer differs from frozen receipt: '+script+(' -O' if optimized else ''))
                observed.append(digest(output))
            finally:
                if mode=='file':(HERE/receipt).write_bytes(baseline)
        need(observed[0]==observed[1],'Fresh normal/-O equality: '+script)
        outcomes.append(dict(script=script,normal_sha256=observed[0],optimized_sha256=observed[1],wall_seconds=round(time.monotonic()-started,3)))
        print('PASS '+script,flush=True)
    return outcomes

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',action='store_true',help='Integration only: require all completed saved evidence, then write sources.json and verification.json')
    parser.add_argument('--rerun',action='store_true',help='Fresh normal/-O execution of routine jobs; full phase counts can take several minutes')
    parser.add_argument('--only',action='append',choices=[j[0] for j in JOBS],help='With --rerun, limit execution to these script paths; repeatable')
    args=parser.parse_args();need(not args.only or args.rerun,'--only requires --rerun')
    required_files();evidence=saved_evidence();external,internal,references=collect_pins();binding=payload_checks();controls=checker_controls()
    generated={p:sha(HERE/p) for p in GENERATED}
    sources=dict(scope='External source inputs verified against existing receipt pins. Internal packet artifacts belong to the separate manifest.',paths_relative_to='repository root',external_file_count=len(external),external_sources=external)
    if args.record:
        (HERE/'sources.json').write_text(json.dumps(sources,indent=2,sort_keys=True)+'\n')
        report=dict(status='PASS',mode='saved-evidence and current-input verification plus bounded trial division; no tournament search or phase count rerun',
            scope='Fixed219-bank counterexample, one named fourth+2 repair, and only its two declared next endpoints. Unlimited switching and complete spectrum OPEN.',
            normal_optimized_evidence=evidence,payload_binding=binding,checker_self_controls=controls,generated_artifact_sha256=generated,
            external_source_files=len(external),internal_pinned_inputs=len(internal),checked_hash_references=references,
            sources_sha256=sha(HERE/'sources.json'),reproducer_sha256=sha(Path(__file__)),
            optional_engine_scope='The optional compiled counter self-test is included as a separate190-query implementation control with its source pins; full count receipts retain their own declared engine provenance.')
        (HERE/'verification.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    else:
        report=read('verification.json')
        need(report['normal_optimized_evidence']==evidence,'Frozen saved execution evidence or producer changed')
        need(read('sources.json')==sources and sha(HERE/'sources.json')==report['sources_sha256'],'Frozen external source index')
        need(report['reproducer_sha256']==sha(Path(__file__)),'Reproducer changed after recording')
        need(report['payload_binding']==binding and report['generated_artifact_sha256']==generated and report['checker_self_controls']==controls,'Frozen mathematical payload binding and integrity controls')
        need(report['checked_hash_references']==references and report['internal_pinned_inputs']==len(internal),'Frozen source-pin counts')
    if args.rerun:
        jobs=[j for j in JOBS if not args.only or j[0] in args.only]
        outcomes=rerun(jobs);saved_evidence();collect_pins();payload_checks()
        print(json.dumps(dict(status='PASS',mode='fresh normal/-O execution of selected jobs',jobs=outcomes),sort_keys=True))
    else:
        print(json.dumps(dict(status='PASS',mode='saved-evidence verification',routine_jobs=len(evidence),external_sources=len(external),internal_pinned_inputs=len(internal),hash_references=references,exact_named_edges=3,total_native_complements=122),sort_keys=True))
if __name__=='__main__':main()
