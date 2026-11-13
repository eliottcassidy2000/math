"""Replay every finite claim engine in this packet, with explicit scopes.

Normal: python -B reproduce_all.py
Optimized: python -B reproduce_all.py --optimized
Use --only id1,id2 only to finish a partial replay; reports merge exact engine
records and state whether the entire declared scope has been completed.
Text hashes normalize CRLF to LF; binary inputs (including gzip) retain bytes.
Each mode retains its first successful console transcript when the script and
receipt hashes match on a later full recomputation; elapsed timing is not
used as mathematical evidence. Stable execution records remain sealed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import runpy
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
TEXT_SUFFIXES={'.py','.md','.json','.jsonl','.txt','.out','.csv','.tsv','.yaml','.yml','.toml'}

def require(ok,message):
    if not ok: raise RuntimeError(message)

def native(path):
    path=Path(path).resolve()
    value=str(path)
    if os.name=='nt' and not value.startswith('\\\\?\\'):
        value='\\\\?\\UNC\\'+value[2:] if value.startswith('\\\\') else '\\\\?\\'+value
    return Path(value)

def read_json(path):
    return json.loads(native(path).read_text(encoding='utf-8-sig'))

def write_json(path,data):
    native(path).parent.mkdir(parents=True,exist_ok=True)
    native(path).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')

def digest(path):
    data=native(path).read_bytes()
    if Path(path).suffix.lower() in TEXT_SUFFIXES:
        data=data.replace(b'\r\n',b'\n')
    return hashlib.sha256(data).hexdigest()

def engine(name,script,receipts,scope,args=(),depends=(),search=False,derived=()):
    return dict(id=name,script=script,receipts=list(receipts),scope=scope,args=list(args),depends=list(depends),bounded_search=search,derived=list(derived))

ENGINES=[
 engine('port-collisions','returns/probe.py',['returns/port_collisions.json'],
        'All7412 inherited ordinary parent classes of orders1..8; every unordered pair of entry ports.'),
 engine('port-advance','returns/follow_port_switch.py',['returns/port_advance.json'],
        'Every retained neutral collision and the explicit birth-following port comparison.',depends=['port-collisions']),
 engine('constructor-language','returns/constructor_language.py',['returns/constructor_language.json'],
        'All inherited parent classes through order5 and indegree-one roots; legal L/T words of cost1..8; physical macro controls.',depends=['port-collisions']),
 engine('inflation-charts','counts/inflation_probe.py',['counts/inflation_probe.json','counts/ready_families.json'],
        'All3615 non-port module charts over the complete collision bank; exact run polynomials and60 literal controls.',depends=['port-collisions']),
 engine('two-block','counts/two_block.py',['counts/two_block.json'],
        'Exact two-block run polynomial for the specified internal-port versus latest-port difference; empty blocks allowed.'),
 engine('count-audit','counts/audit_counts.py',['counts/audit_counts.json'],
        'Selected k1..8, independent parent permutations k1/2, all225 neutral charts at k2, and two-block grid0..3.',depends=['inflation-charts','two-block']),
 engine('count-repair-radius2','counts/repair_search.py',['counts/repair_search.json'],
        'Every arc subset of sizes0..2 of the specified9-vertex parent; all nonmodule ports, target293 and defect1.',search=True),
 engine('six-state','compiler/six_state_audit.py',['compiler/six_state_audit.json'],
        'Exact six-dimensional rank and affine-rank certificates; seven literal controls; all510 E/B programs length1..8;64 fixed tails length0..6.'),
 engine('module-return','compiler/module_return.py',['compiler/module_return.json'],
        'Literal TT1..7 and three nontransitive module controls for the proved universal return.'),
 engine('dynamics','dynamics/verify.py',['dynamics/receipts.json'],
        '7396 rational pairs; integral reservoirs1..200; hexagon radii1..12; all circulants q3/5/7; three Collatz clocks; six20-vertex full subset counts.'),
 engine('strong-wrapper','dynamics/strong_wrapper_control.py',['dynamics/strong_wrapper_receipts.json'],
        'Six16-vertex strong carriers counted in full and from the root; independent local P2 for22-vertex5/17 formula; equal-order H9 hostile.'),
 engine('exchange','dynamics/exchange_control.py',['dynamics/exchange_receipts.json'],
        '1849 positive rational pairs of numerator/denominator1..8;256 odd-input pairs1..31.'),
 engine('compiler-independent','dynamics/audit_compiler.py',['dynamics/compiler_independent_receipts.json'],
        'Independent Fraction determinant computation and direct permutation/weak-cut enumeration on all seven compiler controls.',depends=['six-state']),
 engine('sparse-repair-radius4','returns/sparse_repair.py',['returns/sparse_repair_radius4.json'],
        'Every parent arc subset of sizes0..4 in the specified9-vertex parent; every sole predecessor; target293 and defect<=2; exact module polynomial controls.',
        args=['--radius','4'],search=True,
        derived=[dict(source='returns/sparse_repair_radius4.json',target='returns/sparse_repair.json',radius=3)]),
 engine('module-relocation','returns/relocate_modules.py',['returns/repair_module_relocation.json'],
        'Every retained positive-defect<=2 radius<=4 target and every one-vertex module except final root and predecessor.',depends=['sparse-repair-radius4']),
 engine('module-migration','returns/module_migration.py',['returns/module_migration.json'],
        'Exact all-module polynomial certificate and independent literal transitive k1..5.'),
 engine('cross-audit','counts/cross_audit.py',['counts/cross_audit.json'],
        'Physical L/T macro states on four seeds; exact four-arc root-star switch; two universal ready targets on five modules.'),
 engine('next-repair-radius3','returns/next_repair.py',['returns/next_repair_radius3.json'],
        'All27776 arc subsets of sizes0..3 of the next11-vertex parent; every sole predecessor; target881, defects1/2; all eligible single modules.',
        args=['--radius','3'],search=True,
        derived=[dict(source='returns/next_repair_radius3.json',target='returns/next_repair_radius2.json',radius=2)]),
 engine('fibre-graph','returns/fibre_graph.py',['returns/fibre_graph.json'],
        'Every macro word through depth7; private-bit recovery; intrinsic peeling under original/reversed labels; Johnson x K2 degrees/diameters, every pairwise BFS distance, both distance-generating profiles, adjacent refinement, and weak-source hostile.'),
 engine('growing-successor','returns/growing_successor_probe.py',['returns/growing_successor_probe.json'],
        'All1541 radius0..2 parents of the next source; all strong target881 ports and all defects; all70 eligible one-module lifts on positive hits.',search=True),
 engine('binary-language','returns/binary_language.py',['returns/binary_language.json'],
        'All2046 L/K words through10 births with canonical decoding/relabeling;62 literal state controls; binary/Fibonacci class counts and H histograms;17 all-K closed-form/book-isomorphism controls c0..16.'),
 engine('ancestry-boundary','compiler/ancestry_boundary.py',['compiler/ancestry_boundary.json'],
        'Six module degree/strongness profiles; all7412 ordinary bank classes through order8;279 strong two-root inverse controls and the exact restricted-birth recurrence.'),
]

def clean_local_path(value):
    value=str(value)
    if value.startswith('\\\\?\\UNC\\'): value='\\\\'+value[8:]
    elif value.startswith('\\\\?\\'): value=value[4:]
    return Path(value).resolve()

def child_run(script,capture,extra):
    seen=set()
    def note(value):
        if not isinstance(value,(str,bytes,os.PathLike)): return
        try:
            p=clean_local_path(os.fsdecode(value))
            if p.is_relative_to(ROOT) and not p.is_relative_to(HERE) and p.suffix.lower()!='.pyc' and native(p).is_file():
                seen.add(p.relative_to(ROOT).as_posix())
        except (OSError,ValueError,TypeError):
            pass
    def audit(event,args):
        if event=='open' and args: note(args[0])
    sys.addaudithook(audit)
    path=HERE/script
    sys.path.insert(0,str(path.parent))
    sys.argv=[str(path)]+list(extra)
    try:
        runpy.run_path(str(path),run_name='__main__')
    finally:
        for mod in tuple(sys.modules.values()):
            value=getattr(mod,'__file__',None)
            if value: note(value)
        inputs={name:digest(ROOT/name) for name in sorted(seen)}
        versions={'python':sys.version.split()[0]}
        for name in ('sympy','numpy'):
            mod=sys.modules.get(name)
            if mod is not None: versions[name]=str(getattr(mod,'__version__','unknown'))
        write_json(HERE/capture,{'script':script,'arguments':list(extra),'external_repository_inputs':inputs,'runtime_versions':versions})

def project_radius(spec):
    source=HERE/spec['source']; target=HERE/spec['target']; radius=spec['radius']
    raw=read_json(source)
    require(raw.get('max_radius',-1)>=radius,'projection source radius')
    data=dict(raw)
    data['max_radius']=radius
    data['stats']={k:v for k,v in raw['stats'].items() if int(k)<=radius}
    data['hits']=[row for row in raw['hits'] if row['radius']<=radius]
    before=read_json(target) if native(target).exists() else None
    if before is not None:
        old_derivation=before.pop('derived_receipt',{})
        legacy=old_derivation.get('legacy_metadata')
        if before!=data and spec['target']=='returns/sparse_repair.json' and 'max_radius' not in before:
            require(before.get('status')=='FINITE-EXACT radius-three repair scope plus exact module-polynomial audit','unrecognized legacy radius3 status')
            require(before.get('scope')=='all parent arc subsets of radius <=3 in X=P+z+x (9 vertices), every new sole predecessor, child total H=293, defect<=2; base module vertex1','unrecognized legacy radius3 scope')
            legacy={k:before[k] for k in ('status','scope')}
            before['max_radius']=radius;before['status']=data['status'];before['scope']=data['scope']
        require(before==data,'Maximal-radius projection does not match the previously recorded exact universe: '+spec['target'])
    else: legacy=None
    data['derived_receipt']={'source':spec['source'],'source_lf_sha256':digest(source),'rule':'Keep identical audit/literal/filter metadata, set max_radius, retain stats and hits with radius<=bound. No independent search rerun is claimed.','bound':radius}
    if legacy is not None:data['derived_receipt']['legacy_metadata']=legacy
    write_json(target,data)
    return {'receipt':spec['target'],'sha256':digest(target),'status':'DERIVED, exact radius projection','source':spec['source']}

def replay(args):
    mode='optimized' if args.optimized else 'normal'
    report_path=HERE/('REPLAY_'+mode.upper()+'.json')
    previous=read_json(report_path) if native(report_path).exists() else {}
    completed={row['id']:row for row in previous.get('engines',[])}
    selected={s for s in (args.only or '').split(',') if s}
    known={e['id'] for e in ENGINES}
    require(not(selected-known),'Unknown engine ids: '+str(selected-known))
    chosen=[e for e in ENGINES if not selected or e['id'] in selected]
    for item in chosen:
        for dependency in item['depends']:
            dep=next(e for e in ENGINES if e['id']==dependency)
            require(all(native(HERE/p).is_file() for p in dep['receipts']),'Missing dependency receipt for '+item['id'])
        before={p:digest(HERE/p) if native(HERE/p).exists() else None for p in item['receipts']}
        folder=HERE/'_replay'/mode
        native(folder).mkdir(parents=True,exist_ok=True)
        capture=(folder/(item['id']+'.inputs.json')).relative_to(HERE).as_posix()
        command=[sys.executable,'-B']+(['-O'] if args.optimized else [])+[str(HERE/'reproduce_all.py'),'--engine',item['script'],'--capture-inputs',capture,'--engine-args']+item['args']
        execution_path=folder/(item['id']+'.execution.json')
        if args.finalize_only:
            executed=read_json(execution_path)
            require(executed['exit_code']==0 and executed['script_lf_sha256']==digest(HERE/item['script']),'Cannot finalize changed/failed execution')
            require(executed['receipts']=={p:digest(HERE/p) for p in item['receipts']},'Completed receipt changed before finalization')
            print('FINALIZE '+item['id']+' without rerunning its completed engine',flush=True)
        else:
            print('START '+item['id']+' '+item['scope'],flush=True)
            proc=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,encoding='utf-8',errors='replace',env={**os.environ,'PYTHONIOENCODING':'utf-8','PYTHONUTF8':'1','PYTHONDONTWRITEBYTECODE':'1'})
            fresh={'id':item['id'],'exit_code':proc.returncode,'script_lf_sha256':digest(HERE/item['script']),'receipts':{p:digest(HERE/p) for p in item['receipts'] if native(HERE/p).exists()},'input_capture':capture}
            old_execution=read_json(execution_path) if native(execution_path).exists() else {}
            retain_first=proc.returncode==0 and all(old_execution.get(k)==v for k,v in fresh.items())
            if not retain_first:
                native(folder/(item['id']+'.stdout.txt')).write_text(proc.stdout.replace('\r\n','\n'),encoding='utf-8')
                native(folder/(item['id']+'.stderr.txt')).write_text(proc.stderr.replace('\r\n','\n'),encoding='utf-8')
            require(proc.returncode==0,f"{item['id']} failed ({proc.returncode}); inspect _replay/{mode}/{item['id']}.stderr.txt\n{proc.stderr[-4000:]}")
            if not retain_first:write_json(execution_path,fresh)
        outputs=[]
        for p in item['receipts']:
            after=digest(HERE/p)
            require(before[p] is None or before[p]==after,'Receipt changed during exact replay: '+p+'; inspect before accepting a new baseline')
            outputs.append({'receipt':p,'sha256':after,'status':'RECOMPUTED'})
        outputs.extend(project_radius(spec) for spec in item['derived'])
        completed[item['id']]={**item,'exit_code':0,'outputs':outputs,'input_capture':capture,'script_lf_sha256':digest(HERE/item['script'])}
        write_json(report_path,{'status':('PASS' if known<=completed.keys() else 'PASS so far'),'mode':mode,'hash_basis':'SHA256; text CRLF normalized to LF, binary untouched','declared_engine_ids':[e['id'] for e in ENGINES],'complete':known<=completed.keys(),'engines':[completed[e['id']] for e in ENGINES if e['id'] in completed]})
        print('PASS '+item['id'],flush=True)
    print('REPLAY '+mode+' '+str(len(completed))+'/'+str(len(ENGINES))+' engines complete',flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--optimized',action='store_true')
    ap.add_argument('--finalize-only',action='store_true',help='Finalize saved successful executions after a packaging-only interruption; does not rerun engines')
    ap.add_argument('--only',help='Comma-separated engine ids for resuming a partial replay')
    ap.add_argument('--list',action='store_true')
    ap.add_argument('--engine',help=argparse.SUPPRESS)
    ap.add_argument('--capture-inputs',help=argparse.SUPPRESS)
    ap.add_argument('--engine-args',nargs=argparse.REMAINDER,default=[],help=argparse.SUPPRESS)
    args=ap.parse_args()
    if args.engine:
        require(bool(args.capture_inputs),'child input capture path required')
        child_run(args.engine,args.capture_inputs,args.engine_args)
    elif args.list:
        print(json.dumps(ENGINES,indent=2))
    else:
        replay(args)

if __name__=='__main__': main()
