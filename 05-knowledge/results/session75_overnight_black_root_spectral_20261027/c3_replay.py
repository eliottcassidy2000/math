#!/usr/bin/env python3
"""Portable, fail-closed Session75 C3 replay. Frozen inputs are never modified.

Materialize the checkpoint transport first. All compilation and generated files
live in --work, which must not exist. Python stdlib + g++ are sufficient.
The joint audit's four git-show reads are replaced ONLY in the adapter subprocess
by hash-bound embedded snapshots. This is an OFFLINE SNAPSHOT CHECK, not a new
Git-origin proof; see the separate checkpoint source-origin receipt.
"""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import tempfile
import time

PIN = '96e35642110104cc547516a63f37eebcf930f707'
AUDIT_SHA = 'e0fe4f90de36369aa51450980257e0c566a51ceb52caa517f2139a5d135c611a'
SNAPSHOTS = {
    '05-knowledge/results/paid_root_cycle_budget_20261026/README.md': ('sources/paid_root_README.md', '21b5845d0fea0aa2a5b1f347d208686d05f30b8f410fbf6f894f8d5a0fb16b3c'),
    '05-knowledge/results/paid_root_cycle_budget_20261026/verify.py': ('sources/paid_root_verify.py', 'f6889495e0938a2fa3bece3e536be1a943aa6315a9ee4f3f33f12fbe987c94b2'),
    '05-knowledge/results/root_eisenstein_first_arrival_20261018j.md': ('sources/root_eisenstein_first_arrival_20261018j.md', 'accf35a22c7fc6557e55a9098134eb6e6a74092abc372cc3f24d4f5a41556233'),
    '05-knowledge/results/paid_root_cycle_budget_20261026/frontier.txt': ('frontier.txt', 'df88f03b576b13c98c0c5378fc1cf61685eeb6817073681266d2a781dc329b26'),
}
LANES = {
    'boundary': 'session75-black-boundary',
    'boundary-audit': 'session75-black-boundary-audit',
    'twins': 'session75-black-twins',
    'twins-audit': 'session75-black-twins-audit',
    'root-twins': 'session75-root-twins',
    'joint': 'session75-joint-cutoffs',
    'joint-audit': 'session75-joint-cutoff-audit',
    'tc': 'session75-tc-audit',
    'tc-window': 'session75-tc-window-audit',
}

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')

def snapshot_reader(packet):
    """Allow exactly four declared requests, and no execution of arbitrary git."""
    pins = json.loads((packet/'SOURCE_PINS.json').read_text())
    need(pins['commit'] == PIN, 'offline snapshot pin mismatch')
    actual = {r['repository_path']: (r['local_path'], r['sha256']) for r in pins['sources']}
    need(actual == SNAPSHOTS and len(pins['sources']) == 4, 'offline snapshot allowlist mismatch')
    used = []
    def read(command, *args, **kwargs):
        need(not args and not kwargs, 'unexpected offline snapshot invocation options')
        need(isinstance(command, list) and len(command) == 5 and command[:4] == ['git', '-C', 'OFFLINE_SNAPSHOTS_ONLY', 'show'], 'unexpected offline snapshot command')
        spec = command[4]
        need(isinstance(spec, str) and spec.startswith(PIN+':'), 'unexpected offline snapshot commit')
        name = spec[len(PIN)+1:]
        need(name in SNAPSHOTS, 'unexpected offline snapshot path')
        local, expected = SNAPSHOTS[name]
        data = (packet/local).read_bytes()
        need(hashlib.sha256(data).hexdigest() == expected, 'offline snapshot content hash mismatch')
        used.append(name)
        return data
    return read, used

def offline_joint(script, packet, output):
    need(sha(script) == AUDIT_SHA, 'frozen independent audit source hash mismatch')
    reader, used = snapshot_reader(packet)
    old = subprocess.check_output
    subprocess.check_output = reader
    try:
        sys.argv = [str(script), '--packet', str(packet), '--repo', 'OFFLINE_SNAPSHOTS_ONLY', '--out', str(output)]
        runpy.run_path(str(script), run_name='__main__')
    finally:
        subprocess.check_output = old
    need(sorted(used) == sorted(SNAPSHOTS), 'not all four offline snapshots were requested exactly once')

def self_test():
    with tempfile.TemporaryDirectory(prefix='c3-adapter-hostile-') as td:
        root = Path(td)
        write_json(root/'SOURCE_PINS.json', {'commit': PIN, 'sources': [dict(repository_path=k, local_path=v[0], sha256=v[1]) for k,v in SNAPSHOTS.items()]})
        reader, used = snapshot_reader(root)
        bad = [
            ['git', '-C', 'WRONG_REPO', 'show', PIN+':'+next(iter(SNAPSHOTS))],
            ['git', '-C', 'OFFLINE_SNAPSHOTS_ONLY', 'show', '0'*40+':'+next(iter(SNAPSHOTS))],
            ['git', '-C', 'OFFLINE_SNAPSHOTS_ONLY', 'show', PIN+':../../unlisted'],
            ['git', '-C', 'OFFLINE_SNAPSHOTS_ONLY', 'status'],
        ]
        rejected = []
        for cmd in bad:
            try:
                reader(cmd)
            except RuntimeError as e:
                rejected.append(str(e))
            else:
                raise RuntimeError('hostile adapter request accepted')
        name = next(iter(SNAPSHOTS)); p = root/SNAPSHOTS[name][0]
        p.parent.mkdir(parents=True); p.write_bytes(b'corrupted snapshot\n')
        try:
            reader(['git','-C','OFFLINE_SNAPSHOTS_ONLY','show',PIN+':'+name])
        except RuntimeError as e:
            rejected.append(str(e))
        else:
            raise RuntimeError('corrupted snapshot accepted')
        need(len(rejected) == 5 and not used, 'hostile adapter test incomplete')
        return {'status':'PASS', 'rejected_controls': rejected, 'assertions_required':False}

class Replay:
    def __init__(self, args):
        self.a = args
        self.driver_sha = sha(__file__)
        self.source = args.materialized.resolve()
        self.work = args.work.resolve()
        need(self.source.is_dir(), 'materialized directory missing')
        need(not self.work.exists(), '--work must be a fresh, nonexistent directory')
        need(self.source != self.work and self.source not in self.work.parents, '--work must be outside materialized inputs')
        self.work.mkdir(parents=True)
        self.logs = self.work/'logs'; self.logs.mkdir()
        self.bin = self.work/'bin'; self.bin.mkdir()
        self.jobs = []
        self.results = {}
        self.py = [sys.executable, '-B'] + (['-O'] if args.mode == 'optimized' else [])
        self.flags = ['-std=c++17', '-O2'] if args.mode == 'normal' else ['-std=c++17','-O3','-DNDEBUG']
        if args.mode == 'ubsan':
            self.flags = ['-std=c++17','-O1','-g','-DNDEBUG','-fsanitize=undefined','-fno-sanitize-recover=undefined']
        self.env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
        self.component = None

    def run(self, name, command, cwd, stdin=None, output=None):
        output = output or self.logs/(name+'.out')
        err = self.logs/(name+'.err')
        start = time.monotonic()
        print('RUN '+name, flush=True)
        with output.open('wb') as out, err.open('wb') as e:
            if stdin is None:
                p = subprocess.run(list(map(str,command)), cwd=cwd, env=self.env, stdout=out, stderr=e, timeout=self.a.timeout)
            else:
                with Path(stdin).open('rb') as inp:
                    p = subprocess.run(list(map(str,command)), cwd=cwd, env=self.env, stdin=inp, stdout=out, stderr=e, timeout=self.a.timeout)
        row = {'component': self.component, 'name':name,'command':[str(x) for x in command], 'returncode':p.returncode,'seconds':round(time.monotonic()-start,3),'stdout_sha256':sha(output),'stderr_sha256':sha(err)}
        self.jobs.append(row)
        write_json(self.work/'progress.json', {'status':'RUNNING','mode':self.a.mode,'jobs':self.jobs,'results':self.results})
        need(p.returncode == 0, name+' failed: '+err.read_text()[-4000:])
        need(not err.read_bytes(), name+' wrote stderr (compiler warnings included): '+err.read_text()[-4000:])
        return output

    def compare(self, actual, expected):
        need(Path(actual).read_bytes() == Path(expected).read_bytes(), 'byte mismatch '+str(actual)+' versus '+str(expected))

    def copy_lane(self, component):
        name = LANES[component]; src=self.source/name; dst=self.work/name
        need(src.is_dir(), 'missing materialized lane '+name)
        shutil.copytree(src,dst)
        return src,dst

    def python(self, name, root, script, args=(), expected=None, output=None):
        out=self.run(name,self.py+[root/script]+list(args),root,output=output)
        if expected is not None:self.compare(out,expected)
        return out

    def compile(self, name, root, source):
        output=self.bin/name
        self.run(name+'-compile',['g++']+self.flags+[root/source,'-o',output],root)
        return output

    def boundary(self, audit=False):
        component='boundary-audit' if audit else 'boundary'; src,root=self.copy_lane(component)
        code='audit_universe.cpp' if audit else 'finite/authenticate_universe.cpp'
        script='universe_input.py' if audit else 'finite/universe_input.py'
        exe=self.compile(component+'-universe',root,code)
        data=self.python(component+'-input',root,script)
        out=self.run(component+'-universe',[exe],root,stdin=data)
        self.compare(out,src/('universe-normal.json' if audit else 'finite/universe-normal.json'))
        if self.a.mode != 'ubsan':
            if audit:
                for name,stem in [('families','family'),('survivors','survivors'),('finite','finite'),('controls','controls')]:
                    self.python(component+'-'+name,root,'audit_'+name+'.py',expected=src/(stem+'-normal.json'))
            else:
                for name,script,expected in [('witnesses','finite/verify_witnesses.py','finite/witness-normal.json'),('counterexample','finite/verify_counterexample.py','finite/counterexample-normal.json'),('family','family/check_family.py','family/normal.json'),('odd','family/check_odd.py','family/odd-normal.json'),('survivor','survivors/check_survivors.py','survivors/normal.json')]:
                    self.python(component+'-'+name,root,script,expected=src/expected)
        self.results[component]={'status':'PASS','scope':'all ordinary classes n4..9 authentication'+('; full frozen finite/family/survivor/rejection receipts' if self.a.mode!='ubsan' else '; UBSan C++ authentication only')}

    def twins(self, audit=False):
        component='twins-audit' if audit else 'twins';src,root=self.copy_lane(component)
        script='check_addendum.py' if audit else 'check_twins.py'
        self.python(component,root,script,expected=src/'normal.json')
        self.results[component]={'status':'PASS','scope':'exactly two named n10 sources; all cuts, paths, native routes and cycle polynomials; one pinned incoming recipe'}

    def root_twins(self):
        src,root=self.copy_lane('root-twins')
        self.python('root-twins',root,'check.py',expected=src/'normal.json')
        self.results['root-twins']={'status':'PASS','scope':'independent root controls: exactly two matrices, 1022 named cuts, Hamilton/automorphism/colour/exit counts and two native core routes'}

    def joint(self):
        src,root=self.copy_lane('joint'); manifest=json.loads((src/'MANIFEST.json').read_text());ns=list(map(str,manifest['declared_N']))
        suffix=['--suffix','.optimized'] if self.a.mode=='optimized' else []
        # Delete outputs that this invocation must regenerate; never accept prefilled evidence.
        generated=[]
        for n in ns:
            generated += ['results/native_sources_N'+n+'.json','results/receipt_N'+n+('.optimized' if suffix else '')+'.json','results/literal_receipt_N'+n+('.optimized' if suffix else '')+'.json']
        stems=['mutation_receipt','N4096_clock_border_mechanism','N4096_depth_boundary_stability_hostile','resource_port_diagnostics','pair_birth_receipt_N65536','boundary_invoice_receipt']
        generated += ['results/'+s+('.optimized' if suffix else '')+'.json' for s in stems]
        for name in generated:(root/name).unlink()
        self.python('joint-scan',root,'scan.py',ns+suffix)
        self.python('joint-literal',root,'literal_verify.py',ns+suffix)
        self.python('joint-mutations',root,'mutation_controls.py',suffix)
        self.python('joint-mechanisms',root,'mechanisms.py',suffix)
        self.python('joint-pair',root,'pair_check.py',['65536']+suffix)
        self.python('joint-boundary-invoice',root,'boundary_invoice.py',suffix)
        for name in generated:self.compare(root/name,src/name)
        self.python('joint-summary',root,'summarize.py')
        self.compare(root/'SUMMARY.json',src/'SUMMARY.json')
        self.results['joint']={'status':'PASS','scales':len(ns),'producer_rectangles':245826,'source_records':1000908,'regenerated_files_compared':len(generated)+1,'scope':'producer replay; not an independent audit of all rectangles'}

    def joint_audit(self):
        src=self.source/LANES['joint-audit'];packet=self.source/LANES['joint'];out=self.logs/'joint-independent.json'
        self.run('joint-independent',self.py+[Path(__file__).resolve(),'--offline-joint-adapter',src/'independent_audit.py',packet,out],self.work)
        self.compare(out,src/'independent_receipt.json')
        extra=self.logs/'joint-supplementary.json'
        self.python('joint-supplementary',src,'supplementary_checks.py',['--packet',packet,'--out',extra])
        self.compare(extra,src/'supplementary_receipt.json')
        value=json.loads(out.read_text())
        need(value['total_source_records']==1000908 and value['independently_recomputed_rectangles']==42425 and len(value['direct_cutoff_audits'])==43,'joint independent scope totals')
        self.results['joint-audit']={'status':'PASS','source_records_rebuilt':1000908,'rectangles_structurally_checked':245826,'rectangles_independently_recomputed':42425,'direct_controls':43,'source_authentication':'OFFLINE hash-bound snapshot check only; repository origin is in separate source-origin receipt','adapter_sha256':self.driver_sha,'frozen_audit_sha256':AUDIT_SHA,'snapshots':SNAPSHOTS}

    def tc(self):
        src,root=self.copy_lane('tc')
        for stem in ['analytic','group_checks']:
            self.python('tc-'+stem,root,stem+'_audit.py' if stem=='analytic' else stem+'.py',expected=src/(stem+'.normal.json'))
        exe=self.compile('tc-cayley',root,'cayley_independent.cpp')
        out=self.run('tc-cayley',[exe],root);self.compare(out,src/'cayley_independent.out')
        out=self.run('tc-vt15-independent',[exe,'input'],root,stdin=root/'vt15_independent_inputs.txt');self.compare(out,src/'vt15_independent.out')
        self.python('tc-vt15-source',root,'04-computation/experiments/tc_fibres_n12_20261027/vt15_noncayley.py',expected=src/'vt15.normal.out')
        # This checks frozen receipt cross-bindings, not fresh window computation.
        # The fresh computations above are independently compared byte-for-byte;
        # the window computations have their own mandatory tc-window component.
        window=self.source/LANES['tc-window']
        original=root/'verify_receipts.py'
        need(sha(original)=='cb0e7a6a491611cf6e5ccf0f9040e391329c44eab4aee3c7c48c96ed9d38209f','TC receipt adapter source changed')
        old="W=pathlib.Path('/tmp/session75-tc-window-audit')"
        body=original.read_text();need(body.count(old)==1,'TC sibling adapter occurrence mismatch')
        adapter=root/'verify_receipts.portable.py'
        adapter.write_text(body.replace(old,'W=pathlib.Path('+repr(str(window))+')'))
        self.python('tc-frozen-receipt-consistency',root,adapter.name,expected=src/'reproduce.out')
        self.results['tc']={'status':'PASS','Cayley_inputs':526,'VT15_inputs':16,'greedy_orders':192,'degree9_subgroups':50,'transitive_degree9_subgroups':13,'scope':'analytic/group/exact Cayley acceptance path; does not rerun original large n12 type scan; tc-window replay is separate','frozen_receipt_adapter':{'source_sha256':sha(original),'adapted_sha256':sha(adapter),'replacement_count':1,'purpose':'frozen receipt consistency only; fresh window computation is a separate component'}}

    def tc_window(self):
        src,root=self.copy_lane('tc-window')
        for script,dirs,outname in [('window_configs.py',('jobs-normal','jobs-optimized'),'configs-normal.out'),('window_configs_n12.py',('jobs-n12','jobs-n12-opt'),'configs-n12.out')]:
            for d in dirs:
                shutil.rmtree(root/d)
                out=self.python('tc-window-'+d,root,script,[d]);self.compare(out,src/outname)
                names={p.name for p in (root/d).glob('*.job')};need(names=={p.name for p in (src/d).glob('*.job')},'generated job membership')
                for name in names:self.compare(root/d/name,src/d/name)
        self.python('tc-window-job-audit',root,'audit_jobs.py',expected=src/'audit-jobs-normal.json',output=root/'audit-jobs-normal.json')
        (root/'tc_window_scan_portable.cpp').unlink()
        self.python('tc-window-portable-source',root,'make_portable.py')
        self.compare(root/'tc_window_scan_portable.cpp',src/'tc_window_scan_portable.cpp')
        exe=self.compile('tc-window',root,'tc_window_scan_portable.cpp')
        test=self.compile('tc-aut-control',root,'test_aut.cpp')
        out=self.run('tc-aut-control',[test],root,output=root/'test-aut.out');self.compare(out,src/'test-aut.out')
        if self.a.mode=='ubsan':
            for label,job in [('n22','n22_K0_C5_R0_P1_P.job'),('nonmodule','n16_K0_C3_R0_P1_ker0011.job')]:
                # Match the frozen invocations exactly: the n22 module control
                # has no Held--Karp check (that cross-check is restricted n<=16).
                extra=['check'] if label=='nonmodule' else []
                out=self.run('tc-ubsan-'+label,[exe,'jobs-normal/'+job]+extra,root)
                self.compare(out,src/('ubsan-'+label+'.out'))
            self.results['tc-window']={'status':'PASS','scope':'UBSan: all generated jobs and brute Aut controls; selected order22 module and order16 nonmodule scanner controls only','full_window_replayed':False}
            return
        # Native run_scan.py expects the executable beside itself, still in this disposable work copy.
        shutil.copy2(exe,root/'tc_window_scan_portable')
        for tag in ['normal','n12']:
            (root/('replay-'+tag+'.jsonl')).unlink()
            self.python('tc-window-scan-'+tag,root,'run_scan.py',[tag])
        out=self.python('tc-window-summary',root,'summarize.py')
        value=json.loads(out.read_text());expected=json.loads((src/'SUMMARY.json').read_text())
        for obj in (value,expected):
            for tag in ['normal','n12']:
                del obj[tag]['total_job_seconds'];del obj[tag]['sha256']
        need(value==expected,'TC window semantic summary mismatch (only timing/raw-output hashes ignored)')
        write_json(self.logs/'tc-window.semantic.json',value)
        self.results['tc-window']={'status':'PASS','window_jobs':83,'window_strong':185642,'alternate_n12_jobs':45,'alternate_n12_strong':46968,'projection':'exact whole SUMMARY after removing only per-run timing and raw jsonl hash; each window line suffix checked by frozen summarize.py','full_window_replayed':True,'original_large_n12_scan_replayed':False}

    def main(self):
        components=self.a.component or list(LANES)
        need(len(components)==len(set(components)),'duplicate component')
        for component in components:
            self.component=component
            if component=='boundary':self.boundary()
            elif component=='boundary-audit':self.boundary(True)
            elif component=='twins':self.twins()
            elif component=='twins-audit':self.twins(True)
            elif component=='root-twins':self.root_twins()
            elif component=='joint':self.joint()
            elif component=='joint-audit':self.joint_audit()
            elif component=='tc':self.tc()
            elif component=='tc-window':self.tc_window()
        need(sha(__file__)==self.driver_sha,'driver changed while replay was running')
        result={'status':'PASS','mode':self.a.mode,'adapter_sha256':self.driver_sha,'results':self.results,'jobs':self.jobs,'negative_adapter_controls':self_test(),'immutable_input_policy':'materialized inputs read-only; all compiled/generated artifacts outside them; no source files modified'}
        write_json(self.work/'REPLAY.json',result)
        print(json.dumps({'status':'PASS','mode':self.a.mode,'receipt':str(self.work/'REPLAY.json'),'results':self.results},sort_keys=True),flush=True)

def main():
    if len(sys.argv)>1 and sys.argv[1]=='--offline-joint-adapter':
        need(len(sys.argv)==5,'adapter invocation shape')
        offline_joint(Path(sys.argv[2]),Path(sys.argv[3]),Path(sys.argv[4]));return
    if sys.argv[1:]==['--self-test']:
        print(json.dumps(self_test(),sort_keys=True,indent=2));return
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--materialized',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--mode',choices=['normal','optimized','ubsan'],default='normal')
    p.add_argument('--component',action='append',choices=list(LANES))
    p.add_argument('--timeout',type=int,default=7200,help='per native command seconds; does not silently skip a timed-out job')
    Replay(p.parse_args()).main()

if __name__=='__main__':main()
