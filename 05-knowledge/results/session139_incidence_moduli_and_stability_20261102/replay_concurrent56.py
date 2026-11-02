#!/usr/bin/env python3
"""Portable exact replay and source-integrity controls for the independently audited concurrent n56 frontier."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile

def require(ok,label):
    if not ok:raise RuntimeError(label)

def verify_sources(root):
    manifest=json.loads((root/'CONCURRENT56_SOURCE_FILES.json').read_text())
    for row in manifest['files']:
        f=root/row['path']
        require(f.is_file(),'missing source '+row['path'])
        require(hashlib.sha256(f.read_bytes()).hexdigest()==row['sha256'],'source pin '+row['path'])
    return {'status':'PASS','source_files':len(manifest['files'])}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--verify-sources-only',type=Path)
    args=parser.parse_args()
    if args.verify_sources_only:
        print(json.dumps(verify_sources(args.verify_sources_only),sort_keys=True));return
    source=Path(__file__).resolve().parent
    integrity=verify_sources(source)
    manifest=json.loads((source/'CONCURRENT56_SOURCE_FILES.json').read_text())
    records=[]
    upstream={}
    for pin in manifest['upstream_sources']:
        target=source.parents[2]/pin['path']
        if target.is_file(): data=target.read_bytes()
        else: data=subprocess.check_output(['git','-C',str(source.parents[2]),'show',pin['commit']+':'+pin['path']])
        require(hashlib.sha256(data).hexdigest()==pin['sha256'],'upstream source pin '+pin['path'])
        upstream[Path(pin['path']).name]=data
    with tempfile.TemporaryDirectory(prefix='session139-concurrent56-') as td:
        root=Path(td)/'relocated';root.mkdir();cwd=Path(td)/'wrong-cwd';cwd.mkdir()
        for row in manifest['files']:
            dst=root/row['path'];dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/row['path'],dst)
        shutil.copyfile(source/'CONCURRENT56_SOURCE_FILES.json',root/'CONCURRENT56_SOURCE_FILES.json')
        shutil.copyfile(Path(__file__).resolve(),root/'replay_concurrent56.py')
        (root/'upstream').mkdir();(root/'upstream/verify_stability_56.py').write_bytes(upstream['verify_stability_56.py'])
        def run(label,script,argv,mode,expected_code,expected_output=None,reason=None):
            command=[sys.executable,'-B']+(['-O'] if mode=='optimized' else [])+[str(root/script)]+argv
            cp=subprocess.run(command,cwd=cwd,capture_output=True,timeout=120)
            require(cp.returncode==expected_code,label+' unexpected exit '+str(cp.returncode))
            if expected_output is not None:require(cp.stdout==expected_output,label+' frozen output mismatch')
            if reason is not None:require(reason.encode() in cp.stdout+cp.stderr,label+' expected rejection missing')
            records.append({'label':label,'mode':mode,'returncode':cp.returncode,'stdout_sha256':hashlib.sha256(cp.stdout).hexdigest(),'expected_rejection':expected_code!=0})
        for mode in ('normal','optimized'):
            run('n56-independent','independent_concurrent56/audit_n56.py',[],mode,0,(source/'independent_concurrent56/normal.json').read_bytes())
            run('false-n55','independent_concurrent56/audit_n56.py',['--threshold','55'],mode,1,(source/'independent_concurrent56/negative_threshold55.json').read_bytes(),reason='bootstrap_advance_0')
            run('false-dual','independent_concurrent56/audit_n56.py',['--bad-dual'],mode,1,(source/'independent_concurrent56/negative_dual.json').read_bytes(),reason='dual_remainder_uniform_nonnegative')
            run('source-integrity','replay_concurrent56.py',['--verify-sources-only',str(root)],mode,0)
            run('auditor-curation','independent_concurrent56/verify_curated_audit.py',[],mode,0,(source/'independent_concurrent56/manifest_verification.json').read_bytes())
            run('pinned-upstream-producer','upstream/verify_stability_56.py',[],mode,0)
            require((root/'upstream/stability_56_checks.json').read_bytes()==(source/'independent_concurrent56/producer_normal.json').read_bytes(),'upstream producer frozen report')
        target=root/'independent_concurrent56/AUDIT.md';original=target.read_bytes();require(original.startswith(b'# Independent'),'corruption control precondition');target.write_bytes(original.replace(b'# Independent',b'# Xndependent',1))
        for mode in ('normal','optimized'):
            run('same-size-audit-corruption','replay_concurrent56.py',['--verify-sources-only',str(root)],mode,1,reason='source pin independent_concurrent56/AUDIT.md')
            run('auditor-curation-corruption','independent_concurrent56/verify_curated_audit.py',[],mode,1,reason='pin mismatch')
        target.write_bytes(original)
        require(verify_sources(root)==integrity,'restored source integrity')
    positive=json.loads((source/'independent_concurrent56/normal.json').read_text())
    report={'status':'PASS','total_executions':len(records),'scientific_positive_runs':4,'source_integrity_positive_runs':4,'expected_rejections':8,'active_conditions':len(positive['active_checks']),'producer_gates':4351,'fresh_relocation':True,'wrong_cwd':True,'normal_optimized_both':True,'source_files_verified':integrity['source_files'],'source_manifest_sha256':hashlib.sha256((source/'CONCURRENT56_SOURCE_FILES.json').read_bytes()).hexdigest(),'runs':records,'scope':'Proof controls and exact arithmetic verification, not finite sampling as an all-order proof. False55 control rejects the certificate, not the theorem.'}
    text=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')

if __name__=='__main__':main()
