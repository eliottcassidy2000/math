#!/usr/bin/env python3
"""Relocated exact controls and fail-closed integrity tests for Session140."""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys, tempfile


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def verify_sources(root):
    manifest=json.loads((root/'SOURCE_FILES_CHECKPOINT1.json').read_text())
    for row in manifest['files']:
        f=root/row['path']
        require(f.is_file(), 'missing source '+row['path'])
        require(hashlib.sha256(f.read_bytes()).hexdigest()==row['sha256'], 'source pin '+row['path'])
    return len(manifest['files'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--verify-sources-only', type=Path)
    args=parser.parse_args()
    if args.verify_sources_only:
        print(json.dumps({'status':'PASS','source_files':verify_sources(args.verify_sources_only)},sort_keys=True))
        return
    source=Path(__file__).resolve().parent
    count=verify_sources(source)
    rows=json.loads((source/'SOURCE_FILES_CHECKPOINT1.json').read_text())['files']
    records=[]
    with tempfile.TemporaryDirectory(prefix='session140-checkpoint1-') as td:
        root=Path(td)/'relocated';root.mkdir()
        cwd=Path(td)/'unrelated-cwd';cwd.mkdir()
        for row in rows:
            dst=root/row['path'];dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(source/row['path'],dst)
        for name in ('SOURCE_FILES_CHECKPOINT1.json','replay_checkpoint1.py'):
            shutil.copyfile(source/name,root/name)
        def run(label,script,mode,expected_code=0,output=None,reason=None,argv=()):
            cmd=[sys.executable,'-B']+(['-O'] if mode=='optimized' else [])+[str(root/script),*argv]
            cp=subprocess.run(cmd,cwd=cwd,capture_output=True,timeout=180)
            require(cp.returncode==expected_code,label+' unexpected returncode '+str(cp.returncode)+' '+cp.stderr.decode()[-1000:])
            if output is not None:
                require(cp.stdout==(source/output).read_bytes(),label+' frozen output mismatch')
            if reason is not None:
                require(reason.encode() in cp.stdout+cp.stderr,label+' rejection reason missing')
            records.append({'label':label,'mode':mode,'returncode':cp.returncode,
                'stdout_sha256':hashlib.sha256(cp.stdout).hexdigest(),
                'expected_rejection':expected_code!=0})
        for mode in ('normal','optimized'):
            run('collatz-controls','collatz_hex/verify.py',mode,output='collatz_hex/verification.json')
            run('HN-controls','hadwiger_nelson/verify.py',mode,output='hadwiger_nelson/verification.json')
            run('independent-controls','independent/verify.py',mode,output='independent/normal.json')
            run('source-integrity','replay_checkpoint1.py',mode,argv=('--verify-sources-only',str(root)))
        mutations=[
            ('false-five-step-formula','collatz_hex/verify.py','(16*n-2**(4-k))//3','(15*n-2**(4-k))//3','selection'),
            ('false-mod55-class','hadwiger_nelson/verify.py','[[19,10],[24,10],[29,10],[39,10],[54,2]]','[[19,10],[24,10],[29,10],[38,10],[54,2]]','correct fixed-prime classes'),
        ]
        for label,name,before,after,reason in mutations:
            f=root/name;original=f.read_text()
            require(original.count(before)==1,'mutation target '+label)
            f.write_text(original.replace(before,after))
            for mode in ('normal','optimized'):
                run(label,name,mode,expected_code=1,reason=reason)
            f.write_text(original)
        f=root/'collatz_hex/PROOF.md';original=f.read_bytes()
        require(original.startswith(b'# Collatz'), 'same-size mutation precondition')
        f.write_bytes(original.replace(b'# Collatz',b'# Xollatz',1))
        for mode in ('normal','optimized'):
            run('same-size-proof-corruption','replay_checkpoint1.py',mode,expected_code=1,
                reason='source pin collatz_hex/PROOF.md',argv=('--verify-sources-only',str(root)))
            run('independent-proof-pin-rejection','independent/verify.py',mode,expected_code=1,reason='source pin')
        f.write_bytes(original)
        require(verify_sources(root)==count,'restored-source-integrity')
    report={'status':'PASS','fresh_relocation':True,'unrelated_cwd':True,
        'normal_and_optimized':True,'source_files_verified':count,
        'total_executions':len(records),'scientific_positive_runs':6,
        'source_integrity_positive_runs':2,'expected_rejections':8,
        'source_manifest_sha256':hashlib.sha256((source/'SOURCE_FILES_CHECKPOINT1.json').read_bytes()).hexdigest(),
        'scope':'Finite exact controls and corruption detection; infinite results depend on written proofs and explicitly cited inputs. No Collatz convergence, external HN proof, or Lean build is certified.',
        'runs':records}
    text=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')

if __name__=='__main__':main()
