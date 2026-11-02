#!/usr/bin/env python3
"""Portable exact replay and source-integrity controls for checkpoint five."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile

def require(ok,label):
    if not ok:raise RuntimeError(label)

def verify_sources(root):
    manifest=json.loads((root/'CHECKPOINT5_SOURCE_FILES.json').read_text())
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
    manifest=json.loads((source/'CHECKPOINT5_SOURCE_FILES.json').read_text())
    records=[]
    with tempfile.TemporaryDirectory(prefix='session139-checkpoint5-') as td:
        root=Path(td)/'relocated';root.mkdir();cwd=Path(td)/'wrong-cwd';cwd.mkdir()
        for row in manifest['files']:
            dst=root/row['path'];dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/row['path'],dst)
        shutil.copyfile(source/'CHECKPOINT5_SOURCE_FILES.json',root/'CHECKPOINT5_SOURCE_FILES.json')
        shutil.copyfile(Path(__file__).resolve(),root/'replay_checkpoint5.py')
        def run(label,script,argv,mode,expected_code,expected_output=None,reason=None):
            command=[sys.executable,'-B']+(['-O'] if mode=='optimized' else [])+[str(root/script)]+argv
            cp=subprocess.run(command,cwd=cwd,capture_output=True,timeout=120)
            require(cp.returncode==expected_code,label+' unexpected exit '+str(cp.returncode))
            if expected_output is not None:require(cp.stdout==expected_output,label+' frozen output mismatch')
            if reason is not None:require(reason.encode() in cp.stdout+cp.stderr,label+' expected rejection missing')
            records.append({'label':label,'mode':mode,'returncode':cp.returncode,'stdout_sha256':hashlib.sha256(cp.stdout).hexdigest(),'expected_rejection':expected_code!=0})
        for mode in ('normal','optimized'):
            run('n68-independent','independent_concurrent68/audit_n68.py',[],mode,0,(source/'independent_concurrent68/normal.json').read_bytes())
            run('false-n67','independent_concurrent68/audit_n68.py',['--threshold','67'],mode,1,(source/'independent_concurrent68/negative_threshold67.json').read_bytes(),reason='bootstrap_mass_0')
            run('omit-selected-payment','independent_concurrent68/audit_n68.py',['--omit-selected-payment'],mode,1,(source/'independent_concurrent68/negative_payment.json').read_bytes(),reason='selected_payment_essential')
            run('source-integrity','replay_checkpoint5.py',['--verify-sources-only',str(root)],mode,0)
            run('auditor-curation','independent_concurrent68/verify_curated_checkpoint5.py',[],mode,0,(source/'independent_concurrent68/manifest_verification.json').read_bytes())
        target=root/'concurrent68/PROOF.md';original=target.read_bytes();require(original.startswith(b'# Combining'),'corruption control precondition');target.write_bytes(original.replace(b'# Combining',b'# Xombining',1))
        for mode in ('normal','optimized'):
            run('same-size-proof-corruption','replay_checkpoint5.py',['--verify-sources-only',str(root)],mode,1,reason='source pin concurrent68/PROOF.md')
            run('auditor-curation-corruption','independent_concurrent68/verify_curated_checkpoint5.py',[],mode,1,reason='pin mismatch')
        target.write_bytes(original)
        require(verify_sources(root)==integrity,'restored source integrity')
    positive=json.loads((source/'independent_concurrent68/normal.json').read_text())
    report={'status':'PASS','total_executions':len(records),'scientific_positive_runs':2,'source_integrity_positive_runs':4,'expected_rejections':8,'active_conditions':len(positive['active_checks']),'exact_weighted_graph_controls':int(positive['active_checks']['exact_graph_controls']),'fresh_relocation':True,'wrong_cwd':True,'normal_optimized_both':True,'source_files_verified':integrity['source_files'],'source_manifest_sha256':hashlib.sha256((source/'CHECKPOINT5_SOURCE_FILES.json').read_bytes()).hexdigest(),'runs':records,'scope':'Proof controls and exact arithmetic verification, not finite sampling as an all-order proof. False67 control rejects the certificate, not the theorem.'}
    text=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')

if __name__=='__main__':main()
