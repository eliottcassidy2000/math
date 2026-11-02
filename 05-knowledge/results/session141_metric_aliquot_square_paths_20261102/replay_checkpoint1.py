#!/usr/bin/env python3
"""Portable normal/optimized replay and deliberate fail-closed corruption probes."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile

def require(ok,label):
    if not ok:raise RuntimeError(label)
def digest(data):return hashlib.sha256(data).hexdigest()
def source_rows(root,number):return json.loads((root/('SOURCE_FILES_CHECKPOINT'+number+'.json')).read_text())['files']
def verify_sources(root,number):
    rows=source_rows(root,number)
    for row in rows:
        path=Path(row['path'])
        require(not path.is_absolute() and '..' not in path.parts,'unsafe manifest path')
        f=root/path
        require(f.is_file() and not f.is_symlink(),'missing source '+row['path'])
        require(digest(f.read_bytes())==row['sha256'],'source pin '+row['path'])
    return len(rows)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--verify-sources-only',type=Path)
    args=parser.parse_args()
    source=Path(__file__).resolve().parent
    number=Path(__file__).stem.removeprefix('replay_checkpoint')
    config_name='REPLAY_CONFIG_CHECKPOINT'+number+'.json'
    if args.verify_sources_only:
        print(json.dumps({'status':'PASS','source_files':verify_sources(args.verify_sources_only,number)},sort_keys=True));return
    count=verify_sources(source,number)
    config=json.loads((source/config_name).read_text())
    rows=source_rows(source,number);records=[]
    with tempfile.TemporaryDirectory(prefix='session141-checkpoint'+number+'-') as td:
        root=Path(td)/'relocated';root.mkdir()
        cwd=Path(td)/'unrelated-cwd';cwd.mkdir()
        for row in rows:
            dst=root/row['path'];dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(source/row['path'],dst)
        for name in [Path(__file__).name,config_name,'SOURCE_FILES_CHECKPOINT'+number+'.json',*config.get('extra_root_files',[])]:
            shutil.copyfile(source/name,root/name)
        def run(label,script,mode,expected_code=0,output=None,reason=None,argv=()):
            cmd=[sys.executable,'-B']+(['-O'] if mode=='optimized' else [])+[str(root/script),*argv]
            cp=subprocess.run(cmd,cwd=cwd,capture_output=True,timeout=360)
            require(cp.returncode==expected_code,label+' returncode '+str(cp.returncode)+' '+cp.stderr.decode()[-1600:])
            if output is not None:require(cp.stdout==(source/output).read_bytes(),label+' frozen output mismatch')
            if reason is not None:require(reason.encode() in cp.stdout+cp.stderr,label+' rejection reason absent')
            records.append({'label':label,'mode':mode,'returncode':cp.returncode,'stdout_sha256':digest(cp.stdout),'expected_rejection':expected_code!=0})
        for mode in ('normal','optimized'):
            for test in config['scientific_tests']:
                run(test['label'],test['script'],mode,output=test['output'],argv=test.get('argv',[]))
            run('source-integrity',Path(__file__).name,mode,argv=('--verify-sources-only',str(root)))
        for mutation in config['scientific_mutations']:
            f=root/mutation['script'];original=f.read_text()
            require(original.count(mutation['before'])==1,'mutation target '+mutation['label'])
            f.write_text(original.replace(mutation['before'],mutation['after']))
            for mode in ('normal','optimized'):
                run(mutation['label'],mutation['script'],mode,expected_code=1,reason=mutation['reason'])
            f.write_text(original)
        target=config['proof_corruption'];f=root/target;original=f.read_bytes()
        require(original and original[0]!=0,'proof corruption precondition')
        f.write_bytes(bytes([original[0]^1])+original[1:])
        for mode in ('normal','optimized'):
            run('same-size-proof-corruption',Path(__file__).name,mode,expected_code=1,reason='source pin '+target,argv=('--verify-sources-only',str(root)))
        f.write_bytes(original)
        require(verify_sources(root,number)==count,'restored source integrity')
    report={'status':'PASS','checkpoint':number,'fresh_relocation':True,'unrelated_cwd':True,'normal_and_optimized':True,'source_files_verified':count,'total_executions':len(records),'scientific_positive_runs':2*sum(t.get('kind','scientific')=='scientific' for t in config['scientific_tests']),'lane_integrity_positive_runs':2*sum(t.get('kind')=='integrity' for t in config['scientific_tests']),'source_integrity_positive_runs':2,'expected_rejections':2*(len(config['scientific_mutations'])+1),'source_manifest_sha256':digest((source/('SOURCE_FILES_CHECKPOINT'+number+'.json')).read_bytes()),'scope':config['scope'],'runs':records}
    text=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')
if __name__=='__main__':main()
