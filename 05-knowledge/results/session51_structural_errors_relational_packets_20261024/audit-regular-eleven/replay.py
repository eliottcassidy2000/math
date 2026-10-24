#!/usr/bin/env python3
"""Replay using only physical copies of files explicitly in the audit manifest."""
import argparse,hashlib,json,pathlib,shutil,subprocess,sys,tempfile,time
HERE=pathlib.Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args()
    records=[]
    with tempfile.TemporaryDirectory(prefix='regular-eleven-audit-copy-') as temp:
        copy=pathlib.Path(temp)/'packet';copy.mkdir();out=pathlib.Path(temp)/'outputs';out.mkdir()
        manifest=json.loads((HERE/'MANIFEST.json').read_text())
        for rel in ['MANIFEST.json']+[x['path'] for x in manifest['files']]:
            target=copy/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/rel,target)
        def run(cmd,label):
            start=time.monotonic();c=subprocess.run(list(map(str,cmd)),cwd=out,text=True,capture_output=True)
            if c.returncode:raise RuntimeError(f'{label}: {c.stdout}\n{c.stderr}')
            records.append(dict(step=label,seconds=round(time.monotonic()-start,4),returncode=c.returncode))
        def same(new,old,label):
            if new.read_bytes()!=old.read_bytes():raise RuntimeError(label+' bytes differ')
        for mode,flags in [('normal',[]),('O2',['-O2']),('O3',['-O3','-DNDEBUG'])]:
            binary=out/('independent-'+mode)
            run(['g++','-std=c++17',*flags,copy/'independent.cpp','-o',binary],'compile '+mode)
            for m,label,source in [(5,'seven','parents-5.txt'),(7,'nine','parents-7.txt'),(9,'eleven','parents-nine.txt')]:
                dest=out/(label+'-'+mode+'.json')
                run([binary,m,copy/'producer/sources'/source,copy/'producer'/f'{label}-O2-parents.txt',copy/'producer'/f'{label}-O2-states.tsv',dest],'independent '+label+' '+mode)
                same(dest,copy/f'independent-{label}-{mode}.json','independent '+label+' '+mode)
        for flags,name in [([], 'normal'),(['-O'],'optimized')]:
            run([sys.executable,*flags,copy/'degree_dp.py','--output',out/'degree.json'],'literal degree DP '+name)
            same(out/'degree.json',copy/'degree-normal.json','degree DP')
            run([sys.executable,*flags,copy/'hostiles.py','--producer',copy/'producer','--binary',out/'independent-O3','--output',out/'hostiles.json'],'semantic hostiles '+name)
            same(out/'hostiles.json',copy/'hostile-normal.json','semantic hostiles')
            run([sys.executable,*flags,copy/'verify_audit.py'],'audit frozen validation '+name)
            run([sys.executable,*flags,copy/'producer/verify_packet.py'],'producer frozen validation '+name)
        run([sys.executable,copy/'producer/reproduce.py','--output-directory',out/'producer-replay'],'all producer two-mode replay')
        # Byte-corruption hostility tests against both manifests, restored after each.
        for rel,validator in [('producer/eleven-O2-states.tsv','producer/verify_packet.py'),('independent-eleven-normal.json','verify_audit.py')]:
            target=copy/rel;old=target.read_bytes();target.write_bytes(old+b' ')
            for flags in ([],['-O']):
                c=subprocess.run([sys.executable,*flags,str(copy/validator)],cwd=out,capture_output=True,text=True)
                if c.returncode==0:raise RuntimeError('manifest failed to reject '+rel)
                records.append(dict(step='reject byte mutation '+rel+' '+str(flags),returncode=c.returncode))
            target.write_bytes(old)
    receipt=dict(status='PASS',manifest_sha256=hashlib.sha256((HERE/'MANIFEST.json').read_bytes()).hexdigest(),physical_copy_only=True,steps=records)
    pathlib.Path(args.output).write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
