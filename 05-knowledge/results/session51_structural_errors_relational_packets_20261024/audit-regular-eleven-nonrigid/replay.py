#!/usr/bin/env python3
import argparse,hashlib,json,pathlib,shutil,subprocess,sys,tempfile,time
HERE=pathlib.Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args();records=[];manifest=json.loads((HERE/'MANIFEST.json').read_text())
 with tempfile.TemporaryDirectory(prefix='regular-eleven-nonrigid-audit-copy-') as d:
  root=pathlib.Path(d)/'packet';out=pathlib.Path(d)/'unrelated';root.mkdir();out.mkdir()
  for rel in ['MANIFEST.json']+[f['path'] for f in manifest['files']]:
   dest=root/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/rel,dest)
  def run(cmd,label):
   start=time.monotonic();r=subprocess.run(list(map(str,cmd)),cwd=out,capture_output=True)
   if r.returncode:raise RuntimeError(label+': '+r.stdout.decode()+r.stderr.decode())
   records.append(dict(step=label,seconds=round(time.monotonic()-start,4),returncode=r.returncode));return r.stdout
  def same(a,b,label):
   if a.read_bytes()!=b.read_bytes():raise RuntimeError(label+' differs')
  for mode,flags in [('O2',['-O2']),('O3',['-O3','-DNDEBUG'])]:
   binary=out/('control-'+mode);run(['g++','-std=c++17',*flags,root/'atlas_control.cpp','-o',binary],'compile holdout '+mode)
   run([binary,root/'representatives.txt',root/'sources/parents-ten.txt',root/'sources/states.tsv',out/'control.json'],'complete atlas holdout '+mode);same(out/'control.json',root/f'atlas-control-{mode}.json','holdout')
  for flags,name in [([], 'normal'),(['-O'],'optimized')]:
   run([sys.executable,*flags,root/'independent.py','--output',out/'independent.json','--representatives',out/'representatives.txt'],'independent structural generation '+name);same(out/'independent.json',root/'normal.json','independent result');same(out/'representatives.txt',root/'representatives.txt','independent representatives')
   run([sys.executable,*flags,root/'crosscheck.py','--producer',root/'producer/normal.json','--output',out/'cross.json'],'full producer class bijection '+name);same(out/'cross.json',root/'crosscheck-normal.json','crosscheck')
   run([sys.executable,*flags,root/'hostiles.py','--binary',out/'control-O3','--output',out/'hostiles.json'],'independent semantic controls '+name);same(out/'hostiles.json',root/'hostile-normal.json','hostiles')
   run([sys.executable,*flags,root/'verify_audit.py'],'audit manifest '+name)
   stdout=run([sys.executable,*flags,root/'producer/check.py'],'producer check '+name)
   if stdout!=(root/'producer/normal.json').read_bytes():raise RuntimeError('producer bytes differ')
   run([sys.executable,*flags,root/'producer/verify_manifest.py'],'producer manifest '+name)
  run([sys.executable,root/'producer/replay.py'],'producer own physical-copy and hostile replay')
  for rel,validator in [('normal.json','verify_audit.py'),('producer/check.py','producer/verify_manifest.py')]:
   file=root/rel;old=file.read_bytes();file.write_bytes(old+b' ')
   for flags in ([],['-O']):
    r=subprocess.run([sys.executable,*flags,str(root/validator)],cwd=out,capture_output=True)
    if not r.returncode:raise RuntimeError('failed to reject manifest byte mutation')
    records.append(dict(step='reject manifest mutation '+rel+' '+str(flags),returncode=r.returncode))
   file.write_bytes(old)
 result=dict(status='PASS',physical_manifest_only_copy=True,manifest_sha256=hashlib.sha256((HERE/'MANIFEST.json').read_bytes()).hexdigest(),steps=records);pathlib.Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
