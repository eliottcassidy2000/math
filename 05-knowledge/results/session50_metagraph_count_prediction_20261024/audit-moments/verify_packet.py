#!/usr/bin/env python3
"""Hash verification and real manifest-only normal/optimized reproduction."""
import hashlib,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def require(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    manifest=json.loads((ROOT/'MANIFEST.json').read_text())
    with tempfile.TemporaryDirectory(prefix='session50-moment-audit-replay-') as temp:
        fresh=Path(temp)
        for item in manifest['files']:
            rel=Path(item['path'])
            require(not rel.is_absolute() and '..' not in rel.parts,'unsafe manifest path')
            b=(ROOT/rel).read_bytes()
            require(len(b)==item['bytes'] and hashlib.sha256(b).hexdigest()==item['sha256'],'manifest byte mismatch: '+str(rel))
            target=fresh/rel;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(b)
        env=os.environ.copy();env.pop('PYTHONPATH',None);env['PYTHONDONTWRITEBYTECODE']='1'
        for flag,name in [([], 'fresh-normal.json'),(['-O'],'fresh-optimized.json')]:
            subprocess.run([sys.executable,*flag,'check.py','--output',name],cwd=fresh,env=env,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        normal=(fresh/'fresh-normal.json').read_bytes();optimized=(fresh/'fresh-optimized.json').read_bytes()
        require(normal==optimized,'fresh normal/-O mismatch')
        require(normal==(fresh/'normal.json').read_bytes(),'fresh result mismatch against frozen result')
        report=json.loads(normal)
        if (fresh/'verify_producer.py').exists():
            producer=json.loads(subprocess.check_output([sys.executable,'verify_producer.py'],cwd=fresh,env=env,text=True))
            require(producer['status']=='PASS','producer replay failed')
        else:producer=None
    print(json.dumps({'status':'PASS','actual_fresh_directory_manifest_only_replay':True,'normal_equals_optimized':True,'fresh_equals_frozen':True,'independent_atlas_rows':report['atlas']['rows_checked'],'independent_small_orders':[r['n'] for r in report['small_enumeration']],'prediction':report['atlas']['prediction']['prediction'],'isolated_short_certificate':report['atlas']['short_certificate']['prediction'],'producer_replay':producer},indent=2,sort_keys=True))
if __name__=='__main__':main()
