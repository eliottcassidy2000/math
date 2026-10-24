"""Verify immutable packet bytes and all normal/optimized exact receipts."""
from pathlib import Path
import hashlib,json,subprocess,sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def main():
    manifest=json.loads((HERE/'MANIFEST.json').read_text())
    for name,digest in manifest['sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
            raise RuntimeError('Manifest mismatch: '+name)
    reports=[]
    for name in ('check','native','robust','balanced','receiver'):
        for optimized in (False,True):
            command=[sys.executable]+(['-O'] if optimized else [])+['-B',str(HERE/(name+'.py'))]
            result=subprocess.run(command,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
            suffix='optimized' if optimized else 'normal';expected=(HERE/(name+'.'+suffix+'.json')).read_bytes()
            if result.stdout!=expected:raise RuntimeError('Receipt mismatch: '+name+' '+suffix)
            parsed=json.loads(result.stdout);reports.append({'checker':name,'mode':suffix,'checks':parsed['checks']})
    print(json.dumps({'status':'PASS complete frozen packet replay','manifest_files':len(manifest['sha256']),
                      'runs':reports},indent=2,sort_keys=True))
if __name__=='__main__':main()
