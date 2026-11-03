"""Replay the frozen incoming prime packet, normalizing only stdout newlines."""
from pathlib import Path
import hashlib,json,subprocess,sys

HERE=Path(__file__).resolve().parent
def main():
    target=HERE.parent/'session148_square_sum_20261103'/'prime_alignment'
    manifest=json.loads((target/'MANIFEST.json').read_text(encoding='utf-8'))
    files={};runs={}
    for row in manifest['files']:
        data=(target/row['path']).read_bytes();digest=hashlib.sha256(data).hexdigest()
        if len(data)!=row['bytes'] or digest!=row['sha256']:
            raise RuntimeError('incoming artifact differs: '+row['path'])
        files[row['path']]=digest
    for row in manifest['replay']:
        for flags,label in [([], 'normal'),(['-O'], 'optimized')]:
            result=subprocess.run([sys.executable,'-B',*flags,str(target/row['script'])],
                                  capture_output=True,check=True,cwd=target)
            digest=hashlib.sha256(result.stdout.replace(b'\r\n',b'\n')).hexdigest()
            if digest!=row['stdout_sha256']:raise RuntimeError('incoming replay differs')
            runs[label+':'+row['script']]={'raw_stdout_sha256':hashlib.sha256(result.stdout).hexdigest(),
                                         'lf_stdout_sha256':digest}
    (HERE/'incoming_check.json').write_text(json.dumps({'files':files,'runs':runs},
                         indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('PASS incoming prime packet: exact artifacts; normal/-O stdout agrees after newline normalization')
if __name__=='__main__':main()
