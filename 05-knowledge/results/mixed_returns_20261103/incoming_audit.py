"""Pinned replay of concurrent literal mixed-return controls."""
from pathlib import Path
import hashlib, json, subprocess, sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
    pins=json.loads((HERE/'incoming_pins.json').read_text(encoding='utf8'))
    for path, expected in pins['sha256_lf'].items():
        if digest(ROOT/path)!=expected:
            raise RuntimeError(('changed incoming source',path))
    target=ROOT/'05-knowledge/results/session156_spectrum_20261103/returns/verify_swap2.py'
    cmd=[sys.executable,'-X','utf8','-B']+(['-O'] if sys.flags.optimize else [])+[str(target)]
    run=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,encoding='utf8')
    if run.returncode:
        raise RuntimeError((run.returncode,run.stdout,run.stderr))
    if run.stderr.strip():
        raise RuntimeError(('unexpected diagnostic',run.stderr))
    stdout=run.stdout.replace('\r\n','\n')
    for token in ('H 12784595','H 38353833','H 38353841',
                  'state (231346053, 190887381, 239135436)',
                  'state (231346045, 190887385, 239135436)'):
        if token not in stdout:
            raise RuntimeError(('missing actual control',token))
    result={'status':'PASS','incoming_commit':pins['input_commit'],
            'pinned_inputs':len(pins['sha256_lf']),
            'reproduction_script':str(target.relative_to(ROOT)).replace('\\','/'),
            'stdout_sha256':hashlib.sha256(stdout.encode('utf8')).hexdigest(),
            'stdout':stdout}
    (HERE/'incoming_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Incoming actual mixed-return graphs and ports: PASS')

if __name__=='__main__':
    main()
