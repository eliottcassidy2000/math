"""LF artifact integrity and immutable dependency pins; not a proof checker."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ARTIFACTS=('README.md','PROOFS.md','AUDIT.md','check.py','receipt.json','seal.py')
RUNTIME='05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
SOURCES={
    '5a6d7be9b8':(
        '05-knowledge/results/ramified_symmetry_spectral_20261030/PROOFS.md',
        '05-knowledge/results/ramified_symmetry_spectral_20261030/parent_group.json',
        RUNTIME),
    '3bfd37ac31':(
        '05-knowledge/results/session111_scale_connections_20261030/PORT_INFORMATION.md',
        '05-knowledge/results/session111_scale_connections_20261030/ACTION_CARRIERS.md'),
    '0751bba7c9':(
        '05-knowledge/results/session111_scale_connections_20261030/BANK_RECURRENCE.md',
        '05-knowledge/results/session111_scale_connections_20261030/SIGN_QUOTIENT.md')}


def git(*args):
    return subprocess.check_output(['git',*args],cwd=ROOT)


def sha(data):
    return hashlib.sha256(data.replace(b'\r\n',b'\n')).hexdigest()


def snapshot():
    sources=[]
    for ref,paths in SOURCES.items():
        commit=git('rev-parse',ref).decode().strip()
        for path in paths:
            source=git('show',commit+':'+path)
            if path==RUNTIME and sha(source)!=sha((ROOT/path).read_bytes()):
                raise RuntimeError('Live runtime input differs from pin')
            sources.append({'commit':commit,'path':path,'lf_sha256':sha(source)})
    return {'status':'Integrity, not independent mathematical review',
            'artifacts':{path:sha((HERE/path).read_bytes()) for path in ARTIFACTS},
            'sources':sources,'runtime_inputs_checked':[RUNTIME],
            'primary_literature':{'url':'https://arxiv.org/abs/math/0608469',
                                  'version':'v1','sections_read':['1','3.2']}}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--record',action='store_true')
    args=parser.parse_args()
    current=snapshot()
    path=HERE/'manifest.json'
    if args.record:
        path.write_bytes((json.dumps(current,indent=2,sort_keys=True)+'\n').encode())
    elif current!=json.loads(path.read_text(encoding='utf-8')):
        raise RuntimeError('Manifest mismatch')
    print(json.dumps({'artifacts':len(ARTIFACTS),'source_pins':len(current['sources']),
                      'live_runtime_inputs':1,'match':True}))


if __name__=='__main__':
    main()
