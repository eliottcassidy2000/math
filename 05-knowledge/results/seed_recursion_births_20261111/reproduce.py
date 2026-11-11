"""Replay all packet computations, checking pinned sources and exact outputs."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
JOBS=['head_search','head_completion','port_transport','cones','birth_atoms',
      'marked_orbits','arithmetic_births','verify_new']
OUTPUTS=[name+'.json' for name in JOBS if name!='verify_new']+['verification.json']


def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--optimized',action='store_true')
    args=parser.parse_args()
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for name,expected in pins['sha256_lf'].items():
        if digest(ROOT/name)!=expected:raise RuntimeError('source changed: '+name)
    before={name:digest(HERE/name) for name in OUTPUTS}
    for job in JOBS:
        command=[sys.executable]+(['-O'] if args.optimized else [])+['-u',str(HERE/(job+'.py'))]
        result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,encoding='utf8')
        if result.returncode:
            print(result.stdout);print(result.stderr);raise RuntimeError(job+' failed')
        print(job+': passed',flush=True)
    after={name:digest(HERE/name) for name in OUTPUTS}
    if before!=after:raise RuntimeError('replay changed a certified output')
    report=dict(mode='optimized' if args.optimized else 'normal',jobs=JOBS,
                pinned_sources=len(pins['sha256_lf']),unchanged_outputs_sha256_lf=after)
    name='REPLAY_OPTIMIZED.json' if args.optimized else 'REPLAY.json'
    (HERE/name).write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('All output hashes unchanged.',flush=True)


if __name__=='__main__':main()
