"""Complete labelled eight-core census; graph-polynomial list certificates."""
import argparse,hashlib,json,shutil,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
def sha(raw):return hashlib.sha256(raw).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--order',type=int,choices=(7,8),default=8);args=parser.parse_args()
    pins={name:sha((HERE/name).read_bytes()) for name in ('run.py','census.cpp')}
    compiler=shutil.which('g++')
    if not compiler:raise RuntimeError('existing g++ required')
    with tempfile.TemporaryDirectory(prefix='eight-private-core-') as tmp:
        exe=Path(tmp)/'census.exe';subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'census.cpp'),'-o',str(exe)],check=True,capture_output=True)
        result=json.loads(subprocess.run([str(exe),str(args.order)],check=True,stdout=subprocess.PIPE,text=True,encoding='utf8').stdout)
    if pins!={name:sha((HERE/name).read_bytes()) for name in pins}:raise RuntimeError('source changed')
    if args.order==7:
        if result['min_degree3_graphs']!=236926 or result['metric_filter_graphs']!=8190 or [r['mask'] for r in result['representatives']]!=[120383,120444,120568,122616,378106,380278]:raise RuntimeError('independent inherited order7 control')
    result['source_sha256']=pins;print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
