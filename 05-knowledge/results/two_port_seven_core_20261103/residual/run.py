"""Complete n10 endpoint mod2 / ordered cover mod4 necessary-condition sieve."""
import argparse,gzip,hashlib,json,shutil,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def sha(raw):return hashlib.sha256(raw).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--limit',type=int);args=parser.parse_args()
    if args.limit is not None and not 1<=args.limit<=9733056:raise RuntimeError('limit')
    paths=[HERE/'run.py',HERE/'check.cpp',HERE.parent.parent/'two_port_beyond_nine_20261103/census/census.cpp']
    sources={p.relative_to(ROOT).as_posix():sha(p.read_bytes()) for p in paths}
    bank=ROOT/'05-knowledge/results/consecutive_flags_order10_20261028';manifest=json.loads((bank/'bank_manifest.json').read_text());pins={}
    with tempfile.TemporaryDirectory(prefix='module-congruence-') as tmp:
        tmp=Path(tmp)
        for name in ('parents9.u64','flags10.varint'):
            raw=(bank/(name+'.gz')).read_bytes();data=gzip.decompress(raw);expected=manifest['artifacts'][name+'.gz']
            if sha(raw)!=expected['sha256'] or sha(data)!=expected['uncompressed_sha256']:raise RuntimeError('bank hash')
            pins[name+'.gz']=sha(raw);(tmp/name).write_bytes(data)
        compiler=shutil.which('g++')
        if not compiler:raise RuntimeError('g++ required')
        exe=tmp/'check.exe';subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'check.cpp'),'-o',str(exe)],check=True,capture_output=True)
        command=[str(exe),str(tmp/'parents9.u64'),str(tmp/'flags10.varint')]
        if args.limit is not None:command.append(str(args.limit))
        result=json.loads(subprocess.run(command,check=True,stdout=subprocess.PIPE,text=True,encoding='utf8').stdout)
    if sources!={p.relative_to(ROOT).as_posix():sha(p.read_bytes()) for p in paths}:raise RuntimeError('source changed')
    result.update(source_sha256=sources,bank_sha256=pins,status='FINITE-EXACT',scope='Necessary modular conditions only; a survivor is not a reconstruction collision')
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
