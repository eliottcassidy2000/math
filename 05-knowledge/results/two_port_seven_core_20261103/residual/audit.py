"""Compile the producer parity function with an independent exact-DP audit."""
import hashlib,json,shutil,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def sha(raw):return hashlib.sha256(raw).hexdigest()
def main():
    old=HERE.parent.parent/'two_port_beyond_nine_20261103/census/census.cpp'
    paths=[HERE/'check.cpp',HERE/'audit.cpp',HERE/'audit.py',old]
    pins={p.relative_to(ROOT).as_posix():sha(p.read_bytes()) for p in paths}
    # Use the actual producer function, not an independently transcribed copy.
    prefix=(HERE/'check.cpp').read_text().split('int main(',1)[0]
    prefix=prefix.replace('../../two_port_beyond_nine_20261103/census/census.cpp',old.as_posix())
    with tempfile.TemporaryDirectory(prefix='module-parity-audit-') as tmp:
        tmp=Path(tmp);src=tmp/'check.cpp';exe=tmp/'check.exe'
        src.write_text(prefix+(HERE/'audit.cpp').read_text(),encoding='utf8')
        compiler=shutil.which('g++')
        if not compiler:raise RuntimeError('g++ required')
        subprocess.run([compiler,'-O3','-std=c++17',str(src),'-o',str(exe)],check=True,capture_output=True)
        result=json.loads(subprocess.run([str(exe)],check=True,stdout=subprocess.PIPE,text=True).stdout)
    if result['integer_endpoint_comparisons']!=2123:raise RuntimeError('audit universe')
    if pins!={p.relative_to(ROOT).as_posix():sha(p.read_bytes()) for p in paths}:raise RuntimeError('source changed')
    result['source_sha256']=pins;print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
