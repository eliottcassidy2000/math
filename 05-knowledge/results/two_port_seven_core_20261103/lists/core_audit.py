"""Independent seven-core universe audit, with full relabelling orbit unions."""
import hashlib,json,shutil,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
def main():
    source={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ('core_audit.py','core_audit.cpp')}
    compiler=shutil.which('g++')
    if not compiler:raise RuntimeError('existing g++ required')
    with tempfile.TemporaryDirectory(prefix='unit-seven-core-audit-') as tmp:
        exe=Path(tmp)/'audit.exe'
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'core_audit.cpp'),'-o',str(exe)],check=True,capture_output=True)
        r=subprocess.run([str(exe)],check=True,capture_output=True,text=True,encoding='utf8')
    if source!={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in source}:raise RuntimeError('source changed')
    data=json.loads(r.stdout);data['source_sha256']=source;print(json.dumps(data,indent=2,sort_keys=True))
if __name__=='__main__':main()
