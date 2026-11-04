"""Compile and run exact list controls, without installing software."""
import hashlib,json,shutil,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
def main():
    source={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ('run.py','check.cpp')}
    compiler=shutil.which('g++')
    if not compiler:raise RuntimeError('existing g++ required')
    with tempfile.TemporaryDirectory(prefix='moser-seven-lists-') as tmp:
        exe=Path(tmp)/'lists.exe'
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'check.cpp'),'-o',str(exe)],check=True,capture_output=True)
        result=subprocess.run([str(exe)],check=True,capture_output=True,text=True,encoding='utf8')
        data=json.loads(result.stdout)
    if source!={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in source}:raise RuntimeError('source changed')
    data['source_sha256']=source
    print(json.dumps(data,indent=2,sort_keys=True))
if __name__=='__main__':main()
