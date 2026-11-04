"""Independent deletion-parent enumeration and neighbourhood bipartiteness."""
import hashlib,json,shutil,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
def sha(raw):return hashlib.sha256(raw).hexdigest()
def main():
    paths=['independent.py','independent.cpp','order8.json'];pins={name:sha((HERE/name).read_bytes()) for name in paths}
    records=json.loads((HERE/'order8.json').read_text())['representatives']
    compiler=shutil.which('g++')
    if not compiler:raise RuntimeError('existing g++ required')
    with tempfile.TemporaryDirectory(prefix='eight-core-independent-') as tmp:
        tmp=Path(tmp);exe=tmp/'audit.exe';inp=tmp/'representatives.txt'
        inp.write_text(str(len(records))+'\n'+'\n'.join(f"{r['mask']} {r['orbit']}" for r in records)+'\n',encoding='ascii')
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'independent.cpp'),'-o',str(exe)],check=True,capture_output=True)
        out=json.loads(subprocess.run([str(exe),str(inp)],check=True,stdout=subprocess.PIPE,text=True,encoding='utf8').stdout)
    if pins!={name:sha((HERE/name).read_bytes()) for name in paths}:raise RuntimeError('source changed')
    out['source_sha256']=pins;print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
