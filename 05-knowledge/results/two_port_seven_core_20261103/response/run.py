"""Full named tensors at order11, literal paths versus subset products."""
import hashlib,json,shutil,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
def sha(raw):return hashlib.sha256(raw).hexdigest()
def code(n,rows):return sum(((rows[i]>>j)&1)<<k for k,(i,j) in enumerate(( (i,j) for i in range(n) for j in range(i+1,n))))
def main():
    paths=[HERE/'check.cpp',HERE/'run.py',HERE.parent/'self_converse/census.cpp']
    pins={p.relative_to(HERE.parent).as_posix():sha(p.read_bytes()) for p in paths}
    residues={1,3,4,5,9};paley=[sum(1<<j for j in range(11) if (j-i)%11 in residues) for i in range(11)]
    controls=[(11,(1<<55)-1),(11,code(11,paley))]
    with tempfile.TemporaryDirectory(prefix='two-port-eleven-tensors-') as tmp:
        tmp=Path(tmp);inp=tmp/'input.txt';exe=tmp/'check.exe'
        inp.write_text(str(len(controls))+'\n'+'\n'.join(f'{n} {c}' for n,c in controls)+'\n',encoding='ascii')
        compiler=shutil.which('g++')
        if not compiler:raise RuntimeError('g++ required')
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'check.cpp'),'-o',str(exe)],check=True,capture_output=True)
        result=json.loads(subprocess.run([str(exe),str(inp)],check=True,stdout=subprocess.PIPE,text=True,encoding='utf8').stdout)
    if pins!={p.relative_to(HERE.parent).as_posix():sha(p.read_bytes()) for p in paths}:raise RuntimeError('source changed')
    if [x['sidecar_aut'] for x in result['controls']]!=[1,55]:raise RuntimeError('positive automorphism controls')
    if [sum(x['S']) for x in result['controls']]!=[1,95095]:raise RuntimeError('Hamilton controls')
    result.update(source_sha256=pins,universe='Transitive11 and Paley11; every literal old permutation/cut, every ordered endpoint tensor entry, entire4489-query interpolation battery plus9 high-weight queries each')
    print(json.dumps(result,sort_keys=True,separators=(',',':')))
if __name__=='__main__':main()
