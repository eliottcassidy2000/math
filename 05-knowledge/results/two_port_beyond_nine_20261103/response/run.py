"""Full ordered one/two/three-path tensors and actual two-port query audit."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
BANK=ROOT/'05-knowledge/results/consecutive_flags_order10_20261028'


def sha(body):return hashlib.sha256(body).hexdigest()


def code(n,rows):
    pairs=((i,j) for i in range(n) for j in range(i+1,n))
    return sum(((rows[i]>>j)&1)<<k for k,(i,j) in enumerate(pairs))


def main():
    paths=[HERE/'check.cpp',HERE/'run.py',HERE.parent/'census/census.cpp']
    pins={p.relative_to(HERE.parent).as_posix():sha(p.read_bytes()) for p in paths}
    parents_raw=(BANK/'parents9.u64.gz').read_bytes();flags_raw=(BANK/'flags10.varint.gz').read_bytes()
    parents=struct.unpack('<191536Q',gzip.decompress(parents_raw))
    data=gzip.decompress(flags_raw);manifest=json.loads((BANK/'bank_manifest.json').read_text(encoding='utf8'))
    for name,raw in [('parents9.u64.gz',parents_raw),('flags10.varint.gz',flags_raw)]:
        if sha(raw)!=manifest['artifacts'][name]['sha256']:raise RuntimeError('bank hash')
    selected={0,100000,4866528,9733055};samples=[];address=delta=shift=index=0
    for byte in data:
        delta|=(byte&127)<<shift
        if byte&128:shift+=7;continue
        address+=delta;delta=shift=0
        if index in selected:
            parent=parents[address>>9];mask=address&511;rows=[0]*10
            for k,(i,j) in enumerate(( (i,j) for i in range(9) for j in range(i+1,9))):
                u,v=(i,j) if parent&(1<<k) else (j,i);rows[u]|=1<<v
            for i in range(9):
                u,v=(9,i) if mask&(1<<i) else (i,9);rows[u]|=1<<v
            samples.append({'index':index,'address':address,'code':code(10,rows)})
        index+=1
    if index!=9733056 or shift or delta or len(samples)!=4:raise RuntimeError('bank decoding')
    cyc=code(5,(6,12,24,17,3));controls=[(5,cyc),(5,cyc^((1<<10)-1))]+[(10,s['code']) for s in samples]
    # Positive symmetry control: a universal source above C3[C3,C3,C3].
    wreath=[0]*10
    for i in range(9):
        for j in range(9):
            if i!=j and ((i//3==j//3 and (j-i)%3==1) or (j//3-i//3)%3==1):wreath[i]|=1<<j
    wreath[9]=511
    controls += [(10,(1<<45)-1),(10,code(10,wreath))]
    with tempfile.TemporaryDirectory(prefix='two-port-tensors-') as tmp:
        tmp=Path(tmp);inp=tmp/'controls.txt';binary=tmp/'check.exe'
        inp.write_text(str(len(controls))+'\n'+'\n'.join(f'{n} {c}' for n,c in controls)+'\n',encoding='ascii')
        compiler=shutil.which('g++')
        if not compiler:raise RuntimeError('available g++ required')
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'check.cpp'),'-o',str(binary)],check=True,capture_output=True)
        result=subprocess.run([str(binary),str(inp)],check=True,stdout=subprocess.PIPE,text=True,encoding='utf8')
    artifact=json.loads(result.stdout)
    if pins!={p.relative_to(HERE.parent).as_posix():sha(p.read_bytes()) for p in paths}:raise RuntimeError('source changed')
    artifact.update(source_sha256=pins,bank_sha256={name:sha(raw) for name,raw in [('parents9.u64.gz',parents_raw),('flags10.varint.gz',flags_raw)]},selected_bank_classes=samples)
    if artifact['controls'][0]['battery']==artifact['controls'][1]['battery']:raise RuntimeError('C5 converse hostile missed')
    if artifact['controls'][-2]['sidecar_aut']!=1 or artifact['controls'][-1]['sidecar_aut']!=81:raise RuntimeError('positive symmetry controls')
    print(json.dumps(artifact,sort_keys=True,separators=(',',':')))


if __name__=='__main__':main()
