"""Independent C++ search for the last two binary-head interfaces."""
from pathlib import Path
import os,gzip,json,subprocess,hashlib,shutil
import numpy as np
from census_nine import need,save,tables,defect,ib
HERE=Path(__file__).resolve().parent

def run():
    temp=Path(os.environ['TEMP']);packed=temp/'math-tourn10-mckay.txt.gz';plain=temp/'math-joint-tourn10-mckay.txt'
    if not plain.exists():
        with gzip.open(packed,'rb') as a,plain.open('wb') as b:shutil.copyfileobj(a,b)
    need(plain.stat().st_size==9733056*46,'complete plain catalogue encoding')
    exe=temp/'math-joint-head-target.exe';source=HERE/'head_target_ten.cpp'
    subprocess.run(['g++','-O3','-std=c++17',str(source),'-o',str(exe)],check=True)
    output=temp/'math-joint-head-target.jsonl'
    with output.open('w') as f:subprocess.run([str(exe),str(plain)],stdout=f,check=True)
    records=[json.loads(s) for s in output.read_text().splitlines()];hits=[r for r in records if r['kind']=='hit'];summary=records[-1]
    for r in hits:
        rows=r['rows'];O=[i for i in range(10) if r['selected']>>i&1];h,u,w,_=tables(np.array(rows,dtype=np.int64))
        need(int(h)==r['H'] and defect(h,w,O)==r['defect'] and sum(int(u[i]) for i in O)==r['U'],'independent Python endpoint replay')
        ih,iu,iw,*_=ib.tc.path_tables(rows);need(ih==int(h) and iu==u.tolist() and iw==w.tolist(),'second independent Python whole interface')
    result={'scope':'Deterministic prefix of McKay ordinary order10 catalogue, stopping only after finding both specified endpoint interface types or exhausting the catalogue. Every hit replayed by two independent Python endpoint engines. No absence claim for an unfinished prefix.',
            'source_gzip_sha256':hashlib.sha256(packed.read_bytes()).hexdigest(),'cpp_source_lf_sha256':hashlib.sha256(source.read_text().encode()).hexdigest(),'summary':summary,'records':hits}
    save('head_target_ten.json.gz',result);print(json.dumps(result),flush=True)
if __name__=='__main__':run()
