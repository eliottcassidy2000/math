"""Independent per-record H, delta and U replay of all order10 unit seeds."""
from pathlib import Path
import os,subprocess,json,hashlib
from census_nine import read,save,need
HERE=Path(__file__).resolve().parent

def run():
    d=read('unit_ten.json.gz');temp=Path(os.environ['TEMP']);inp=temp/'math-joint-unit10-interfaces.txt';exe=temp/'math-joint-unit10-interfaces.exe';src=HERE/'audit_unit_slopes.cpp'
    with inp.open('w') as f:
        for r in d['records']:
            values=[r['H'],*r['rows'],len(r['marks'])]
            for m in r['marks']:values.extend([m['selected'],m['U']])
            f.write(' '.join(map(str,values))+'\n')
    subprocess.run(['g++','-O3','-std=c++17',str(src),'-o',str(exe)],check=True)
    result=json.loads(subprocess.check_output([str(exe),str(inp)],text=True));need(result['sources']==465023 and result['marked_interfaces']==478684,'every source and marked orbit replayed')
    need(result['sum_U']==sum(m['U'] for r in d['records'] for m in r['marks']),'whole-bank slope checksum')
    result.update(scope='Every stored ordinary strong order10 unit source and every marked orbit: individual H, delta and U values recounted independently in C++ using selected boundary vectors, not the Python full W-matrix engine.',input_sha256=hashlib.sha256((HERE/'unit_ten.json.gz').read_bytes()).hexdigest(),cpp_source_lf_sha256=hashlib.sha256(src.read_text().encode()).hexdigest())
    save('audit_unit_slopes.json.gz',result);print(json.dumps(result),flush=True)
if __name__=='__main__':run()
