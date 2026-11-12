"""Run the independent catalogue executable and compare every aggregate."""
from pathlib import Path
import os,subprocess,gzip,json,hashlib
from census_nine import read,save,need
HERE=Path(__file__).resolve().parent

def run():
    exe=Path(os.environ['TEMP'])/'math-joint-audit-catalogue.exe'
    source=HERE/'audit_catalogue.cpp'
    subprocess.run(['g++','-O3','-std=c++17',str(source),'-o',str(exe)],check=True)
    inp=Path(os.environ['TEMP'])/'math-tourn10-mckay.txt.gz'
    output=Path(os.environ['TEMP'])/'math-joint-audit-catalogue.json'
    with output.open('wb') as f:
        process=subprocess.Popen([str(exe),'10'],stdin=subprocess.PIPE,stdout=f)
        with gzip.open(inp,'rb') as g:
            while chunk:=g.read(1<<20):process.stdin.write(chunk)
        process.stdin.close();need(process.wait()==0,'independent executable completion')
    d=json.loads(output.read_text());baseline=read('unit_ten.json.gz')
    mapping={'classes':'ordinary_classes','strong_classes':'strong_classes','unit_sources':'all_unit_sources','unit_masks':'all_unit_masks','strong_unit_sources':'strong_unit_sources'}
    for a,b in mapping.items():need(d[a]==baseline[b],'independent aggregate '+a)
    need(d['H_histogram']==baseline['H_histogram'],'complete independent H distribution')
    need(d['strong_unit_masks']==sum(r['Aut']*len(r['marks']) for r in baseline['records']),'independent marked presentation total')
    need(d['classes']==9733056 and d['strong_classes']==9355949 and len(d['H_histogram'])==7663 and max(map(int,d['H_histogram']))==15745,'inherited complete census controls')
    d['scope']='Independent C++ replay of every ordinary order-ten catalogue class, full H distribution, strong connectivity and every selected unit mask. Python marking orbits checked separately against actual automorphisms.'
    d['input_sha256']=hashlib.sha256(inp.read_bytes()).hexdigest();d['cpp_source_lf_sha256']=hashlib.sha256(source.read_text().replace('\r\n','\n').encode()).hexdigest()
    save('audit_catalogue.json.gz',d);print(json.dumps({k:v for k,v in d.items() if k!='H_histogram'}),flush=True)
if __name__=='__main__':run()
