"""Complete order9 H221/H223 module-response screen from pinned input."""
from pathlib import Path
from collections import Counter
import csv,gzip,hashlib,json,subprocess,tempfile
from classify_modules import H,strong,converse,two_birth,need
HERE=Path(__file__).resolve().parent
def main():
    source=HERE.parent.parent/'blue_quotient_20261020b/n9_classes.tsv.gz'
    with tempfile.TemporaryDirectory(prefix='uniform-module9-') as td:
        td=Path(td);inp=td/'input.txt';out=td/'hits.txt';exe=td/'filter9.exe'
        count=0
        with gzip.open(source,'rt') as f,inp.open('w') as w:
            for rec in csv.DictReader(f,delimiter='\t'):
                code,dual=int(rec['code']),int(rec['converse_code'])
                if code>dual:continue
                rows=[0]*9;b=35
                for i in range(9):
                    for j in range(i+1,9):
                        if code>>b&1:rows[i]|=1<<j
                        else:rows[j]|=1<<i
                        b-=1
                w.write(' '.join(map(str,[code,dual,int(rec['Aut']),int(rec['SC_H'])]+rows))+'\n');count+=1
        need(count==97144,'complete inherited merged universe')
        subprocess.run(['g++','-O3','-std=c++17',str(HERE/'filter9.cpp'),'-o',str(exe)],check=True)
        with inp.open() as f,out.open('w') as w:subprocess.run([str(exe)],stdin=f,stdout=w,check=True)
        candidates=[];hits=[]
        for line in out.read_text().splitlines():
            code,dual,h,*rr=map(int,line.split());rows=tuple(rr)
            need(H(rows)==h,'fresh independent selected H')
            if not strong(rows):continue
            candidates.append({'code':code,'converse':dual,'h':h,'rows':rows})
            for v in range(9):
                co=two_birth.module_runs(rows,v)
                if co[2]!=102 or any(co[3:]):continue
                hits.append({'code':code,'converse':dual,'h':h,'rows':rows,'module':v,'co':co})
        report={'status':'PASS','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                'scope':'All97144 converse-merged order9 representatives. Order10 remains outside this screen.',
                'strong_candidate_counts':dict(Counter(r['h'] for r in candidates)),
                'strong_candidates':candidates,'templates':hits}
        (HERE/'scan9.json').write_text(json.dumps(report,indent=2)+'\n')
        print('template hits',[(r['h'],r['code'],r['module']) for r in hits])
if __name__=='__main__':main()
