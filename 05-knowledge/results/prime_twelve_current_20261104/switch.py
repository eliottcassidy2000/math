"""Complete specified skew-conference switching orbit, with exact endpoint counts."""
import argparse,json,subprocess,tempfile
from pathlib import Path
from hashlib import sha256
from collections import Counter
HERE=Path(__file__).resolve().parent
def need(ok,what):
 if not ok:raise RuntimeError(what)
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'switch.json');args=p.parse_args()
 pins={f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['switch.py','switch.cpp']}
 base=[sum(1<<v for v in range(11) if (v-u)%11 in [1,3,4,5,9]) for u in range(11)]+[2047]
 s=[[0 if u==v else 1 if base[u]>>v&1 else -1 for v in range(12)] for u in range(12)]
 need(all(sum(s[u][k]*s[k][v] for k in range(12))==(-11 if u==v else 0) for u in range(12) for v in range(12)), 'conference identity')
 with tempfile.TemporaryDirectory(prefix='conference-orbit-') as t:
  exe=Path(t)/'switch.exe';subprocess.run(['g++','-O3','-std=c++17',str(HERE/'switch.cpp'),'-o',str(exe)],check=True)
  out=subprocess.run([str(exe)],text=True,capture_output=True,check=True).stdout
 records=[]
 for line in out.splitlines():
  x=list(map(int,line.split()));m,h,sp,ep,prime,strong=x[:6]
  need(len(x)==42,'record width')
  records.append(dict(mask=m,H=h,S_parity=sp,E_parity=ep,prime=bool(prime),strong=bool(strong),S=x[6:30:2],E=x[7:30:2],simple_paths_by_vertices=x[30:]))
 need([x['mask'] for x in records]==list(range(2048)),'orbit universe')
 need(all(x['S_parity']^x['E_parity']==4095 for x in records),'coordinatewise endpoint imbalance')
 need(all(sum(x['S'])==sum(x['E'])==x['H'] for x in records),'endpoint sums')
 need(all(x['simple_paths_by_vertices'][:2]==[12,66] and x['simple_paths_by_vertices'][-1]==x['H'] for x in records),'simple path profile')
 result=dict(status='FINITE-EXACT entire declared switching orbit; all cores excluded by endpoint parity',
             base_rows=base,order=12,presentations=2048,prime_presentations=sum(x['prime'] for x in records),
             strong_presentations=sum(x['strong'] for x in records),H_histogram=dict(Counter(x['H'] for x in records)),
             common_Seidel_characteristic_polynomial='(t^2+11)^6',records=records,source_sha256=pins)
 need(pins=={f:sha256((HERE/f).read_bytes()).hexdigest() for f in pins},'source changed')
 args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
 print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
if __name__=='__main__':main()

