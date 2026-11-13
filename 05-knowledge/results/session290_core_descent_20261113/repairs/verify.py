from pathlib import Path
import argparse,gzip,csv,subprocess,tempfile,json,hashlib,sys
HERE=Path(__file__).resolve().parent

def need(x,m):
 if not x:raise RuntimeError(m)
def endpoint(rows):
 n=len(rows);d=[[0]*n for _ in range(1<<n)]
 incoming=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)]
 for mask in range(1,1<<n):
  todo=mask
  while todo:
   bit=todo&-todo;todo-=bit;v=bit.bit_length()-1;rest=mask^bit
   if not rest:d[mask][v]=1;continue
   inc=incoming[v]&rest;s=0
   while inc:
    bb=inc&-inc;inc-=bb;s+=d[rest][bb.bit_length()-1]
   d[mask][v]=s
 return d

def inflate(a):
 base=[2,4,9,19,7];labs=[0]*a+[1,2,3,4]
 rows=[sum(1<<j for j,w in enumerate(labs) if (v==w and i<j) or(v!=w and base[v]>>w&1)) for i,v in enumerate(labs)]
 return labs,rows

def check(a):
 labs,rows=inflate(a);n=len(rows);full=(1<<n)-1
 ends=endpoint(rows);rev=[full^(1<<v)^rows[v] for v in range(n)];starts=endpoint(rev)
 O=[labs.index(3),labs.index(4)];H=sum(ends[full]);E=sum(ends[full][v] for v in O)
 U=sum(starts[full][v] for v in O);W=0
 for mask in range(1,full):
  rest=full^mask;U+=sum(ends[mask])*sum(starts[rest][v] for v in O)
  W+=ends[mask][O[0]]*starts[rest][O[1]]+ends[mask][O[1]]*starts[rest][O[0]]
 out=sum(1<<v for v in O);Q=[r|(0 if i in O else 1<<n) for i,r in enumerate(rows)]+[out]
 qend=endpoint(Q)[-1];R=-3+12*2**a+8*3**a
 need(H==3*2**a+3,'source H');need(E==2,'selected endpoint invoice');need([ends[full][v] for v in O]==[1,1],'selected end vector')
 need(U==-4+12*2**a+8*3**a,'full selected slope');need(H-W==1,'unit defect');need(sum(qend)==R-2,'literal precursor H');need(qend[-1]==H-2,'new root endpoint sidecar')
 if a==1:need(qend==[9,8,6,6,7,7],'precursor full endpoint boundary')
 return {'a':a,'source_order':n,'H_source':H,'U':U,'delta':H-W,'source_selected_ends':[ends[full][v] for v in O],'precursor_selected_ends':[qend[v] for v in O],'precursor_H':sum(qend),'target':R,'root_end':qend[-1]}

def main():
 p=argparse.ArgumentParser();p.add_argument('--census');p.add_argument('--fullscan',action='store_true');args=p.parse_args()
 controls=[check(a) for a in range(1,8)]
 scan=(HERE/'SCAN.txt').read_text();need(scan.rstrip().endswith('TOTAL 191536 sources 68 tests 14892 secondhits 24 thirdhits 0'),'declared finite scan summary')
 lines=[s for s in scan.splitlines() if s and not s.startswith('TOTAL')];need(len(lines)==24,'complete retained second-parameter hits')
 for s in lines:
  z=s.split();need(z[1]=='factor' and z[3]=='slopes' and z[13]=='third','typed near-match record');G=int(z[2]);ss=list(map(int,z[4:13]));need(G in[1,3] and 1<=sum(ss)<=3,'finite domain');need(int(z[14])*G!=6597,'third parameter rejects every retained chart')
 full=False;census_sha=None
 if args.fullscan:
  census=Path(args.census) if args.census else HERE.parents[1]/'blue_quotient_20261020b/n9_classes.tsv.gz'
  need(census.exists(),'supply inherited census using --census');data=census.read_bytes();census_sha=hashlib.sha256(data).hexdigest()
  codes=[r['code'] for r in csv.DictReader(gzip.decompress(data).decode().splitlines(),delimiter='\t')];need(len(codes)==191536,'inherited class count')
  with tempfile.TemporaryDirectory() as td:
   binary=Path(td)/'scan';subprocess.run(['g++','-O3','-std=c++17',str(HERE/'scan.cpp'),'-o',str(binary)],check=True)
   r=subprocess.run([str(binary)],input='\n'.join(codes)+'\n',capture_output=True,text=True,check=True)
   need(r.stdout+r.stderr==scan,'entire independent replay matches saved finite certificate');full=True
 out={'status':'PROVED endpoint corollary; FINITE-EXACT declared chart scan; OPEN repair','source':[2,4,9,19,7],'selected':[3,4],'clock':0,'controls':controls,'full_census_replayed':full,'census_sha256':census_sha,'scan_sha256':hashlib.sha256(scan.encode()).hexdigest(),'source_count':68,'slope_presentations':14892,'second_hits':24,'third_hits':0,'unclosed_family_ledger_unchanged':473}
 name='RESULTS.optimized.json' if not __debug__ else 'RESULTS.json';(HERE/name).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
