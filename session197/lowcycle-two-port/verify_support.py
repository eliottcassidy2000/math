"""Independent Python control checks and C++ closure certificate reproduction."""
from control import *
from pathlib import Path
import hashlib,subprocess,tempfile
P=Path(__file__).resolve().parent


def literal_cycles(r):
 return sum(all(r[u]>>v&1 for u,v in zip((0,)+p,p+(0,))) for p in permutations(range(1,len(r)))) if len(r)>=3 else 0

def main():
 data=machine();need(data==json.loads((P/'control.json').read_text()),'finite control certificate')
 # Check the serialized C++ transition table against the audited combinatorial generator.
 seed=brute([2,0],[0,1]);code='const int hc[]={'+','.join(str(i) for i,a in enumerate(shapes(2)) if len(a)==1)+'};\n'
 code+='const State seed={{'+','.join(map(str,seed))+'}};\nconst vector<vector<pair<int,int>>> tr[12]={\n'
 for S in range(4):
  for d in range(3):code+='{'+','.join('{'+','.join('{'+str(j)+','+str(m)+'}' for j,m in row)+'}' for row in transitions(2,S,d))+'},\n'
 code+='};\n';need(code==(P/'automaton.inc').read_text(),'serialized physical automaton')
 physical_checks=0
 for n in range(2,5):
  pairs=list(combinations(range(n),2))
  for bits in range(1<<len(pairs)):
   r=[0]*n
   for i,(u,v) in enumerate(pairs):
    if bits>>i&1:u,v=v,u
    r[u]|=1<<v
   for ports in permutations(range(n),2):
    v=brute(r,ports)
    for S in range(4):
     out,_=physical(r,list(ports),S,2);C=literal_cycles(out)
     need(cycle_readout(v,S)==C,'literal cycle readout')
     need(min(3,cycle_readout(tuple(min(3,x) for x in v),S))==min(3,C),'cycle cap')
     physical_checks+=1
 receipts=[];witnesses=0
 with tempfile.TemporaryDirectory() as tmp:
  exe=Path(tmp)/'support'
  subprocess.run(['g++','-std=c++17','-O2',str(P/'support.cpp'),'-o',str(exe)],check=True)
  for X in (21,173,501):
   got=json.loads(subprocess.check_output([str(exe),str(X),'5000000'],text=True))
   if X==21:
    py=bounded_python(X)
    need(all(got[k]==py[k] for k in ('states','admissible_edges','cycle_rejected_edges')),'independent Python closure')
    need(sorted(map(int,got['witnesses']))==py['support'],'Python support')
   else:need(got==json.loads((P/f'support{X}.json').read_text()),'frozen closure receipt')
   need(got['complete'] and got['processed']==got['states'],'incomplete closure')
   for h,w in got['witnesses'].items():
    v=seed;rows=[2,0];ports=[0,1];C=0
    for S,d in w['word']:
     C=cycle_readout(v,S);need(C<=2,'witness forbidden prefix')
     v=step(v,2,S,d);rows,ports=physical(rows,ports,S,d)
    need(H(v,2)==int(h) and C==w['C'] and C in (1,2),'uncapped witness readout');witnesses+=1
   receipts.append({k:v for k,v in got.items() if k!='witnesses'})
 return dict(status='PASS',raw_control_states=len(data['raw_states']),edge_control_classes=len(data['quotient']),
             refinement_rounds=data['refinement_rounds'],literal_cycle_readout_checks=physical_checks,
             uncapped_physical_witness_replays=witnesses,closures=receipts)

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
