from pathlib import Path
import json,subprocess
from explore import op,aug
D=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def main():
 Q=[480,449,387,263,15,30,60,120,240];queries=[(0,7),(8,264)];cases=[aug(Q if not side else op(Q),a,c) for a,c in queries for side in range(2)]
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'even_runs.cpp'),'-o',str(D/'even_runs')],check=True);raw=subprocess.check_output([str(D/'even_runs')],input=('4\n'+''.join(' '.join(map(str,T))+'\n' for T in cases)).encode()).decode();(D/'even_values.txt').write_text(raw);values=[list(map(int,l.split())) for l in raw.splitlines()];need(len(values)==880,'all220 run-excess shapes per side')
 bank=[]
 for qi,query in enumerate(queries):
  S=values[qi*440:qi*440+220];U=values[qi*440+220:qi*440+440];need([r[:9] for r in S]==[r[:9] for r in U],'aligned shape counts');delta=[s[9]-u[9] for s,u in zip(S,U)];g=sum((-2)**sum(k-1 for k in r[:9])*d for r,d in zip(S,delta))%16;bits=[]
  for i in range(9):
   k=[1]*9;k[i]=2;ind=next(j for j,r in enumerate(S) if r[:9]==k);bits.append(delta[ind]%2)
  need(g==8*qi and bits==[0,0,0,0,1,0,0,1,0],'two affine even-profile local identities')
  bank.append(dict(query=list(query),source_augmented_rows=cases[2*qi],target_augmented_rows=cases[2*qi+1],shapes=[r[:9] for r in S],source_counts=[r[9] for r in S],target_counts=[r[9] for r in U],geometric_contrast_mod16=g,q2_derivative_parities=bits))
 (D/'even_bank.json').write_text(json.dumps(dict(rows=Q,queries=bank),separators=(',',':'))+'\n');print(json.dumps(dict(status='PASS',complete_run_shapes=220,exact_counts=880,even_profile_identities=['8(e4+e7)','8(1+e4+e7)'])))
if __name__=='__main__':main()
