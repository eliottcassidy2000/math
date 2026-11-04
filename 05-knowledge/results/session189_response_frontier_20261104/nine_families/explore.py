from pathlib import Path
from itertools import combinations,permutations
import json,subprocess
D=Path(__file__).resolve().parent

def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def mods(Q):
 n=len(Q);out=[]
 for s in range(1,(1<<n)-1):
  if s.bit_count()<2:continue
  if all(not(s>>i&1) and (Q[i]&s) in (0,s) or s>>i&1 for i in range(n)):out.append(s)
 return out
def endpoints(Q):
 n=len(Q);joint=[[0]*n for _ in range(n)]
 for start in range(n):
  dp=[[0]*n for _ in range(1<<n)];dp[1<<start][start]=1
  for m in range(1,1<<n):
   for last in range(n):
    value=dp[m][last]
    if not value:continue
    todo=Q[last]&~m
    while todo:
     bit=todo&-todo;todo-=bit;j=bit.bit_length()-1;dp[m|bit][j]+=value
  joint[start]=dp[-1]
 return joint

def main():
 models=[]
 for mask in range(16):
  steps=[d if mask>>(d-1)&1 else 9-d for d in range(1,5)]
  rows=[sum(1<<((i+d)%9) for d in steps) for i in range(9)]
  models.append(dict(name=f'circulant_{mask:02d}',construction=dict(step_choice_mask=mask,steps=steps),rows=rows))
 T=[sum(1<<((i+d)%9) for d in range(1,5)) for i in range(9)]
 U=T[:]
 for i in range(4):
  for j in range(i+1,4):U[i]^=1<<j;U[j]^=1<<i
 W=[sum(1<<j for j in range(i+1,8))|((i%2)<<8) for i in range(8)]+[sum(1<<i for i in range(0,8,2))]
 models += [dict(name='critical_U9',construction='Converse the first four consecutive vertices of carousel T9.',rows=U),dict(name='critical_W9',construction='Transitive vertices0..7; vertex8 dominates even positions and is dominated by odd positions.',rows=W)]
 for q in models:
  Q=q['rows'];q['nontrivial_modules']=mods(Q);q['prime']=not q['nontrivial_modules'];J=endpoints(Q);q['joint_endpoints']=J;q['S']=list(map(sum,J));q['E']=[sum(J[i][j] for i in range(9)) for j in range(9)];q['p']=[(a+b)%2 for a,b in zip(q['S'],q['E'])]
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'anchors.cpp'),'-o',str(D/'anchors')],check=True)
 inp=str(len(models))+'\n'+''.join(' '.join(map(str,q['rows']))+'\n' for q in models)
 raw=subprocess.check_output([str(D/'anchors')],input=inp.encode()).decode();(D/'anchor_values.txt').write_text(raw)
 for q,line in zip(models,raw.splitlines()):
  values=list(map(int,line.split()));q['start_run_skews']=[values[9*i:9*(i+1)] for i in range(9)];q['one_port_base']=values[81:90];q['one_port_double']=[values[90+9*i:90+9*(i+1)] for i in range(9)];q['patterns']=[]
  for V in range(512):
   d=[(q['one_port_base'][i]+2*sum(q['one_port_double'][i][j] for j in range(9) if V>>j&1))%4 for i in range(9)]
   t=[(q['S'][i]-q['E'][i]+2*sum(q['start_run_skews'][i][j] for j in range(9) if V>>j&1))%4 for i in range(9)]
   a=[(2*t[i]+2*d[i])%8 if not q['p'][i] else None for i in range(9)]
   q['patterns'].append(dict(nontrivial_child_mask=V,first=d,second=a))
  q['anchor_blind_patterns']=[p['nontrivial_child_mask'] for p in q['patterns'] if not any(q['p']) and not any(p['first']) and not any(p['second'])]
  q['sigma']=[(q['S'][i]-q['E'][i]+2*sum(q['start_run_skews'][i]))%4 for i in range(9)]
 (D/'families.json').write_text(json.dumps(models,separators=(',',':'))+'\n')
 print(json.dumps([dict(name=q['name'],rows=q['rows'],prime=q['prime'],p=q['p'],sigma=q['sigma'],blind_count=len(q['anchor_blind_patterns']),blind=q['anchor_blind_patterns']) for q in models],indent=2))
if __name__=='__main__':main()
