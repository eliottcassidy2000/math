from pathlib import Path
import json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
Q=[24,17,3,6,12]
def build(Q):
 T=[0]*22
 for i in range(5):
  for u in range(4):T[4*i+u]=((15^((1<<(u+1))-1))<<(4*i))|sum(15<<(4*j) for j in range(5) if Q[i]>>j&1)|(1<<20)
 T[20]=1<<21;T[21]=(1<<20)-1
 return T
def phys(a):return sum(15<<(4*i) for i in range(5) if a>>i&1)|(1<<21)
records=[];lines=[]
for rev in [False,True]:
 R=[31^(1<<i)^v for i,v in enumerate(Q)] if rev else Q;T=build(R);queries=[(phys(2),phys(18))];lines.extend([f'22 1',' '.join(map(str,T))]+[f'{a} {c}' for a,c in queries]);records.append({'reverse':rev,'old_rows':T,'physical_query':queries[0],'block_query':[2,18],'outer_quotient':[2,4,1],'root':0})
with tempfile.TemporaryDirectory() as td:
 exe=Path(td)/'physical';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'physical.cpp'),'-o',str(exe)],check=True);vv=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
for r,v in zip(records,vv):r['value_mod16']=v
print(json.dumps({'status':'PASS','physical_order':24,'records':records,'contrast_mod16':(vv[0]-vv[1])%16},sort_keys=True,indent=2))
