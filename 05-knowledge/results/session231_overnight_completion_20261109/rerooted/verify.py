from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import json
from explore import *
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def strong(rows):
 full=(1<<len(rows))-1
 for a in range(len(rows)):
  seen=1<<a
  while True:
   z=seen
   for u in range(len(rows)):
    if seen>>u&1:z|=rows[u]
   if z==seen:break
   seen=z
  if seen!=full:return False
 return True
def induced(rows,vs):return [sum(((rows[u]>>v)&1)<<i for i,v in enumerate(vs)) for u in vs]
out=[]
bare=[]
for k in range(4):
 R=4096**k;J=(R-1)//4095;A=7516*R+343260*J+45
 D=(245051264*R+1602510*k+1118023)//24843;F=D-A-45
 for j in range(k+1):
  rows,blocks,names=family(k,j);r=names.index('r');v=names.index('v');tail=names.index('trail')
  bs=[[u for u in b if u!=v] for b in blocks]
  z=count(rows,bs,[(None,r),(r,None)])
  need(z==F-2*j,'deleted-v scalar')
  state=(z,count(rows,bs,[(None,r)]),count(rows,bs,[(r,None),(None,None)]),count(rows,bs,[(None,None)]))
  need(state==(F-2*j,0,A+45,45),'deleted-v full sidecar')
  reduced=[u for u in range(len(rows)) if u not in [v,tail]]
  need(strong(induced(rows,reduced)),'reduced strongness')
  bs=[[u for u in b if u not in [v,tail]] for b in blocks]
  need(count(rows,bs,[(None,r),(r,None)])==z,'remove sink')
  excluded=[names.index(x) for x in ['lead','zero','trail','v']]
  bs=[[u for u in b if u not in excluded] for b in blocks]
  B=F-2*(1440*R-577-2383200*J)
  need(24843*B==13617152*R+1602510*k+717259,'bare closed formula')
  Ab=878*R+171630*J
  barestate=(count(rows,bs,[(None,r),(r,None)]),count(rows,bs,[(None,r)]),count(rows,bs,[(r,None),(None,None)]),count(rows,bs,[(None,None)]))
  need(barestate==(B-2*j,45-(14 if j==k else 0),Ab+45,45),'bare full state')
  need(strong(induced(rows,[u for u in range(len(rows)) if u not in excluded])),'bare strongness')
 bare.append({'k':k,'upper':F-2*(1440*R-577-2383200*J),'order':12*k+9})
 out.append({'k':k,'F':F,'A':A,'D':D,'strong_order':12*k+11,'interval':[F-2*k,F]})
# Stored literal controls are regenerated, not merely trusted.
import literal
ns={'state':literal.state}
rows,_,names=family(0,0);s=states(0,0)
for a,name in enumerate(names):need(s[name]==ns['state'](rows,a),'literal rooted state '+name)
vs=[u for u,name in enumerate(names) if name!='v'];g=induced(rows,vs)
need(ns['state'](g,len(g)-1)==(2303,0,7606,45),'literal deleted-v state')
s0,s1=states(1,0),states(1,1)
diffs={}
for name,x in s0.items():
 delta=tuple(v-u for u,v in zip(x,s1[name]));diffs[name]=delta
 dD,dE,dU,dd=delta;dJ=dD-dE+dU;dV=dU-dd
 if name in ['r','v']:need(dJ==-2 and dV==0,'closing port gap')
 elif name=='trail':need((dJ,dV)==(8,-2),'trailing port parity')
 else:need(dJ < -2 and dV<0,'all-gadget gap obstruction '+name)
print(json.dumps({'status':'PASS','reduced_intervals':out,'bare_intervals':bare,'literal_rooted_states':s,'k1_root_differences':diffs},indent=2))
