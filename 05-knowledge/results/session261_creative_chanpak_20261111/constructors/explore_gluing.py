from pathlib import Path
import importlib.util
P=Path(__file__).resolve().parents[2]/'session260_quasiperfect_chanpak_20261111/constructors/module_ladders.py'
s=importlib.util.spec_from_file_location('ladders',P);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def glue(T,O,S,P):
 a=[i for i in range(len(T)) if O>>i&1];b=[i for i in range(len(S)) if P>>i&1]
 if not T[a[0]]>>a[1]&1:a.reverse()
 if not S[b[0]]>>b[1]&1:b.reverse()
 ea=[i for i in range(len(T)) if i not in a];eb=[i for i in range(len(S)) if i not in b]
 n=2+len(ea)+len(eb);out=[0]*n
 for Z,order in [(T,a+ea),(S,b+eb)]:
  tar=list(range(2))+list(range(2,2+len(ea))) if Z is T else list(range(2))+list(range(2+len(ea),n))
  for i,x in enumerate(order):
   for j,y in enumerate(order):
    if Z[x]>>y&1:out[tar[i]]|=1<<tar[j]
 for i in range(2,2+len(ea)):
  for j in range(2+len(ea),n):out[i]|=1<<j
 return tuple(out),3
A=(20,1,2,7,14);B=(24,1,3,6,14)
seeds=[(A,6),(A,24),(B,6),(B,24)]
for i,(T,O) in enumerate(seeds):
 for j,(S,P) in enumerate(seeds):
  Z,Q=glue(T,O,tuple(list(S)),P); h,w,u=m.data(Z,Q)
  print(i,j,(h,h-w,u))
for T,O in seeds:
 R=tuple(sum(1<<i for i in range(len(T)) if T[i]>>j&1) for j in range(len(T)))
 print('REV',T,O,m.data(T,O),m.data(R,O))
