"""Literal central-bank overlap, with explicit unbounded linear parameter family."""
import json,subprocess,tempfile
from pathlib import Path
from central import *
from replay import L,conv,fusion,F
from radix import State
from carrier import profile

def run():
 U=circulant(7,0);C=apex(U,21);Cp=apex(U,11);p=profile(C);pp=profile(Cp)
 expected=dict(H=641,h=175,Q=1988,k=0,U=952,V=1211,d=100,e=75)
 if p!=expected or pp!={**expected,'H':643}:raise RuntimeError('core')
 checks=0;examples=[]
 for t in list(range(100))+[10**30,10**100]:
  B=7+25*t
  for A in [170+568*t,172+568*t]:
   N=F(p,A,B)
   if F(p,A+1,B)!=N+2 or F(pp,A,B)!=N+2 or F(pp,A+1,B)!=N+4:raise RuntimeError('overlap')
   checks+=1
   if t<3:examples.append(dict(t=t,A=A,B=B,N=N,interval=[N,N+2,N+4]))
 literal=[]
 for core in [C,Cp]:
  for A in [170,171]:
   T=fusion(fusion(L(A),core),conv(L(7)));literal.append(dict(A=A,coreH=profile(core)['H'],T=T,expected=F(profile(core),A,7)))
 source=Path(__file__).with_name('literal_dp.cpp')
 with tempfile.TemporaryDirectory(prefix='s276-overlap-') as d:
  exe=Path(d)/'dp';subprocess.run(['g++','-O2','-std=c++17',str(source),'-o',str(exe)],check=True)
  inputs=''.join(str(len(z['T']))+' '+' '.join(map(str,z['T']))+'\n' for z in literal)
  values=list(map(int,subprocess.run([str(exe)],input=inputs,text=True,capture_output=True,check=True).stdout.split()))
  intersections=[tuple(x&y for x,y in zip(literal[j]['T'],literal[j+1]['T'])) for j in [0,2]]
  raw=''.join(str(len(T))+' '+' '.join(map(str,T))+'\n' for T in intersections)
  common=[]
  for T,line in zip(intersections,subprocess.run([str(exe),'--path'],input=raw,text=True,capture_output=True,check=True).stdout.splitlines()):
   num,path=line.split(':');path=list(map(int,path.split()));num=int(num)
   if num<=0 or sorted(path)!=list(range(len(T))) or not all(T[a]>>b&1 for a,b in zip(path,path[1:])):raise RuntimeError('common frame')
   common.append(dict(count=num,path=path));checks+=1
 for z,v in zip(literal,values):
  if v!=z['expected']:raise RuntimeError(('literal path DP',z,v))
  z['literalH']=v;checks+=1
 print(json.dumps(dict(status='PASS',checks=checks,core641=C,core643=Cp,examples=examples,literal=literal,common_frames=common,scope='all t theorem plus finite algebra checks, four independent literal nineteen-vertex path DPs'),indent=2))
if __name__=='__main__':run()
