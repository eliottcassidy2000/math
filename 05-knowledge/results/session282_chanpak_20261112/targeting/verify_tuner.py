#!/usr/bin/env python3
"""Scalar-neutral physical receiver tuning and its exact sandwich price."""
from pathlib import Path
import runpy,contextlib,io,json
with contextlib.redirect_stdout(io.StringIO()):v=runpy.run_path(str(Path(__file__).with_name('verify.py')))
ck=v['ck'];profile=v['profile'];conv=v['converse'];join=v['join'];tt=v['tt'];tape=v['tape'];mul=v['mul'];C3=v['C3']
def matrix(p):
 h=p['h'];return [[h,p['U']-h,p['d']-h,p['H']-h],[0,p['Q'],0,p['V']-h],[0,0,p['k'],p['e']-h],[0,0,0,h]]
def read(M):
 h=M[0][0];return dict(H=h+M[0][3],h=h,Q=M[1][1],U=h+M[0][1],V=h+M[1][3],d=h+M[0][2],e=h+M[2][3],k=M[2][2])
def fuse(A,B):
 n=len(A);m=len(B);out=[0]*(n+m-1)
 for T,ids in((A,list(range(n))),(B,[0]+list(range(n,n+m-1)))):
  for i in range(len(T)):
   for j in range(len(T)):
    if T[i]>>j&1:out[ids[i]]|=1<<ids[j]
 for i in range(1,n):
  for j in range(n,n+m-1):out[i]|=1<<j
 return tuple(out)
X=[join((0,),tt(6),tape((0,1)*3)),join((0,),C3,C3,tape((0,0,0,1,1,1)))]
contexts=[(0,),(2,0),(0,1),C3];ps=[profile(T)for T in contexts]
a=9;delta=-3072;controls=0;literal=[]
for converse in(False,True):
 twins=[conv(T)if converse else T for T in X];px=[profile(T)for T in twins]
 for i,(L,pL)in enumerate(zip(contexts,ps)):
  for j,(R,pR)in enumerate(zip(contexts,ps)):
   got=[read(mul(mul(matrix(pL),matrix(p)),matrix(pR)))for p in px]
   DH=got[1]['H']-got[0]['H'];DV=got[1]['V']-got[0]['V']
   wantH=delta*(pL['U']-pL['h'])*pR['V'] if converse else delta*pL['U']*(pR['V']-pR['h'])
   wantV=delta*pL['Q']*pR['V'] if converse else delta*pL['Q']*(pR['V']-pR['h'])
   ck((DH,DV)==(wantH,wantV),'sandwich difference law')
   ck(got[1]['h']==got[0]['h'] and got[1]['e']==got[0]['e'],'remaining receiver differences')
   neutral=pL['U']==pL['h'] if converse else pR['V']==pR['h']
   ck((DH==0)==neutral,'physical neutrality iff')
   if neutral:
    wanted=pL['h']*a*pR['H'] if converse else pL['H']*a*pR['h']
    ck(got[0]['H']==got[1]['H']==wanted and wanted%9==0,'neutral content price')
   if(i,j)in((0,0),(1,0),(0,2),(3,0),(0,3),(2,1)):
    actual=[profile(fuse(fuse(L,T),R))for T in twins]
    ck(actual==got,'literal sandwich full carrier')
    literal.append(dict(converse=converse,left=i,right=j,H=[p['H']for p in actual],V=[p['V']for p in actual]))
   controls+=1
# Direct positive tuner with B=C3.
for p in [read(mul(matrix(profile(conv(T))),matrix(profile(C3))))for T in X]:
 ck(p['H']==27 and p['h']==9 and p['e']==0,'positive tuner common coordinates')
print(json.dumps(dict(status='PASS',exact_matrix_sandwiches=controls,literal_sandwiches=literal,delta=delta,content_divisor=9,scope='common fusion contexts only'),indent=2,sort_keys=True))
