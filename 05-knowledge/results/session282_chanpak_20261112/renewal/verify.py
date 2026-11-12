"""Standalone verification of the receiver inverse and literal renewable eight-bank."""
import json
from pathlib import Path
from itertools import combinations,permutations
from receiver_fiber import *

def brute(D,V):
 C=cores(D);h=C[0]['h'];Q=C[0]['Q'];out=set()
 for a in range(max(1,V.bit_length()+1)):
  for p in range(1<<a):
   v0=h*(1+p)
   if v0>V:continue
   for c in C:
    rem=V-v0-(1<<a)*(c['V']-h)
    if rem<0 or rem%((1<<a)*Q):continue
    q=rem//((1<<a)*Q);out.add(fuse(fuse(unit(p,a),c),unit(q))['H'])
 return out

def core_census():
 out=[]
 for n in range(1,6):
  pairs=list(combinations(range(n),2));strongs=five=0
  for code in range(1<<len(pairs)):
   T=[0]*n
   for j,(u,v) in enumerate(pairs):
    if code>>j&1:T[u]|=1<<v
    else:T[v]|=1<<u
   reach=list(T)
   for k in range(n):
    for i in range(n):
     if reach[i]>>k&1:reach[i]|=reach[k]
   if not(n==1 or all((reach[i]|(1<<i))==(1<<n)-1 for i in range(n))):continue
   strongs+=1
   if endpoints(T)[0][-1]!=5:continue
   five+=1;require(n==4,'strong count5 order')
   require(any(all(((T[i]>>j)&1)==((S4[p[i]]>>p[j])&1) for i in range(n) for j in range(n)) for p in permutations(range(n))),'unique strong count5 isomorphism type')
  out.append(dict(n=n,strong_labelled=strongs,H5_strong=five))
 require(out[3]['H5_strong']==24,'four-vertex orbit')
 return out

def outctx(P):
 n=P.bit_length();q=sum((1-((P>>(n-1-i))&1))<<i for i in range(n));return unit(q,n)

def main():
 census=core_census();fibers=0
 for D in (TRI,S4):
  for V in range(1,161):
   require({w['profile']['H'] for w in fiber(D,V)}==brute(D,V),'inverse brute equality');fibers+=1
 for r in range(4096):
  got={w['profile']['H'] for w in fiber(TRI,3*(r+1))}
  require(got=={3*unit(r)['H']+d for d in triangle_deltas(r)},'triangle defect formula')
  require(not any(h+2 in got and h+4 in got for h in got),'triangle no three')
 require({w['profile']['H'] for w in fiber(TRI,6)}=={3,5},'sharp h3 boundary')
 require(fiber(S4,4)==[],'below-h boundary');require({w['profile']['H'] for w in fiber(S4,5)}=={5},'stationary source')
 candidates=fiber(S4,1906);actual=[]
 for w in candidates:
  T=witness(S4,w);require(profile(T)==w['profile'],'literal h5 recovered packet')
  require(profile(graph_fuse(T,unit_graph(0,2)))['H']==w['profile']['H'],'neutral suffix zeros')
  actual.append(dict(**w,rows=T,order=len(T)))
 require({w['profile']['H'] for w in candidates}=={195,741,819,821,823,825},'complete h5 selected fiber')
 selected=[next(w for w in candidates if w['profile']['H']==H) for H in range(819,826,2)]
 C2=unit(6,3);completed=[]
 for w in selected:
  pr=fuse(C2,w['profile']);require(profile(graph_fuse(unit_graph(6,3),witness(S4,w)))==pr,'literal h5 completion');completed.append(pr)
 require([p['H'] for p in completed]==[8443,8445,8447,8449],'renewable four')
 require({(p['h'],p['V'],p['e']) for p in completed}=={(5,15243,5)},'four receiver')
 levels=[];bank=completed
 for k in range(41):
  V=bank[0]['V'];require(V>13 and (V-13)%10==0,'carry phase');P=2*((V-13)//5)
  outs=[fuse(outctx(p),w) for p in (P,P+1) for w in bank];vals=[w['H'] for w in outs]
  require(vals==list(range(vals[0],vals[0]+16,2)),'eight output interval')
  require(outs[4]['V']-outs[0]['V']==-5*(1<<(P.bit_length()-1)),'receiver debt retained')
  levels.append(dict(k=k,V=str(V),P=str(P),start=str(vals[0])));bank=[fuse(unit(14,4),w) for w in bank]
 require(levels[0]['start']=='92862629','h5 eight baseline')
 literal8=[];renewable8=[];C1=unit(2,2)
 for D,w in nine_witnesses():
  T=witness(D,w);p=profile(T);require(p==w['profile'],'literal new nine packet')
  require(any(c['profile']==p for c in fiber(D,139878)),'nine witness recovered inverse')
  literal8.append(dict(core=D,**w,rows=T,order=len(T)));renewable8.append(fuse(C1,p))
  print('literal9',p['H'],'order',len(T),flush=True)
 require([p['H'] for p in renewable8]==list(range(334121,334138,2)),'renewable nine exact counts')
 require({(p['h'],p['V'],p['e']) for p in renewable8}=={(15,559497,15)},'renewable nine receiver')
 cheap=[fuse(E,w['profile']) for D,w in nine_witnesses()]
 require([p['H'] for p in cheap]==list(range(54365,54382,2)),'count-neutral source completion')
 require({(p['h'],p['V'],p['e']) for p in cheap}=={(15,279756,0)},'count-neutral common receiver')
 transport_checks=0
 for n in range(6):
  for x in range(1<<n):
   bank=[fuse(unit(x,n),p) for p in renewable8]
   require([p['H']-bank[0]['H'] for p in bank]==list(range(0,18,2)),'all finite unit controls')
   require(len({(p['h'],p['V'],p['e']) for p in bank})==1,'receiver still common');transport_checks+=1
 orbit={139878*2**j%15 for j in range(4)};require(orbit=={3,6,12,9} and 8 not in orbit and 9 in orbit,'ninth state repairs carry phase')
 eighteen=[];bank=[fuse(unit(14,4),w['profile']) for D,w in nine_witnesses()]
 require([p['H'] for p in bank]==list(range(1173389,1173406,2)),'prepared nine counts')
 require({(p['h'],p['V'],p['e']) for p in bank}=={(15,2238033,15)},'prepared nine receiver')
 for k in range(41):
  V=bank[0]['V'];require(V>33 and (V-33)%30==0,'eighteen phase');P=2*((V-33)//15)
  outs=[fuse(outctx(p),w) for p in (P,P+1) for w in bank];vals=[p['H'] for p in outs]
  require(vals==list(range(vals[0],vals[0]+36,2)),'eighteen exact interval')
  require(outs[9]['V']-outs[0]['V']==-15*(1<<(P.bit_length()-1)),'eighteen receiver defect')
  require(V==2238032*16**k+1,'receiver closed form')
  require(bank[0]['H']==1173389+17904256*(16**k-1)//15+8*k,'baseline closed form')
  eighteen.append(dict(k=k,V=str(V),P=str(P),start=str(vals[0])));bank=[fuse(unit(14,4),w) for w in bank]
 require(eighteen[0]['start']=='667830121109','first eighteen baseline')
 result=dict(status='PASS',strong_core_census=census,inverse_brute_fibers=fibers,triangle_receiver_codes=4096,exact_h5_fiber=actual,renewable_four=completed,h5_eight_levels=levels,literal_nine=literal8,cheap_renewable_nine=cheap,renewable_nine=renewable8,eighteen_levels=eighteen,unit_context_controls=transport_checks,scope='All-order fixed-core inverse; GLOBAL h3 no-three and h5 inverse via SCC/Camion+finite core certificate; renewable h15 nine-bank and infinite eighteen-wide family. No Chan-Pak completeness or repeated width doubling.')
 Path(__file__).with_name('RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
 print('PASS',fibers,4096,len(literal8),transport_checks)
if __name__=='__main__':main()
