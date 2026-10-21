#!/usr/bin/env python3
"""Bounded exact checks for native affine phase lifts. Standard library only."""
import argparse,hashlib,json
from collections import defaultdict
from pathlib import Path
HERE=Path(__file__).resolve().parent
CHECKS=0
L=((1,2),(-2,-1),(1,-1)); U=((0,1),(-1,-1),(1,0))
def check(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise ValueError(msg)
def enc(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def source(s):return type(s) is int and s>0 and s%2==1 and s%8!=5
def Q(v):a,b=v;return a*a-a*b+b*b
def R(v):a,b=v;return -b,a-b
def rot(v,c):
 for _ in range(c%3):v=R(v)
 return v
def add(v,w):return v[0]+w[0],v[1]+w[1]
def sub(v,w):return v[0]-w[0],v[1]-w[1]
def branch(s):check(source(s),'primitive source');return 2 if s%8==1 else 1
def lift(s):return 64*s+21+7*2**branch(s)
def edge(s):
 a=branch(s);x=(3*s+1)//2**a;r=0
 while x%8==5:x=(x-1)//4;r+=1
 check(source(x),'primitive target')
 return x,a+2*r,r

def charge(N,s):
 a=b=0
 while s<4*N:
  p,q=(L if s<2*N else U)[s%3];a+=p;b+=q;s=4*s+1
 return a,b

def independent_charge(N,s):
 K=(0,-3,3);T=(-1,-1,2);c0=c1=0
 while s<4*N:
  row=K if s<2*N else T;c0+=row[s%3];c1+=row[(s+1)%3];s=4*s+1
 check((c0-c1)%3==0 and (c0+2*c1)%3==0,'charge coordinates')
 return (c0-c1)//3,-(c0+2*c1)//3

def height(N):
 h=0;p=4
 while p<=12*N:h+=1;p*=4
 return h

def exception_events(N):
 events=[];p=1;b=0
 for j in range(height(N)):
  for K in (2,4):
   s=(K*N-b-1)//p
   if s%2==0:s-=1
   if source(s) and s<4*N:
    a=branch(s);c=21+7*2**a;m=p*s+b;ml=64*m+(c-21)*p+21
    if m<K*N and 64*K*N<=ml:
     check(j>=1 and N<m<4*N,'boundary localization')
     events.append([s,j,K,c,m,ml])
  p*=4;b=4*b+1
 check(len(events)<=2*(height(N)-1),'sparse boundary count')
 check(len({e[0] for e in events})==len(events),'one event per source')
 return events

def defect_from_events(N,s,events):
 c=21+7*2**branch(s);v=(0,0)
 for z,j,K,cc,m,ml in events:
  if z!=s:continue
  w=U[m%3] if K==4 else sub(L[m%3],U[m%3])
  v=sub(v,rot(w,c))
 return v

def route(s,V):
 seen=set();clock=depth=0
 while s not in V:
  check(s not in seen and depth<10000,'finite native route')
  seen.add(s);s,e,r=edge(s);clock+=e;depth+=1
 return s,clock,depth

def source_row(N,s,V):
 u,t,d=route(s,V);v=charge(N,s);r15=s%15;r35=s%35
 check(r15%5==r35%5 and (r15+15*((r35-r15)%7))%105==s%105,'joint 105 register')
 return [s,r15,r35,u,t,d,v[0],v[1],Q(v)]

def inverse_base(u,a):
 check(source(u) and a in (1,2),'inverse branch')
 r=((2 if a==1 else 1)-u)%3;x=u
 for _ in range(r):x=4*x+1
 check((2**a*x-1)%3==0,'inverse integrality')
 return (2**a*x-1)//3,r

def run():
 data=json.loads((HERE/'frontier.json').read_text());V=set(data['frontier']);ranks=dict(data['frontier_ranks'])
 check(len(V)==627 and max(V)==63665,'inherited frontier')
 first=next(s for s in range(1,max(V)+2,2) if source(s) and s not in V)
 check(first==1009 and 64*first+35>max(V),'outside-frontier lift guard')
 for s in V:
  check(source(s),'frontier primitive')
  if s!=1:check(edge(s)[0] in V and ranks[edge(s)[0]]<ranks[s],'frontier grounding')
 # Complete direct small source test, including zero-charge sources.
 phase_rows=[];phase_digest=hashlib.sha256();max_exception_count=0
 for N in range(1,257):
  events=exception_events(N);bad=[];max_exception_count=max(max_exception_count,len(events))
  for s in range(1,4*N,2):
   if not source(s):continue
   ls=lift(s);a=branch(s);c=ls-64*s;u,t,r=edge(s);ul,tl,rl=edge(ls)
   check(source(ls) and branch(ls)==a and ls<256*N,'lift admission/branch')
   check((ul,tl,rl)==(u,t+6,r+3),'exact native edge lift')
   v=charge(N,s);w=charge(64*N,ls)
   check(v==independent_charge(N,s) and w==independent_charge(64*N,ls),'independent shell charges')
   defect=sub(w,rot(v,c));expected=defect_from_events(N,s,events)
   check(defect==expected,'exact boundary defect formula')
   check(Q(defect) in (0,1),'unit or zero defect')
   if defect!=(0,0):bad.append(s)
   phase_digest.update(enc([N,s,ls,c,list(v),list(w),list(defect)]))
  check(sorted(bad)==sorted(e[0] for e in events),'complete exception enumeration')
  phase_rows.append([N,len(events),height(N)])
 # Sparse-event enumeration at large integer scales without large source sweeps.
 sparse=[]
 for N in (2019,4096,16384,10**40+7):
  ev=exception_events(N)
  for s,j,K,c,m,ml in ev:
   check(sub(charge(64*N,lift(s)),rot(charge(N,s),c))==defect_from_events(N,s,ev),'large exact exception')
  sparse.append({'N':N,'events':ev})
 # Exhaust all modest native sources and audit their unique inverse-chain address.
 inverse_digest=hashlib.sha256()
 for s in range(1,65536,2):
  if not source(s):continue
  a=branch(s);u,t,r=edge(s);base,r0=inverse_base(u,a)
  check(r>=r0 and (r-r0)%3==0,'inverse residue precision')
  z=base
  for k in range((r-r0)//3):z=lift(z)
  check(z==s and branch(base)==a,'unique inverse-chain reconstruction')
  inverse_digest.update(enc([s,u,a,r,base,r0,(r-r0)//3]))
 # Explicit outside-frontier simultaneous packets at scale 4096*256.
 packets=[]
 for s in (1009,1015):
  N=256;M=4096*N;ss=[s,lift(s),lift(lift(s))];rows=[source_row(M,z,V) for z in ss]
  check(all(z not in V for z in ss),'packet outside frontier')
  check(len(set(z[3] for z in rows))==1 and len(set(z[5] for z in rows))==1,'actual common terminal and depth')
  check([z[4]-rows[0][4] for z in rows]==[0,6,12],'full integer clock shifts')
  signed=[((1-2*(z[4]%2))*z[6],(1-2*(z[4]%2))*z[7]) for z in rows]
  check(tuple(map(sum,zip(*signed)))==(0,0),'source-authenticated packet zero')
  check(len(set(edge(z)[0] for z in ss))==1,'immediate target equality')
  first_signed=[tuple((1-2*(edge(z)[1]%2))*x for x in charge(M,z)) for z in ss]
  check(tuple(map(sum,zip(*first_signed)))==(0,0),'one-edge moving cancellation')
  packets.append({'base_N':N,'scale_M':M,'rows':rows,'immediate_target':edge(s)[0],'edge_clocks':[edge(z)[1] for z in ss]})
 # Exact missing-boundary hostile, using an actual outside-frontier source.
 N=2019;s=1009;c=49
 bad={'N':N,'source':s,'lift':lift(s),'old_charge':charge(N,s),'rotated_old':rot(charge(N,s),c),'lifted_charge':charge(64*N,lift(s)),'events':[e for e in exception_events(N) if e[0]==s]}
 check(bad['rotated_old']==(2,1) and bad['lifted_charge']==(1,0),'literal phase hostile')
 badrows=[source_row(4096*N,z,V) for z in (s,lift(s),lift(lift(s)))];badtotal=tuple(map(sum,zip(*[(z[6],z[7]) for z in badrows])))
 check(badtotal==(0,-1),'false unguarded triple hostile');bad['triple_rows']=badrows;bad['triple_sum']=badtotal
 # Frontier guard hostile: the edge identity is not an arrival identity here.
 ss=(1,lift(1),lift(lift(1)));frontrows=[source_row(4096,z,V) for z in ss]
 check(ss==(1,113,7281) and len(set(z[3] for z in frontrows))>1,'frontier-stop hostile')
 # Pair debt and polynomial-factor algebra on the entire signed alphabet.
 alphabet={(a,b) for a in range(-2,3) for b in range(-2,3) if Q((a,b)) in (1,3)}
 check(len(alphabet)==12,'signed alphabet')
 for v in alphabet:
  for c in (1,2):
   rv=rot(v,c);r2=rot(v,2*c)
   check(add(add(v,rv),r2)==(0,0),'three-phase zero')
   check(Q(add(v,rv))==Q(v),'pair diagonal exact payment')
   # (x-1)(R^(2c)x-1) has coefficients 1,R^c,R^(2c).
   check(sub((0,0),add(v,r2))==rv,'full clock polynomial factor')
 return {'status':'PROVED scoped lift/circuit identities; finite checks only; universal annulus bound and entry OPEN','small_N_range':[1,256],'phase_rows_schema':['N','exception_count','h_N'],'phase_rows':phase_rows,'phase_source_rows_sha256':phase_digest.hexdigest(),'max_exception_count':max_exception_count,'sparse_large_controls':sparse,'inverse_source_universe':'every primitive source1<=s<65536','inverse_address_sha256':inverse_digest.hexdigest(),'row_schema':['s','mod15','mod35','actual_terminal','full_integer_clock','P_depth','a','b','Q'],'packets':packets,'boundary_hostile':bad,'frontier_hostile_rows':frontrows,'checks':CHECKS,'frontier_sha256':hashlib.sha256((HERE/'frontier.json').read_bytes()).hexdigest()}

def main():
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();result=run();blob=enc(result);out=HERE/'phase_receipt.json'
 if args.write:out.write_bytes(blob)
 else:check(out.read_bytes()==blob,'saved receipt equality')
 print(json.dumps({'status':'PASS','checks':result['checks'],'receipt_sha256':hashlib.sha256(blob).hexdigest(),'max_small_exception_count':result['max_exception_count'],'packet_sources':[[z[0] for z in p['rows']] for p in result['packets']]},sort_keys=True))
if __name__=='__main__':main()
