#!/usr/bin/env python3
"""All-order binary-return coefficient certificate. No target census is run."""
import hashlib,itertools,json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent

def check(ok,why):
 if not ok:raise RuntimeError(why)


def dp(a):
 n=len(a);N=1<<n;H=[0]*N;H[0]=1;E=[[0]*n for _ in range(N)];S=[[0]*n for _ in range(N)]
 for i in range(n):
  check(0<=a[i]<N and not(a[i]>>i&1),'invalid row')
  for j in range(i):check(((a[i]>>j)^(a[j]>>i))&1,'invalid pair')
 for mask in range(1,N):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v)
    if not rest:E[mask][v]=S[mask][v]=1
    else:
     E[mask][v]=sum(E[rest][u] for u in range(n) if a[u]>>v&1)
     S[mask][v]=sum(S[rest][u] for u in range(n) if a[v]>>u&1)
  H[mask]=sum(E[mask]);check(H[mask]==sum(S[mask]),'start/end mismatch')
 full=N-1;out=[]
 for r in range(n):
  deleted=full^(1<<r)
  L=sum(E[u][r]*H[full^u] for u in range(N) if u>>r&1)
  R=sum(S[u][r]*H[full^u] for u in range(N) if u>>r&1)
  B=sum(H[u]*H[deleted^u] for u in range(N) if u&deleted==u)
  out.append((H[-1],S[-1][r],L,E[-1][r],R,B,int(n==1),H[deleted]))
 return out


def literal(a,r):
 n=len(a);h=s=e=L=R=B=d=0
 def path(p):return all(a[u]>>v&1 for u,v in zip(p,p[1:]))
 for p in itertools.permutations(range(n)):
  if path(p):h+=1;s+=p[0]==r;e+=p[-1]==r
  for k in range(n+1):
   if path(p[:k]) and path(p[k:]):
    L+=bool(k and p[k-1]==r);R+=bool(k<n and p[k]==r)
 for p in itertools.permutations([i for i in range(n) if i!=r]):
  d+=path(p)
  for k in range(n):B+=path(p[:k]) and path(p[k:])
 return h,s,L,e,R,B,int(n==1),d


def state(p):return p[0],p[3],p[4],p[7]


def cone(x):
 h,e,R,d=x;return e,d-e,h-d,R+e-h-d


def tr(x,op,bank):
 h,e,R,d=x
 if op=='A':return h+R-e,h-e,h+R,h
 if op=='K':return h+R-e,d,2*R-d,d
 if op=='D':return h,0,2*R,d
 g,s,L,E,Q,B,delta,dG=bank[op]['sig']
 return g*h-s*e+L*R,E*h-delta*e,Q*h+B*R,h*dG


def build(rows,r,op,bank):
 rows=list(rows);n=len(rows)
 if op=='D':return [x|(1<<n) for x in rows]+[0],r
 if op in ('A','K'):
  return [x|((i!=r)<<n) for i,x in enumerate(rows)]+[1<<r],n if op=='A' else r
 y=bank[op];G=y['rows'];b=y['root'];mask=((1<<len(G))-1)<<n
 out=[x|mask for x in rows]+[x<<n for x in G]
 out[r]^=1<<(n+b);out[n+b]|=1<<r
 return out,n+b


def main():
 raw=json.loads((HERE/'inventory67.json').read_text());bank=[]
 for y in raw:
  sig=tuple(y['signature'])+(y['H_deleted_root'],)
  check(dp(y['rows'])[y['root']]==sig==literal(y['rows'],y['root']),'physical gadget profile')
  g,s,L,E,Q,B,delta,dG=sig;U=Q+E-g-dG;V=B-L
  a=y['rows'];r=y['root'];vertices=[v for v in range(len(a)) if v!=r]
  deleted=[sum((1<<j) for j,v in enumerate(vertices) if a[u]>>v&1) for u in vertices]
  ps=dp(deleted)
  missing=sum(ps[j][2] for j,u in enumerate(vertices) if a[r]>>u&1)
  check(V==missing and U>=0 and L>=g>=E,'universal missing-bin identity')
  check(2*Q+2*B-g-L-2*dG==(g-E)+(L-E)+2*U+2*V>=0,'universal DK coefficient decomposition')
  if a[r]:check(V>=dG,'non-sink gadget lower bound')
  else:check(U==V==0 and g==E==Q==dG and s==delta,'sink gadget identities')
  bank.append({'sig':sig,'rows':y['rows'],'root':y['root']})
 check(len(bank)==67 and len({y['sig'] for y in bank})==67,'67 distinct physical matrices')
 cert=[]
 for i,y in enumerate(bank):
  g,s,L,E,Q,B,delta,dG=y['sig']
  forms=[('K',(Q+B-g-L+E-2*dG,B-L+s-delta,B-L)),
         ('DK',(2*Q+2*B-g-L-2*dG,2*B-L+s,2*B-L))]
  if not delta and E==dG:forms.insert(0,('direct',(Q+B-g-L-dG,B-L+s,B-L)))
  for name,(A,B0,C) in forms:
   T=3*A+B0
   if A>=0 and T>=0 and 3*T-abs(C)>=2:
    direction='lower';bound=3*T-abs(C)
   elif A<0 and T<=0 and 3*T+abs(C)<=-2:
    direction='upper';bound=3*T+abs(C)
   else:raise RuntimeError(('unproved coefficient cone',i,name,A,B0,C))
   cert.append({'gadget_index':i,'test':name,'coefficients':[A,B0,C],
                'bound_direction':direction,'universal_bound':bound})
 (HERE/'coefficient_certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 # Literal physical positive control demonstrating d>=3 is necessary.
 a=[4,1,2];r=0;x=state(dp(a)[r]);check(x==(3,1,3,1),'C3 state')
 trace=[x]
 for op in ['A','D','K']:
  a,r=build(a,r,op,bank);x=tr(x,op,bank);check(state(dp(a)[r])==x,'physical small-deletion reset');trace.append(x)
 check(x==(17,3,21,3),'positive d=1 reset')
 control={'rows':a,'root':r,'state':x,'trace':trace}
 # Root switches followed by short padding on an actual d=3 ready graph.
 actual_checks=0
 for op in [0,1,65,66]:
  z,t=build(a,r,op,bank);y=tr(x,op,bank)
  check(state(dp(z)[t])==y,'actual gadget graph')
  for suffix in ['', 'K','DK']:
   zz=z;tt=t;yy=y
   for step in suffix:zz,tt=build(zz,tt,step,bank);yy=tr(yy,step,bank)
   check(state(dp(zz)[tt])==yy,'actual padded graph')
   hh,ee,RR,dd=yy
   check(not(ee==dd and abs(RR-hh-dd)==1),'unexpected physical ready return');actual_checks+=1
 for k in range(6):
  check(state(dp(a)[r])==x and x[0]==18*2**k-1,'unbounded-family finite control')
  a,r=build(a,r,'K',bank);x=tr(x,'K',bank)
 summary={}
 for name in ('direct','K','DK'):
  group=[z for z in cert if z['test']==name]
  pos=[z['universal_bound'] for z in group if z['bound_direction']=='lower']
  neg=[z['universal_bound'] for z in group if z['bound_direction']=='upper']
  summary[name]={'tested_forms':len(group),'positive_forms':len(pos),'negative_forms':len(neg),
                 'smallest_positive_lower_bound':min(pos) if pos else None,
                 'largest_negative_upper_bound':max(neg) if neg else None}
 print(json.dumps({'status':'PASS','physical_gadgets':len(bank),'universal_identity_controls':len(bank),'coefficient_forms':len(cert),
                   'cone':'All h>=3d, d>=3, sigma in {-1,+1}, with no upper bound',
                   'coefficient_summary':summary,'actual_padded_graph_controls':actual_checks,
                   'small_deletion_positive_control':control,
                   'inventory_sha256':hashlib.sha256((HERE/'inventory67.json').read_bytes()).hexdigest(),
                   'scope':'Controls for an arbitrary-gadget structural theorem: one root-switch followed by arbitrary K/D padding cannot reset a binary-ready input with d>=3. Multiple root switches are not excluded.'},indent=2))


if __name__=='__main__':main()
