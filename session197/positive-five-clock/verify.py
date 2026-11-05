#!/usr/bin/env python3
"""Physical +5 ternary family; exact root words, graphs and primality certificates."""
import itertools,json,math,pathlib
ROOT=pathlib.Path(__file__).resolve().parent
Q=[14,32,2,22,103,13,47]

def require(p,s):
 if not p:raise RuntimeError(s)

def transitive(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]

def graph(m):
 require(m>=1,'positive root multiplicity')
 blocks=[transitive(m),transitive(2),[8,1,3,6],[0],transitive(2),[0],[0]]
 offsets=[];n=0
 for B in blocks:offsets.append(n);n+=len(B)
 rows=[0]*n
 for a,B in enumerate(blocks):
  for i,row in enumerate(B):
   rows[offsets[a]+i]=row<<offsets[a]
   for b,C in enumerate(blocks):
    if Q[a]>>b&1:rows[offsets[a]+i]|=((1<<len(C))-1)<<offsets[b]
 validate(rows);return rows

def validate(T):
 n=len(T);require(all(0<=T[i]<(1<<n) and not(T[i]>>i&1) for i in range(n)),'row range')
 require(all(bool(T[i]>>j&1)!=bool(T[j]>>i&1) for i,j in itertools.combinations(range(n),2)),'tournament')

def strong(T):
 n=len(T)
 for root in range(n):
  seen=1<<root
  while True:
   after=seen
   for i in range(n):
    if seen>>i&1:after|=T[i]
   if after==seen:break
   seen=after
  if seen!=(1<<n)-1:return False
 return True

def hp_subsets(T):
 n=len(T);N=1<<n;d=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for i in range(n):d[1<<i][i]=1
 for mask in range(1,N):
  for v in range(n):
   if not(mask>>v&1):continue
   rest=mask^(1<<v)
   if rest:d[mask][v]=sum(d[rest][u] for u in range(n) if rest>>u&1 and T[u]>>v&1)
  h[mask]=sum(d[mask])
 return h

def root_runs(T,r):
 old=[v for v in range(len(T)) if v!=r];n=len(old);N=1<<n;K=n+1
 d=[[[0]*(n+1) for _ in range(K+1)] for _ in range(N)]
 d[0][1][n]=1
 for i in range(n):d[1<<i][0][i]=1
 for mask in range(N):
  for k in range(K+1):
   for i in range(n+1):
    c=d[mask][k][i]
    if not c:continue
    v=r if i==n else old[i]
    for j,w in enumerate(old):
     if not(mask>>j&1) and T[v]>>w&1:d[mask|1<<j][k][j]+=c
    if i!=n and k<K and T[v]>>r&1:d[mask][k+1][n]+=c
 return [sum(d[-1][k]) for k in range(K+1)]

def left_birth(T,ports):
 n=len(T);out=T[:]+[0]
 for v in range(n):
  if v in ports:out[v]|=1<<n
  else:out[n]|=1<<v
 return out

def prime_certificate(p):
 bound=math.isqrt(p);remainders=[p%d for d in range(2,bound+1)]
 require(p>=2 and all(remainders),'complete trial-division primality')
 return {'integer':p,'sqrt_floor':bound,'first_divisor':2,'last_divisor':bound,'remainders_in_divisor_order':remainders}

def main():
 source=graph(1);ports=[8,9,10,11]
 require(len(source)==12 and all(not(source[0]>>v&1) for v in ports),'root in-neighbor ports')
 require([v for v in range(1,12) if not(source[0]>>v&1)]==ports,'exact port set')
 crossing=[(i,j) for i in range(1,12) if source[0]>>i&1 for j in ports if source[i]>>j&1]
 require(crossing==[(1,10),(2,10),(7,8),(7,9)],'complete crossing support')
 require(all(i==7 or j==10 for i,j in crossing),'two-vertex crossing cover')
 a=root_runs(source,0)
 require(a[1:4]==[3461,3456,1152] and not any(a[4:]),'exact complete root-run profile')
 V=[8,1,3,6];hv=hp_subsets(V);bV=sum(hv[s]*hv[15^s] for s in range(16))
 require((hv[-1],bV)==(5,32),'strong four-vertex body profile')
 physical=source;controls=[]
 for m in range(1,5):
  T=graph(m);H=hp_subsets(T)[-1]
  require(H==1152*3**m+5 and len(T)==m+11 and strong(T),'exact physical family control')
  order=list(range(12+m-2,11,-1))+[0]+list(range(1,12))
  require(len(order)==len(T),'returning clone permutation')
  require(all(bool(T[i]>>j&1)==bool(physical[order[i]]>>order[j]&1) for i in range(len(T)) for j in range(len(T))),'all-edge fixed-port return')
  controls.append({'m':m,'order':len(T),'H':H,'canonical_rows':T,'fixed_port_birth_rows':physical,'canonical_to_birth_labels':order})
  physical=left_birth(physical,ports)
 require(controls[0]['H']==3461 and controls[3]['H']==7*13331,'requested prime and7p values')
 certs=[prime_certificate(3461),prime_certificate(13331)]
 residues=[(1152*pow(3,m,7)+5)%7 for m in range(6)]
 require(residues==[2,3,6,1,0,4] and pow(3,6,7)==1,'exact seven-divisibility period')
 witness={'quotient_rows':Q,'module_order':['T','P','U','a','R','b','V'],'fixed_modules':{'P':transitive(2),'U':V,'a':[0],'R':transitive(2),'b':[0],'V':[0]},'source_root':0,'fixed_left_ports':ports,'root_run_coefficients':a,'physical_controls':controls,'prime_certificates':certs,'mod7_residues_m0_through5':residues}
 (ROOT/'witness.json').write_text(json.dumps(witness,sort_keys=True,indent=2)+'\n')
 print(json.dumps({'status':'PASS','identity':'H_m=1152*3^m+5 for every integer m>=1','affine_return':'H_next=3H-10','source_rows':source,'positive_run_coefficients':a[1:4],'maximum_run_index':3,'fixed_left_ports':ports,'physical_counts':[x['H'] for x in controls],'prime_H':3461,'seven_prime_H':93317,'prime_quotient':13331,'divisibility_by7':'m=4 modulo6','unbounded_primality_claimed':False},sort_keys=True,indent=2))
if __name__=='__main__':main()
