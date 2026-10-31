#!/usr/bin/env python3
"""Inert pinned Gerbicz tables: edge-edit statistics, no external code execution."""
import argparse,hashlib,itertools,json,math,re
from pathlib import Path
from functools import lru_cache
from fractions import Fraction
PIN='602c1919ee51101519a2f19f2c1817407006824f6f09c85a3b772a818a0b736e'
def need(ok,msg):
 if not ok:raise ValueError(msg)
def parse(path):
 raw=Path(path).read_bytes();need(hashlib.sha256(raw).hexdigest()==PIN,'source digest')
 t=raw.decode('ascii');out={}
 for name in ('S','data'):
  q=re.search(r'int\s+'+name+r'\[[^;=]+?=\s*(\{.*?\});',t,re.S)
  need(q is not None,'array');s=q.group(1);need(re.fullmatch(r'[\s\d{},\-]+',s)!=None,'integer grammar')
  out[name]=json.loads(s.replace('{','[').replace('}',']'))
 return out

def edges(p):return {tuple(sorted(e)) for e in zip(p,p[1:])}
def dist(a,b):return len(edges(a)-edges(b))
def cycle_edges(a):return edges(a)|{tuple(sorted((a[0],a[-1])))}
def cycle_dist(a,b):return len(cycle_edges(a)-cycle_edges(b))
def phase(p):return {x:j%2 for j,x in enumerate(p)}
def validate_pair(n,a,b):
 for k,p in enumerate((a,b)):
  need(len(p)==n+k and set(p)==set(range(1,n+k+1)), 'expanded coverage')
  need(p[0]==1 and p[-1]==(3 if (n+k)%2 else 8), 'expanded ports')
  need(all(math.isqrt(u+v)**2==u+v for u,v in zip(p,p[1:]+p[:1])), 'expanded square cycle')
 pa=phase(a);pb=phase(b);need(all(pa[v]==pb[v] for v in a),'expanded phase')
def atom_ends(row):
 out=[];_,d,p,ex=row[0]
 for c,o,e in row[1:]:
  if o==0:out.append((c,c))
  else:
   ell=(p+e)%2;a=49+c;b=49*(3 if ell else 8)+(-1)**(ell+1)*c
   out.append((a,b) if o==1 else (b,a))
 return out

def connectors(row):
 a=atom_ends(row);return {tuple(sorted((x[1],y[0]))) for x,y in zip(a,a[1:])}

def run(path):
 D=parse(path);words=D['data'];S=D['S']
 @lru_cache(maxsize=16)
 def pair(n):
  if 41<=n<=2032:return tuple(tuple(p) for p in S[2*(n-41):2*(n-41)+2])
  need(n>=2033,'order');m,d=divmod(n-24,49);d+=24;a,b=pair(m);out=[]
  for ex in (0,1):
   w=[]
   for c,o,e in words[4*(d-24)+2*(m%2)+ex][1:]:
    if o==0:w.append(c)
    else:
     block=[49*x+(c if j%2==0 else -c) for j,x in enumerate((a,b)[e])]
     w.extend(block if o==1 else reversed(block))
   out.append(tuple(w))
  return tuple(out)
 base=[]
 for n in range(41,2033):
  a,b=pair(n);base.append((n,dist(a,b),cycle_dist(a,b)))
 changes=[];K=[];J=[];Hword=[]
 for p in (0,1):
  for d in range(24,73):
   r0=words[4*(d-24)+2*p];r1=words[4*(d-24)+2*p+1]
   a={c:e for c,o,e in r0[1:] if o};b={c:e for c,o,e in r1[1:] if o}
   cs=[c for c in a if a[c]!=b[c]];need(len(cs)==1,'one changed block')
   need(cs==[(-1)**p*(d+1-49)],'threshold offset')
   changes.append((d,p,cs[0]));K.append((d,p,len(connectors(r0)-connectors(r1))))
   if d<72:
    nxt=words[4*(d+1-24)+2*p]
    J.append((d,p,len(connectors(r0)-connectors(nxt))))
    Hword.append((d,p,len(connectors(r1)-connectors(nxt))))
 carry_K={p:len(connectors(words[4*48+2*(1-p)+1])-connectors(words[2*p])) for p in (0,1)}
 Hbase=[(r,dist(pair(r-1)[1],pair(r)[0])) for r in range(42,2033)]
 samples=[]
 for n in [2033,2034,2080,2081,2082,100000]:
  a,b=pair(n);m=(n-24)//49;d=n-49*m;k=next(z for dd,pp,z in K if (dd,pp)==(d,m%2))
  validate_pair(n,a,b);dn=dist(a,b);need(dn==dist(*pair(m))+k,'exact internal+connector recurrence')
  samples.append(dict(n=n,parent=m,d=d,path_deletions=dn,cycle_deletions=cycle_dist(a,b),connector_deletions=k))
 canonical=[]
 for n in [2033,2034,2080,2081,2082,98170,102041]:
  a,b=pair(n);c=pair(n+1)[0];h=dist(b,c);dn=dist(a,b);dc=cycle_dist(a,b);cn=dist(a,c);cc=cycle_dist(a,c)
  need(h-dn-1<=cn<=h+dn,'canonical path triangle bounds')
  need(h-dc-1<=cc<=h+dc,'canonical cycle triangle bounds')
  canonical.append(dict(n=n,path_deletions=cn,cycle_deletions=cc,representation_gap=h,within_pair_path_deletions=dn))
 hs=[]
 for r in [42,43,49*42+24,49*2003+24,49*(49*42+24)+24]:
  x,y=pair(r-1)[1],pair(r)[0];px=phase(x);py=phase(y)
  if r>=2082:
   parent=(r-24)//49;need(dist(x,y)==49*dist(pair(parent-1)[1],pair(parent)[0])+carry_K[parent%2],'exact carry recurrence')
  hs.append(dict(r=r,path_deletions=dist(x,y),phase_disagreements=sum(px[v]!=py[v] for v in range(1,r+1))))
 roots=Hbase+[(2033,dist(pair(2032)[1],pair(2033)[0]))]
 ranked=sorted(((Fraction(48*h+carry_K[r%2],24*(2*r+1)),r,h) for r,h in roots),reverse=True)
 maxratio,seed,seed_h=ranked[0];outside=Fraction(48*72+max(carry_K.values()),24*(2*2034+1))
 need(maxratio>outside,'finite roots dominate all larger noncarry roots')
 global_limsup=dict(ratio=str(maxratio),maximizing_seed=seed,seed_gap=seed_h,mixed_boundary_gap=roots[-1][1],outside_root_upper_bound=str(outside),top_ten=[dict(ratio=str(q),seed=r,gap=h) for q,r,h in ranked[:10]])
 selector=[]
 for carry in [2082,102042]:
  orders=list(range(carry-2,carry+3));banks=[[cycle_edges(pair(n)[0]),cycle_edges(pair(n-1)[1])] for n in orders]
  costs=[[[len(banks[i][j]-banks[i+1][k]) for k in (0,1)] for j in (0,1)] for i in range(4)]
  dp=[(0,[]),(0,[])]
  for matrix in costs:
   dp=[min((max(dp[j][0],matrix[j][k]),dp[j][1]+[j]) for j in (0,1)) for k in (0,1)]
  value,choice,k=min((v,ch,k) for k,(v,ch) in enumerate(dp));choice=choice+[k]
  brute=min(max(costs[i][seq[i]][seq[i+1]] for i in range(4)) for seq in itertools.product((0,1),repeat=5))
  need(value==brute,'selector dynamic programming vs all32 choices')
  h=dist(pair(carry-1)[1],pair(carry)[0]);bp=cycle_dist(*pair(carry-1));bn=cycle_dist(*pair(carry));lower=(h-bp-bn-146)/2
  need(value>=lower,'selector theorem lower bound')
  selector.append(dict(carry=carry,orders=orders,bank_order=['short_n','long_n_minus_1'],transition_deletions=costs,minimax_deletions=value,one_optimal_choice=choice,theorem_lower_bound=lower))
 out=dict(status='PASS',source_sha256=PIN,external_code_executed=False,within_pair_base_range=[41,2032],base_path_max=max(base,key=lambda x:x[1]),base_cycle_max=max(base,key=lambda x:x[2]),base_path_min=min(base,key=lambda x:x[1]),K_max=max(K,key=lambda x:x[2]),K_min=min(K,key=lambda x:x[2]),noncarry_canonical_connector_max=max(J,key=lambda x:x[2]),noncarry_representative_change_max=max(Hword,key=lambda x:x[2]),H_base_min=min(Hbase,key=lambda x:x[1]),H_base_max=max(Hbase,key=lambda x:x[1]),H_base_first=Hbase[:12],transition_K=K,carry_connector_deletions=carry_K,expanded_pair_samples=samples,carry_samples=hs,canonical_samples=canonical,selector_window_controls=selector,global_pinned_bank_limsup=global_limsup)
 return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--output');a=ap.parse_args();out=run(a.source);text=json.dumps(out,indent=2)+'\n'
 if a.output:Path(a.output).write_text(text)
 print(text)
