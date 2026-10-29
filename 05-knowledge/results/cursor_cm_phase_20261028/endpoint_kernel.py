"""Constructive endpoint-response kernel for supplied-run quotient-walk partners."""
from pathlib import Path
from itertools import permutations
import json
import run_probe as r
c=r.c
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(msg)
def reverse_forest(t,F):
 keep={frozenset((x,y)) for R in F for x,y in zip(R,R[1:])}
 return tuple(sum((c.c.edge(t,i,j) if frozenset((i,j)) in keep else 1-c.c.edge(t,i,j))<<j
                  for j in range(len(t)) if i!=j) for i in range(len(t)))
def feedback(n):
 t=list(r.trans(n))
 t[0]^=1<<(n-1);t[-1]|=1
 return tuple(t)
def rail(a,F):
 n=len(a);k=len(F)
 if k<2:raise ValueError('k>=2')
 L=n+3;NB=1+k*L;NC=k+1;children=(a,feedback(NB),feedback(NC))
 owners=[(i,j) for i,t in enumerate(children) for j in range(len(t))]
 src=tuple(sum((c.c.edge(children[i],x,y) if i==j else ((j-i)%3==1))<<v
               for v,(j,y) in enumerate(owners) if (i,x)!=(j,y)) for i,x in owners)
 B=[[n]]+[list(range(n+1+j*L,n+1+(j+1)*L)) for j in range(k)]
 C=list(range(n+NB,n+NB+NC))
 path=B[0]+[C[0]]
 for j,R in enumerate(F):path+=list(R)+B[j+1]+[C[j+1]]
 need(len(path)==len(src)==len(set(path)),'spanning path')
 need(all(c.c.edge(src,x,y) for x,y in zip(path,path[1:])),'legal quotient walk')
 edges={frozenset((x,y)) for x,y in zip(path,path[1:])}
 out=tuple(sum((c.c.edge(src,i,j) if frozenset((i,j)) in edges else 1-c.c.edge(src,i,j))<<j
               for j in range(len(src)) if i!=j) for i in range(len(src)))
 return out
def decode(t,n,k):
 L=n+3;NB=1+k*L;NC=k+1;degree=[row.bit_count() for row in t]
 tags=[]
 def unique(values,why):
  need(len(values)==1,why);return values[0]
 for j in range(k):
  start=1+j*L;end=(j+1)*L
  anchor=NC if j==0 else start+1
  v=unique([i for i,d in enumerate(degree) if d==n+anchor],'unique degree anchor')
  run={anchor:v}
  for pos in range(anchor,start,-1):
   run[pos-1]=unique([i for i,d in enumerate(degree) if d==n+pos-1 and c.c.edge(t,i,run[pos])],'backward rail')
  for pos in range(anchor,end):
   expected=n+min(pos+1,NB-2)
   run[pos+1]=unique([i for i,d in enumerate(degree) if d==expected and c.c.edge(t,run[pos],i)],'forward rail')
  tags.append([run[pos] for pos in range(start,end+1)])
 known={v for tag in tags for v in tag}
 Cs=[i for i,d in enumerate(degree) if i not in known and d>=NB]
 need(len(Cs)==NC,'auxiliary C degree band')
 low=[i for i in range(len(t)) if i not in known and i not in Cs]
 b0=unique([i for i in low if sum(c.c.edge(t,i,j) for j in Cs)==1],'initial B singleton')
 core=[i for i in low if i!=b0]
 need(len(core)==n,'core size')
 c0=unique([j for j in Cs if c.c.edge(t,b0,j)],'initial C port')
 cmarks=[c0]+[unique([v for v in Cs if c.c.edge(t,tag[-1],v)],'C port after rail') for tag in tags]
 need(len(set(cmarks))==NC,'all C ports')
 starts=[unique([v for v in core if c.c.edge(t,cmarks[j],v)],'run start') for j in range(k)]
 ends=[unique([v for v in core if c.c.edge(t,v,tags[j][0])],'run end') for j in range(k)]
 inv={v:i for i,v in enumerate(core)}
 return c.c.induced(t,core),tuple((inv[s],inv[e]) for s,e in zip(starts,ends)),core
def check(a,F,shift=0):
 t=rail(a,F);N=len(t);order=tuple(range(shift%N,N))+tuple(range(shift%N))
 b,p,verts=decode(c.c.relabel(t,order),len(a),len(F))
 inv={order[v]:i for i,v in enumerate(verts)}
 expected=reverse_forest(a,F)
 need(all(c.c.edge(expected,i,j)==c.c.edge(b,inv[i],inv[j]) for i in range(len(a)) for j in range(len(a))),'core relation decoded')
 need(p==tuple((inv[R[0]],inv[R[-1]]) for R in F),'ordered endpoint pairs decoded')
 return N
def endpoint_group(b,F,lengths=False):
 pairs={(R[0],R[-1]):i for i,R in enumerate(F)}
 H=[];images=set()
 for g in c.o.isomorphisms(b,b,True):
  moved=[(g[R[0]],g[R[-1]]) for R in F]
  if all(p in pairs for p in moved) and (not lengths or all(len(F[i])==len(F[pairs[p]]) for i,p in enumerate(moved))):
   H.append(g);images.add(tuple(pairs[p] for p in moved))
 return H,images
def coupled(m):
 n=2*m+2;B=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if i==0 or j==n-1:bit=1
   elif (i-1)//2==(j-1)//2:bit=1
   else:bit=1<=(((j-1)//2-(i-1)//2)%m)<=(m-1)//2
   B[i if bit else j]|=1<<(j if bit else i)
 B=tuple(B);F=(tuple([0]+list(range(1,n-1,2))+[n-1]),)+tuple((v,) for v in range(2,n-1,2))
 return reverse_forest(B,F),B,F
def short_coupled(m):
 n=2*m+1;B=[0]*n;v=2*m
 for i in range(n):
  for j in range(i+1,n):
   if j==v:bit=i%2==0
   elif i//2==j//2:bit=1
   else:bit=1<=((j//2-i//2)%m)<=(m-1)//2
   B[i if bit else j]|=1<<(j if bit else i)
 B=tuple(B);F=((0,v,1),)+tuple((2*j,2*j+1) for j in range(1,m))
 return reverse_forest(B,F),B,F

def equal_coupled(m):
 # Three blocks X,V,Y in a directed cycle. X,Y are cyclic; V is rigid.
 children=[]
 for block in range(3):
  children.append(tuple(sum((int(i<j) if block==1 else int(1<=((j-i)%m)<=(m-1)//2))<<j
                            for j in range(m) if i!=j) for i in range(m)))
 B=tuple(sum((c.c.edge(children[i//m],i%m,j%m) if i//m==j//m else int((j//m-i//m)%3==1))<<j
             for j in range(3*m) if i!=j) for i in range(3*m))
 F=tuple((i,m+i,2*m+i) for i in range(m))
 return reverse_forest(B,F),B,F
def main():
 rows=[];mergers=[]
 for n in range(2,6):
  count=0;orders_count=0
  for a in r.representatives(n):
   aut=c.o.isomorphisms(a,a,True);seen=set()
   for F in r.covers(a):
    key=min(tuple(sorted(r.act_cover(F,g))) for g in aut)
    if key in seen or len(F)<2:continue
    seen.add(key);B=reverse_forest(a,F)
    G=[g for g in aut if tuple(sorted(r.act_cover(F,g)))==F]
    H,images=endpoint_group(B,F)
    need(len(images)>=len(G),'faithful forest action injects into endpoint action')
    Ps=sorted({min(r.act_cover(P,g) for g in G) for P in permutations(F)})
    for i,P in enumerate(Ps):check(a,P,i);orders_count+=1
    if len(images)>len(G):mergers.append({'A':a,'F':F,'joint_group':len(G),'endpoint_image':len(images)})
    count+=1
  rows.append({'n':n,'forests':count,'order_orbits_decoded':orders_count})
  print(json.dumps(rows[-1]),flush=True)
 family=[]
 for kind,m in [(s,m) for s in ('class_only','length_aware','equal_length') for m in ((3,) if s=='equal_length' else (3,5,7))]:
  a,b,F={'class_only':short_coupled,'length_aware':coupled,'equal_length':equal_coupled}[kind](m)
  autA=c.o.isomorphisms(a,a,True);autB=c.o.isomorphisms(b,b,True)
  H,images=endpoint_group(b,F,lengths=kind!='class_only')
  G=[g for g in autB if tuple(sorted(r.act_cover(F,g)))==tuple(sorted(F))]
  need(len(G)==1 and len(images)==m,'coupled odd rotation family')
  need(len(autB)==(m*m if kind=='equal_length' else m),'target automorphism group')
  if kind=='class_only':need(all(row.bit_count()==m for row in a),'regular strong source family')
  Q=(F[:1]+F[2:]+F[1:2]) if kind=='length_aware' else F[1:]+F[:1]
  if kind!='class_only':need(tuple(map(len,F))==tuple(map(len,Q)),'identical ordered lengths')
  need(not any(r.act_cover(F,g)==Q for g in G),'different faithful component orders')
  p=next(g for g in H if [(g[R[0]],g[R[-1]]) for R in F]==[(R[0],R[-1]) for R in Q])
  left,right=rail(a,F),rail(a,Q);extension=tuple(p)+tuple(range(len(a),len(left)))
  need(c.c.relabel(right,extension)==left,'literal universal endpoint transporter realizes rail collision')
  check(a,F);check(a,Q,3)
  family.append({'kind':kind,'m':m,'A':a,'B':b,'F':F,'second_order':Q,'AutA':len(autA),'AutB':len(autB),
                 'joint_group':len(G),'endpoint_image':len(images),'transporter':p,'context_order':len(left)})
 hostile=c.o.hostile()
 for a in (tuple(hostile['rows_AAA']),tuple(hostile['rows_BAA'])):
  F=tuple((i,) for i in range(len(a)));B=reverse_forest(a,F)
  H,im=endpoint_group(B,F)
  need(len(im)==len(c.o.isomorphisms(a,a,True)),'81/27 keeps full source action on singleton ports')
  check(a,F,7)
 out={'status':'FINITE-EXACT','small_universe':rows,'small_endpoint_mergers':mergers,'coupled_family':family,
      'hostile_81_27':hostile,'checks':CHECKS}
 (HERE/'endpoint_kernel.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':'PASS','checks':CHECKS,'small_mergers':len(mergers)}))
if __name__=='__main__':main()
