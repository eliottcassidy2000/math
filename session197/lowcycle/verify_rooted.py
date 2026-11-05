#!/usr/bin/env python3
"""Physical independent tests for ROOTED.md. Standard library, explicit checks."""
from itertools import permutations
from pathlib import Path
import json
from verify import require, all_labelled, count_paths, strong, transitive_order, extend_classes
from unique_start import gadgets, word, transfer, bounded

def induced(rows,vs):
 return [sum(1<<j for j,v in enumerate(vs) if rows[u]>>v&1) for u in vs]

def profile(rows,root):
 n=len(rows);N=1<<n;full=N-1
 dp=[[0]*n for _ in range(N)];rp=[[0]*n for _ in range(N)]
 for u in range(n):dp[1<<u][u]=1
 rp[1<<root][root]=1
 for S in range(1,N):
  for u in range(n):
   if not (dp[S][u] or rp[S][u]):continue
   todo=rows[u]&~S
   while todo:
    bit=todo&-todo;todo-=bit;v=bit.bit_length()-1
    dp[S|bit][v]+=dp[S][u];rp[S|bit][v]+=rp[S][u]
 H=[sum(x) for x in dp];H[0]=1
 starts=[sum(x) for x in rp]
 R=sum(starts[S]*H[full^S] for S in range(N) if S>>root&1)
 C=sum(rp[full][u] for u in range(n) if rows[u]>>root&1) if n>=3 else 0
 rest=full^(1<<root)
 B=sum(H[S]*H[rest^S] for S in range(N) if not S&~rest)
 return (H[full],dp[full][root],R,starts[full],C,H[rest],B)

def converse(rows):
 return [sum(1<<v for v in range(len(rows)) if rows[v]>>u&1) for u in range(len(rows))]

def gadget_graph(L):
 bits=word(L);m=len(bits);r=[0]*(m+1)
 for i,b in enumerate(bits,1):
  if b:r[0]|=1<<i
  else:r[i]|=1
  for j in range(i+1,m+1):r[i]|=1<<j
 return r

def glue(A,a,G,b):
 m=len(A);n=len(G)
 r=[x|(((1<<n)-1)<<m) for x in A]+[x<<m for x in G]
 r[a]^=1<<(m+b);r[m+b]|=1<<a
 return r,m+b

def components(rows):
 n=len(rows);reach=[r|1<<u for u,r in enumerate(rows)]
 for k in range(n):
  for u in range(n):
   if reach[u]>>k&1:reach[u]|=reach[k]
 unseen=set(range(n));out=[]
 while unseen:
  u=min(unseen);c=[v for v in unseen if reach[u]>>v&1 and reach[v]>>u&1]
  unseen.difference_update(c);out.append(sorted(c))
 return sorted(out,key=lambda c:-sum(bool(rows[c[0]]>>d[0]&1) for d in out if c!=d))

def decompose_start_one(rows,r):
 if len(rows)==1:return
 require(strong(rows),'decomposition input not strong')
 vs=[v for v in range(len(rows)) if v!=r];small=induced(rows,vs)
 cs=[[vs[u] for u in c] for c in components(small)]
 U=cs[0];W=sum(cs[1:],[])
 require(all(len(c)==1 for c in cs[1:]),'nontransitive tail')
 aa=[u for u in U if rows[r]>>u&1];require(len(aa)==1,'root has wrong port count')
 a=aa[0];A=induced(rows,U);ar=U.index(a)
 require(profile(A,ar)[3]==1,'recursive root-start mismatch')
 require(transitive_order(rows,W)==W,'tail order')
 G=induced(rows,[r]+W)
 require(profile(G,0)[1]==1,'gadget ending count')
 actual,newroot=glue(A,ar,G,0)
 require(actual==induced(rows,U+[r]+W),'physical reconstruction')
 decompose_start_one(A,ar)

def main():
 bank=gadgets(255);gadget_checks=0
 for L in range(1,256):
  G=gadget_graph(L);p=profile(G,0);q=profile(converse(G),0)
  g,Q,B,delta=bank[L]
  require((p[0],p[1],p[2],p[3],p[5],p[6],q[2])==(g,1,Q,L%2,1,B,L),f'compiled gadget {L}')
  gadget_checks+=1
 inputs=[(r,a) for n in range(1,4) for r in all_labelled(n) for a in range(n)]
 glue_checks=0
 for A,a in inputs:
  h,e,R,s,C,d,B=profile(A,a)
  for G,b in inputs:
   g,E,Q,t,CG,dG,BG=profile(G,b);L=profile(converse(G),b)[2]
   T,r=glue(A,a,G,b);p=profile(T,r);delta=int(len(G)==1)
   require(p[:3]==(g*h-t*e+L*R,E*h-delta*e,Q*h+BG*R),'inherited response')
   require(p[3]==dG*s,'root-start product')
   if len(T)>=3:require(p[4]==s*E,'cycle product')
   glue_checks+=1
 # Add gadgets with nontransitive root-deletion to test a genuinely nonunit dG.
 for G in all_labelled(4):
  for b in range(4):
   for A,a in [([2,4,1],0),([8,1,3,6],1),([0],0)]:
    h,e,R,s,C,d,B=profile(A,a);g,E,Q,t,CG,dG,BG=profile(G,b);L=profile(converse(G),b)[2]
    T,r=glue(A,a,G,b);p=profile(T,r)
    require(p[:3]==(g*h-t*e+L*R,E*h,Q*h+BG*R),'order-four response')
    require(p[3]==dG*s and p[4]==s*E,'nonunit root-start/cycle response')
    glue_checks+=1
 ts={(0,)};deletion_checks=0
 for n in range(1,7):
  if n>1:ts=extend_classes(ts)
  for rows in ts:
   if n>1 and not strong(rows):continue
   for r in range(n):
    if profile(rows,r)[3]==1:decompose_start_one(rows,r);deletion_checks+=1
 family=[];rows=[0];root=0;marks=[];blocks=[];x=(1,1,1)
 for t in range(1,5):
  offset=len(rows);G=gadget_graph(3);rows,root=glue(rows,root,G,0);marks.append(root);blocks.append(list(range(offset,offset+3)))
  x=transfer(x,3,bank[3]);p=profile(rows,root)
  require(p[:3]==x and p[3:5]==(1,1),'triangle family')
  require(transitive_order(rows,[v for v in range(len(rows)) if v not in marks]) is not None,'feedback upper certificate')
  require(all(strong(induced(rows,b)) for b in blocks),'feedback lower certificate')
  family.append({'t':t,'order':len(rows),'rows':rows[:],'root':root,'H':p[0],'C':1,'feedback_vertex_number':t,'feedback_set':marks[:]})
 out501=bounded(501);out5001=bounded(5001)
 require((out501['states'],out501['edges'])==(2727,3900),'image501')
 require((out5001['states'],out5001['edges'])==(81585,117980),'image5001')
 require(67 in out5001['missing_odd'],'67 root-start exclusion')
 receipt={'compiled_gadget_graph_checks':gadget_checks,'physical_gluing_checks':glue_checks,'rooted_recursive_decomposition_checks':deletion_checks,'family':family,'image501':{k:v for k,v in out501.items() if k!='values'},'image5001':{k:v for k,v in out5001.items() if k!='values'}}
 text=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(Path(__file__).resolve().parent/'rooted_receipt.json').write_text(text);print(text,end='')

if __name__=='__main__':main()
