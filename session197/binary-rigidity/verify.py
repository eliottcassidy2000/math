#!/usr/bin/env python3
import itertools,json,math

def need(p,s):
 if not p:raise RuntimeError(s)
def h(T):
 n=len(T)
 if not n:return 1
 D=[[0]*n for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for mask in range(1,1<<n):
  for v,c in enumerate(D[mask]):
   if c:
    C=T[v]&~mask
    while C:
     bit=C&-C;C-=bit;D[mask|bit][bit.bit_length()-1]+=c
 return sum(D[-1])
def induced(T,S):return tuple(sum(1<<j for j,v in enumerate(S)if T[u]>>v&1)for u in S)
def response(S,a):
 n=len(S);q=0
 for k in range(n):
  for seq in itertools.permutations([v for v in range(n)if v!=a],k):
   P=seq+(a,)
   if all(S[u]>>v&1 for u,v in zip(P,P[1:])):q+=h(induced(S,[v for v in range(n)if v not in P]))
 return q

def build(S,a,m):
 t=len(S);T=[0]*(m+2*t);O=[m+i for i in range(t)];I=[m+t+i for i in range(t)]
 for r in range(m):
  T[r]|=sum(1<<s for s in range(r+1,m))+sum(1<<v for v in O)
  for v in I:T[v]|=1<<r
 for i in range(t):
  for j in range(t):
   if S[i]>>j&1:T[O[i]]|=1<<O[j];T[I[j]]|=1<<I[i]
   if i==j==a:T[O[i]]|=1<<I[j]
   else:T[I[j]]|=1<<O[i]
 tau=list(range(m+2*t))
 for r in range(m):tau[r]=m-1-r
 for i in range(t):tau[O[i]]=I[i];tau[I[i]]=O[i]
 for u in range(len(T)):
  for v in range(len(T)):
   need((T[u]>>v&1)==(T[tau[v]]>>tau[u]&1),'anti-involution')
 return T

def matching(C):
 used={0}
 for row in C:
  nxt=set(used)
  for mask in used:
   X=row&~mask
   while X:
    bit=X&-X;X-=bit;nxt.add(mask|bit)
  used=nxt
 return max(map(int.bit_count,used))

def main():
 matrices=rankone=0
 for t in range(6):
  entries=[(i,j)for i in range(t)for j in range(i,t)]
  for mask in range(1<<len(entries)):
   C=[0]*t
   for e,(i,j)in enumerate(entries):
    if mask>>e&1:C[i]|=1<<j;C[j]|=1<<i
   nu=matching(C);ones=[(i,j)for i in range(t)for j in range(t)if C[i]>>j&1]
   need((nu==1)==(len(ones)==1 and ones[0][0]==ones[0][1]),'symmetric matching-one classification')
   matrices+=1;rankone+=nu==1
 sources=clones=units=0;amplitudes={};maxorder=0
 for t in range(1,5):
  pairs=list(itertools.combinations(range(t),2))
  for mask in range(1<<len(pairs)):
   S=[0]*t
   for e,(i,j)in enumerate(pairs):
    if mask>>e&1:S[i]|=1<<j
    else:S[j]|=1<<i
   hs=h(S)
   for a in range(t):
    q=response(S,a);sources+=1
    for m in (1,2,3):
     T=build(S,a,m);H=h(T);need(H==hs*hs+q*q*2**m,'literal square formula');clones+=1;maxorder=max(maxorder,len(T))
     if H%7==0:need(hs%7==q%7==0 and H%49==0 and H>=147,'general rank-one no-7p obstruction')
     if hs==1:need(H%7!=0,'unit no-seven obstruction')
    if hs==1:
     units+=1;order=sorted(range(t),key=lambda v:-S[v].bit_count())
     need(all(S[order[i]]>>order[j]&1 for i in range(t)for j in range(i+1,t)),'H1 transitive')
     j=order.index(a);need(q==2**j,'unit amplitude classification')
     amplitudes[q*q]=amplitudes.get(q*q,0)+1
     T=build(S,a,1);O=[1+i for i in order];I=[1+t+i for i in order]
     core=I[:j+1][::-1]+[0]+O[:j+1];pre=I[j+1:][::-1];post=O[j+1:]
     total=pre+core+post
     for u in range(len(total)):
      for v in range(u+1,len(total)):
       back=(total[u]==core[0] and total[v]==core[-1])
       need(bool(T[total[u]]>>total[v]&1)==(not back),'complete core and padding description')
 residues=0
 for hs in range(1,100,2):
  for q in range(1,100):
   for m in range(1,7):
    H=hs*hs+q*q*2**m
    if H%7==0:need(hs%7==q%7==0 and H%49==0 and H>=147,'universal arithmetic rank-one no-7p')
    residues+=1
 print(json.dumps({'rank_one_norm_controls':residues,'status':'PASS','symmetric_matrices_through_size5':matrices,'rank_one_matrices':rankone,'all_labelled_side_sources_through4':sources,'physical_clone_graphs':clones,'maximum_physical_order':maxorder,'unit_offset_sources':units,'unit_amplitude_histogram':amplitudes},sort_keys=True,indent=2))
if __name__=='__main__':main()
