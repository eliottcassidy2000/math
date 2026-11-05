#!/usr/bin/env python3
import itertools,json,math
from fractions import Fraction

def need(p,s):
 if not p:raise RuntimeError(s)
def permcode(T,p):return tuple(sum(1<<j for j in range(len(T))if T[p[i]]>>p[j]&1)for i in range(len(T)))
def canonical(T,r=None):
 P=itertools.permutations(range(len(T)))if r is None else ((r,)+q for q in itertools.permutations([v for v in range(len(T))if v!=r]))
 return min(permcode(T,p)for p in P)
def autom(T,p,anti=False):return all(bool(T[u]>>v&1)==bool(T[p[v] if anti else p[u]]>>(p[u] if anti else p[v])&1)for u in range(len(T))for v in range(len(T)))
def h(T):return sum(all(T[u]>>v&1 for u,v in zip(p,p[1:]))for p in itertools.permutations(range(len(T))))
def build(S,a):
 t=len(S);O=list(range(1,t+1));I=list(range(t+1,2*t+1));T=[0]*(2*t+1)
 T[0]=sum(1<<v for v in O)
 for v in I:T[v]|=1
 for i in range(t):
  for j in range(t):
   if S[i]>>j&1:T[O[i]]|=1<<O[j];T[I[j]]|=1<<I[i]
   if i==j==a:T[O[i]]|=1<<I[j]
   else:T[I[j]]|=1<<O[i]
 return T,O,I

def closure(T,u):
 s=1<<u
 while True:
  q=s
  for v in range(len(T)):
   if s>>v&1:q|=T[v]
  if q==s:return s
  s=q

def main():
 A=[1];R=[0];C=[0];L=[0];CL=[0];data=[];allhalf=allaut=allanti=0
 for t in range(1,5):
  pairs=list(itertools.combinations(range(t),2));plain=set();rooted=set()
  for bits in range(1<<len(pairs)):
   S=[0]*t
   for e,(u,v)in enumerate(pairs):
    if bits>>e&1:S[u]|=1<<v
    else:S[v]|=1<<u
   plain.add(canonical(S))
   for a in range(t):rooted.add(canonical(S,a))
  A.append(len(plain));R.append(len(rooted));C.append(0);L.append(0);CL.append(0);outputcodes=set()
  for S in sorted(rooted):
   a=0;hs=h(S);G=[p for p in itertools.permutations(range(t))if p[a]==a and autom(S,p)];g=len(G)
   T,O,I=build(S,a);n=len(T);aut=0;codes=[]
   for po in itertools.permutations(range(t)):
    for pi in itertools.permutations(range(t)):
     p=[0]+[O[j]for j in po]+[I[j]for j in pi]
     aut+=autom(T,p);codes.append(permcode(T,p));allaut+=1
   need(aut==g*g,'literal root stabilizer product');outputcodes.add(min(codes))
   anti=0;cent=0
   for po in itertools.permutations(range(t)):
    inv=[0]*t
    for j,k in enumerate(po):inv[k]=j
    tau=[0]+[I[j]for j in po]+[O[j]for j in inv]
    anti+=autom(T,tau,True);allanti+=1
   need(anti==g,'literal anti-involution count')
   tau=[0]+I+O
   for po in itertools.permutations(range(t)):
    for pi in itertools.permutations(range(t)):
     p=[0]+[O[j]for j in po]+[I[j]for j in pi]
     if autom(T,p)and all(p[tau[v]]==tau[p[v]]for v in range(n)):cent+=1
   need(cent==g,'literal centralizer')
   fixed=0
   for bits in range(1<<t):
    choices=[I[i]if bits>>i&1 else O[i]for i in range(t)]
    for half in itertools.permutations(choices):
     word=half+(0,)+tuple(tau[v]for v in reversed(half))
     fixed+=all(T[u]>>v&1 for u,v in zip(word,word[1:]));allhalf+=1
   need(fixed==hs,'literal fixed Hamilton words')
   strong=all(closure(T,v)==(1<<n)-1 for v in range(n));terminal=all(closure(S,v)>>a&1 for v in range(t))
   need(strong==terminal,'strong iff marked terminal component')
   need(hs%g==0,'integral fixed mass');mass=hs//g;L[t]+=mass
   if strong:C[t]+=1;CL[t]+=mass
   data.append({'side':S,'h':hs,'g':g,'fixed':fixed,'root_stabilizer':aut,'anti_involutions_fixing_root':anti,'centralizer':cent,'strong':strong,'mass':mass})
  need(len(outputcodes)==len(rooted),'rooted physical isomorphism injection')
  need(L[t]==t*2**((t-1)*(t-2)//2),'all rank-one mass')
  need(R[t]==sum(A[i]*C[t-i]for i in range(t)),'ordinary A000568 convolution')
  M=[1]+[2**((j-1)*(j-2)//2)for j in range(1,t+1)]
  need(t*M[t]==sum(M[i]*CL[t-i]for i in range(t)),'strong weighted convolution')
 T,_,_=build((0,),0);fullaut=sum(autom(T,p)for p in itertools.permutations(range(3)));need(fullaut==3,'C3 full group differs from root stabilizer')
 need(A==[1,1,1,2,4]and R==[0,1,2,4,12]and C==[0,1,1,2,7]and CL==[0,1,1,4,25],'literal initial inventories')
 print(json.dumps({'status':'PASS','A000568_through4':A,'rooted_classes':R,'strong_rank_one_classes':C,'full_rank_one_mass':L,'strong_rank_one_mass':CL,'distinct_source_classes':len(data),'root_stabilizer_permutations_checked':allaut,'anti_permutations_checked':allanti,'symmetric_half_words_checked':allhalf,'class_certificates':data},sort_keys=True,indent=2))
if __name__=='__main__':main()
