from core import inventories,birth,require,tournament
from pathlib import Path
import json
P=Path(__file__).parent

def module(U,r,S):
 n=len(U);T=list(U)+[x<<n for x in S]
 for a in range(n):
  for b in range(len(S)):
   if a==r:T[n+b]|=1<<r
   else:T[a]|=1<<(n+b)
 return T

def sub(T,vs):return [sum(((T[a]>>b)&1)<<j for j,b in enumerate(vs))for a in vs]
def terminal_scc(T,vs):
 reach={a:{a}|{b for b in vs if T[a]>>b&1}for a in vs}
 for z in vs:
  for a in vs:
   if z in reach[a]:reach[a]|=reach[z]
 groups=[];left=set(vs)
 while left:
  a=min(left);C={b for b in left if b in reach[a] and a in reach[b]};groups.append(C);left-=C
 sink=[C for C in groups if not any(T[a]>>b&1 for a in C for b in set(vs)-C)]
 require(len(sink)==1,'unique sink SCC');return sorted(sink[0])
def main():
 banks=json.loads((P/'controls.json').read_text());modchecks=0;statechecks=0;unitchecks=0;found=[]
 for n in range(1,6):
  for z in banks[str(n)]:
   U=z['rows'];p=inventories(U);F=(1<<n)-1
   for r in range(n):
    h,e,R,d=p['H'],p['e'][r],p['R'][r],p['allH'][F^(1<<r)]
    for s in range(1,5):
     for zz in banks[str(s)]:
      S=zz['rows'];q=inventories(S);G=q['H'];B=sum(q['allH'][M]*q['allH'][((1<<s)-1)^M]for M in range(1<<s))
      T=module(U,r,S);v=inventories(T);D=v['allH'][((1<<len(T))-1)^(1<<r)]
      require((v['H'],v['e'][r],v['R'][r],D)==((R-d)*B+(h-e-R+2*d)*G,d*G,(R-d)*B+d*G,d*G),'arbitrary module profile')
      require(v['R'][r]-v['H']-D==(R-h+e-2*d)*G,'offset factor');modchecks+=1
 for n in range(2,8):
  for z in banks[str(n)]:
   T=z['rows'];p=inventories(T);F=(1<<n)-1
   for r in range(n):
    d=p['allH'][F^(1<<r)];e=p['e'][r]
    if e!=d:continue
    S=terminal_scc(T,[v for v in range(n)if v!=r]);require(all(T[v]>>r&1 for v in S),'all terminal vertices beat root')
    G=inventories(sub(T,S))['H'];beta=p['R'][r]-p['H']-d;require(beta%G==0,'sink factor');statechecks+=1
    if abs(beta)==1:
     require(len(S)==1 and T[S[0]]==1<<r,'unit implies outdegree-one vertex');unitchecks+=1
     U=T[:]
     for t in range(1,4):
      U,_=birth(U,[r],1,1);q=inventories(U)
      require(q['H']==(p['H']+beta)*2**t-beta,'returning scalar clock')
     if p['H'] in [1,3,5,9,11,13,23,67]and n<=7:found.append([T,r,beta,p['H'],S[0]])
 # A nontrivial sink component gives a nonunit offset; no stronger deletion claim is made.
 T=module([2,4,1],0,[8,1,3,6]);p=inventories(T);require((p['H'],p['e'][0],p['R'][0])==(69,5,69),'nonunit H69 hostile')
 T=[8,1,3,6];r=2;vals=[5]
 for _ in range(3):
  T,_=birth(T,[r],1,1);vals.append(inventories(T)['H'])
 require(vals==[5,13,25,49],'unit tail does not remove first transient')
 T=[2,12,9,17,7];r=1;p=inventories(T);require((p['H'],p['e'][r],p['R'][r],p['allH'][29])==(13,3,15,3),'H13 opposite sign control')
 out=dict(status='PASS',transient_counts=vals,arbitrary_module_checks=modchecks,physical_e_equals_d_states=statechecks,unit_states=unitchecks,nonunit_hostile=dict(rows=module([2,4,1],0,[8,1,3,6]),root=0,H=69,e=5,R=69,d=5,beta=-5),unit_examples=found)
 print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__':main()
