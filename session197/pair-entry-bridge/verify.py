from pathlib import Path
import json,itertools,hashlib
from core import *
P=Path(__file__).resolve().parent

def main():
 bank=json.loads((P/'control_graphs.json').read_text());counts=dict(pair_grafts=0,general_pair_grafts=0,adjacent_profiles=0,nonadjacent_maps=0,tail_moves=0,literal_hamilton_counts=0)
 for n in range(1,7):
  for entry in bank[str(n)]:
   T=entry['rows'];p=inventories(T);require(p['H']==entry['H'],'control graph H')
   if n<=5:
    for a in range(n):
     h,e,R,d=p['H'],p['e'][a],p['R'][a],p['allH'][((1<<n)-1)^(1<<a)]
     U,r=glue(T,{a},[0]);q=inventories(U)
     require((q['H'],q['e'][r],q['R'][r],q['allH'][((1<<len(U))-1)^(1<<r)])==(h-e+R,h-e,h+R,h),'J')
     require((q['H'],q['e'][a],q['R'][a],q['allH'][((1<<len(U))-1)^(1<<a)])==(h-e+R,d,2*R-d,d),'Z')
     V,_=birth(T,[a],0,1);q=inventories(V)
     require((q['H'],q['e'][a],q['R'][a],q['allH'][((1<<len(V))-1)^(1<<a)])==(h,0,2*R,d),'O')
     counts['tail_moves']+=3
   if 2<=n<=5:
    for a,b in itertools.combinations(range(n),2):
     response=pair_response(T,a,b)
     for L in (1,2,3,4,7):
      U,r=glue(T,{a,b},gadget(L));require(root_profile(U,r)==pair_transfer(response,L),'general pair transfer');counts['general_pair_grafts']+=1
   if cycles(T)!=1:continue
   ones=[a for a,x in enumerate(p['s']) if x==1]
   for a,b in itertools.combinations(ones,2):
    response=pair_response(T,a,b)
    for L in range(1,9):
     G=gadget(L);g,Q,B=constants(L);gp=inventories(G)
     require((gp['H'],gp['R'][0],gp['L'][0],gp['s'][0],gp['e'][0])==(g,Q,L,L%2,1),'physical gadget constants')
     U,r=glue(T,{a,b},G);q=inventories(U);want=pair_transfer(response,L)
     require((q['H'],q['e'][r],q['R'][r])==want and q['s'][r]==2,'pair transfer')
     C=cycles(U);require(C in (1,2) and (L==1 or C==2),'pair cycle count')
     if L>=2:require(want[0]>=response[0]+(L-1)*response[2],'pair inverse bound')
     # Literal two-port construction of every gadget incidence.
     V,ports=birth(T,[a,b],3,0);root=n
     for w in range(1,len(G)):
      mask=0 if G[0]>>w&1 else 1<<ports.index(root)
      V,ports=birth(V,ports,mask,2)
     require(V==U and root in ports,'physical pair gadget word')
     if len(U)<=7:require(literal_h(U)==q['H'],'literal Hamilton count');counts['literal_hamilton_counts']+=1
     counts['pair_grafts']+=1
   for a in range(n):
    if T[a].bit_count()!=1:continue
    b=T[a].bit_length()-1
    require(a in ones and b in ones,'adjacent start-one roots')
    response=pair_response(T,a,b);require(response[3]==p['R'][a]-p['e'][a],'adjacent K identity')
    U,ports=birth(T,[a,b],3,1);f=n;q=inventories(U);hs=p['H']-p['e'][b]+p['R'][b]
    require(ports==[a,f] and root_profile(U,f)==(hs,p['H']-p['e'][a]-p['e'][b],p['H']+p['R'][a]+p['R'][b]),'first clone')
    require(root_profile(U,a)==(hs,p['H']-p['e'][b],p['H']+p['R'][b]),'last clone')
    require(q['s'][f]==2 and q['s'][a]==1 and cycles(U)==1,'clone starts and cycles')
    V,new=glue(T,{b},[0]);perm=list(range(n+1));perm[a],perm[new]=new,a
    require(adjacency_map(V,U,perm),'same graph first/last isomorphism')
    for L in (1,2,3,4,7):
     first,_=glue(U,{f},gadget(L));last,_=glue(U,{a},gadget(L))
     require(inventories(first)['H']-inventories(last)['H']==(L%2)*p['e'][a]+L*p['R'][a],'physical next-gadget difference')
    counts['adjacent_profiles']+=1
    # Nonadjacent canonical graph and precise retained-LAST-clone relabeling.
    F,ports=birth(U,ports,2,1);require(ports==[a,n+1],'nonadjacent retained last clone')
    A,u=glue(T,{a},[0]);V,r=glue(A,{u,b},[0]);perm=list(range(n+2));perm[a]=n;perm[u]=n+1;perm[r]=a
    require(adjacency_map(V,F,perm),'nonadjacent graph map')
    pf=inventories(F);pv=inventories(V)
    require((pf['H'],pf['e'][a],pf['R'][a],pf['s'][a])==(pv['H'],pv['e'][r],pv['R'][r],2),'nonadjacent marked root')
    require(pf['H']==p['H']+p['e'][a]+p['R'][a]+2*p['R'][b] and cycles(F)==2,'nonadjacent scalar')
    counts['nonadjacent_maps']+=1
 T=[2,4,1];U,p=birth(T,[0,1],3,1);F,p=birth(U,p,2,1)
 require(F==[18,20,25,3,8] and root_profile(F,0)==(13,2,18) and root_profile(F,3)==(13,3,15),'marked-root hostile')
 U,r=glue([0],{0},gadget(3));require(root_profile(U,r)==(5,1,7),'physical inverse hostile')
 h,e,R=(5,1,7)
 # Inverse M1: old h=R_new-h_new+e_new; old e=R_new-h_new; old R=h_new-e_new.
 old=(R-h+e,R-h,h-e);require(old==(3,2,4),'formal inverse hostile')
 out={'status':'PASS','checks':counts,'nonadjacent_example':{'rows':F,'ports':p,'desired_root':0,'desired_profile':[13,2,18],'other_root':3,'other_profile':[13,3,15]},'formal_inverse_hostile':{'physical':[5,1,7],'nonphysical_inverse':list(old)},'scope':'All-order identities proved by bijections; finite physical controls are not a coverage extrapolation.'}
 return out
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
