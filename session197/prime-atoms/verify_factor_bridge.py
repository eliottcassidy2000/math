#!/usr/bin/env python3
from verify import *
from collections import Counter


def converse(rows):
 return [sum(1<<j for j in range(len(rows)) if rows[j]>>i&1) for i in range(len(rows))]


def value2(x):
 require(x>0,'positive valuation input')
 return (x & -x).bit_length()-1


def kbound(f):return (1+math.isqrt(1+8*f))//2

def ell(f):return max(1,f).bit_length()-1

def source_bound(U,beta):return 1+(2*U+2*beta-1)//3+kbound(U-1)+ell(U-1)+value2(U)

def partner_bound(L):return 1 if L==1 else L+kbound(L-1)+2*ell(L-1)


def star(X,r,Y,s,m):
 a=[i for i in range(len(X)) if i!=r];b=[i for i in range(len(Y)) if i!=s]
 A=list(range(1,1+len(a)));B=list(range(1+len(a),1+len(a)+len(b)))
 Z=list(range(1+len(a)+len(b),1+len(a)+len(b)+m));n=1+len(a)+len(b)+m
 rows=[0]*n;mx={r:0,**dict(zip(a,A))};my={s:0,**dict(zip(b,B))}
 for old,mp in [(X,mx),(Y,my)]:
  for i in range(len(old)):
   for j in range(len(old)):
    if old[i]>>j&1:rows[mp[i]]|=1<<mp[j]
 for x in A:
  for y in B+Z:rows[x]|=1<<y
 for x in Z:
  rows[x]|=1
  for y in B:rows[x]|=1<<y
 for i,x in enumerate(Z):
  for y in Z[i+1:]:rows[x]|=1<<y
 check_tournament(rows)
 return rows,Z[0]


def normalize_source(rows,r):
 C,c,j,_=normalize(rows,r)
 return k_graph(C,c,j),c


def normalize_partner(rows,r,L):
 if L==1:return [0],0
 C,c,j,_=normalize(converse(rows),r)
 return converse(k_graph(C,c,j)),c


def prepend_sources(rows,r,k):
 n=len(rows)
 out=[sum(1<<j for j in range(i+1,n+k)) for i in range(k)]
 out.extend(x<<k for x in rows)
 return out,r+k


def append_sinks(rows,r,k):
 out=rows[:]
 for _ in range(k):
  n=len(out);out=[x|(1<<n) for x in out]+[0]
 return out,r


def check_pair(X,r,Y,s,actual=True):
 h,e,R,d=profiles(X)[0][r];U=R-d;f=R+e-h-d
 g,starts,L,t=profiles(converse(Y))[0][s]
 beta=d*g-f*t
 require(U>0 and beta in [-1,1],'unit pair')
 require(d+2*f<=2*U-1,'general source endpoint budget')
 require(3*d<=2*U+2*beta-1,'general factor budget')
 require(f<=U-1,'bounded slack')
 if beta==-1:require(f>d,'negative determinant source')
 if beta==1 and f<d:require(g==t==1 and d-f==1,'positive source dichotomy')
 NX,nr=normalize_source(X,r);NY,ns=normalize_partner(Y,s,L)
 require(profiles(NX)[0][nr]==(h,e,R,d),'source profile preservation')
 ng,_,nL,nt=profiles(converse(NY))[0][ns]
 require((ng,nL,nt)==(g,L,t),'partner profile preservation')
 require(len(NX)<=source_bound(U,beta),'source order budget')
 require(len(NY)<=partner_bound(L),'partner order budget')
 if actual:
  for m in [1,2]:
   T,z=star(NX,nr,NY,ns,m)
   require(profiles(T)[0][0][0]==U*L*(1<<m)+beta,'physical star count')
 return U*L,beta


def main():
 sources={};partners={};source_mult=Counter();partner_mult=Counter()
 for n in range(1,5):
  for rows in labelled(n):
   ps,_=profiles(rows);qs,_=profiles(converse(rows))
   for r,z in enumerate(ps):
    h,e,R,d=z
    if R>d:
     sources.setdefault(z,(rows,r));source_mult[z]+=1
    gg,ee,L,t=qs[r];key=(gg,t,L)
    partners.setdefault(key,(rows,r));partner_mult[key]+=1
 pairs=0;labelled_pairs=0;actual=0;small_rays=set();chosen=[]
 for z,(X,r) in sources.items():
  h,e,R,d=z;U=R-d;f=R+e-h-d
  for key,(Y,s) in partners.items():
   g,t,L=key
   if d*g-f*t not in [-1,1]:continue
   pairs+=1;labelled_pairs+=source_mult[z]*partner_mult[key]
   small_rays.add(check_pair(X,r,Y,s));actual+=2
   chosen.append((X,r,Y,s))
 padded=0
 for X,r,Y,s in chosen[:24]:
  PX,pr=prepend_sources(k_graph(X,r,2),r,2)
  PY,ps=prepend_sources(Y,s,2);PY,ps=append_sinks(PY,ps,2)
  check_pair(PX,pr,PY,ps,actual=False);padded+=1
 examined,witnesses=classify([4,1,2],3,6)
 require(sorted(witnesses)==[9,15,27,36,60],'deletion3 slack6 table')
 for U,w in witnesses.items():
  z=tuple(w['profile'])
  for m in [1,2]:
   T=k_graph(w['rows'],w['root'],m)
   require(profiles(T)[0][w['root']]==k_profile(z,m),'slack6 K extensions')
 require(profiles([10,12,9,16,7])[0][0][0]==15,'H15 exists')
 # Explicit all-factor inverse budget for alpha=8, beta=-1.
 candidates=[]
 for U in [1,2,4,8]:
  for d in range(1,max(0,(2*U-3)//3)+1,2):
   candidates.append((U,d))
 require(candidates==[(4,1),(8,1),(8,3)],'all alpha8 factor budgets')
 require(all(U* (1<<j)!=8 for U in [9,15,27,36,60] for j in range(4)),
         'slack6 power-of-two exclusion')
 print(json.dumps({'status':'PASS','unique_unit_profile_pairs':pairs,
                   'represented_labelled_unit_pairs_through_order4':labelled_pairs,
                   'actual_star_graph_checks':actual,'larger_padded_pair_checks':padded,
                   'small_physical_rays':sorted(small_rays),
                   'deletion3_slack6':{'examined_reduced_cores':examined,'U_base':sorted(witnesses),
                                        'witnesses':witnesses},
                   'alpha8_negative_factor_candidates':candidates,
                   'scope':'Exact unit binary-amplifier bridge; full spectrum remains open.'},
                  indent=2,sort_keys=True))


if __name__=='__main__':main()
