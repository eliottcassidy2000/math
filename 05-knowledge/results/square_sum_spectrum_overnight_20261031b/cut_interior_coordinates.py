#!/usr/bin/env python3
"""Literal-cover audit of positive cut/interior coordinates for marked blocks.

No incoming cut-response implementation is imported. The proof uses its audited
canonical weighted-cut law; these controls reconstruct every cover literally.
"""
from itertools import permutations
import json
from diamond_profiles import profiles,apply,single
from audit_feedback import require,subset_dp


def literal(rows):
 n=len(rows)
 def path(p):return all(rows[u]>>v&1 for u,v in zip(p,p[1:]))
 starts=[0]*n;ends=[0]*n;W=[[0]*n for _ in rows];H=0
 for word in permutations(range(n)):
  if path(word):
   H+=1;starts[word[0]]+=1;ends[word[-1]]+=1
  for k in range(1,n):
   if path(word[:k]) and path(word[k:]):W[word[k-1]][word[k]]+=1
 weights=[[0]*n for _ in rows]
 for u in range(n):
  for v in range(u):
   require((W[u][v]+W[v][u])%2==0,'integral symmetric cover weight')
   weights[u][v]=weights[v][u]=(W[u][v]+W[v][u])//2
 return H,starts,ends,W,weights


def convert(state):
 F,q,g,h=state
 require((F-h)%2==0 and (q-F+h-g)%2==0,'coordinate integrality')
 return ((F-h)//2,(q-F+h-g)//2,g,h)


def coefficients(b):
 H,E,p,s=b['H'],b['E'],b['p'],b['s']
 numerators=(p+H-s,E-H-p+s,b['F']-H,b['q']-b['F']+b['t'])
 require(all(x>=0 and x%2==0 for x in numerators),'positive integral block coefficients')
 alpha,beta,a,ell=(x//2 for x in numerators)
 require(p<=E-H and alpha+beta==E//2,'positive complementary coefficients')
 return alpha,beta,a,ell


def transport(state,b):
 C,L,g,h=state;H,E,p=b['H'],b['E'],b['p']
 alpha,beta,a,ell=coefficients(b)
 return ((H+p)*C+p*L+alpha*g+a*h,
         (E-H-p)*C+(E-p)*L+beta*g+ell*h,
         (H-b['t'])*h,H*h)


def main():
 blocks=profiles();parents=cuts=transports=coefficient_checks=0
 for n in range(1,6):
  edges=[(u,v) for u in range(n) for v in range(u+1,n)]
  for code in range(1<<len(edges)):
   rows=[0]*n
   for i,(u,v) in enumerate(edges):
    if code>>i&1:rows[u]|=1<<v
    else:rows[v]|=1<<u
   H,starts,ends,W,weights=literal(rows);parents+=1
   E=2*H+sum(map(sum,W))
   for mask in range(1<<n):
    inside=[v for v in range(n) if mask>>v&1]
    C=sum(weights[u][v] for u in inside for v in range(n)
          if not(mask>>v&1) and rows[u]>>v&1)
    L=sum(weights[u][v] for u in inside for v in inside if u<v)
    g=sum(ends[v] for v in inside)
    q=sum(starts[v]+sum(W[u][v] for u in range(n)) for v in inside)
    child=list(rows)+[mask]
    for v in range(n):
     if not(mask>>v&1):child[v]|=1<<n
    F=subset_dp(child)
    require(convert((F,q,g,H))==(C,L,g,H),'literal cut/interior interpretation')
    cuts+=1
    s=sum(starts[v] for v in inside);t=H-g
    p=sum(ends[u]+sum(W[u][v] for v in range(n)) for u in range(n) if not(mask>>u&1))
    b={'H':H,'E':E,'p':p,'s':s,'t':t,'q':q,'F':F}
    coefficients(b);coefficient_checks+=1
    state=(F,q,g,H);new=(C,L,g,H)
    require(convert(single(state,0))==(2*C+L+g,L,0,H),'singleton zero split')
    require(convert(single(state,1))==(C,C+2*L+g,H,H),'singleton one split')
    for block in blocks:
     require(convert(apply(state,block))==transport(new,block),'exact positive conjugacy')
     transports+=1
 print(json.dumps({'status':'passed','universe':'all labelled parent tournaments of orders1 through5, every marked cut',
                   'parents':parents,'cuts_with_literal_child_DP':cuts,
                   'positive_profile_checks':coefficient_checks,'diamond_transport_checks':transports,
                   'coordinates':['C=(F-h)/2','L=(q-F+h-g)/2','g','h'],
                   'necessary_guard':'q>=F-h+g and q-F+h-g even',
                   'scope':'The cut/interior identity is universal for nonempty parents; positivity alone is not actual-state reachability.'},indent=2,sort_keys=True))


if __name__=='__main__':main()
