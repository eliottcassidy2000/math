#!/usr/bin/env python3
"""Checks actual core routes; frozen helpers are co-delivered under inputs/."""
import sys,json,pathlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent/'inputs'))
from check_family import build,need,strong,anti_failure,apply,converse,literal_k,route_for
from check_odd import odd_build

def check(a,odd):
 t,u,r=(odd_build(a) if odd else build(a));m=a+4;w=2*m;n=len(t);v=a+2
 p,_=route_for(a)
 if odd:p=p[:-2]+(w,)+p[-2:]
 q=(r,a+1)+((w,) if odd else ())+tuple(range(2*m-1,m-1,-1))+tuple(range(a))+(a,a+2)
 curr=t;moves=[];need(strong(curr),'source not strong')
 for j,word in enumerate((p,q)):
  need(anti_failure(curr,word) is not None,'blue path');nxt=apply(curr,word)
  need(anti_failure(nxt,word[::-1]) is not None,'blue reversed path')
  ks=curr if not j else converse(curr);kw=word if not j else word[::-1];kn=literal_k(ks,kw)
  need(kn==(converse(nxt) if not j else nxt),'literal K');need(anti_failure(ks,kw) is not None and anti_failure(kn,kw) is not None,'K blue')
  moves.append({'R_word':word,'R_source':curr,'R_target':nxt,'failed_anti':anti_failure(curr,word),'K_word':kw,'K_source':ks,'K_target':kn})
  curr=nxt
  if not j:need(strong(curr),'intermediate nonstrong')
 need(curr[v].bit_count()==n-1,'named source not reached')
 need(strong(curr,[x for x in range(n) if x!=v]),'deleted core not strong')
 need(tuple((curr[m+x]>>m)&((1<<m)-1) for x in range(m))==u,'U not restored')
 need(strong(curr,[x for x in range(m) if x!=v]),'V minus source not strong')
 need(min(min(row.bit_count(),n-1-row.bit_count()) for row in t)==2,'delta changed')
 return {'a':a,'odd':odd,'n':n,'initial_rows':t,'source_label':v,'delta':2,'black_core_distance':2,'route':moves,'terminal_core_strong':True}

def main():
 print(json.dumps({'status':'PASS','scope':'literal two-black core routes on even and odd infinite-family controls a1..20','cases':[check(a,odd) for a in range(1,21) for odd in (False,True)]},sort_keys=True,indent=2))
if __name__=='__main__':main()
