#!/usr/bin/env python3
"""Odd middle-vertex extension; imports only the co-delivered stdlib verifier."""
import json
from check_family import build,arc,need,paths,anti_failure,apply,strong,converse,literal_k,cuts,mincut_score_certificate,route_for

def odd_build(a):
 t,u,r=build(a);m=a+4;n=2*m+1;w=2*m;t=list(t)+[0]
 # Undo direct bridge Uroot -> Vroot to restore all V -> U.
 t[r]^=1<<(m+r);t[m+r]^=1<<r
 # V -> w -> U, except Uroot -> w -> Vroot.
 for v in range(m):t[v]|=1<<w
 for v in range(m,2*m):t[w]|=1<<v
 t[r]^=1<<w;t[w]^=1<<r;t[w]^=1<<(m+r);t[m+r]^=1<<w
 return tuple(t),u,r

def check(a,exhaustive=False):
 t,u,r=odd_build(a);m=a+4;w=2*m;n=w+1;shore=((1<<m)-1)<<m;shore2=shore|(1<<w)
 need(strong(t),'initial not strong');ds=[x.bit_count() for x in t];need(min(min(d,n-1-d) for d in ds)==2,'delta')
 mc=mincut_score_certificate(t);need(mc['lambda']==1 and mc['size_minima']==[m,m+1],'minimum-cut sizes')
 order=sorted(range(n),key=lambda v:ds[v]);need(sum(1<<v for v in order[:m])==shore and sum(1<<v for v in order[:m+1])==shore2,'minshore labels');need(ds[order[m-1]]<ds[order[m]]<ds[order[m+1]],'minshore uniqueness')
 endings=[tuple(range(a))+tuple(a+(i+j)%3 for j in range(3))+(r,) for i in range(3)]
 dec=[tuple(x+m for x in p)+(w,)+q[::-1] for p in endings for q in endings]
 for p in dec:apply(t,p);need(anti_failure(t,p) is None,'decreasing word black')
 pe,qe=route_for(a);p=pe[:-2]+(w,)+pe[-2:];q=qe+(w,);curr=t;moves=[];states=[t]
 for j,word in enumerate((p,q)):
  need(anti_failure(curr,word) is not None,'route is blue');nxt=apply(curr,word);need(anti_failure(nxt,word[::-1]) is not None,'reverse route blue')
  ks=curr if j%2==0 else converse(curr);kw=word if j%2==0 else word[::-1];kn=literal_k(ks,kw);need(kn==(converse(nxt) if j%2==0 else nxt),'literal K')
  need(anti_failure(ks,kw) is not None and anti_failure(kn,kw) is not None,'K blue')
  moves.append({'R_word':word,'R_source':curr,'R_target':nxt,'failed_anti':anti_failure(curr,word),'K_word':kw,'K_source':ks,'K_target':kn});curr=nxt;states.append(curr)
 need(strong(states[1]) and not strong(states[2]),'strongness')
 need(sum((curr[v]&(((1<<n)-1)^shore)).bit_count() for v in range(m,2*m))==0,'terminal cut')
 endpoint_sum=states[1][q[0]].bit_count()+states[1][q[-1]].bit_count();need(endpoint_sum==n-4,'second endpoint sum')
 report={'a':a,'n':n,'initial_rows':t,'scores':ds,'delta':2,'minimum_cut':mc,'minimum_shores':[shore,shore2],'decreasing_words':dec,'all_decreasing_words_blue':True,'route':moves,'intermediate_strong':True,'terminal_strong':False,'Q_endpoint_score_sum':endpoint_sum,'black_nonstrong_distance':2}
 if exhaustive:
  pp=list(paths(t));cone=[p for p in pp if any((s>>p[0]&1) and not(s>>p[-1]&1) for s in (shore,shore2))];need(set(cone)==set(dec),'complete cone')
  bp=[p for p in pp if anti_failure(t,p) is not None];need(all(strong(apply(t,p)) for p in bp),'black one-step exit')
  cc=cuts(t);need([s for s,v in cc.items() if v==1]==[shore,shore2],'explicit cuts')
  report['exhaustive']={'all_global_hamilton_paths':len(pp),'black_paths':len(bp),'all_black_neighbors_strong':True,'decreasing_paths':len(cone),'direct_cuts_checked':len(cc)}
 return report

def main():
 print(json.dumps({'status':'PASS','scope':'odd all-order proof separately; a1..20 n11..49 controls, complete n11 Hamilton universe','cases':[check(a,a==1) for a in range(1,21)]},sort_keys=True,indent=2))
if __name__=='__main__':main()
