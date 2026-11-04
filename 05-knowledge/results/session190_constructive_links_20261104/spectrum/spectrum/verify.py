#!/usr/bin/env python3
"""Independent permutation/cut audit of the flat-response incidence patch."""
import hashlib,itertools,json,pathlib
HERE=pathlib.Path(__file__).resolve().parent

def need(test,msg):
 if not test:raise RuntimeError(msg)
def valid(a):
 n=len(a)
 need(all(0<=a[i]<1<<n and not a[i]>>i&1 for i in range(n)),'row bounds')
 need(all(((a[i]>>j)&1)+((a[j]>>i)&1)==1 for i in range(n)for j in range(i)),'tournament')
def strong(a):
 n=len(a)
 for x in range(n):
  S=1<<x
  while True:
   Q=S
   for y in range(n):
    if S>>y&1:Q|=a[y]
   if Q==S:break
   S=Q
  if S!=(1<<n)-1:return False
 return True

def permutation_profile(a,v):
 n=len(a);h=s=e=F=q=0
 def path(p):return all(a[x]>>y&1 for x,y in zip(p,p[1:]))
 for p in itertools.permutations(range(n)):
  good=path(p)
  h+=good;s+=good and p[0]==v;e+=good and p[-1]==v
  cut=p.index(v)+1
  F+=path(p[:cut])and path(p[cut:])
 q=sum(path(p)for p in itertools.permutations([i for i in range(n)if i!=v]))
 return dict(h=h,s=s,e=e,F=F,q=q)

def dp(a):
 n=len(a);N=1<<n;d=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for i in range(n):d[1<<i][i]=1
 # Predecessor recurrence, unlike the C++ forward extension.
 for M in range(1,N):
  if M&(M-1):
   for last in range(n):
    if M>>last&1:d[M][last]=sum(d[M^(1<<last)][p]for p in range(n)if M>>p&1 and a[p]>>last&1)
  h[M]=sum(d[M])
 return h,d

def profile(a,v):
 h,e=dp(a);op=[((1<<len(a))-1)^(1<<i)^a[i]for i in range(len(a))];_,s=dp(op);full=len(h)-1
 return dict(h=h[-1],s=s[-1][v],e=e[-1][v],q=h[full^(1<<v)],F=sum(e[M][v]*h[full^M]for M in range(full+1)))
def endpoint_cover_profile(a,v):
 n=len(a);h,e=dp(a);full=len(h)-1;op=[full^(1<<i)^a[i]for i in range(n)];_,s=dp(op)
 exact=dict(P2=sum(h[M]*h[full^M]for M in range(full+1)),S2=sum(s[M][v]*h[full^M]for M in range(full+1)),F_by_first_bin_size=[sum(e[M][v]*h[full^M]for M in range(full+1)if M.bit_count()==k)for k in range(n+1)])
 direct=dict(P2=0,S2=0,F_by_first_bin_size=[0]*(n+1))
 def path(p):return all(a[x]>>y&1 for x,y in zip(p,p[1:]))
 for p in itertools.permutations(range(n)):
  for cut in range(n+1):
   if path(p[:cut])and path(p[cut:]):
    direct['P2']+=1
    direct['S2']+=cut>0 and p[0]==v
    if cut>0 and p[cut-1]==v:direct['F_by_first_bin_size'][cut]+=1
 need(exact==direct,'independent full endpoint two-bin checks')
 return exact
def relabel(a,p):return[sum(((a[i]>>j)&1)<<k for k,j in enumerate(p))for i in p]
def wrapper(a,v,k):
 n=len(a);r=a+[0]*(k+1);w=n+k
 for z in range(n,n+k):
  r[z]=(1<<n)-1
  for zz in range(z+1,n+k):r[z]|=1<<zz
  r[w]|=1<<z
 for u in range(n):
  if u==v:r[u]|=1<<w
  else:r[w]|=1<<u
 return r

def main():
 A=[6,28,8,1,109,15,47];B=[4,25,2,37,109,7,47];v=2
 states=[A,[6,24,10,1,109,15,47],[4,25,10,1,109,15,47],[4,25,2,5,109,15,47],B]
 expect=[(29,1,157),(41,2,134),(29,1,89),(25,1,121),(31,1,157)]
 checked=[]
 for i,a in enumerate(states):
  valid(a);need(strong(a),'strong patch state');p=permutation_profile(a,v);need(p==profile(a,v),'independent permutation/DP');need((p['h'],p['s'],p['F'])==expect[i],'specified response');checked.append(dict(rows=a,profile=p))
  if i:need(sum((x^y).bit_count()for x,y in zip(a,states[i-1]))==2,'exactly one flipped arc')
 endpoint_covers=[endpoint_cover_profile(a,v)for a in(A,B)]
 need(endpoint_covers==[dict(P2=580,S2=34,F_by_first_bin_size=[0,15,27,27,27,27,24,10]),dict(P2=644,S2=31,F_by_first_bin_size=[0,9,25,31,28,30,25,9])],'distinct full endpoint responses')
 wrappers=[]
 for k in range(1,5):
  vals=[]
  for a,p in zip(states,checked):
   r=wrapper(a,v,k);valid(r);need(strong(r),'strong wrapper');h,_=dp(r);pred=p['profile']['h']+(1<<k)*p['profile']['F'];need(h[-1]==pred,'literal wrapper identity');vals.append(h[-1])
  need(vals[-1]-vals[0]==2 and max(vals[1:-1])<vals[0],'nonovershooting sibling patch');wrappers.append(dict(k=k,counts=vals))
 # Exact rooted comparison: 6! possible identifications; no converse quotient.
 min_dist=99;align=[]
 for pp in itertools.permutations([x for x in range(7)if x!=v]):
  p=list(pp);p.insert(v,v);r=relabel(B,p);ed=[(i,j)for i in range(7)for j in range(i+1,7)if(A[i]^r[i])>>j&1]
  if len(ed)<min_dist:min_dist=len(ed);align=[]
  if len(ed)==min_dist:align.append((p,r,ed))
 need(min_dist==4 and len(align)==5,'minimum rooted incidence edit')
 strong_paths=unit_paths=0;min_peak=99;matrix=[]
 for p,r,edges in align:
  rec=[]
  for M in range(16):
   z=A[:]
   for k,(i,j)in enumerate(edges):
    if M>>k&1:z[i]^=1<<j;z[j]^=1<<i
   a=permutation_profile(z,v);rec.append(dict(mask=M,rows=z,strong=strong(z),profile=a))
  for order in itertools.permutations(range(4)):
   m=0;masks=[0]
   for k in order:m|=1<<k;masks.append(m)
   if all(rec[m]['strong']for m in masks):
    strong_paths+=1;peak=max(rec[m]['profile']['s']for m in masks);min_peak=min(min_peak,peak);unit_paths+=peak==1
  matrix.append(dict(permutation=p,edges=edges,states=rec))
 need(strong_paths==38 and unit_paths==0 and min_peak==2,'sharp temporary start budget')
 twins=[]
 for a,b,v,h,F in [([22,12,8,33,14,1367,159,63,735,255,991],[6,12,8,49,7,1111,927,63,191,447,991],3,575,3819),([6,12,8,849,39,15,183,63,247,1527,511],[22,12,8,161,14,87,31,1399,639,255,895],2,605,15909)]:
  for r,hh in[(a,h),(b,h+2)]:valid(r);need(strong(r),'large twin strong');p=profile(r,v);need((p['h'],p['s'],p['F'])==(hh,1,F),'large twin profile')
  twins.append(dict(A=a,B=b,root=v,H=h,F=F,necessary_parent_H=(h-2)//3,necessary_parent_F=F//3))
 receipt=dict(status='PASS',patch=checked,named_two_bin_endpoint_profiles=endpoint_covers,wrapper_controls=wrappers,rooted_identifications=720,min_arc_edits=min_dist,shortest_alignments=5,all_shortest_strong_paths=strong_paths,minimum_peak_root_start=min_peak,unit_only_shortest_paths=unit_paths,larger_twins=twins,complete_shortest_patch_data=matrix)
 (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print('PASS: independent seven-vertex permutations/cuts; five literal strong states; distinct full endpoint two-bin responses; 20 wrapper DPs; all720 rooted identifications; all80 intermediate relations; all120 shortest orders; two larger twin controls.')
if __name__=='__main__':main()
