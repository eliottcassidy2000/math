#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations,permutations
import ast,json,hashlib,subprocess,random
P=Path(__file__).resolve().parent

def need(x,s):
 if not x:raise RuntimeError(s)
def graph(n,z):
 rows=[0]*n
 for b,(i,j)in enumerate(combinations(range(n),2)):
  if z>>b&1:rows[i]|=1<<j
  else:rows[j]|=1<<i
 return rows
def valid(r):return all(not(x>>i&1)and x<(1<<len(r))and all(((x>>j)&1)+((r[j]>>i)&1)==1 for j in range(len(r))if j!=i)for i,x in enumerate(r))
def birth(r,t):
 n=len(r);out=r[:]+[1<<t]
 for i in range(n):
  if i!=t:out[i]|=1<<n
 return out
def nearjoin(r,t,k):
 n=len(r);out=[x|(((1<<k)-1)<<n)for x in r]+[sum(1<<j for j in range(n+i+1,n+k))for i in range(k)]
 return birth(out,t)
def literal(r,t):
 n=len(r);h=e=B=0
 for p in permutations(range(n)):
  bad={j for j in range(n-1)if not(r[p[j]]>>p[j+1]&1)}
  if not bad:h+=1;e+=p[-1]==t
  pos=p.index(t)
  if bad<={pos-1} if pos else not bad:B+=1
 return h,e,B

def main():
 requests=[];towers=[];near=[];root_records=[]
 for n in range(1,5):
  for z in range(1<<(n*(n-1)//2)):
   r=graph(n,z)
   for t in range(n):root_records.append((r,t,5))
 rng=random.Random(4865)
 for n in range(5,8):
  for _ in range(12):root_records.append((graph(n,rng.getrandbits(n*(n-1)//2)),rng.randrange(n),3))
 for r,t,height in root_records:
  inds=[];rows=r[:];root=t
  for k in range(height+1):
   need(valid(rows),'tower actual tournament');inds.append(len(requests));requests.append((rows,root));rows=birth(rows,root);root=len(rows)-1
  towers.append({'base':r,'root':t,'indices':inds})
  for k in range(1,4):
   rows=nearjoin(r,t,k);need(valid(rows),'nearjoin actual tournament');near.append({'base_index':inds[0],'k':k,'query':len(requests)});requests.append((rows,len(rows)-1))
 # Parse only the literal data assignment from the incoming verifier, not execute its code.
 tree=ast.parse((P/'input/verify_memberships.py').read_text());node=next(x for x in tree.body if isinstance(x,ast.Assign)and any(isinstance(y,ast.Name)and y.id=='WITNESSES'for y in x.targets));witnesses=ast.literal_eval(node.value);members=[]
 for name,(rows,digit)in witnesses.items():
  n=len(rows);need(valid(rows)and rows[0].bit_count()==1,'actual named readiness graph');root=rows[0].bit_length()-1;deleted=[x>>1 for x in rows[1:]];entry={'name':name,'rows':rows,'wanted_H':int(name.split()[0]),'digit':digit,'root':root,'whole':len(requests),'deleted':len(requests)+1};requests.extend([(rows,0),(deleted,root-1)])
  if name.startswith(('12875 in R1(14)','12877 in R1(14)')):
   child=birth(birth(rows,0),n);entry['C2_index']=len(requests);entry['C2_rows']=child;requests.append((child,len(child)-1))
  members.append(entry)
 # Rerooted tower is a constructive digit-one family, independent of any marked bijection.
 reroot=[]
 for n in range(4,14):
  rows=[0];root=0
  for _ in range(n-2):rows=birth(rows,root);root=len(rows)-1
  new=birth(rows,0);reroot.append({'n':n,'core':len(requests),'ready':len(requests)+1});requests.extend([(rows,0),(new,len(new)-1)])
 text=''.join(str(len(r))+' '+str(t)+' '+' '.join(map(str,r))+'\n'for r,t in requests);(P/'actual_graphs.in').write_text(text);subprocess.run(['g++','-O3','-std=c++17',str(P/'ham.cpp'),'-o',str(P/'ham')],check=True);data=subprocess.check_output([str(P/'ham')],input=text,text=True);(P/'actual_values.txt').write_text(data);values=[tuple(map(int,l.split()))for l in data.splitlines()];need(len(values)==len(requests),'all actual graph counts')
 for item in towers:
  base=values[item['indices'][0]][:3]
  if len(item['base'])<=4:need(base==literal(item['base'],item['root']),'literal rooted h/e/B')
  prev=None;hs=[]
  for idx in item['indices']:
   h,e,B,st=values[idx];hs.append(h)
   if prev is not None:
    H,E,b=prev;need((h,e,B)==(H-E+b,H-E,b+H),'tower three-state exact law')
   prev=h,e,B
  for k in range(3,len(hs)):need(hs[k]==sum(hs[k-3:k]),'actual tribonacci')
  if values[item['indices'][0]][3]and len(item['base'])>=3:need(all(values[i][3]for i in item['indices']),'strong tower extension')
 for item in near:
  h,e,B,st=values[item['base_index']];H,_,_,ST=values[item['query']];need(H==h+(1<<item['k'])*B,'actual nearjoin formula');need(not st or len(requests[item['base_index']][0])<3 or ST,'strong nearjoin')
 for item in members:
  h,e,B,st=values[item['whole']];digit=values[item['deleted']][1];need((h,digit,st,B-h)==(item['wanted_H'],item['digit'],1,item['digit']),'actual strong ready membership')
  item['verified_H']=h;item['verified_digit']=digit
  if 'C2_index'in item:need(values[item['C2_index']][0]==3*h+2*digit and values[item['C2_index']][3],'actual order16 +2 child');item['C2_H']=values[item['C2_index']][0]
 seq=[]
 for item in reroot:
  need(values[item['core']][1]==1 and values[item['ready']][3],'rerooted family digit-one strong');seq.append(values[item['ready']][0])
 need(seq==[5,13,31,77,193,489,1255,3253,8493,22293],'explicit reroot sequence')
 M=[[1,-1,1],[1,-1,0],[1,0,1]]
 def mm(A,B):return[[sum(A[i][k]*B[k][j]for k in range(3))for j in range(3)]for i in range(3)]
 M2=mm(M,M);M3=mm(M2,M);need(M3==[[M2[i][j]+M[i][j]+int(i==j)for j in range(3)]for i in range(3)],'all-order matrix polynomial')
 (P/'membership_controls.json').write_text(json.dumps(members,sort_keys=True,indent=2)+'\n');(P/'tower_controls.json').write_text(json.dumps({'towers':towers,'nearjoin':near,'rerooted':reroot,'rerooted_counts':seq},separators=(',',':'))+'\n');out={'status':'PASS','rooted_base_controls':len(towers),'actual_graph_counts':len(values),'nearjoin_controls':len(near),'readiness_memberships':len(members),'source_sha256':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['controls.py','ham.cpp']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['actual_graphs.in','actual_values.txt','membership_controls.json','tower_controls.json']}};(P/'controls_receipt.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
