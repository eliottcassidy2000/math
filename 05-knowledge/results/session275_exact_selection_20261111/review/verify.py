"""Independent small controls. No imports from producer scripts."""
import itertools,json,pathlib,collections,hashlib
HERE=pathlib.Path(__file__).resolve().parent;PACKET=HERE.parent
checks=0
def need(v):
 global checks
 checks+=1
 if not v:raise ValueError('failed check '+str(checks))
def valid(a):
 n=len(a)
 return all(not(a[i]>>i&1) for i in range(n)) and all(((a[i]>>j)&1)+((a[j]>>i)&1)==1 for i in range(n) for j in range(i))
def paths(a,vertices=None):
 vertices=tuple(range(len(a))) if vertices is None else tuple(vertices)
 return [p for p in itertools.permutations(vertices) if all(a[x]>>y&1 for x,y in zip(p,p[1:]))]
def packet(a):
 n=len(a); sets={mask:paths(a,[i for i in range(n) if mask>>i&1]) for mask in range(1<<n)};full=(1<<n)-1
 H=len(sets[full]); R=full^1;h=len(sets[R]);Q=0;U=V=d=e=0
 for mask,P in sets.items():
  comp=len(sets[full^mask])
  if mask&1:
   start=sum(bool(p) and p[0]==0 for p in P);end=sum(bool(p) and p[-1]==0 for p in P)
   sign=(-1)**(mask.bit_count()+1)
   U+=start*comp;V+=end*comp;d+=sign*start*comp;e+=sign*end*comp
  else:Q+=len(P)*len(sets[R^mask])
 return dict(H=H,h=h,Q=Q,U=U,V=V,d=d,e=e)
def two(a):
 n=len(a);return sum(len(paths(a,[i for i in range(n) if mask>>i&1]))*len(paths(a,[i for i in range(n) if not(mask>>i&1)])) for mask in range(1<<n))
def canon(a):
 return min(tuple((a[p[i]]>>p[j])&1 for i in range(len(a)) for j in range(i+1,len(a))) for p in itertools.permutations(range(len(a))))
def small(n):
 E=[(i,j) for i in range(n) for j in range(i+2,n)];m=len(E);rows=[]
 for mask in range(1<<m):
  a=[0]*n
  for i in range(n-1):a[i]|=1<<(i+1)
  for k,(i,j) in enumerate(E):
   x,y=(j,i) if mask>>k&1 else (i,j);a[x]|=1<<y
  rows.append(a)
 values=[len(paths(a)) for a in rows];codes=[canon(a) for a in rows]
 moves={1<<i for i in range(m)}|{(1<<i)^(1<<j) for i in range(m) for j in range(i)}|{(1<<m)-1}
 for subset in range(1<<n):moves.add(sum(1<<k for k,(i,j) in enumerate(E) if ((subset>>i)^(subset>>j))&1))
 adjacency=[set() for _ in rows]
 for i in range(len(rows)):
  for j in range(i):
   if codes[i]==codes[j] or (i^j in moves and abs(values[i]-values[j])<=2):adjacency[i].add(j);adjacency[j].add(i)
 unseen=set(range(len(rows)));comps=[]
 while unseen:
  todo=[unseen.pop()];seen=set(todo)
  while todo:
   i=todo.pop()
   for j in adjacency[i]&unseen:unseen.remove(j);seen.add(j);todo.append(j)
  vs=set(values[i] for i in seen);need(vs==set(range(min(vs),max(vs)+1,2)));comps.append((min(vs),max(vs),len(seen)))
 return dict(n=n,classes=len(set(codes)),components=len(comps),intervals=sorted(comps))
known=json.loads((PACKET/'landscape/RESULTS.json').read_text());small_results=[]
for n in range(3,6):
 s=small(n);r=known[n-3];need(s['classes']==r['isomorphism_classes']);need(s['components']==r['graphs']['reframed_full']['components']);small_results.append(s)
for r in known:
 runs=[]
 for h in r['spectrum']:
  if runs and h==runs[-1][1]+2:runs[-1][1]=h
  else:runs.append([h,h])
 need(len(runs)==r['graphs']['reframed_full']['components'])
 need(sorted((x['min'],x['max']) for x in r['graphs']['reframed_full']['largest_intervals'])==[tuple(x) for x in runs])
seed=(14,60,56,48,33,1);need(valid(seed));need(len(paths(seed))==23)
letters={'A':(6,16,2,7,13,31),'B':(4,37,0,7,15,29),'X':(46,8,2,4,15,30),'Y':(46,20,0,6,13,30)}
profiles={k:packet(a) for k,a in letters.items()}
expected={'A':dict(H=9,h=5,Q=64,U=16,V=32,d=2,e=0),'B':dict(H=5,h=5,Q=64,U=10,V=32,d=0,e=0),'X':dict(H=9,h=3,Q=48,U=36,V=9,d=0,e=3),'Y':dict(H=9,h=3,Q=48,U=42,V=9,d=0,e=3)}
for k,a in letters.items():need(valid(a));need(profiles[k]==expected[k])
need(two(letters['A'])==168);need(two(letters['B'])==128)
# Independent recursive endpoint counter for literal composed controls.
from functools import lru_cache
def count_dp(a):
 @lru_cache(None)
 def end(mask,v):
  if mask==1<<v:return 1
  before=mask^(1<<v)
  return sum(end(before,u) for u in range(len(a)) if before>>u&1 and a[u]>>v&1)
 return sum(end((1<<len(a))-1,v) for v in range(len(a)))
def compose(a,b):
 n=len(a); out=list(a)+[0]*(len(b)-1)
 mapping=[0]+list(range(n,n+len(b)-1))
 for i in range(len(b)):
  for j in range(len(b)):
   if b[i]>>j&1:out[mapping[i]]|=1<<mapping[j]
 for i in range(1,n):
  for j in range(n,len(out)):out[i]|=1<<j
 return tuple(out)
paircounts=[]
for x,y in itertools.product('BA',repeat=2):
 a=compose(letters[x],letters[y]);need(valid(a)); h=count_dp(a)
 need(h==185+(172 if x=='A' else 0)+(20 if y=='A' else 0));paircounts.append(h)
for spacer in [(0,),(2,0),(0,1),(2,4,1)]:
 cube=[]
 for x,y in itertools.product('BA',repeat=2):cube.append(count_dp(compose(compose(letters[x],spacer),letters[y])))
 need(cube[0]-cube[1]-cube[2]+cube[3]==0)
B=(108,41,2,4,47,12,62);T=(2,508,433,165,9,17,189,49,249)
need(valid(B));need(valid(T));need(len(paths(B))==29);need(two(B)==580);need(len(paths(T))==609)
# Three literal unit-response moves stitch the two order-seven corridors.
portal=json.loads((PACKET/'landscape/PORTAL63.json').read_text())
bridge=json.loads((PACKET/'landscape/BRIDGE63.json').read_text())
portal_receipt=[]
for item in portal['chain']:
 a=item['before'];b=item['after_raw'];c=item['after_canonical']
 changed=[(i,j) for i in range(8) for j in range(i+1,8) if ((a[i]^b[i])>>j)&1]
 need(changed==[tuple(x) for x in item['flips']]);need(valid(a) and valid(b) and valid(c))
 need(len(paths(a))==item['before_H']);need(len(paths(b))==item['after_H'])
 pp=paths([a[i]&b[i] for i in range(8)]);need(bool(pp))
 phi=next(q for q in itertools.permutations(range(8)) if all(((b[q[i]]>>q[j])&1)==((c[i]>>j)&1) for i in range(8) for j in range(8)))
 portal_receipt.append(dict(counts=[item['before_H'],item['after_H']],common_path=pp[0],canonical_map=phi))
a=bridge['seed'];b=bridge['neighbors']['65'][0]['rows']
need(len(paths(a))==63);need(len(paths(b))==65);pp=paths([a[i]&b[i] for i in range(8)]);need(bool(pp));need(b[6]==191)
need(portal['found'][7]==127)
portal_receipt.append(dict(counts=[63,65],common_path=pp[0],source=6))
# New bifan invariant and moving-root witnesses.
def seam_packet(a,root):
 old=[i for i in range(len(a)) if i!=root];allset=set(old);q=R=L=0
 outgoing={i for i in old if a[root]>>i&1}
 for mask in range(1<<len(old)):
  S=[old[j] for j in range(len(old)) if mask>>j&1];P=paths(a,S);hc=len(paths(a,allset-set(S)))
  q+=len(P)*hc
  R+=sum(bool(p) and p[0] in outgoing for p in P)*hc
  L+=sum(bool(p) and p[-1] in outgoing for p in P)*hc
 full=paths(a,old);h=len(full);s=sum(bool(p) and p[0] in outgoing for p in full);e=sum(bool(p) and p[-1] in outgoing for p in full)
 return (q,-L-h,R-q,len(paths(a))+L+2*h-s-R+e)
def bifan(a,root,m,n):
 k=len(a);X=list(range(k,k+m));Y=list(range(k+m,k+m+n));b=list(a)+[0]*(m+n)
 for x in X:
  b[root]|=1<<x
  for u in range(k):
   if u!=root:b[x]|=1<<u
 for y in Y:
  b[y]|=1<<root
  for u in range(k):
   if u!=root:b[u]|=1<<y
 for block in [X,Y]:
  for i,x in enumerate(block):
   for y in block[i+1:]:b[x]|=1<<y
 for x in X:
  for y in Y:b[x]|=1<<y
 return tuple(b)
seam_examples=[((96,97,59,115,35,0,52),0,(168,-104,-137,109)),((96,97,59,115,35,0,52),4,(168,-134,-107,107)),((58,36,17,6,10,28),1,(84,-31,-50,41)),((38,16,34,7,45,10),5,(84,-34,-50,43)),((44,25,34,36,45,2),2,(84,-33,-60,41)),((16,5,17,7,74,31,175,63),4,(128,-26,-52,13)),((16,5,17,23,66,31,175,63),4,(128,-18,-60,11))]
for a,r,p in seam_examples:
 need(valid(a));need(seam_packet(a,r)==p)
 for m,n in [(1,1),(1,2),(2,1),(2,2)]:
  out=bifan(a,r,m,n);need(valid(out));need(count_dp(out)==p[0]*2**(m+n)+p[1]*2**m+p[2]*2**n+p[3])
# Entire vertex-transitive triangle apex-cut cube, not an assumed common frame.
for cut in range(8):
 a=list((2,4,1))+[cut]
 for v in range(3):
  if not(cut>>v&1):a[v]|=8
 a=tuple(a);p=seam_packet(a,3);need(p[0]==12 and p[1]+p[2]==-15)
 need(p[3]==count_dp(a)+6);need(count_dp(bifan(a,3,1,1))==count_dp(a)+24)
# New mixed-gcd bank and quadratic gates.
N=(26,32,3,6,14,29);P=(52,41,2,5,14,28)
need(packet(N)==dict(H=37,h=9,Q=84,U=59,V=43,d=5,e=7));need(packet(P)==dict(H=37,h=9,Q=84,U=62,V=43,d=4,e=7))
def dual(a):return tuple(sum(1<<j for j in range(len(a)) if a[j]>>i&1) for i in range(len(a)))
for spacer,expected in [((0,),(75,115,115,195)),((2,0),(125,225,205,377)),((2,4,1),(275,475,475,819))]:
 got=tuple(count_dp(compose(compose(letters[x],spacer),dual(letters[y]))) for x,y in itertools.product('BA',repeat=2));need(got==expected)
# The literal portal's single flip stays outside the shared frame after fan insertion.
low=(16,5,17,23,66,31,175,63);high=(16,5,17,7,74,31,175,63)
for m in (1,2,3):
 lo=bifan(low,4,m,m);hi=bifan(high,4,m,m)
 need([(i,j) for i in range(len(lo)) for j in range(i+1,len(lo)) if ((lo[i]^hi[i])>>j)&1]==[(3,4)])
 need(count_dp(hi)-count_dp(lo)==2)
# Re-marked content escape. Move original vertex7 to index0 for packet notation.
p=(7,0,1,2,3,4,5,6,8);marked=tuple(sum(1<<j for j in range(9) if T[p[i]]>>p[j]&1) for i in range(9));pr=packet(marked)
need((pr['H'],pr['U'],pr['h'])==(609,1066,241))
# Here weighted V is not endpoint E; directly count full paths ending at root.
need(sum(x[-1]==0 for x in paths(marked))==108)
for k in range(1,4):need(count_dp(bifan(marked,0,0,k))==825*2**k-83)
need(pow(2,28,29)==1);need(all(pow(2,d,29)!=1 for d in [1,2,4,7,14]));need((2**28-1)%29==0 and (2**28-1)%29**2!=0)
receipt=dict(status='PASS',checks=checks,landscape_independent_orders=small_results,landscape_higher_order_scope='Inspected producer algorithm and compared reported maximal intervals; no independent full n6/n7 replay here',letters=profiles,cofactor=dict(H=29,Q=580,graftH=609),portal63=portal_receipt,producer_hashes={str(p.relative_to(PACKET)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PACKET/'landscape/experiment.py',PACKET/'switches/PROOF.md',PACKET/'descent/PROOF.md',PACKET/'descent/ESCAPE.md',PACKET/'switches/EXPANSION.md',PACKET/'seams/PROOF.md']})
(HERE/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
