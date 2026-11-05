import collections, itertools, json, math

def need(x,s):
 if not x:raise RuntimeError(s)
def tournament(n,code):
 rows=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   u,v=(i,j)if code>>k&1 else(j,i);rows[u]|=1<<v;k+=1
 return tuple(rows)
def dp(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)];H=[0]*(1<<n);H[0]=1
 for i in range(n):d[1<<i][i]=1
 for S in range(1,1<<n):
  H[S]=sum(d[S])
  for i,z in enumerate(d[S]):
   if z:
    todo=a[i]&~S
    while todo:
     b=todo&-todo;todo-=b;d[S|b][b.bit_length()-1]+=z
 return H,d
def profiles(a):
 n=len(a);H,d=dp(a);op=tuple(sum(((a[j]>>i)&1)<<j for j in range(n))for i in range(n));_,start=dp(op);full=(1<<n)-1
 return [(H[-1],d[-1][r],sum(start[S][r]*H[full^S]for S in range(1,full+1)))for r in range(n)]
def gadget(L):
 s=(L-1).bit_length();a=[sum(1<<j for j in range(i+1,s))for i in range(s)]+[0]
 for i in range(s):
  if (L-1)>>i&1:a[i]|=1<<s
  else:a[s]|=1<<i
 return tuple(a)
def seam(a,r,b):
 n=len(a);m=len(b);out=[x|(((1<<m)-1)<<n)for x in a]+[x<<n for x in b];out[r]^=1<<(n+m-1);out[n+m-1]|=1<<r;return tuple(out)
def exact_g(L):
 n=0;f=1
 for bit in bin(L-1)[2:]:
  if bit=='0':f+=n+n%2;n*=2
  else:n=2*n+1
 return f
def automaton(M):
 initial=(0,1);prev={initial:None};queue=collections.deque([initial])
 while queue:
  n,f=queue.popleft()
  for bit,state in [(0,(2*n%M,(f+n+n%2)%M)),(1,((2*n+1)%M,f))]:
   if state not in prev:prev[state]=((n,f),bit);queue.append(state)
 need(len(prev)==M*M//2 and all(f%2 for n,f in prev),'complete physical binary-word automaton')
 return prev
def witness(state,prev):
 bits=[]
 while prev[state]is not None:state,b=prev[state];bits.append(str(b))
 return int(''.join(reversed(bits))or'0',2)+1

def run():
 P=set();literal=0;small=[]
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   a=tournament(n,code);v=profiles(a)
   for r,p in enumerate(v):
    h,e,R=p;need(h%2 and R>=h>=1 and 0<=e<=h and (R-e)%2==0,'all physical prefix invariants');P.add(p);literal+=1
    if n<=3:small.append((a,r,p))
 # A literal named prefix from L32: transitive order with its unique long arc reversed.
 a=[sum(1<<j for j in range(i+1,7))for i in range(7)];a[0]^=1<<6;a[6]|=1;a=tuple(a);p=profiles(a)[6];need(p==(33,1,33),'literal composite-branch prefix');P.add(p)
 gadget_checks=seam_checks=0
 for L in range(1,129):
  b=gadget(L);H,ending=dp(b);full=(1<<len(b))-1;anchored=sum(ending[S][-1]*H[full^S]for S in range(1,full+1));need((H[-1],ending[-1][-1],anchored)==(exact_g(L),1,L),'literal one-apex gadget profile');gadget_checks+=1
 for a,r,p in small+[(a,6,p)]:
  for L in range(1,17 if len(a)<=3 else 9):
   H=dp(seam(a,r,gadget(L)))[0][-1];h,e,R=p;need(H==h*exact_g(L)+R*L-e*(L%2),'literal seam count');seam_checks+=1
 P=sorted(P);images=actual=states=0;maxbits=0;moduli=list(range(1,41))+[45,49,63,64,65,81,99,101]
 chosen=[(1,1,1),(3,1,3),(9,3,9),(33,1,33),(13,2,18)]
 for p in chosen:need(p in P,'named example is an actual rooted profile')
 for m in moduli:
  M=2*m;prev=automaton(M);states+=len(prev)
  for h,e,R in P:
   image={ (h*f+R*(n+1)-e*((n+1)%2))%m for n,f in prev}
   D=math.gcd(math.gcd(h,R),m);expected={r for r in range(m)if(m%2 or r%2)and(r%D==0 or(r+e)%D==0)}
   need(image==expected,'actual binary-reachable residue image')
   size=(m//D if m%2 else m//(2*D))*(1 if e%D==0 else 2);need(len(image)==size,'projected cardinality');images+=1
  if m<=24:
   for h,e,R in chosen:
    selected={}
    for state in prev:
     n,f=state;r=(h*f+R*(n+1)-e*((n+1)%2))%m
     if r not in selected:selected[r]=state
    for r,state in selected.items():
     L=witness(state,prev);g=exact_g(L);need(((L-1)%M,g%M)==state,'literal integer lifts BFS state');need((h*g+R*L-e*(L%2))%m==r,'actual shortest-word local witness');actual+=1;maxbits=max(maxbits,(L-1).bit_length())
 need(143%3==2 and 143%11==0 and 143%33==11,'mixed-prime branch hostile')
 need(all(143%p in(0,p-1)for p in(3,11))and 143%33 not in(0,32),'same branch required')
 print(json.dumps({'status':'PASS','literal_rooted_profiles':literal,'distinct_physical_profiles':len(P),'literal_gadget_profiles':gadget_checks,'literal_seam_counts':seam_checks,'moduli':moduli,'reachable_binary_states':states,'exact_image_and_cardinality_checks':images,'actual_BFS_integer_witnesses':actual,'max_witness_bit_length':maxbits,'composite_branch_hostile':[143%3,143%11,143%33]},sort_keys=True,indent=2))
if __name__=='__main__':run()
