import contextlib,io,itertools,json
with contextlib.redirect_stdout(io.StringIO()):
 from replay import Field,H,matrices

def switch(A,bits):
 n=len(A);B=[0]*n
 for i in range(n):
  for j in range(n):
   if i!=j and bool(A[i]>>j&1) ^ bool((bits>>i&1)^(bits>>j&1)):B[i]|=1<<j
 return B

def power(f,x,n):
 ans=1
 while n:
  if n&1:ans=f.mul(ans,x)
  x=f.mul(x,x);n//=2
 return ans
f=Field(9);A=[14,12,8,0];xs=f.embed(A);bank={}
for a in range(f.q):
 if a in xs:continue
 bits=sum((f.chi[f.sub(x,a)]==-1)<<i for i,x in enumerate(xs))
 bits=min(bits,bits^15)
 if bits in bank:continue
 ys=[f.sub(0,power(f,f.sub(x,a),f.q-2)) for x in xs]
 B=[sum(f.edge(x,y)<<j for j,y in enumerate(ys) if i!=j) for i,x in enumerate(ys)]
 assert B==switch(A,bits)
 bank[bits]={'pole':a,'image':ys,'matrix':B,'H':H(B)}
 if len(bank)==8:break
assert len(bank)==8
checks=0;profiles={}
for A in matrices(5):
 hist={}
 for b in range(16):
  h=H(switch(A,b));hist[h]=hist.get(h,0)+1
 assert sum(h*c for h,c in hist.items())==120
 profiles[tuple(sorted(hist.items()))]=A;checks+=1
print(json.dumps({'q':f.q,'polynomial':f.mod+[1],'tuple':xs,'switch_bank':bank,'five_vertex_sum_checks':checks,'five_vertex_switch_histograms':[{'histogram':dict(h),'example':A} for h,A in profiles.items()]},indent=2))
