V=[(a,b) for b in range(3) for a in range(7)]
D={1,2,4}; S=[{1,2,3},{0,1,2,4},D]
def A(u,v):
 a,b=u;c,d=v
 return (pow(2,-b,7)*(c-a))%7 in S[(d-b)%3]
def need(ok):
 if not ok: raise ValueError('control failed')
for u in V:
 for v in V:
  if u==v:continue
  need(A(u,v)!=A(v,u))
  common=sum(A(u,w) and A(v,w) for w in V)
  need(common==sum(A(w,u) and A(w,v) for w in V))
  need((common==7)==(u[0]==v[0]))
  if u[0]!=v[0]:need(common in {3,4,5})
for a in range(7):
 for d in D:
  c=(a+d)%7; z=[1,2,4].index(d)
  exceptions=[(b,e) for b in range(3) for e in range(3) if not A((a,b),(c,e))]
  need(exceptions==[((z-2)%3,(z-2)%3)])
print('Fano triangle and unique dissent port verified')
def twisted(u,v,mask):
 return A(u,v)^bool(u[0]==v[0] and u!=v and ((mask>>u[0])&1))
bad=[]
for mask in range(128):
 wrong=[(u,v) for u in V for v in V if u[0]!=v[0] and sum(twisted(u,w,mask) and twisted(v,w,mask) for w in V)==7]
 if wrong:bad.append((mask,wrong[0]))
print('twist masks breaking canonical triangles',len(bad),bad[:1])
def orbit(m):
 return {sum(((m>>a)&1)<<((t+pow(2,b,7)*a)%7) for a in range(7)) for t in range(7) for b in range(3)}
print('representatives up to complement size:',[(m,len(orbit(m))) for m in range(128) if m.bit_count()<4 and m==min(orbit(m))])
for b in range(3):
 vals=[]
 for a in range(1,7):
  u=(0,0);v=(a,b)
  vals.append(max(sum(twisted(u,w,m) and twisted(v,w,m) for w in V) for m in (0,1,1<<a,1+(1<<a))))
 print('cross maxima b=',b,vals)
need(not bad)
reps=[(m,len(orbit(m))) for m in range(128) if m.bit_count()<4 and m==min(orbit(m))]
need(reps==[(0,1),(1,7),(3,21),(7,21),(11,7),(13,7)])
