M=50000
# state h,R,e,d. Begin two-vertex edge rooted source.
cur={(1,2,0,1)}; support={1}; states=0
n=2
while cur:
 support.update(s[0] for s in cur)
 print(n,len(cur),len(support),[i for i in range(23,min(M,300)+1,2) if i not in support][:25],flush=True)
 nxt=set()
 for h,r,e,d in cur:
  # Any next non-O step has count >= h+r-e; O only raises r.
  if h+r-e<=M:
   nxt.add((h+r-e,2*r-d,d,d))
   nxt.add((h+r-e,h+r,h-e,h))
  if h+2*r<=M:
   nxt.add((h,2*r,0,d))
 cur=nxt
 n+=1
print('MISSING', [i for i in range(1,M+1,2) if i not in support])
