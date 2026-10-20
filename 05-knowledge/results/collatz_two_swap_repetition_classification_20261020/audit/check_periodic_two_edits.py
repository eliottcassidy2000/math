import itertools,json

def carry(w):
    # Independently use positional sum, rather than edit perturbations.
    p=sum(w); seen=0; out=0
    for i,t in enumerate(w):
        if t:
            seen+=1; out+=(1<<i)*3**(p-seen)
    return out

def swap(w,i):
    j=(i+1)%len(w)
    if w[i]==w[j]: return None
    v=list(w); v[i],v[j]=v[j],v[i]
    return tuple(v)
counts={'bases':0,'sources':0,'paths':0,'divisible':0,'undo':0,'mirrored':0}
for q0 in range(2,8):
  for v in itertools.product((0,1),repeat=q0):
    p0=sum(v)
    if p0 in (0,q0):continue
    counts['bases']+=1
    A=2**q0; B=3**p0
    for d in (2,3,4):
      w=v*d; q=len(w); H=sum(A**(d-1-j)*B**j for j in range(d))
      if carry(w)%H:raise RuntimeError('source factor')
      counts['sources']+=1
      for i in range(q):
        u=swap(w,i)
        if u is None:continue
        for j in range(q):
          z=swap(u,j)
          if z is None:continue
          counts['paths']+=1
          div=carry(z)%H==0
          expected=(z==w) or (d==2 and (j-i)%q==q0)
          if div!=expected:raise RuntimeError((q0,v,d,i,j,H,carry(z),expected))
          if div:
            counts['divisible']+=1
            counts['undo' if z==w else 'mirrored']+=1
            if z!=w and z[:q0]!=z[q0:]:raise RuntimeError('not repeated')
print(json.dumps(counts,sort_keys=True,indent=2))
