from core import birth,require

def induced(rows,vertices):
 return [sum(((rows[u]>>v)&1)<<j for j,v in enumerate(vertices))for u in vertices]

def prune(seed,ports,word,keep):
 n0=len(seed);k=len(ports)
 require(set(range(n0))<=keep and keep<=set(range(n0+len(word))),'valid retained set containing whole initial seed')
 old=seed[:];op=ports[:];new=seed[:];np=ports[:]
 labels=list(range(n0));ix={v:v for v in labels};shadow={p:p for p in op};outword=[]
 for side,mask,drop in word:
  z=len(old);before=op[:]
  if z in keep:
   wins={shadow[p]for i,p in enumerate(before)if (mask>>i&1)and p in keep}
   sm=sum(1<<i for i,p in enumerate(np)if p in wins)
   sd=k if drop==k else np.index(shadow[before[drop]])
   zz=len(new);new,np=birth(new,np,side,sm,sd);outword.append((side,sm,sd))
   labels.append(z);ix[z]=zz
   if drop<k:
    del shadow[before[drop]];shadow[z]=zz
  else:
   if drop<k:shadow[z]=shadow.pop(before[drop])
  old,op=birth(old,op,side,mask,drop)
  require(set(shadow)==set(op) and set(shadow.values())==set(np),'port bijection')
  require(len(np)==k and len(set(np))==k,'distinct physical ports')
  for p in op:
   if p in keep:require(shadow[p]==ix[p],'surviving original port fixed')
   else:require(labels[shadow[p]] not in op,'dummy is original nonport')
  require(new==induced(old,labels),'exact induced graph after every step')
 return new,np,outword,shadow
