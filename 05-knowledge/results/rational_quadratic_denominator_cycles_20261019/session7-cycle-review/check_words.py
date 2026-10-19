import itertools,json

def primitive(w):
 n=len(w)
 return not any(n%d==0 and all(w[i]==w[i%d] for i in range(n)) for d in range(1,n))
def exponent(w,i,j):
 n=len(w)
 if w[i%n]==w[j%n]:return 0
 return next(k for k in range(1,n) if w[(i+k)%n]!=w[(j+k)%n])+0
# Opposite branches at time zero: k=1+L(next), i.e. first mismatch at positive time.
def table(N):
 return {(0,0,0)}|{(0,1,k) for k in range(2,N)}|{(0,k,1) for k in range(2,N)}|{(k,k-1,0) for k in range(2,N)}|{(1,0,k) for k in range(1,N)}
result=[]
for n in range(3,19):
 vals=set();count=0
 for w in itertools.product((0,1),repeat=n-1):
  w=(0,)+w
  if not primitive(w):continue
  count+=1
  vals.add((exponent(w,0,1),exponent(w,0,2),exponent(w,1,2)))
 missing=vals-table(n)
 result.append(dict(n=n,primitive_words_up_to_complement=count,actual_signature_count=len(vals),superset_size=len(table(n)),missed=sorted(missing),extra=sorted(table(n)-vals)))
 assert not missing
print(json.dumps(result,indent=2))
