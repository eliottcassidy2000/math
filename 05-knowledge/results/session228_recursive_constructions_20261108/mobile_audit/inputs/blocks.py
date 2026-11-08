from inherited import profile

def graph(bits,pos):
 # expanded binary word, replace positions pos,pos+1 by triangle; root outdegree into triangle is binary digit 2b+b
 blocks=[]; out=[];i=0
 while i<len(bits):
  if i==pos:
   k=2*bits[i]+bits[i+1]; blocks.append([2,4,1]);out.append(set(range(k)));i+=2
  else:blocks.append([0]);out.append({0} if bits[i] else set());i+=1
 rows=[0];off=1
 for b,o in zip(blocks,out):
  for j,v in enumerate(b):
   rows.append((v<<off)|sum(1<<q for q in range(off+len(b),len(bits)+2))|(1 if j not in o else 0))
   if j in o:rows[0]|=1<<(off+j)
  off+=len(b)
 return rows
if __name__=='__main__':
 for m in range(3,12):
  bits=[i%2 for i in range(m-1)]+[1]
  vals=[profile(graph(bits,p),0) for p in range(m-2)]
  print(bits, vals)
