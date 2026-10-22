"""Independent canonical-label enumeration of all supplied almostregular8 classes."""
import json,itertools,math
from pathlib import Path
HERE=Path(__file__).resolve().parent
bank=json.loads((HERE.parent/'session30-almost-regular/census8_bank.json').read_text())
seen=set();mass=0
for row in bank['almostregular8']:
 n=8;s=[[0]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):s[i][j]=1 if row['bits'][k]=='1'else-1;s[j][i]=-s[i][j];k+=1
 low=[i for i in range(n)if sum(s[i])==-1];high=[i for i in range(n)if sum(s[i])==1]
 if len(low)!=4 or len(high)!=4:raise RuntimeError('wrong score vector')
 least=None;aut=0
 for p0 in itertools.permutations(low):
  for p1 in itertools.permutations(high):
   p=p0+p1;code=0
   for i in range(n):
    for j in range(i+1,n):code=2*code+(s[p[i]][p[j]]==1)
   if least is None or code<least:least=code;aut=1
   elif code==least:aut+=1
 if least in seen:raise RuntimeError('isomorphic duplicate')
 if aut!=row['aut']:raise RuntimeError('wrong automorphism size')
 seen.add(least);mass+=math.factorial(8)//aut
if mass!=3230080 or len(seen)!=85:raise RuntimeError('incorrect mass or universe')
result={'classes':len(seen),'mass':mass,'method':'All 4!4! score-respecting relabellings per supplied order8 graph; distinct canonical integers and exact stabilizer counts. Degree-count recurrence inspected separately.','status':'PASS'}
(HERE/'census8-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
