import sys,csv
from collections import Counter
from math import factorial
from verify import cover,rigid
def main(path):
 M={49:Counter(),63:Counter(),343:Counter()}
 for row in csv.DictReader(open(path),delimiter='\t'):
  h=int(row['h'])
  if h not in M:continue
  t=tuple(map(int,row['outgoing_masks'].split(',')))
  q=cover(t);d=[factorial(i)*v for i,v in enumerate(q)]
  first=next((2*i for i,a in enumerate(d) if a%7),None)
  M[h][first]+=1
  if h in (49,63):assert rigid(t)
 assert M=={49:Counter({2:5,None:2}),63:Counter({2:2}),343:Counter({2:10,4:2})}
 print(M)
if __name__=='__main__':main(sys.argv[1])
