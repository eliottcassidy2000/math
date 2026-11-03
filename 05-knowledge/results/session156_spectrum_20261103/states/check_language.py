from itertools import product
from bank import rows,strong
count=0
for m in range(0,8):
 for word in product(range(4),repeat=m):
  # Apex 0 beats 1. The empty word contains two vertices and is not strong.
  regular_non=(m>0 and (word[0]==3 or word[-1]==0)) or any(all(x in (2,3) for x in word[:p]) and all(x in (0,2) for x in word[p:]) for p in range(m+1))
  actual=not strong(rows(word,1))
  if regular_non!=actual:raise SystemExit(f'counterexample {word}')
  count+=1
print('verified exact cut-language against connectivity for',count,'words of lengths 0..7')
