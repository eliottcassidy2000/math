from itertools import permutations
E4={(1,2),(1,3),(4,1),(2,3),(2,4),(4,2),(3,4),(4,3)}
E5={(1,2),(1,3),(1,4),(5,1),(2,3),(4,2),(3,4),(2,5),(5,2),(3,5),(5,3),(4,5),(5,4)}
for n,E,want in [(4,E4,7),(5,E5,21)]:
 got=sum(all((a,b) in E for a,b in zip(p,p[1:])) for p in permutations(range(1,n+1)))
 if got!=want: raise ValueError((n,got,want))
 print(n,got)
