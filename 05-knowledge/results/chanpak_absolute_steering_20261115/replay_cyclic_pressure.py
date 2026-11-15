from itertools import product

E=((1,1,-1),(1,1,0),(1,0,-1))
A=((1,2,0),(1,2,0),(1,0,0))
B=((1,1,-1),(2,2,0),(0,0,0))
C=((5,9,-1),(6,12,0),(2,0,0))
D=((5,6,-2),(9,12,0),(1,0,0))
for name,cyclic in [('EABC',C),('EABD',D)]:
    for s,d in product((1,4),(0,1)):
        initial=(1,s,d)
        rows=[(initial,-1)]
        for _ in range(8):
            rows=[(tuple(sum(a*b for a,b in zip(row,v)) for row in mat),i)
                  for v,last in rows for i,mat in enumerate((E,A,B,cyclic))
                  if (last,i) not in {(2,0),(1,1),(2,2)}]
        if len(rows)!=16897:
            raise AssertionError(len(rows))
        unit=(1+s)*10**9
        ub=sum((unit+w[0]+w[1]-1)//(w[0]+w[1]) for w,_ in rows)
        if 20*ub>=19*10**9:
            raise AssertionError(ub)
        print(name,(1,s,d),len(rows),ub,'over 1000000000')
if not 21**8*19**1000<20**1000:
    raise AssertionError("rational exponent validation")
