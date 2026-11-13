from replay_switch import *
for k, expected in zip((4,5,6),((1416,1440),(3636,3732),(9624,9936))):
    vals=[]
    for w,B in (('2',2**k+3),('30',2**(k-1)+1)):
        D=subtable(rows(tuple(map(int,bin(B)[2:]+w))))
        H=[1]+list(map(sum,D[1:])); N=len(H)-1
        vals.append(sum(H[s]*H[N^s] for s in range(N+1)))
    if tuple(vals)!=expected or vals[0]-vals[1]!=12-4*3**(k-2):
        raise ValueError((k,vals))
