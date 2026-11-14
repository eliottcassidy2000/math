from replay_source_switch import state,packing,physical_rows,literal,prime,matrix

PAIRS=[
 (26235,(1,0,1,0,2,1),(165,330),(1,0,0,2,3),(53,424)),
 (1031563,(1,0,1,0,4,0,1,0,1),(5701,13636),(4,0,0,2,0,1,0,0,1),(5589,13730)),
 (3338481,(1,0,3,0,0,1,0,1,0,0,0,1),(18687,43932),(3,0,0,1,0,3,0,0,0,0,1),(18575,44026)),
 (6491169,(4,0,0,3,0,1,0,0,0,0,1),(35935,85754),(1,0,3,0,0,3,0,0,0,0,0,1),(35823,85848)),
]
assert literal((1,0,1,0,2))==165 and literal((1,0,0,2))==53
assert literal(PAIRS[0][1])==495 and literal(PAIRS[0][3])==477
for U,A,z,B,t in PAIRS:
    assert state(A)==z and state(B)==t
    assert tuple(x-y for x,y in zip(z,t))==(112,-94)
    assert 47*z[0]+56*z[1]==47*t[0]+56*t[1]==U
    for q,d in ((10,-2),(13,0),(14,2)):
        for r in (1,2,3,5,17):
            assert packing(A+(0,)*q+(r,))-packing(B+(0,)*q+(r,))==d
    print('equal-coordinate baseline',U,'physical prefix orders',len(physical_rows(A)),len(physical_rows(B)))
_,A,_,B,_=PAIRS[1]
N=packing(A+(0,)*10+(1,))
assert N==26360017 and prime(N) and packing(B+(0,)*10+(1,))==N+2
print('new ordinary prime witness',N,'at order',len(physical_rows(A+(0,)*10+(1,))))
