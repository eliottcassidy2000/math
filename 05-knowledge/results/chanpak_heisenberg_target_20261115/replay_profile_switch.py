from heisenberg import dp, literal
for w,h,q in [('ABBBABB',61,1744),('ABABBAB',61,1864)]:
    got=dp(literal(w))
    if got!=(h,q):raise AssertionError((w,got))
    print(w,got)
for j in range(1,6):
    m=5**j
    r=1864*pow(1744,-1,m)%m
    u={pow(r,i,m) for i in range(m//5)}
    if len(u)!=m//5 or {i%5 for i in u}!={1}:raise AssertionError(j)
    print('fixed-h coset',m,len(u))
