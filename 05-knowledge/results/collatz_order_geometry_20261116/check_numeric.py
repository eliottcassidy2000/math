from check_tree_orders import H,direct_tri

def A(n):
    g=[[0]*n for _ in range(n)]
    for a in range(1,n+1):
        for b in range(a+1,n+1):
            if a>1 and a%2 and (3*a+1)//2==b:g[a-1][b-1]=1
            else:g[b-1][a-1]=1
    return g

def reach(g,i):
    seen={i};todo=[i]
    while todo:
        j=todo.pop()
        for k,e in enumerate(g[j]):
            if e and k not in seen:seen.add(k);todo.append(k)
    return seen

for n in range(2,81):
    j=(n-2)//3
    g=A(n)
    assert direct_tri(g)==j*(j+1)//2+(n+1)//9
    reaches=[reach(g,i) for i in range(n)]
    nontriv={i+1 for i in range(n) if any(k!=i and i in reaches[k] for k in reaches[i])}
    assert nontriv==set(range(3,3*j+3)) if j else not nontriv
counts={n:H(A(n)) for n in range(2,17)}
for n in range(3,17):
    assert counts[n]>counts[n-1] if n>=5 and n%3==2 else counts[n]==counts[n-1]
print('triangle and SCC formulas through 80; literal Hamilton counts:',counts)
