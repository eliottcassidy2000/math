from itertools import product, permutations
cross=((1,3),(1,4),(2,3),(2,4))
fixed={(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)}
heights=[]
for word in product((0,1),repeat=4):
    edges=fixed|{(v,u) if b else (u,v) for (u,v),b in zip(cross,word)}
    h=sum(all((u,v) in edges for u,v in zip(p,p[1:])) for p in permutations(range(5)))
    heights.append(h)
assert set(heights)=={9,11,13,15}
assert min(heights)==9
print('all 16 one-root triangle completions certify sharp addition toll >=4')
