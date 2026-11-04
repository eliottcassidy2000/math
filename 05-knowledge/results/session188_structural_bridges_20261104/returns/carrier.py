"""Literal inflation of the S185 address24 mixed homogeneous pair.

Source rows use original names0,...,6. In output rows, clone0 is old3,
clone(m-1) is old0 for m>=2; intermediate clones are successive fresh births.
The five exterior vertices keep original names1,2,4,5,6.
"""
base=(2,20,9,3,45,79,31)
outer=[1,2,4,5,6]
v=3

def physical_names(m):
 if m<1:raise ValueError('positive module size required')
 if m==1:return ['contracted_pair_proxy']+['old'+str(u) for u in outer]
 return ['old3']+['fresh'+str(i) for i in range(1,m-1)]+['old0']+['old'+str(u) for u in outer]

def inflate(m):
 if m<1:raise ValueError('positive module size required')
 rows=[sum(1<<j for j in range(m) if j>i)+sum(1<<(m+j) for j,u in enumerate(outer) if base[v]>>u&1) for i in range(m)]
 rows += [sum(1<<(m+j) for j,w in enumerate(outer) if base[u]>>w&1)+(((1<<m)-1) if base[u]>>v&1 else 0) for u in outer]
 return rows,{m+j for j,u in enumerate(outer) if 42>>u&1}
