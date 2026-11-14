import sympy as sp
def step(s,m):
    H,E,L,h=s
    if m=='O': return (H,0,2*L,h)
    if m=='Z': return (H+L-E,h,2*L-h,h)
    return (H+L-E,H-E,H+L,H)
def blockmat(w):
    cols=[]
    for e in ((1,0,0),(0,1,0),(0,0,1)):
        s=step((e[0],e[1],e[2],0),'J')
        for ch in w: s=step(s,ch)
        cols.append(s[:3])
    return sp.Matrix([[cols[c][r] for c in range(3)] for r in range(3)])
W11=['','O','OO','OZ','OOO','OOZ','OZO','OZZ','OOOZ','OZOZ','OZZZ']
for w in W11:
    M=blockmat(w)
    print("J%-5s det=%4s  rows (H,E,lam):"%(w,M.det()), M.tolist())
