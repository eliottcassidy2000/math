from replay_busch_inverse import data, need
M=[238,220,216,16,225,206,8,72]
A=[[(r>>j)&1 for j in range(8)] for r in M]
pairs=[(i,j) for i in range(8) for j in range(i+1,8) if not ((M[i]^M[j])&~((1<<i)|(1<<j)))]
need(pairs==[(0,5),(1,2),(6,7)])
counts=[]
for _,b in pairs:
    I=[i for i in range(8) if i!=b]
    counts.append(data([[A[i][j] for j in I] for i in I],0)[0])
need(counts==[25,33,37])
counts=[]
for v in range(8):
    B=[r+[r[v]] for r in A]+[A[v][:]+[0]]
    B[v][-1]=1; B[-1][v]=0
    counts.append(data(B,0)[0])
need(counts==[139,123,123,165,539,139,115,115])
print("63 rail boundary: all three twin contractions and eight duplications checked")
