from itertools import product

def negative(w,width):
    d=len(w); edges=[]
    for i in range(d):
        for k in range(1,width+d+1):
            j=(i+k)%d
            if k<=width and w[j]>>(k-1)&1: edges.append((j,i,-k))
            else: edges.append((i,j,k))
    v=[0]*d
    for _ in range(d):
        change=False
        for a,b,s in edges:
            if v[b]>v[a]+s:
                v[b]=v[a]+s
                change=True
        if not change:return False
    return True

def literal(w,width,limit=48):
    d=len(w)
    queue=[0]; visited={0}
    for i in queue:
        if i==-d:return True
        for j in range(-limit,limit+1):
            if j in visited:continue
            gap=abs(j-i)
            if i<j: arc=gap>width or not(w[j%d]>>(gap-1)&1)
            else: arc=gap<=width and w[i%d]>>(gap-1)&1
            if arc: visited.add(j);queue.append(j)
    return False

def main():
    assert 5**12<=1000**3<5**13
    for d,n in enumerate((7,59,490,3767),1):
        true=0
        for w in product(range(8),repeat=d):
            a=negative(w,3)
            assert a==literal(w,3),(w,a)
            true+=a
        assert true==n
        print(d,true)

if __name__=="__main__":main()
