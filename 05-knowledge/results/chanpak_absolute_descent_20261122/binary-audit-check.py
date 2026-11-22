def flipped(rows, i, j):
    out = list(rows)
    out[i] ^= 1 << j
    out[j] ^= 1 << i
    return tuple(out)

def star_rows(core, tags, x, y):
    """Ordinal left flank, marked core, right flank, and last vertex as apex."""
    left = bin(x)[2:] if x else ''
    right = bin(y)[2:][::-1].translate(str.maketrans('01', '10')) if y else ''
    tape = left + ''.join(str((tags >> i) & 1) for i in range(len(core))) + right
    n = len(tape)
    mask = (1 << n) - 1
    rows = [mask ^ ((1 << (i + 1)) - 1) for i in range(n)]
    width = len(core)
    core_mask = ((1 << width) - 1) << len(left)
    for i, row in enumerate(core):
        index = len(left) + i
        rows[index] = (rows[index] & ~core_mask) | (row << len(left))
    apex = 0
    for i, bit in enumerate(tape):
        if bit == '1':
            apex |= 1 << i
        else:
            rows[i] |= 1 << n
    return tuple(rows + [apex])

def subpaths(rows, starts=None):
    n = len(rows)
    size = 1 << n
    end = [[0] * n for _ in range(size)]
    for i in range(n):
        if starts is None or starts >> i & 1:
            end[1 << i][i] = 1
    for subset in range(1, size):
        for i, count in enumerate(end[subset]):
            if count:
                choices = rows[i] & ~subset
                while choices:
                    bit = choices & -choices
                    end[subset | bit][bit.bit_length() - 1] += count
                    choices ^= bit
    return [sum(row) for row in end]

def profile(core, tags):
    n = len(core)
    full = (1 << n) - 1
    paths = subpaths(core)
    paths[0] = 1
    h = paths[-1]
    q = sum(paths[s] * paths[full ^ s] for s in range(full + 1))

    def eo(graph, marked):
        rooted = subpaths(graph, marked)
        return tuple(sum(rooted[s] * paths[full ^ s] for s in range(1, full + 1)
                         if s.bit_count() % 2 == parity) for parity in (0, 1))

    reverse = tuple(row ^ full ^ (1 << i) for i, row in enumerate(core))
    k = subpaths(star_rows(core, tags, 0, 0))[-1]
    return h, q, eo(core, tags), eo(reverse, full ^ tags), k

def strong(rows):
    target=(1<<len(rows))-1
    for u in range(len(rows)):
        mask=1<<u
        while True:
            new=mask
            for i,r in enumerate(rows):
                if mask>>i&1:new|=r
            if new==mask:break
            mask=new
        if mask!=target:return False
    return True

from itertools import product,combinations
P=(22,44,88,49,98,69,11)
B=flipped(P,0,1)
assert strong(P) and strong(B)
assert profile(P,1)==(189,2016,(117,144),(702,864),423)
def add(rows,S):
 n=len(rows)
 return tuple(r|((1-((S>>i)&1))<<n) for i,r in enumerate(rows))+(S,)
def H(rows): return subpaths(rows)[-1]
def singleton(i,t): return (1<<i) if t==0 else 127^(1<<i)
pairs=[]
for t,u in product((0,1),repeat=2):
 for v in (0,1,3):
  for e in (0,128):pairs.append(H(add(add(P,singleton(0,t)),singleton(v,u)|e)))
assert len(pairs)==24 and min(pairs)==855
A=[(T,profile(B,T)) for T in range(1,127)]
assert {x[-1] for _,x in A if x[-1]<=379}=={315,353,373,375,379}
assert min(2*min(x[2][1],x[3][1]) for _,x in A if x[-1]<=379)==222
all_data=[]
for e in combinations(range(7),2):
 if e==(0,1):continue
 C=flipped(B,*e)
 assert strong(C)
 for t in range(1,127):
  h,q,eo,eop,k=profile(C,t)
  all_data.append((h,q,k,k-h,2*min(eo[1],eop[1])))
assert tuple(map(min,zip(*all_data)))==(111,1428,189,78,144)
assert all(377%x[0] for x in all_data)
Z=flipped(flipped(B,0,2),0,3)
assert H(Z)==129 and H(add(Z,59))==377
assert H(add(B,8))==379
repair=((753,flipped(B,0,2),32,177),(755,B,1,247),(757,B,4,191),(759,B,4,239))
for target,C,s,t in repair:
 D=add(add(C,s),t)
 assert H(D)==target
 print(target,'strong',strong(D))
print('Paley cut controls, repaired host bounds, and all six target witnesses replayed.')
