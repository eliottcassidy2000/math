from collections import deque

def orbit(m):
    toda = deque([1 % m])
    seen = set(toda)
    while toda:
        x=toda.popleft()
        for y in (2*x % m, (3*x+3) % m):
            if y not in seen:
                seen.add(y)
                toda.append(y)
    return seen

for p in (2,3,5,7,11,13,17,19):
    q=p
    while q<100000:
        if len(orbit(q))!=q:
            raise RuntimeError(q)
        q*=p
if orbit(6)!={0,1,2,3,4}:
    raise RuntimeError(orbit(6))
print("Prime-power orbit controls to 100000; mixed modulus 6 omits 5.")
