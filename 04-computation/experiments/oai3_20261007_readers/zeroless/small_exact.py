# Exact check for n <= 956 (where the 288-digit window may exceed the length of 2^n) and the A007377 list.
zl = [n for n in range(0, 957) if '0' not in str(2**n)]
print('zeroless 2^n, 0<=n<=956:', zl)
print('count', len(zl), ' max', max(zl))
A007377 = [0,1,2,3,4,5,6,7,8,9,13,14,15,16,18,19,24,25,27,28,31,32,33,34,35,36,37,39,49,51,67,72,76,77,81,86]
print('matches A007377:', zl == A007377)
# near misses: exactly one zero digit, n > 86
one = [(n, len(str(2**n))) for n in range(87, 957) if str(2**n).count('0') == 1]
print('n in [87,956] with exactly one 0 digit:', one)
print('2^86 =', 2**86)
