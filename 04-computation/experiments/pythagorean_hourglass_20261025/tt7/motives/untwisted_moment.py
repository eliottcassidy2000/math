# Untwisted fourth Legendre moment sum_{l != 0,1} a_l^4 against span{q^i, q^i a(8.4.a.a), a(8.6.a.a), a(4.6.a.a)}.
src = open('legendre_moments.py').read()
exec(src.split('rows = []')[0])
exec('def in_span' + src.split('def in_span')[1].split('cols = {')[0])
exec('cols = {' + src.split('cols = {')[1].split('for label, names')[0])
def target(q, M): return M['u4']
rows2 = []
for q in primerange(5, P):
    chi = np.full(q, -1, dtype=np.int64); chi[0] = 0; chi[[(x * x) % q for x in range(1, q)]] = 1
    t = np.arange(q); base_ = chi[(t * (t - 1)) % q]
    a = np.array([-int((base_ * chi[(t - l) % q]).sum()) for l in range(2, q)])
    rows2.append((q, {'u4': int((a ** 4).sum()), 'u2': int((a ** 2).sum())}))
rows = rows2
for lab, names in [('poly + a84', ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84']), ('+ a46', ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84', 'a46']),
                   ('+ a86 + a46', ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84', 'a86', 'a46']), ('+ chi_-4 polys', ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84', 'a86', 'a46', 'chi-4', 'q chi-4', 'q^2 chi-4'])]:
    s = in_span(cols, names, target)
    print(f"untwisted sum a_l^4 in span({lab}) over {len(rows)} primes:", ('EXACT ' + str({k: str(v) for k, v in s.items() if v})) if s else 'NOT')
s2 = in_span(cols, ['1', 'q', 'q^2', 'chi-4', 'q chi-4'], lambda q, M: M['u2'])
print("untwisted sum a_l^2:", ('EXACT ' + str({k: str(v) for k, v in s2.items() if v})) if s2 else 'NOT')
