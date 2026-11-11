"""a8_thm448b.py (AUDIT) - CD_5 versus the THM-448 tower T_31: triple-intersection spectra (isomorphism invariant)."""
import sys, os, collections, itertools
import a8_thm448 as m
out = []
def P(s):
    print(s); out.append(s)
def triple_spectrum(T):
    n = len(T)
    c = collections.Counter()
    for u, v, w in itertools.combinations(range(n), 3):
        c[len(T[u] & T[v] & T[w])] += 1
    return sorted(c.items())
for k in (4, 5):
    Tc = m.cd_tournament(k)
    Tt = m.core(m.tower_matrix(k))
    P("n = %d: triple out-intersection spectrum of CD_n: %s" % (k, triple_spectrum(Tc)))
    P("n = %d: triple out-intersection spectrum of the THM-448 tower: %s" % (k, triple_spectrum(Tt)))
    P("        (THM-448 (c) records for T_31: {0:28, 1:84, 2:252, 3:3192, 4:812, 5:84, 6:28, 7:15})")
open(os.path.join(m.HERE, "a8_thm448b.out"), "w").write("\n".join(out) + "\n")
