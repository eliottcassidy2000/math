"""Compare the auditor's pair censuses (pairsdp_6_40.out by DP; pairsx_*.out exhaustive) with the session's
pairs_check.out / pairs_check_insample.out, and test the three-family gap law of HYP-9590 (A) / THM-5020 (6).
usage: python compare_pairs.py
"""
import re
import ast

def parse(path, tag_re=r"n\s*=\s*(\d+),?\s*k\s*=\s*(\d+)"):
    out = {}
    for line in open(path):
        m = re.search(tag_re, line)
        if not m:
            continue
        n = int(m.group(1))
        bl = re.search(r"black (?:pairs )?missing (\[[^\]]*\])", line)
        bu = re.search(r"blue \(non-tight\) (?:pairs )?missing (\[[^\]]*\])", line)
        adm = re.search(r"admissible words (\d+)", line)
        pal = re.search(r"(?:admissible palindromes|palindromic) (\d+)", line)
        rec = {'black': set(ast.literal_eval(bl.group(1))), 'blue': set(ast.literal_eval(bu.group(1)))}
        if adm:
            rec['adm'] = int(adm.group(1))
        if pal:
            rec['pal'] = int(pal.group(1))
        allb = re.search(r"ALL words: black (\[[^\]]*\]), blue (\[[^\]]*\])", line)
        if allb:
            rec['all_black'] = set(ast.literal_eval(allb.group(1)))
            rec['all_blue'] = set(ast.literal_eval(allb.group(2)))
        out[n] = rec
    return out

def families(n):
    k = n // 3
    black = {(k, k)} if n % 3 == 0 else set()
    blue = set()
    if n % 3 in (0, 1):
        blue.add((k, k))
    if n % 2 == 0:
        blue |= {(i, j) for i in range(1, k + 1) for j in range(i, k + 1) if i + j > k and i % 2 and j % 2}
    return black, blue

mine = parse('pairsdp_6_40.out')
sess = parse('../pairs_check.out')
sess.update(parse('../pairs_check_insample.out'))
ex = {}
for f in ('pairsx_6_21.out', 'pairsx_22_28.out'):
    ex.update(parse(f))

agree = all(mine[n]['black'] == sess[n]['black'] and mine[n]['blue'] == sess[n]['blue'] for n in range(6, 41))
counts = all(mine[n]['adm'] == sess[n]['adm'] and mine[n]['pal'] == sess[n]['pal'] for n in range(6, 41))
print(f"auditor DP census = session census (missing pairs), n = 6..40: {agree}")
print(f"admissible-word and admissible-palindrome counts identical, n = 6..40: {counts}")
exag = all(ex[n]['black'] == mine[n]['black'] and ex[n]['blue'] == mine[n]['blue'] for n in range(6, 29))
print(f"exhaustive census = DP census, n = 6..28: {exag}")
allsame = all(ex[n]['all_black'] == ex[n]['black'] and ex[n]['all_blue'] == ex[n]['blue'] for n in range(6, 29))
print(f"pairs missing for admissible words = pairs missing for ALL words (admissibility costs nothing), n = 6..28: {allsame}")
law = []
for n in range(6, 41):
    fb, fu = families(n)
    if mine[n]['black'] != fb or mine[n]['blue'] != fu:
        law.append(n)
print(f"gap law (exactly the three families) fails at n in {law} (range 6..40)")
for n in law:
    fb, fu = families(n)
    print(f"  n = {n}: missing black {sorted(mine[n]['black'])} vs families {sorted(fb)}; missing blue "
          f"{sorted(mine[n]['blue'])} vs families {sorted(fu)}")
