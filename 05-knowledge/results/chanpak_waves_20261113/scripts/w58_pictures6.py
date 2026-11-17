"""Wave 58: prime-slot matching pictures at order 6 (the 15 prime 6-tournaments of THM-1960): distinct count forms."""
import io, contextlib
src = open('w58_pictures.py').read().replace("if n in (3, 5) and not modules(n, out)", "if n == 6 and not modules(n, out)").replace("for _ in range(12):", "for _ in range(3):")
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(src)
lines = buf.getvalue().splitlines()
forms = sorted({l.split('->')[1].split('(checked')[0].strip() for l in lines if '->' in l})
print(lines[0].replace('order 3 and 5', 'order 6'))
print("slot types:", sum(1 for l in lines if '->' in l), "| distinct count forms:", len(forms), "| all checked on random cores:", all('checked: True' in l for l in lines if '->' in l))
print("forms with no chi(3) term (no third free run):", [f for f in forms if 'chi(3)' not in f])
print("largest chi(3) coefficient:", max(int(f.split('*chi(3)')[0].split()[-1]) for f in forms if 'chi(3)' in f))
