#!/usr/bin/env python3
"""Extract plain text from a (downloaded, untrusted) PDF for literature checks (S43 audit).
Usage: python -I pdf_to_text.py <pdf> <out.txt>
"""
import sys
import fitz  # PyMuPDF

doc = fitz.open(sys.argv[1])
with open(sys.argv[2], 'w', encoding='utf-8') as f:
    for i, page in enumerate(doc):
        f.write('\n=== PAGE %d ===\n' % (i + 1))
        f.write(page.get_text())
print('pages', len(doc))
