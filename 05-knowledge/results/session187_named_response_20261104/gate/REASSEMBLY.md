# Complete gate data reassembly

The complete 2,149,227-byte named-response JSON is stored losslessly in gzip/base64 chunks to keep publication requests bounded. No coefficient has been omitted. Run `python reassemble.py` in this directory before running the unchanged author or audit verifier. The script checks every chunk and the exact original SHA-256, then writes `full_response.json`. All tools used are Python standard library.
