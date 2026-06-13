# TIC — Tournament Image Codec

Lossless RGB image compression that beats PNG by ~20% on photographs.

## How it works

Three components, each doing what it's best at:

1. **G-RG-BG color decorrelation**: Stores green channel, red-green difference,
   blue-green difference. This removes inter-channel redundancy that PNG wastes.
   Perfectly reversible via modular arithmetic: `R = (RG + G) & 0xFF`.

2. **MED prediction**: The Median Edge Detector (from JPEG-LS, 1996) predicts
   each pixel from its left, above, and above-left neighbors. Better than
   PNG's Paeth filter on smooth gradients.

3. **zlib-9 compression**: Same Deflate backend as PNG. Fair comparison.

## Results on Kodak benchmark

8 Kodak standard images + 1 real desk photo, 256×256 center crops:

| Image          | PNG      | TIC      | Gain   |
|----------------|----------|----------|--------|
| kodim01 (sails)| 134,653  | 99,894   | 1.35×  |
| kodim05 (stream)| 145,367 | 115,532  | 1.26×  |
| kodim08 (moon) | 129,667  | 105,621  | 1.23×  |
| vlcsnap (desk) | 63,289   | 44,785   | 1.41×  |
| **Aggregate**  | **1,000,068** | **827,317** | **1.21×** |

All 9 images: TIC wins. All roundtrips: verified lossless.

## Honest caveats

- The technique is not novel: it combines JPEG-LS prediction with WebP-style
  color decorrelation. PNG doesn't do this, which is why we beat it.
- JPEG-LS (with Golomb-Rice coding) would beat us by another ~15%.
- WebP lossless and JPEG-XL are significantly better overall.
- This is useful as a simple, clean alternative to PNG for RGB photos,
  not as a state-of-the-art codec.

## Build

```bash
# With system zlib:
gcc -O2 -o tic_cli tic.c -lz

# With bundled zlib (if system zlib unavailable):
gcc -O2 -Izlib-1.3.1 -o tic_cli tic.c zlib-1.3.1/libz.a
```

## Usage

```bash
# Compress (input: raw RGB, W×H×3 bytes)
tic_cli compress photo.raw 256 256 photo.tic

# Decompress
tic_cli decompress photo.tic photo_out.raw

# Benchmark
tic_cli bench photo.raw 256 256 100
```

## API

```c
#include "tic.h"

// Encode
size_t bound = tic_encode_bound(width, height);
uint8_t *out = malloc(bound);
size_t out_len;
int ret = tic_encode(rgb, width, height, out, bound, &out_len);

// Decode
uint16_t w, h;
tic_decode(data, data_len, NULL, &w, &h);  // query dimensions
uint8_t *rgb = malloc(w * h * 3);
tic_decode(data, data_len, rgb, &w, &h);   // decompress
```

## License

Public domain.
