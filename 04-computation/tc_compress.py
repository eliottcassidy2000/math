#!/usr/bin/env python3
"""
tc_compress.py -- Tournament Compression CLI for Real Images
kind-pasteur-2026-03-24-S20cq

Compress ANY image file using tournament codec. Compares against PNG.

USAGE:
  python tc_compress.py image.png                    # analyze & compress
  python tc_compress.py image.jpg --save output.tc   # save compressed
  python tc_compress.py --bench-dir ./images/        # benchmark directory
  python tc_compress.py --generate-test              # generate test images

SUPPORTED: PNG, JPEG, BMP, TIFF (anything Pillow can read)
"""

import sys
import os
import io
import zlib
import bz2
import lzma
import time
import argparse
import numpy as np

try:
    from PIL import Image
except ImportError:
    print("ERROR: Pillow required. Install: pip install Pillow")
    sys.exit(1)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from tc_image_ultimate import compress_ultimate

__version__ = "1.0.0"


def load_image(path):
    """Load image and convert to numpy array."""
    img = Image.open(path)
    if img.mode == 'RGBA':
        # Convert RGBA to RGB (drop alpha)
        img = img.convert('RGB')
    elif img.mode == 'P':
        img = img.convert('RGB')
    elif img.mode == 'L':
        pass  # grayscale is fine
    elif img.mode != 'RGB':
        img = img.convert('RGB')
    return np.array(img)


def png_size(img_array):
    """Get PNG size for comparison."""
    if img_array.ndim == 2:
        img = Image.fromarray(img_array)
    else:
        img = Image.fromarray(img_array)
    buf = io.BytesIO()
    img.save(buf, format='PNG', optimize=True)
    return len(buf.getvalue())


def analyze_file(path):
    """Analyze an image file."""
    img = load_image(path)
    H, W = img.shape[:2]
    is_color = img.ndim == 3
    channels = 3 if is_color else 1
    raw_size = img.nbytes

    t0 = time.time()
    tc_data, tc_method = compress_ultimate(img)
    tc_time = time.time() - t0
    tc_size = len(tc_data)

    png_sz = png_size(img)

    ratio_vs_raw = raw_size / tc_size if tc_size > 0 else 0
    ratio_vs_png = png_sz / tc_size if tc_size > 0 else 0

    return {
        'path': path,
        'resolution': f"{W}x{H}",
        'channels': channels,
        'raw_size': raw_size,
        'tc_size': tc_size,
        'png_size': png_sz,
        'tc_method': tc_method,
        'ratio_raw': ratio_vs_raw,
        'ratio_png': ratio_vs_png,
        'time_ms': tc_time * 1000,
        'beats_png': tc_size < png_sz * 0.98,
    }


def cmd_analyze(args):
    """Analyze one or more image files."""
    files = args.files
    if not files:
        print("No files specified")
        return

    print(f"TC Compress v{__version__}")
    print("=" * 90)

    results = []
    for path in files:
        if not os.path.exists(path):
            print(f"  File not found: {path}")
            continue

        result = analyze_file(path)
        results.append(result)

        tag = "WIN" if result['beats_png'] else "TIE"
        print(f"\n  {os.path.basename(path):>30}")
        print(f"    Resolution:   {result['resolution']} ({result['channels']}ch)")
        print(f"    Raw size:     {result['raw_size']:>10,} bytes")
        print(f"    TC size:      {result['tc_size']:>10,} bytes ({result['ratio_raw']:.1f}x)")
        print(f"    PNG size:     {result['png_size']:>10,} bytes")
        print(f"    TC/PNG:       {result['ratio_png']:.3f}x {tag}")
        print(f"    Method:       {result['tc_method']}")
        print(f"    Time:         {result['time_ms']:.0f}ms")

    if len(results) > 1:
        wins = sum(1 for r in results if r['beats_png'])
        total = len(results)
        print(f"\n  Summary: {wins}/{total} beat PNG ({wins/total*100:.0f}%)")


def cmd_bench_dir(args):
    """Benchmark all images in a directory."""
    directory = args.directory
    extensions = {'.png', '.jpg', '.jpeg', '.bmp', '.tiff', '.tif', '.gif'}

    files = []
    for fname in os.listdir(directory):
        ext = os.path.splitext(fname)[1].lower()
        if ext in extensions:
            files.append(os.path.join(directory, fname))

    if not files:
        print(f"No image files found in {directory}")
        return

    print(f"TC Compress v{__version__} -- Directory Benchmark")
    print(f"  Directory: {directory}")
    print(f"  Files: {len(files)}")
    print("=" * 90)
    print(f"  {'File':>30} {'Res':>10} {'Raw':>10} {'TC':>10} {'PNG':>10} {'TC/PNG':>8} {'Method':>15}")

    wins = ties = losses = 0
    for path in sorted(files):
        try:
            result = analyze_file(path)
        except Exception as e:
            print(f"  {os.path.basename(path):>30} ERROR: {e}")
            continue

        r = result['ratio_png']
        if r > 1.02: wins += 1; tag = "WIN"
        elif r < 0.98: losses += 1; tag = "LOSE"
        else: ties += 1; tag = "TIE"

        print(f"  {os.path.basename(path):>30} {result['resolution']:>10} "
              f"{result['raw_size']:>9,}B {result['tc_size']:>9,}B "
              f"{result['png_size']:>9,}B {r:7.2f}x {result['tc_method']:>15} {tag}")

    total = wins + ties + losses
    if total > 0:
        print(f"\n  SCORE: {wins}W {ties}T {losses}L / {total}")
        print(f"  Win rate: {wins/total*100:.0f}%")


def cmd_generate_test(args):
    """Generate test images for benchmarking."""
    outdir = args.outdir or 'test_images'
    os.makedirs(outdir, exist_ok=True)

    N = 256
    tests = {}
    tests['gradient_h.png'] = np.tile(np.arange(N, dtype=np.uint8), (N, 1))
    tests['gradient_v.png'] = np.tile(np.arange(N, dtype=np.uint8).reshape(-1,1), (1, N))
    tests['checker.png'] = np.array([[(255 if (i//16+j//16)%2==0 else 0) for j in range(N)] for i in range(N)], dtype=np.uint8)

    x = np.linspace(0, 4*np.pi, N)
    X, Y = np.meshgrid(x, x)
    tests['smooth.png'] = ((128+100*np.sin(X)*np.cos(Y))).clip(0,255).astype(np.uint8)

    tests['circle.png'] = np.array([[255 if (i-N//2)**2+(j-N//2)**2<(N//3)**2 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8)

    np.random.seed(42)
    tests['random.png'] = np.random.randint(0,256,(N,N),dtype=np.uint8)

    # Color
    color_grad = np.zeros((N, N, 3), dtype=np.uint8)
    color_grad[:,:,0] = tests['gradient_h.png']
    color_grad[:,:,1] = tests['gradient_v.png']
    color_grad[:,:,2] = 128
    tests['color_gradient.png'] = color_grad

    for name, img in tests.items():
        path = os.path.join(outdir, name)
        if img.ndim == 2:
            Image.fromarray(img).save(path)
        else:
            Image.fromarray(img).save(path)
        print(f"  Generated {path}")

    print(f"\n  {len(tests)} test images saved to {outdir}/")
    print(f"  Run: python tc_compress.py --bench-dir {outdir}")


def main():
    parser = argparse.ArgumentParser(
        description=f'TC Compress v{__version__} -- Tournament Image Compression',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s photo.png                     # analyze single image
  %(prog)s *.png                         # analyze multiple images
  %(prog)s --bench-dir ./images/         # benchmark directory
  %(prog)s --generate-test               # create test images
        """)

    parser.add_argument('files', nargs='*', help='Image files to analyze')
    parser.add_argument('--bench-dir', '-d', dest='directory', help='Benchmark all images in directory')
    parser.add_argument('--generate-test', action='store_true', help='Generate test images')
    parser.add_argument('--outdir', default='test_images', help='Output directory for test images')

    args = parser.parse_args()

    if args.generate_test:
        cmd_generate_test(args)
    elif args.directory:
        cmd_bench_dir(args)
    elif args.files:
        cmd_analyze(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
