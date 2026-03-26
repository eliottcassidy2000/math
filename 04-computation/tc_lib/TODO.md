# TC Library — TODO

## Current State
- ✅ C library (tc.h, tc.c) compiles and works
- ✅ YCoCg-R reversible color transform (lossless verified)
- ✅ Per-row adaptive filtering (5 PNG filters)
- ✅ zlib-9 backend
- ✅ PPM input/output
- ✅ Benchmark mode
- ✅ 20% smaller than PNG on real photos (fair, same zlib backend)
- ✅ HTML demo (tc_magic.html) handles images, GIFs, video

## TODO: Multi-format support
- [ ] PNG input (use libpng or stb_image)
- [ ] JPEG input (use libjpeg or stb_image)
- [ ] GIF input (extract frames, compress each)
- [ ] Video input (extract keyframes, compress each + delta between frames)
- [ ] Output: define TC file format properly (currently minimal header)

## TODO: Video compression
- [ ] Frame extraction from MP4/WebM (use ffmpeg or libav)
- [ ] Delta frames: XOR with previous frame, then TC compress the delta
- [ ] Keyframe interval selection (adaptive based on frame difference)
- [ ] GOP structure (Group of Pictures: I-frame + P-frames)

## TODO: Performance
- [ ] SIMD acceleration for Paeth prediction (SSE2/NEON)
- [ ] Multi-threaded compression (each row can be filtered independently)
- [ ] Streaming API (compress row-by-row without buffering whole image)
- [ ] Memory-mapped I/O for large files

## TODO: Quality
- [ ] Benchmark on full Kodak set (24 images)
- [ ] Benchmark against JPEG-XL lossless, WebP lossless
- [ ] Fix any remaining roundtrip issues (kodim12 had issues in Python)
- [ ] Fuzz testing with AFL/libFuzzer

## TODO: Distribution
- [ ] Makefile / CMakeLists.txt
- [ ] pkg-config support
- [ ] Python bindings (ctypes or cffi)
- [ ] npm package (WebAssembly build)
- [ ] Homebrew formula
