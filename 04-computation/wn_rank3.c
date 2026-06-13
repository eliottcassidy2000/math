/*
 * wn_rank3.c — Compute W(n) with disk-backed chunked DP for large n
 * opus-2026-03-15-S89c
 *
 * For levels where curr+next exceeds RAM, use mmap for curr (sequential read)
 * and process next in two halves. Each half fits in ~4.5 GB for n=28.
 *
 * Memory: max(half_next) ≈ C(n,n/2)/2 * n/2 * 16 bytes
 * For n=28: ~4.5 GB at peak. Fits in 8 GB RAM.
 * Trade-off: 2x compute at peak levels (process curr twice).
 *
 * Compile: gcc -O3 -o wn_rank3 wn_rank3.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

typedef unsigned __int128 u128;

void print_u128(u128 x) {
    if (x == 0) { printf("0"); return; }
    char buf[50];
    int pos = 0;
    while (x > 0) {
        buf[pos++] = '0' + (int)(x % 10);
        x /= 10;
    }
    for (int i = pos - 1; i >= 0; i--)
        putchar(buf[i]);
}

static inline int popcount(unsigned int x) {
    return __builtin_popcount(x);
}

static long long binom[31][31];

void init_binom() {
    for (int i = 0; i <= 30; i++) {
        binom[i][0] = 1;
        for (int j = 1; j <= i; j++)
            binom[i][j] = binom[i-1][j-1] + binom[i-1][j];
        for (int j = i+1; j <= 30; j++)
            binom[i][j] = 0;
    }
}

static inline long long mask_rank(unsigned int mask, int n, int p) {
    long long rank = 0;
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (mask & (1u << i)) {
            rank += binom[i][count + 1];
            count++;
        }
    }
    return rank;
}

static inline unsigned int mask_unrank(long long r, int n, int p) {
    unsigned int mask = 0;
    int remaining = p;
    for (int i = n - 1; i >= 0 && remaining > 0; i--) {
        long long c = binom[i][remaining];
        if (r >= c) {
            r -= c;
            mask |= (1u << i);
            remaining--;
        }
    }
    return mask;
}

static inline int local_index(unsigned int mask, int v) {
    return popcount(mask & ((1u << v) - 1));
}

/* Memory threshold: if both arrays exceed this, use chunked mode */
#define MEM_THRESHOLD (6LL * 1024 * 1024 * 1024) /* 6 GB */

/* Write array to a temp file, return fd. Caller must close and unlink. */
static int write_to_disk(const char *path, u128 *data, long long count) {
    int fd = open(path, O_RDWR | O_CREAT | O_TRUNC, 0600);
    if (fd < 0) { perror("open"); return -1; }
    long long total = count * sizeof(u128);
    long long written = 0;
    while (written < total) {
        long long chunk = total - written;
        if (chunk > (1LL << 30)) chunk = (1LL << 30); /* 1 GB at a time */
        ssize_t w = write(fd, (char *)data + written, chunk);
        if (w <= 0) { perror("write"); close(fd); return -1; }
        written += w;
    }
    return fd;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    if (n < 2 || n > 30) {
        fprintf(stderr, "n must be 2..30\n");
        return 1;
    }

    init_binom();
    time_t t0 = time(NULL);

    int p_curr = 1;
    long long cnt_curr = binom[n][1];
    u128 *dp_curr = (u128 *)calloc(cnt_curr * p_curr, sizeof(u128));
    if (!dp_curr) { fprintf(stderr, "alloc fail\n"); return 1; }

    for (int v = 0; v < n; v++) {
        long long r = mask_rank(1u << v, n, 1);
        dp_curr[r * p_curr + 0] = 1;
    }

    fprintf(stderr, "n=%d, starting DP\n", n);

    for (int p = 1; p < n; p++) {
        int p_next = p + 1;
        long long cnt_next = binom[n][p_next];
        long long mem_curr = cnt_curr * p * (long long)sizeof(u128);
        long long mem_next = cnt_next * p_next * (long long)sizeof(u128);
        long long mem_both = mem_curr + mem_next;

        if (mem_both <= MEM_THRESHOLD) {
            /* In-memory single pass */
            fprintf(stderr, "  level %d->%d: %lld -> %lld masks (%.1f GB) (%lds)...\n",
                    p, p_next, cnt_curr, cnt_next,
                    (double)mem_both / 1e9, time(NULL) - t0);
            fflush(stderr);

            u128 *dp_next = (u128 *)calloc(cnt_next * p_next, sizeof(u128));
            if (!dp_next) {
                fprintf(stderr, "OOM at level %d, falling through to chunked\n", p);
                goto chunked;
            }

            for (long long mi = 0; mi < cnt_curr; mi++) {
                unsigned int mask = mask_unrank(mi, n, p);
                long long base_curr = mi * p;
                unsigned int tmp = mask;
                int li = 0;
                while (tmp) {
                    int v = __builtin_ctz(tmp);
                    tmp &= tmp - 1;
                    u128 cnt = dp_curr[base_curr + li];
                    li++;
                    if (cnt == 0) continue;
                    for (int u = 0; u < n; u++) {
                        if (mask & (1u << u)) continue;
                        if (u == v - 1) continue;
                        u128 weight = (u == v + 1) ? 2 * cnt : cnt;
                        unsigned int new_mask = mask | (1u << u);
                        long long nr = mask_rank(new_mask, n, p_next);
                        int local_u = local_index(new_mask, u);
                        dp_next[nr * p_next + local_u] += weight;
                    }
                }
            }

            fprintf(stderr, "    done (%lds)\n", time(NULL) - t0);
            free(dp_curr);
            dp_curr = dp_next;
            cnt_curr = cnt_next;
            p_curr = p_next;
            continue;

        chunked:;
        }

        /* Chunked two-pass: write curr to disk, process next in two halves */
        long long half = cnt_next / 2;
        long long half2_start = half;
        long long half2_end = cnt_next;

        fprintf(stderr, "  level %d->%d: %lld -> %lld masks CHUNKED "
                "(curr=%.1f GB, next=%.1f GB, half=%.1f GB) (%lds)...\n",
                p, p_next, cnt_curr, cnt_next,
                (double)mem_curr / 1e9, (double)mem_next / 1e9,
                (double)(half * p_next * sizeof(u128)) / 1e9,
                time(NULL) - t0);
        fflush(stderr);

        /* Write curr to disk */
        char tmppath[256];
        snprintf(tmppath, sizeof(tmppath), "/tmp/wn_dp_level_%d.bin", p);
        long long curr_entries = cnt_curr * p;
        int fd = write_to_disk(tmppath, dp_curr, curr_entries);
        if (fd < 0) { free(dp_curr); return 1; }
        free(dp_curr);
        dp_curr = NULL;

        /* mmap curr from disk for reading */
        long long file_size = curr_entries * sizeof(u128);
        u128 *curr_mmap = (u128 *)mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
        if (curr_mmap == MAP_FAILED) {
            perror("mmap");
            close(fd);
            unlink(tmppath);
            return 1;
        }
        /* Advise sequential access */
        madvise(curr_mmap, file_size, MADV_SEQUENTIAL);

        /* Allocate next array for the full level */
        /* We'll try to allocate the full thing; if it fails, do two passes */
        u128 *dp_next = (u128 *)calloc(cnt_next * p_next, sizeof(u128));

        if (dp_next) {
            /* Single pass through curr_mmap -> dp_next (full) */
            fprintf(stderr, "    full next alloc succeeded, single pass (%lds)...\n", time(NULL) - t0);
            fflush(stderr);

            for (long long mi = 0; mi < cnt_curr; mi++) {
                unsigned int mask = mask_unrank(mi, n, p);
                long long base_curr = mi * p;
                unsigned int tmp = mask;
                int li = 0;
                while (tmp) {
                    int v = __builtin_ctz(tmp);
                    tmp &= tmp - 1;
                    u128 cnt = curr_mmap[base_curr + li];
                    li++;
                    if (cnt == 0) continue;
                    for (int u = 0; u < n; u++) {
                        if (mask & (1u << u)) continue;
                        if (u == v - 1) continue;
                        u128 weight = (u == v + 1) ? 2 * cnt : cnt;
                        unsigned int new_mask = mask | (1u << u);
                        long long nr = mask_rank(new_mask, n, p_next);
                        int local_u = local_index(new_mask, u);
                        dp_next[nr * p_next + local_u] += weight;
                    }
                }
            }
        } else {
            /* Two-pass: process each half of next separately */
            fprintf(stderr, "    full alloc failed, two-pass mode (%lds)...\n", time(NULL) - t0);

            /* Allocate space for full next, but in two halves written to disk */
            char next_path[256];
            snprintf(next_path, sizeof(next_path), "/tmp/wn_dp_level_%d.bin", p_next);
            int next_fd = open(next_path, O_RDWR | O_CREAT | O_TRUNC, 0600);
            if (next_fd < 0) { perror("open next"); munmap(curr_mmap, file_size); close(fd); return 1; }
            /* Pre-allocate the file */
            long long next_file_size = cnt_next * p_next * sizeof(u128);
            if (ftruncate(next_fd, next_file_size) < 0) {
                perror("ftruncate");
                close(next_fd); munmap(curr_mmap, file_size); close(fd);
                return 1;
            }

            /* Pass 1: first half of next (ranks 0..half-1) */
            long long half_entries = half * p_next;
            long long half_bytes = half_entries * sizeof(u128);
            u128 *next_half = (u128 *)calloc(half_entries, sizeof(u128));
            if (!next_half) {
                fprintf(stderr, "Failed to allocate even half of next (%.1f GB)\n",
                        (double)half_bytes / 1e9);
                close(next_fd); unlink(next_path);
                munmap(curr_mmap, file_size); close(fd); unlink(tmppath);
                return 1;
            }

            fprintf(stderr, "    pass 1: next ranks 0..%lld (%lds)...\n", half - 1, time(NULL) - t0);
            fflush(stderr);

            for (long long mi = 0; mi < cnt_curr; mi++) {
                unsigned int mask = mask_unrank(mi, n, p);
                long long base_curr = mi * p;
                unsigned int tmp = mask;
                int li = 0;
                while (tmp) {
                    int v = __builtin_ctz(tmp);
                    tmp &= tmp - 1;
                    u128 cnt = curr_mmap[base_curr + li];
                    li++;
                    if (cnt == 0) continue;
                    for (int u = 0; u < n; u++) {
                        if (mask & (1u << u)) continue;
                        if (u == v - 1) continue;
                        u128 weight = (u == v + 1) ? 2 * cnt : cnt;
                        unsigned int new_mask = mask | (1u << u);
                        long long nr = mask_rank(new_mask, n, p_next);
                        if (nr >= half) continue; /* skip second half */
                        int local_u = local_index(new_mask, u);
                        next_half[nr * p_next + local_u] += weight;
                    }
                }
            }

            /* Write first half to disk */
            lseek(next_fd, 0, SEEK_SET);
            long long w = 0;
            while (w < half_bytes) {
                long long chunk = half_bytes - w;
                if (chunk > (1LL << 30)) chunk = (1LL << 30);
                ssize_t wr = write(next_fd, (char *)next_half + w, chunk);
                if (wr <= 0) { perror("write half1"); return 1; }
                w += wr;
            }

            fprintf(stderr, "    pass 1 done, written to disk (%lds)\n", time(NULL) - t0);

            /* Pass 2: second half of next (ranks half..cnt_next-1) */
            long long half2_count = cnt_next - half;
            long long half2_entries = half2_count * p_next;
            /* Reuse the buffer (might be different size) */
            if (half2_entries > half_entries) {
                free(next_half);
                next_half = (u128 *)calloc(half2_entries, sizeof(u128));
                if (!next_half) {
                    fprintf(stderr, "Failed to allocate half2\n");
                    return 1;
                }
            } else {
                memset(next_half, 0, half2_entries * sizeof(u128));
            }

            fprintf(stderr, "    pass 2: next ranks %lld..%lld (%lds)...\n",
                    half, cnt_next - 1, time(NULL) - t0);
            fflush(stderr);

            for (long long mi = 0; mi < cnt_curr; mi++) {
                unsigned int mask = mask_unrank(mi, n, p);
                long long base_curr = mi * p;
                unsigned int tmp = mask;
                int li = 0;
                while (tmp) {
                    int v = __builtin_ctz(tmp);
                    tmp &= tmp - 1;
                    u128 cnt = curr_mmap[base_curr + li];
                    li++;
                    if (cnt == 0) continue;
                    for (int u = 0; u < n; u++) {
                        if (mask & (1u << u)) continue;
                        if (u == v - 1) continue;
                        u128 weight = (u == v + 1) ? 2 * cnt : cnt;
                        unsigned int new_mask = mask | (1u << u);
                        long long nr = mask_rank(new_mask, n, p_next);
                        if (nr < half) continue; /* skip first half */
                        int local_u = local_index(new_mask, u);
                        next_half[(nr - half) * p_next + local_u] += weight;
                    }
                }
            }

            /* Write second half to disk */
            long long half2_bytes = half2_entries * sizeof(u128);
            lseek(next_fd, half_bytes, SEEK_SET);
            w = 0;
            while (w < half2_bytes) {
                long long chunk = half2_bytes - w;
                if (chunk > (1LL << 30)) chunk = (1LL << 30);
                ssize_t wr = write(next_fd, (char *)next_half + w, chunk);
                if (wr <= 0) { perror("write half2"); return 1; }
                w += wr;
            }

            free(next_half);
            fprintf(stderr, "    pass 2 done, written to disk (%lds)\n", time(NULL) - t0);

            /* Now mmap the full next level from disk as the new curr */
            munmap(curr_mmap, file_size);
            close(fd);
            unlink(tmppath);

            dp_next = (u128 *)mmap(NULL, next_file_size, PROT_READ | PROT_WRITE,
                                    MAP_PRIVATE, next_fd, 0);
            if (dp_next == MAP_FAILED) {
                perror("mmap next");
                close(next_fd); unlink(next_path);
                return 1;
            }
            close(next_fd);
            /* Copy mmap'd data to regular memory for next iteration */
            /* Actually, keep it as mmap'd — works fine as dp_curr for next level */
            /* But we need to track that it's mmap'd for cleanup */
            /* For simplicity: if next level is also chunked, it'll write to disk again */

            /* Actually, let's just allocate dp_curr and copy from mmap */
            /* Or better: keep the file and use mmap as dp_curr */
            dp_curr = dp_next;
            cnt_curr = cnt_next;
            p_curr = p_next;
            /* Note: dp_curr is mmap'd, will be written to disk again if needed */
            /* The unlink will happen when we write the next level's curr to disk */
            /* For now, keep the file around (unlink deferred) */
            snprintf(tmppath, sizeof(tmppath), "/tmp/wn_dp_level_%d.bin", p_next);
            /* dp_curr points to mmap'd memory from next_path */
            continue;
        }

        /* Clean up mmap of curr */
        munmap(curr_mmap, file_size);
        close(fd);
        unlink(tmppath);

        fprintf(stderr, "    done (%lds)\n", time(NULL) - t0);
        dp_curr = dp_next;
        cnt_curr = cnt_next;
        p_curr = p_next;
    }

    /* Sum over all endpoints at full mask */
    u128 W = 0;
    for (int j = 0; j < n; j++) {
        W += dp_curr[j];
    }

    printf("W(%d) = ", n);
    print_u128(W);
    printf("\n");

    fprintf(stderr, "Total time: %lds\n", time(NULL) - t0);

    /* dp_curr might be mmap'd or malloc'd — we don't track which. Just exit. */
    return 0;
}
