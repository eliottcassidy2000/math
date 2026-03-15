/*
 * wn_rank4.c — Compute W(n) with multi-pass disk-backed DP
 * opus-2026-03-15-S89c
 *
 * Key improvement over wn_rank3: adaptive number of passes.
 * Each pass processes a chunk of next-level ranks.
 * Chunk size chosen to fit in MAX_CHUNK_GB of RAM.
 * Curr level is always mmap'd from disk with MADV_SEQUENTIAL.
 *
 * Memory: max(chunk_size) ≈ MAX_CHUNK_GB
 * For n=28 with MAX_CHUNK_GB=2: about 4-5 passes at peak levels.
 *
 * Compile: gcc -O3 -o wn_rank4 wn_rank4.c
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

/* Max chunk size in bytes for next-level allocation */
#define MAX_CHUNK_BYTES (2LL * 1024 * 1024 * 1024) /* 2 GB */

/* In-memory threshold: if both arrays fit, do single pass in memory */
#define MEM_THRESHOLD (5LL * 1024 * 1024 * 1024) /* 5 GB */

/* Write array to a file. Returns 0 on success. */
static int write_array(const char *path, u128 *data, long long num_elements) {
    int fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0600);
    if (fd < 0) { perror("open for write"); return -1; }
    long long total = num_elements * sizeof(u128);
    long long written = 0;
    while (written < total) {
        long long chunk = total - written;
        if (chunk > (1LL << 30)) chunk = (1LL << 30);
        ssize_t w = write(fd, (char *)data + written, chunk);
        if (w <= 0) { perror("write"); close(fd); return -1; }
        written += w;
    }
    close(fd);
    return 0;
}

/* Write a chunk from buffer to a specific offset in an already-open file */
static int write_chunk_to_fd(int fd, long long offset, u128 *data, long long num_elements) {
    lseek(fd, offset, SEEK_SET);
    long long total = num_elements * sizeof(u128);
    long long written = 0;
    while (written < total) {
        long long chunk = total - written;
        if (chunk > (1LL << 30)) chunk = (1LL << 30);
        ssize_t w = write(fd, (char *)data + written, chunk);
        if (w <= 0) { perror("write chunk"); return -1; }
        written += w;
    }
    return 0;
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
    fflush(stderr);

    /* Track whether dp_curr is mmap'd or malloc'd */
    int curr_is_mmap = 0;
    long long curr_mmap_size = 0;

    for (int p = 1; p < n; p++) {
        int p_next = p + 1;
        long long cnt_next = binom[n][p_next];
        long long mem_curr = cnt_curr * p * (long long)sizeof(u128);
        long long mem_next = cnt_next * p_next * (long long)sizeof(u128);
        long long mem_both = mem_curr + mem_next;

        if (mem_both <= MEM_THRESHOLD && !curr_is_mmap) {
            /* In-memory single pass */
            fprintf(stderr, "  level %d->%d: %lld -> %lld masks (%.1f GB) (%lds)...\n",
                    p, p_next, cnt_curr, cnt_next,
                    (double)mem_both / 1e9, time(NULL) - t0);
            fflush(stderr);

            u128 *dp_next = (u128 *)calloc(cnt_next * p_next, sizeof(u128));
            if (!dp_next) {
                fprintf(stderr, "OOM in single-pass, switching to disk mode\n");
                goto disk_mode;
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
            fflush(stderr);
            free(dp_curr);
            dp_curr = dp_next;
            cnt_curr = cnt_next;
            p_curr = p_next;
            continue;
        }

    disk_mode:;
        /* Disk-backed multi-pass mode */

        /* Step 1: ensure curr is on disk */
        char curr_path[256];
        snprintf(curr_path, sizeof(curr_path), "/tmp/wn4_level_%d.bin", p);

        if (!curr_is_mmap) {
            /* Write curr to disk */
            if (write_array(curr_path, dp_curr, cnt_curr * p) < 0) return 1;
            free(dp_curr);
            dp_curr = NULL;
        }
        /* else: curr is already on disk from previous chunked iteration */

        /* mmap curr for reading */
        long long curr_file_size = cnt_curr * p * (long long)sizeof(u128);
        int curr_fd = open(curr_path, O_RDONLY);
        if (curr_fd < 0) { perror("open curr"); return 1; }
        u128 *curr_mmap = (u128 *)mmap(NULL, curr_file_size, PROT_READ, MAP_PRIVATE, curr_fd, 0);
        if (curr_mmap == MAP_FAILED) { perror("mmap curr"); close(curr_fd); return 1; }
        madvise(curr_mmap, curr_file_size, MADV_SEQUENTIAL);
        close(curr_fd);

        /* Step 2: determine number of passes */
        long long bytes_per_mask = p_next * (long long)sizeof(u128);
        long long masks_per_chunk = MAX_CHUNK_BYTES / bytes_per_mask;
        if (masks_per_chunk < 1) masks_per_chunk = 1;
        int num_passes = (int)((cnt_next + masks_per_chunk - 1) / masks_per_chunk);
        if (num_passes < 1) num_passes = 1;

        fprintf(stderr, "  level %d->%d: %lld -> %lld masks DISK "
                "(curr=%.1f GB, next=%.1f GB, %d passes, chunk=%.1f GB) (%lds)...\n",
                p, p_next, cnt_curr, cnt_next,
                (double)curr_file_size / 1e9, (double)mem_next / 1e9,
                num_passes, (double)(masks_per_chunk * bytes_per_mask) / 1e9,
                time(NULL) - t0);
        fflush(stderr);

        /* Pre-create next level file */
        char next_path[256];
        snprintf(next_path, sizeof(next_path), "/tmp/wn4_level_%d.bin", p_next);
        int next_fd = open(next_path, O_RDWR | O_CREAT | O_TRUNC, 0600);
        if (next_fd < 0) { perror("open next"); return 1; }
        long long next_file_size = cnt_next * p_next * (long long)sizeof(u128);
        if (ftruncate(next_fd, next_file_size) < 0) {
            perror("ftruncate next");
            close(next_fd);
            munmap(curr_mmap, curr_file_size);
            return 1;
        }

        /* Step 3: process each chunk */
        for (int pass = 0; pass < num_passes; pass++) {
            long long rank_lo = (long long)pass * masks_per_chunk;
            long long rank_hi = rank_lo + masks_per_chunk;
            if (rank_hi > cnt_next) rank_hi = cnt_next;
            long long chunk_masks = rank_hi - rank_lo;
            long long chunk_entries = chunk_masks * p_next;

            u128 *chunk_buf = (u128 *)calloc(chunk_entries, sizeof(u128));
            if (!chunk_buf) {
                fprintf(stderr, "Failed to allocate chunk %d (%.1f GB)\n",
                        pass, (double)(chunk_entries * sizeof(u128)) / 1e9);
                /* Try halving the chunk */
                munmap(curr_mmap, curr_file_size);
                close(next_fd); unlink(next_path); unlink(curr_path);
                return 1;
            }

            fprintf(stderr, "    pass %d/%d: next ranks %lld..%lld (%lds)...\n",
                    pass + 1, num_passes, rank_lo, rank_hi - 1, time(NULL) - t0);
            fflush(stderr);

            /* Re-advise sequential on curr for each pass */
            madvise(curr_mmap, curr_file_size, MADV_SEQUENTIAL);

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
                        unsigned int new_mask = mask | (1u << u);
                        long long nr = mask_rank(new_mask, n, p_next);
                        if (nr < rank_lo || nr >= rank_hi) continue;
                        u128 weight = (u == v + 1) ? 2 * cnt : cnt;
                        int local_u = local_index(new_mask, u);
                        chunk_buf[(nr - rank_lo) * p_next + local_u] += weight;
                    }
                }

                /* Periodically release read pages */
                if ((mi & 0xFFFFF) == 0 && mi > 0) {
                    long long done_bytes = mi * p * sizeof(u128);
                    /* Release pages we've already read */
                    long long release_end = done_bytes & ~4095LL; /* page align */
                    if (release_end > 0) {
                        madvise(curr_mmap, release_end, MADV_DONTNEED);
                    }
                }
            }

            /* Write chunk to next file */
            long long file_offset = rank_lo * p_next * sizeof(u128);
            if (write_chunk_to_fd(next_fd, file_offset, chunk_buf, chunk_entries) < 0) {
                free(chunk_buf);
                munmap(curr_mmap, curr_file_size);
                close(next_fd);
                return 1;
            }

            free(chunk_buf);
            fprintf(stderr, "    pass %d done (%lds)\n", pass + 1, time(NULL) - t0);
            fflush(stderr);
        }

        /* Clean up curr */
        munmap(curr_mmap, curr_file_size);
        unlink(curr_path);
        close(next_fd);

        /* Next level's data is on disk at next_path */
        /* For in-memory levels later (decreasing side), try to load into RAM */
        if (mem_next <= MEM_THRESHOLD) {
            /* Load next from disk into memory */
            int nfd = open(next_path, O_RDONLY);
            if (nfd < 0) { perror("open next for load"); return 1; }
            dp_curr = (u128 *)malloc(next_file_size);
            if (dp_curr) {
                long long rd = 0;
                while (rd < next_file_size) {
                    long long chunk = next_file_size - rd;
                    if (chunk > (1LL << 30)) chunk = (1LL << 30);
                    ssize_t r = read(nfd, (char *)dp_curr + rd, chunk);
                    if (r <= 0) { perror("read"); close(nfd); return 1; }
                    rd += r;
                }
                close(nfd);
                unlink(next_path);
                curr_is_mmap = 0;
            } else {
                /* Keep as file, mmap on next iteration */
                close(nfd);
                dp_curr = NULL;
                curr_is_mmap = 1;
            }
        } else {
            dp_curr = NULL;
            curr_is_mmap = 1;
        }

        cnt_curr = cnt_next;
        p_curr = p_next;
    }

    /* Sum over all endpoints at full mask */
    /* dp_curr might be NULL (on disk) or in memory */
    u128 W = 0;
    if (dp_curr) {
        for (int j = 0; j < n; j++) {
            W += dp_curr[j];
        }
    } else {
        /* Load from disk */
        char final_path[256];
        snprintf(final_path, sizeof(final_path), "/tmp/wn4_level_%d.bin", n);
        int ffd = open(final_path, O_RDONLY);
        if (ffd < 0) { perror("open final"); return 1; }
        u128 vals[32];
        read(ffd, vals, n * sizeof(u128));
        close(ffd);
        unlink(final_path);
        for (int j = 0; j < n; j++) W += vals[j];
    }

    printf("W(%d) = ", n);
    print_u128(W);
    printf("\n");

    fprintf(stderr, "Total time: %lds\n", time(NULL) - t0);
    return 0;
}
