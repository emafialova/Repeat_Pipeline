#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <sys/stat.h>
#include <time.h>
#include <errno.h>
#include <ctype.h>
#include <stdbool.h>
#include "./zhash-c/src/zhash.h"
#include <pthread.h>
#include <stdatomic.h> 
#include <stddef.h>

#define MAX_LENGTH 2000
#define MIN_COUNT 7
const char NT_letter_arr[5] = {'A','T','C','G','N'};
typedef uint64_t word_t;  // 64-bit unsigned integer

int length = 5;
/*
int compare_indexed_doubles_desc(const void *a, const void *b, void *arr_ptr) {
    const double *arr = arr_ptr;
    size_t ia = *(const size_t *)a;
    size_t ib = *(const size_t *)b;
    return (arr[ia] > arr[ib]) ? -1 : (arr[ia] < arr[ib]) ? 1 : 0;
}

void sort_indices_4(const double arr[4], size_t indices[4]) {
    for (size_t i = 0; i < 4; i++) indices[i] = i;
    qsort_r(indices, 4, sizeof(size_t), compare_indexed_doubles_desc, (void*)arr);
}*/

// THE FOLLOWING WORKS FOR MAC OS, the above for linux
static int compare_indexed_doubles_desc(void *arr_ptr, const void *a, const void *b) {
    const double *arr = arr_ptr;
    size_t ia = *(const size_t *)a;
    size_t ib = *(const size_t *)b;

    if (arr[ia] > arr[ib]) return -1;  // reverse order (descending)
    if (arr[ia] < arr[ib]) return 1;
    return 0;
}


 //* Returns the indices of the 4 elements in descending order (largest first).
 //* Works correctly with qsort_r on macOS/BSD.
 
void sort_indices_4(const double arr[4], size_t indices[4]) {
    for (size_t i = 0; i < 4; i++)
        indices[i] = i;

    // macOS/BSD qsort_r: context (arr) comes BEFORE comparator
    qsort_r(indices, 4, sizeof(size_t), (void *)arr, compare_indexed_doubles_desc);
}


typedef struct {
    size_t frequency;
} Base_Freq;

typedef struct {
    double A;
    double C;
    double G;
    double T;
    double N;
} RefFreq;

static inline double now_seconds_monotonic(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

static inline bool per_core_time_exceeded(double start_time, double max_seconds) {
    if (max_seconds <= 0.0) return false; // non-positive -> no limit
    double now = now_seconds_monotonic();
    return (now - start_time) >= max_seconds;
}

typedef struct {
    size_t size;
    size_t n_words;
    uint64_t *words;
} EncKmer;

// used to write one base as 2 bits (not the whole kmer, just one base)
static inline uint8_t base_to_bits(char base) {
    switch(toupper((unsigned char)base)) {
        case 'A': return 0; // 00
        case 'T': return 1; // 01
        case 'C': return 2; // 10
        case 'G': return 3; // 11
        case 'N': return -1; // unknown base
        default: return -1;
    }
};
// used to shift words by X bits to the left - makes space for the new base to be entered, X is 2 usually
void left_shift_words(uint64_t *words, size_t n_words, unsigned shift) {
    if (shift == 0) return;
    unsigned rshift = 64u - shift;
    uint64_t carry = 0;
    for (size_t i = 0; i < n_words; ++i) {
        uint64_t cur = words[i];
        uint64_t new_low = (cur << shift) | carry;
        carry = (cur >> rshift);
        words[i] = new_low;
    }
}

void right_shift_words(uint64_t *words, size_t n_words, unsigned shift) {
    if (shift == 0) return;
    unsigned lshift = 64u - shift;
    uint64_t carry = 0;
    for (size_t i = n_words; i-- > 0; ) {
        uint64_t cur = words[i];
        uint64_t new_high = (cur >> shift) | carry;
        carry = cur << lshift;
        words[i] = new_high;
    }
}

// used to encode kmer from string to bits - first creates an istance on EncKmer, fills the properties
EncKmer *encode_kmer(char *kmer, size_t k) {
    if (!kmer || k == 0) return NULL;
    size_t n_words = (2*k + 63)/64;
    EncKmer *encoded_kmer = malloc(sizeof(EncKmer));
    encoded_kmer->size = k;
    encoded_kmer->n_words = n_words;
    encoded_kmer->words = calloc(n_words, sizeof(uint64_t));
    if (!encoded_kmer->words) { free(encoded_kmer); return NULL; }

    //here it creates the words themselves - loops through kmer and encodes each base as 2 bits, adds to word and shifts
    for (size_t i = 0; i < k; ++i) {
        left_shift_words(encoded_kmer->words, n_words, 2u);
        uint8_t bits = base_to_bits(kmer[i]) & 3u;
        encoded_kmer->words[0] |= (uint64_t)bits; /* OR into least-significant bits */
    }
    size_t used_bits_last_word = (2*k) - (n_words - 1) * 64;
    if (used_bits_last_word < 64) {
        uint64_t mask = (used_bits_last_word == 0) ? 0ULL : ((1ULL << used_bits_last_word) - 1ULL);
        encoded_kmer->words[n_words - 1] &= mask;
    }
    return encoded_kmer;
}
//used to free kmer space
void free_kmer(EncKmer *kmer) {
    free(kmer->words);
    free(kmer);
}

char *decode_kmer(const EncKmer *ek) {
    if (!ek) return NULL;
    char *out = malloc(ek->size + 1);
    if (!out) return NULL;

    uint64_t *tmp = malloc(sizeof(uint64_t) * ek->n_words);
    if (!tmp) { free(out); return NULL; }
    memcpy(tmp, ek->words, ek->n_words * sizeof(uint64_t));

    for (size_t i = 0; i < ek->size; ++i) {
        /* extract least-significant 2 bits (this is the last appended base) */
        uint8_t b = (uint8_t)(tmp[0] & 3ULL);
        char c;
        switch (b) {
            case 0: c = 'A'; break;
            case 1: c = 'T'; break;
            case 2: c = 'C'; break;
            case 3: c = 'G'; break;
            //default: c = 'G'; break;
        }
        out[ek->size - 1 - i] = c; /* fill from the right (reverse order) */

        /* shift right the whole tmp by 2 bits to get next base into LSB */
        uint64_t carry = 0;
        for (size_t w = ek->n_words; w-- > 0; ) {
            uint64_t cur = tmp[w];
            uint64_t newlow = (cur >> 2) | (carry << 62);
            carry = cur & 3ULL; /* keep lowest 2 bits for next more-significant word */
            tmp[w] = newlow;
        }
    }
    out[ek->size] = '\0';
    free(tmp);
    return out;
}

static inline int hamming_encoded(const EncKmer *a, const EncKmer *b) {
    if (!a || !b || a->size != b->size) return -1;
    const uint64_t ONES_PAIR_MASK = 0x5555555555555555ULL; /* pattern 01 01 01 ... */
    int diff = 0;
    for (size_t i = 0; i < a->n_words; ++i) {
        uint64_t x = a->words[i] ^ b->words[i];
        uint64_t t = (x | (x >> 1)) & ONES_PAIR_MASK;
        diff += __builtin_popcountll(t);
    }
    return diff; /* number of bases that differ */
}

struct Kmer_values {
    uint8_t left_shift;
    uint8_t right_shift;
};

typedef struct {
    size_t size;
    size_t capacity;
    char **data; 
} StrArray;

void strarray_initialize(StrArray *arr) {
    arr->size = 0;
    arr->capacity = 10;
    arr->data = malloc(arr->capacity * sizeof(char*));
    if (!arr->data) { perror("malloc"); exit(1); }
}

void strarray_push(StrArray *arr, const char *s) {
    if (arr->size >= arr->capacity) {
        arr->capacity *= 2;
        arr->data = realloc(arr->data, arr->capacity * sizeof(char*));
        if (!arr->data) { perror("realloc"); exit(1); }
    }
    arr->data[arr->size] = strdup(s);
    arr->size++;
    //printf("something\n");
}

typedef struct {
    size_t size;
    size_t capacity;
    EncKmer **data; 
} EncKmerArray;

void enc_kmer_array_init(EncKmerArray *arr) {
    arr->size = 0;
    arr->capacity = 10;
    arr->data = malloc(arr->capacity * sizeof(EncKmer*));
    if (!arr->data) { perror("malloc"); exit(1); }
}

void enc_kmer_array_push(EncKmerArray *arr, EncKmer *kmer) {
    if (arr->size >= arr->capacity) {
        arr->capacity *= 2;
        //printf("We will reallocate\n");
        arr->data = realloc(arr->data, arr->capacity * sizeof(EncKmer*));
        if (!arr->data) { perror("realloc"); exit(1); }
    }
    EncKmer *kmer_copy = malloc(sizeof(EncKmer));
    kmer_copy->size = kmer->size;
    kmer_copy->n_words = kmer->n_words;
    kmer_copy->words = malloc(kmer->n_words * sizeof(uint64_t));
    memcpy(kmer_copy->words, kmer->words, kmer->n_words * sizeof(uint64_t));
    arr->data[arr->size++] = kmer_copy;
    //printf("We reallocated\n");
}

typedef struct {
    size_t *positions;
    size_t size;
    size_t capacity;
} PositionList;

void position_list_init(PositionList *plist, size_t initial_capacity) {
    plist->size = 0;
    plist->capacity = initial_capacity;
    plist->positions = malloc(initial_capacity * sizeof(size_t));
    if (!plist->positions) { perror("malloc"); exit(1); }
}

void positions_list_free(PositionList *plist) {
    plist->size = 0;
    plist->capacity = 0;
    free(plist->positions);
    plist->positions = NULL;
}

void position_list_push(PositionList *plist, size_t pos) {
    if (plist->size >= plist->capacity) {
        plist->capacity *= 2;
        plist->positions = realloc(plist->positions, plist->capacity * sizeof(size_t));
        if (!plist->positions) { perror("realloc"); exit(1); }
    }
    plist->positions[plist->size++] = pos;
}

void position_list_set2(PositionList *plist, size_t index, size_t pos) {        // if out of bounds, resize capacity of the list
    if (index >= plist->capacity) {
        plist->capacity = index + 10;
        plist->positions = realloc(plist->positions, plist->capacity * sizeof(size_t));
        if (!plist->positions) { perror("realloc"); exit(1); }
    }
    plist->positions[index] = pos;
}

void position_list_set(PositionList *plist, size_t index, size_t pos) {
    if (index >= plist->capacity) {
        // grow capacity to accommodate index
        size_t newcap = index + 10;
        plist->positions = realloc(plist->positions, newcap * sizeof(size_t));
        if (!plist->positions) { perror("realloc"); exit(1); }
        plist->capacity = newcap;
    }
    plist->positions[index] = pos;
    // update logical size if we wrote past current size
    if (index >= plist->size) plist->size = index + 1;
}

int filter_overlapping_positions(PositionList *plist, int kmer_length, double overlap_percentage) {
    if (plist->size <= 1) return plist->size;

    int overlap_distance = (int)ceil((double)kmer_length * overlap_percentage);
    if (overlap_distance <= 0) overlap_distance = 0;
    //size_t *positions;
    //size_t size;
    //size_t capacity;
    size_t new_size = 0;
    plist->positions[new_size++] = plist->positions[0]; // Always keep the first position
    for(size_t i = 1; i < plist->size; ++i) {
        //if(plist->positions[i] >= plist->positions[new_size-1] + overlap_distance) {
        //if(overlap_distance >= plist->positions[new_size-1]+kmer_length-plist->positions[i]) {
        if ((plist->positions[i] - plist->positions[new_size-1]) >= (kmer_length - overlap_distance)){
            //printf("Yay, saving this! %zu\n", plist->positions[i]);
            if (i != new_size) {
                    plist->positions[new_size] = plist->positions[i];
                }
                new_size++;
        }
    }
    plist->size = new_size;
    //("This is the new size: %zu\n", new_size);
    return new_size;
};

typedef struct {
    size_t words;
    size_t pattern_length;
    word_t *A;
    word_t *C;
    word_t *G;
    word_t *T;
} BitMasks;

static word_t *calloc_words(size_t words) {
    word_t *p = (word_t*)calloc(words, sizeof(word_t));
    if (!p) { perror("calloc"); exit(1); }
    return p;
};

BitMasks build_pattern_masks(char *pattern) {
    size_t n = strlen(pattern);
    size_t words_needed = (n+63)/64;
    BitMasks pattern_mask;
    pattern_mask.words = words_needed;
    pattern_mask.pattern_length = n;
    pattern_mask.A = calloc_words(words_needed);
    pattern_mask.T = calloc_words(words_needed);
    pattern_mask.C = calloc_words(words_needed);
    pattern_mask.G = calloc_words(words_needed);

    for (size_t i = 0; i < n; ++i) {
        size_t w = i / 64;
        unsigned shift = i % 64;
        char ch = pattern[i];
        if (ch == 'A') pattern_mask.A[w] |= (1ULL << shift);
        else if (ch == 'C') pattern_mask.C[w] |= (1ULL << shift);
        else if (ch == 'G') pattern_mask.G[w] |= (1ULL << shift);
        else if (ch == 'T') pattern_mask.T[w] |= (1ULL << shift);
    }
    return pattern_mask;
};

BitMasks build_text_masks(char *text) {
    size_t n = strlen(text);
    size_t words_needed = (n+63)/64;
    BitMasks text_masks;
    text_masks.words = words_needed;
    text_masks.pattern_length = n;
    text_masks.A = calloc_words(words_needed);
    text_masks.T = calloc_words(words_needed);
    text_masks.C = calloc_words(words_needed);
    text_masks.G = calloc_words(words_needed);

    for (size_t i = 0; i < n; ++i) {
        size_t w = i / 64;
        unsigned shift = i % 64;
        char ch = text[i];
        if (ch == 'A') text_masks.A[w] |= (1ULL << shift);
        else if (ch == 'C') text_masks.C[w] |= (1ULL << shift);
        else if (ch == 'G') text_masks.G[w] |= (1ULL << shift);
        else if (ch == 'T') text_masks.T[w] |= (1ULL << shift);
    }
    return text_masks;
};

static void free_bitmasks(BitMasks *mask) {
    if (!mask) return;
    free(mask->A);
    free(mask->C);
    free(mask->G);
    free(mask->T);
    mask->A = mask->C = mask->G = mask->T = NULL;
}

// popcount for 64-bit word
static inline int popcount_word(word_t x) {
    return __builtin_popcountll(x);
}

// Shift pattern mask left by s bits into 'out', handling multiple 64-bit words
static void shift_pattern_to_target(const word_t *pat_words, size_t pat_words_len,
                                    word_t *out, size_t out_words, size_t s) {
    if (s == 0) {
        for (size_t i = 0; i < pat_words_len && i < out_words; ++i)
            out[i] |= pat_words[i];
        return;
    }
    size_t word_shift = s / 64;
    unsigned bit_shift = s % 64;
    for (size_t i = 0; i < pat_words_len; ++i) {
        size_t target = i + word_shift;
        if (target >= out_words) break;
        word_t v = pat_words[i];
        out[target] |= v << bit_shift;
        if (bit_shift != 0 && target + 1 < out_words) {
            out[target + 1] |= v >> (64 - bit_shift);
        }
    }
}

// Count matches for a single base
static int count_matches_for_base(const word_t *text_mask, size_t text_words,
                                  const word_t *pat_mask, size_t pat_words, size_t shift_s) {
    word_t *shifted = calloc_words(text_words);
    shift_pattern_to_target(pat_mask, pat_words, shifted, text_words, shift_s);
    int sum = 0;
    for (size_t w = 0; w < text_words; ++w)
        sum += popcount_word(text_mask[w] & shifted[w]);
    free(shifted);
    return sum;
}

// Scan either the whole text (position_list == NULL) or only the given starts.
int hamming_bitparallel_search(
    const char *pattern,
    const char *text,
    int max_HD,
    BitMasks *text_masks,            // may be NULL -> built locally
    PositionList *positions,               // output buffer (optional, may be NULL)
    size_t max_positions,            // capacity of positions[]
    const size_t *position_list,     // candidate starts (NULL -> full scan)
    size_t len_position_list         // number of candidates in position_list
) {
    BitMasks local_text_masks;
    bool using_local = false;

    if (text_masks == NULL) {
        local_text_masks = build_text_masks((char*)text);
        text_masks = &local_text_masks;
        using_local = true;
    }

    BitMasks pattern_masks = build_pattern_masks((char*)pattern);
    size_t m = strlen(pattern);
    size_t n = strlen(text);
    if (m == 0 || m > n) {
        if (using_local) free_bitmasks(&local_text_masks);
        free_bitmasks(&pattern_masks);
        return 0;
    }

    size_t start_count = position_list ? len_position_list : (n - m + 1);
    int nmatches = 0;

    for (size_t i = 0; i < start_count; ++i) {
        size_t s = position_list ? position_list[i] : i;
        if (s > n - m) continue; // bounds check

        int matches = 0;

        // A, C, G, T in that order
        for (int b = 0; b < 4; ++b) {
            const word_t *t_mask = NULL;
            const word_t *p_mask = NULL;

            switch (b) {
                case 0: t_mask = text_masks->A; p_mask = pattern_masks.A; break;
                case 1: t_mask = text_masks->C; p_mask = pattern_masks.C; break;
                case 2: t_mask = text_masks->G; p_mask = pattern_masks.G; break;
                case 3: t_mask = text_masks->T; p_mask = pattern_masks.T; break;
            }

            size_t word_shift = s / 64;
            unsigned bit_shift = s % 64;

            /*for (size_t w = 0; w < pattern_masks.words; ++w) {
                if (word_shift + w >= text_masks->words) break;

                word_t shifted = p_mask[w] << bit_shift;

                if (bit_shift != 0 &&
                    word_shift + w + 1 < text_masks->words &&
                    w + 1 < pattern_masks.words) {
                    shifted |= p_mask[w + 1] >> (64 - bit_shift);
                }
            
                word_t common = t_mask[word_shift + w] & shifted;
                matches += __builtin_popcountll(common);
            }*/

            for (size_t w = 0; w < pattern_masks.words; ++w) {
                if (word_shift + w >= text_masks->words) break;

                word_t v = p_mask[w];

                // main aligned part
                word_t shifted_low = v << bit_shift;
                word_t common_low = t_mask[word_shift + w] & shifted_low;
                matches += __builtin_popcountll(common_low);

                // carry-over into next word
                if (bit_shift != 0 && word_shift + w + 1 < text_masks->words) {
                    word_t shifted_high = v >> (64 - bit_shift);
                    word_t common_high = t_mask[word_shift + w + 1] & shifted_high;
                    matches += __builtin_popcountll(common_high);
                }
            }

        }

        int mismatches = (int)m - matches;
        if (mismatches <= max_HD) {
            if (positions && (size_t)nmatches < max_positions) {
                position_list_set(positions, nmatches, s);
            }
            nmatches++;
        }
    }

    free_bitmasks(&pattern_masks);
    if (using_local) free_bitmasks(&local_text_masks);
    return nmatches;
}

// bool present_in_starting_kmers(struct ZHashTable *set, char *kmer) {
//     //char *encoded_kmer = decode_kmer(kmer);
//     bool exists = zhash_exists(set, kmer);
//     //("We managed to check the list and the asnwer is: \n");
//     //free(encoded_kmer);
//     return exists;
// }

int hamming_distance(int length, char kmer_1[length], char kmer_2[length]) {
    int count = 0;
    for (int i = 0; i < length; i++) {
        if (kmer_1[i] != kmer_2[i]) {
            count++;
        }
    }
    return count;
};

bool is_similar_to_known_kmer(EncKmerArray *all_kmers, EncKmer *new_kmer, int max_HD) {
    for (size_t i = 0; i < all_kmers->size; ++i) {
        EncKmer *known_kmer = all_kmers->data[i];
        if (known_kmer->size == new_kmer->size) {
            int HD = hamming_encoded(known_kmer, new_kmer);
            if (HD <= max_HD) return true;
        }
    }
    return false;
};

void simple_extend_kmers(EncKmer *new_kmer, EncKmerArray*all_kmers, /*kmer_dict*/ char *region_seq, 
                    int hd_percentage, EncKmer *old_kmer, FILE *abc_file, FILE *tsv_file, 
                    int left_shift, int right_shift, struct ZHashTable *starting_kmers,
                    PositionList *prev_positions, size_t prev_count, double start_time, double max_seconds,
                    bool *early_stopping_flag, RefFreq *ref_freq);

void simple_extending(char *new_kmer, EncKmerArray *all_kmers, /*kmer_dict*/ char *region_seq, 
                    int hd_percentage, char *old_kmer, FILE *abc_file, FILE *tsv_file, 
                    int left_shift, int right_shift, struct ZHashTable *starting_kmers,
                    size_t *prev_positions, size_t prev_count, double start_time, double max_seconds,
                    bool *early_stopping_flag, RefFreq *ref_freq, double base_freq) {
    
    if (per_core_time_exceeded(start_time, max_seconds)) {
        // optional: fprintf(stderr, "Timeout: stopping extension of %s\n", old_kmer);
        *early_stopping_flag = true;
        free(new_kmer);
        return;
    }
    int max_HD = (int)ceil(strlen(new_kmer) * (hd_percentage / 100.0));
    EncKmer *ek_new_kmer = encode_kmer(new_kmer, strlen(new_kmer));
    if (is_similar_to_known_kmer(all_kmers, ek_new_kmer, max_HD)) {
        fprintf(abc_file, "%s\tSIMILAR\t%s\n", old_kmer, new_kmer);
        free(new_kmer);
        //abc zapis
    } else {
        PositionList positions;
        //printf("we even got here\n");
        position_list_init(&positions, prev_count > 0 ? prev_count : 16);
        //printf("I will use hamming bitparalel with this kmer: %s\n", new_kmer);
        int occurences = hamming_bitparallel_search(new_kmer, region_seq, max_HD, 
                                            NULL, &positions, positions.capacity, prev_positions, prev_count);
        if (occurences > 1) {
            double overlap_percentage = 0.0;
            //printf("This is the occurences before: %d\n", occurences);
            occurences = filter_overlapping_positions(&positions, strlen(new_kmer), overlap_percentage);
            //printf("These are now: %d\n", occurences);
        }
        //printf("We want to check the occurences: %d\n", occurences);
        if (occurences >= MIN_COUNT) {
            int starting_kmer_flag = (int)zhash_exists(starting_kmers, new_kmer);
            fprintf(tsv_file, "%s\t%d\t%zu\t1\t[", new_kmer, occurences, strlen(new_kmer));
            for (int pi = 0; pi < occurences; ++pi) {
                fprintf(tsv_file, "%zu%s", positions.positions[pi], (pi + 1 < occurences) ? "," : "");
            }
            fprintf(tsv_file, "]\t%d\t%f\n", starting_kmer_flag, base_freq);
            fprintf(abc_file, "%s\t1\t%s\n", old_kmer, new_kmer);

            /*struct Kmer_values* added_kmer = malloc(sizeof(struct Kmer_values));
            added_kmer->left_shift = left_shift;
            added_kmer->right_shift = right_shift;
            //printf("This is the left shift: %d\n", right_shift);
            zhash_set(kmer_dict, new_kmer, added_kmer);*/
            //printf("We got here.\n");
            //strarray_push(all_kmers, new_kmer);
            //printf("The error is after this: \n");
            enc_kmer_array_push(all_kmers, ek_new_kmer);
            //printf("The error is after this: \n");
            EncKmer *ek_old_kmer = encode_kmer(old_kmer, strlen(old_kmer));
            //printf("We got here.\n");
            //printf("We also got here!!\n");
            //printf("occurences n2: %d\n", occurences);
            simple_extend_kmers(ek_new_kmer, all_kmers, /*kmer_dict*/ region_seq, hd_percentage, ek_old_kmer,
                    abc_file, tsv_file, left_shift, right_shift, starting_kmers,
                    &positions, occurences, start_time, max_seconds, early_stopping_flag, ref_freq);
            //printf("We have and issue here!!!!");
            free_kmer(ek_old_kmer);
            free(new_kmer);
        } else {
            free(new_kmer);
        }
        positions_list_free(&positions);
    }
    free_kmer(ek_new_kmer);
};

void simple_extend_kmers(EncKmer *ek_new_kmer, EncKmerArray *all_kmers, /*struct ZHashTable *kmer_dict,*/ char *region_seq, 
                    int hd_percentage, EncKmer *ek_old_kmer, FILE *abc_file, FILE *tsv_file, 
                    int left_shift, int right_shift, struct ZHashTable *starting_kmers,
                    PositionList *prev_positions, size_t prev_count, double start_time, double max_seconds, 
                    bool *early_stopping_flag, RefFreq *ref_freq) {
    //printf("We got here, inside simple extend!!\n");
    if ((int)ek_new_kmer->size > MAX_LENGTH) return;
    PositionList right_positions, left_positions;
    position_list_init(&right_positions, prev_count > 0 ? prev_count : 16);
    position_list_init(&left_positions, prev_count > 0 ? prev_count : 16);
    //printf("We got a positions list!!\n");
    if (prev_positions != NULL) {
    for (size_t j = 0; j < prev_count; ++j) {
        position_list_push(&right_positions, prev_positions->positions[j]);
        if (prev_positions->positions[j] > 0) {
            position_list_push(&left_positions, prev_positions->positions[j] - 1);
        }
    }
    }
   //printf("This is the positions list: %i", prev_positions);
    size_t next_counts[5] = {0,0,0,0,0};  // A, T, C, G
    size_t prev_counts[5] = {0,0,0,0,0};
    size_t total_next = 0;
    size_t total_prev = 0;
    
    size_t k = ek_new_kmer->size;
    
    for (size_t j = 0; j < prev_count; ++j) {
        size_t pos = prev_positions->positions[j];
    
        // Right extension base (after k-mer)
        if (pos + k < strlen(region_seq)) {
            char nb = region_seq[pos + k];
            switch (nb) {
                case 'A': next_counts[0]++; break;
                case 'T': next_counts[1]++; break;
                case 'C': next_counts[2]++; break;
                case 'G': next_counts[3]++; break;
                case 'N': next_counts[4]++; break;
            }
            total_next++;
        }
    
        // Left extension base (before k-mer)
        if (pos > 0) {
            char pb = region_seq[pos - 1];
            switch (pb) {
                case 'A': prev_counts[0]++; break;
                case 'T': prev_counts[1]++; break;
                case 'C': prev_counts[2]++; break;
                case 'G': prev_counts[3]++; break;
                case 'N': prev_counts[4]++; break;
            }
            total_prev++;
        }
    }
    
    // Compute relative frequencies
    double next_freqs[5] = {0,0,0,0,0};
    double prev_freqs[5] = {0,0,0,0,0};
    for (int b = 0; b < 5; ++b) {
        if (total_next > 0) next_freqs[b] = (double)next_counts[b] / (double)total_next;
        if (total_prev > 0) prev_freqs[b] = (double)prev_counts[b] / (double)total_prev;
    }
    bool extend_next[4] = {false,false,false,false};
    bool extend_prev[4] = {false,false,false,false};
    size_t freq_array_next[4] = {0,0,0,0};
    size_t freq_array_prev[4] = {0,0,0,0};
    sort_indices_4(next_freqs, freq_array_next);
    sort_indices_4(prev_freqs, freq_array_prev);
    //printf("next %d\n", freq_array_next);
    //printf("next %d\n", freq_array_prev);
    for (int b = 0; b < 4; ++b) {
        switch (b) {
            case 0: extend_next[b] = (next_freqs[b] >= ref_freq->A+0.05); extend_prev[b] = (prev_freqs[b] >= ref_freq->A+0.05); break;
            case 1: extend_next[b] = (next_freqs[b] >= ref_freq->T+0.05); extend_prev[b] = (prev_freqs[b] >= ref_freq->T+0.05); break;
            case 2: extend_next[b] = (next_freqs[b] >= ref_freq->C+0.05); extend_prev[b] = (prev_freqs[b] >= ref_freq->C+0.05); break;
            case 3: extend_next[b] = (next_freqs[b] >= ref_freq->G+0.05); extend_prev[b] = (prev_freqs[b] >= ref_freq->G+0.05); break;
        }
    }
    char *decoded_new_kmer = decode_kmer(ek_new_kmer);
    for (size_t j = 0; j < 4; ++j) {
        size_t i = freq_array_next[j];
        //printf("base: %zu\n", i);
        if (extend_next[i]) {
            uint8_t base_bits = base_to_bits(NT_letter_arr[i]);
            // --- Extend to the right ---
            EncKmer *ek_right = malloc(sizeof(EncKmer));
            ek_right->size = ek_new_kmer->size + 1;
            ek_right->n_words = (size_t)ceil(2.0 * ek_right->size / 64.0);
            ek_right->words = calloc(ek_right->n_words, sizeof(uint64_t));
            // Copy existing kmer bits
            for (size_t w = 0; w < ek_new_kmer->n_words; ++w) {
                ek_right->words[w] = ek_new_kmer->words[w];
            }
            // Shift left 2 bits to append new base
            left_shift_words(ek_right->words, ek_right->n_words, 2);
            ek_right->words[0] |= (uint64_t)base_bits;
            //printf("This is the kmer before extension: %s\n", decode_kmer(ek_new_kmer));
            //printf("We managed extension to right!\n");
            //printf("This is the decoded kmer: %s\n", decode_kmer(ek_right));
            // Recursive extension
            char *right_kmer_str = decode_kmer(ek_right);
            simple_extending(right_kmer_str, all_kmers,/*kmer_dict,*/region_seq,hd_percentage,
                            decoded_new_kmer,abc_file,tsv_file,left_shift,right_shift + 1,
                            starting_kmers,(right_positions.size > 0 ? right_positions.positions : NULL),
                            right_positions.size, start_time, max_seconds, early_stopping_flag, ref_freq, next_freqs[i]);
            free_kmer(ek_right);
            //free(right_kmer_str);
        }
        i = freq_array_prev[j];
        if (extend_prev[i]) {
            uint8_t base_bits = base_to_bits(NT_letter_arr[i]);
        // --- Extend to the left ---
            EncKmer *ek_left = malloc(sizeof(EncKmer));
            //printf("I got to do some left extension too!\n");
            //printf("this is the new_kmer that will be edxtended: %s\n", decode_kmer(ek_new_kmer));
            ek_left->size = ek_new_kmer->size + 1;
            ek_left->n_words = (size_t)ceil(2.0 * ek_left->size / 64.0);
            ek_left->words = calloc(ek_left->n_words, sizeof(uint64_t));

            // Copy existing bits shifted right 2 to make room for new base at MSB
            for (size_t w = 0; w < ek_new_kmer->n_words; ++w) {
                ek_left->words[w] = ek_new_kmer->words[w];
            }
            //right_shift_words(ek_left->words, ek_left->n_words, 2);
            // Set new base in the MSB of leftmost word
            size_t msb_bit_index = 2 * (ek_left->size - 1); // index of the new first base
            size_t word_idx = msb_bit_index / 64;
            unsigned bit_offset = msb_bit_index % 64;
            ek_left->words[word_idx] |= ((uint64_t)base_bits << bit_offset);
            char *left_kmer_str = decode_kmer(ek_left);
            //printf("THis is the left extended kmer: %s\n", decode_kmer(ek_left));
            simple_extending(left_kmer_str, all_kmers,/*kmer_dict*/region_seq,hd_percentage,
                            decoded_new_kmer,abc_file,tsv_file,left_shift + 1,right_shift,
                            starting_kmers,(left_positions.size > 0 ? left_positions.positions : NULL),
                            left_positions.size, start_time, max_seconds, early_stopping_flag, ref_freq, prev_freqs[i]);
            free_kmer(ek_left);
            //free(left_kmer_str);
        }      
    }
    positions_list_free(&right_positions);
    positions_list_free(&left_positions);
    free(decoded_new_kmer);
};

void process_row(int start, int end, int region_id, int hd_percentage, char* sequence, 
                char *output_prefix, const char *gene_name, char **kmer_list, size_t kmer_list_length,
                double per_core_timeout, FILE *summary_file) {
    //clock_t region_start_time = clock();
    double region_start = now_seconds_monotonic();
    //printf("Time starting\n");
    int region_size = end - start + 1;
    char *region = malloc(region_size +1);
    strncpy(region, sequence + start, region_size);
    region[region_size] = '\0';
    //printf("This is the region size: %d\n", region_size);
    //printf("Calling build_text_masks...\n");
    BitMasks text_masks = build_text_masks((char*)region);
    //printf("Back from build_text_masks...\n");
    //struct ZHashTable *kmer_dict = zcreate_hash_table(); 
    EncKmerArray all_kmers = {0, 10, NULL};
    enc_kmer_array_init(&all_kmers);
    //StrArray all_kmers = {0, 10, NULL};
    //strarray_initialize(&all_kmers);

    Base_Freq A;
    Base_Freq C;
    Base_Freq T;
    Base_Freq G;
    Base_Freq N;
    A.frequency = 0;
    C.frequency = 0;
    T.frequency = 0;
    G.frequency = 0;
    N.frequency = 0;
    for(size_t i = 0; i < region_size; ++i) {
        char base = region[i];
        if (base == 'A') {
            A.frequency += 1;
        } else if (base == 'C') {
            C.frequency += 1;
        } else if (base == 'T') {
            T.frequency += 1;
        } else if (base == 'G') {
            G.frequency += 1;
        } else if (base == 'N') {
            N.frequency += 1;
        }
    }
    RefFreq ref_freq;
    ref_freq.A = (double)A.frequency / (double)region_size;
    ref_freq.C = (double)C.frequency / (double)region_size;
    ref_freq.G = (double)G.frequency / (double)region_size;
    ref_freq.T = (double)T.frequency / (double)region_size;
    ref_freq.N = (double)N.frequency / (double)region_size;
    //printf(&ref_freq);

    struct ZHashTable *starting_kmers = zcreate_hash_table();
    for (size_t i = 0; i < kmer_list_length; ++i) {
        char* kmer = kmer_list[i];
        zhash_set(starting_kmers, kmer, NULL);
    }
    //printf("Back from creating dictionaries too!\n");
    char fname_tsv[128];
    char fname_abc[128];
    snprintf(fname_tsv, sizeof(fname_tsv), "%s.tsv", output_prefix);
    snprintf(fname_abc, sizeof(fname_abc), "%s.abc", output_prefix);
    FILE *file_tsv;
    FILE *file_abc;
    file_tsv = fopen(fname_tsv, "w");
    if (!file_tsv) {
    perror("fopen TSV");
    exit(1);
    }
    file_abc = fopen(fname_abc, "w");

    fprintf(file_tsv, "%s|%d|%d|%d\n", gene_name, start, end, region_id);
    //printf("Created files??\n");
    fprintf(file_tsv, "kmer\tcount\tlength\tcompactness\tpositions\tstarting_kmer_flag\tbase_rel.freq\n");
    //printf("Created files??\n");
    fprintf(file_abc, "%s|%d|%d|%d\n", gene_name, start, end, region_id);
    fprintf(file_abc, "kmer\tconnection\tkmer_2\n");    
    
    bool early_stopping_flag = false;

    for (int i = 0; i < kmer_list_length; ++i) {
        double kmer_start = now_seconds_monotonic();
        char *kmer = kmer_list[i];
        EncKmer *ek_kmer = encode_kmer(kmer, strlen(kmer));
        int right_shift = 0;
        int left_shift = 0;
        int max_HD = (int)ceil(strlen(kmer) * (hd_percentage / 100.0));
        //size_t positions[20000];
        PositionList positions;
        position_list_init(&positions, 10000);
        int occurences = hamming_bitparallel_search(kmer, region, max_HD, &text_masks, &positions, 1000000000, NULL, 0);
        fprintf(file_tsv, "%s\t%d\t%zu\t1\t[", kmer, occurences, strlen(kmer));
        for (int pi = 0; pi < occurences; ++pi) {
            fprintf(file_tsv, "%zu%s", positions.positions[pi],
                    (pi + 1 < occurences) ? "," : "");
        }
        fprintf(file_tsv, "]\t1\t1\n");
        simple_extend_kmers(ek_kmer, &all_kmers, /*kmer_dict,*/ region, hd_percentage, ek_kmer, file_abc, file_tsv,
             left_shift, right_shift, starting_kmers, &positions, occurences, kmer_start, per_core_timeout, &early_stopping_flag
            , &ref_freq);
        free_kmer(ek_kmer);
        positions_list_free(&positions);
    };
    
    fprintf(summary_file, "%s|%d|%d|%d|%d\n", gene_name, start, end, region_id, early_stopping_flag);
    fclose(file_abc);
    fclose(file_tsv);
    free_bitmasks(&text_masks);
    for (size_t i = 0; i < all_kmers.size; i++) free_kmer(all_kmers.data[i]);
    free(all_kmers.data);

    /*for (size_t i = 0; i < kmer_dict->entry_count; i++) {
        free_kmer(kmer_dict->entries[i]->val);
        free(kmer_dict->entries[i]->key);
    }*/
    zfree_hash_table(starting_kmers);
    //zfree_hash_table(kmer_dict);
    free(region);

    //clock_t region_end_time = clock();
    //double time_spent = (double)(region_end_time - region_start_time) / CLOCKS_PER_SEC;
    double region_end = now_seconds_monotonic();
    double time_spent = region_end - region_start;
    printf("Finished processing %s in %.2f seconds.\n", output_prefix, time_spent);
}

char *read_fasta_single(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) { fprintf(stderr, "Error opening FASTA: %s\n", filename); perror("fopen"); exit(1); }
    size_t cap = 1024*1024;
    char *seq = (char*)malloc(cap);
    seq[0] = '\0';
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    size_t current_len = 0;
    // Read lines; skip header lines starting with '>'
    while ((read = getline(&line, &len, f)) != -1) {
        if (read <= 0) continue;
        if (line[0] == '>') continue;
        // strip newline
        while (read > 0 && (line[read-1] == '\n' || line[read-1] == '\r')) { line[--read] = 0; }
        // ensure capacity
        //size_t need = strlen(seq) + read + 1;
        size_t need = current_len + read + 1;
        if (need > cap) {
            while (cap < need) cap *= 2;
            seq = (char*)realloc(seq, cap);
        }
        //strcat(seq, line);
        for(int i=0; i<read; ++i) {
            char c = toupper((unsigned char)line[i]);
            if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
                seq[current_len++] = 'N'; 
            } else {
                seq[current_len++] = c;
            }
        }
        seq[current_len] = '\0';
    }
    free(line);
    fclose(f);
    return seq;
}

typedef struct {
    int start, end;
    char **kmers;
    int kmer_count;
} RegionTask;

RegionTask *load_tsv_tasks(const char *filename, int *out_n) {
    FILE *f = fopen(filename, "r");
    if (!f) { fprintf(stderr, "Error opening TSV: %s\n", filename); perror("fopen"); exit(1); }
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    // skip two header lines
    getline(&line, &len, f);
    getline(&line, &len, f);

    RegionTask *tasks = NULL;
    int capacity = 0;
    int count = 0;
    while ((read = getline(&line, &len, f)) != -1) {
        if (read <= 1) continue;
        // strip newline
        while (read > 0 && (line[read-1] == '\n' || line[read-1] == '\r')) { line[--read] = 0; }
        // parse first two ints and the rest as kmers string
        int start = 0, end = 0;
        char *kmers_str = NULL;
        // We use strtok-ish parse: split by tabs
        char *copy = strdup(line);
        char *p = copy;
        char *tok1 = strsep(&p, "\t");
        char *tok2 = p ? strsep(&p, "\t") : NULL;
        char *tok3 = p ? p : NULL;
        if (!tok1 || !tok2 || !tok3) { free(copy); continue; }
        start = atoi(tok1);
        end = atoi(tok2);
        kmers_str = strdup(tok3);
        // Now parse kmers_str like "['ATG','CGA']" or with double quotes
        // We'll extract contiguous alpha sequences of A/T/C/G letters.
        char **kmers = NULL;
        int kcap = 0, kcount = 0;
        char *q = kmers_str;
        while (*q) {
            // skip non-ATCG letters
            while (*q && !( (*q=='A')||(*q=='T')||(*q=='C')||(*q=='G')||(*q=='N')||
                            (*q=='a')||(*q=='t')||(*q=='c')||(*q=='g')||(*q=='n') )) q++;
            if (!*q) break;
            char tmp[4096]; int ti = 0;
            while (*q && ( (*q=='A')||(*q=='T')||(*q=='C')||(*q=='G')||(*q=='N')||
                           (*q=='a')||(*q=='t')||(*q=='c')||(*q=='g')||(*q=='n')  )) {
                tmp[ti++] = toupper((unsigned char)*q);
                q++;
            }
            tmp[ti] = 0;
            if (ti > 0) {
                if (kcount + 1 > kcap) { kcap = kcap ? kcap * 2 : 8; kmers = realloc(kmers, sizeof(char*) * kcap); }
                kmers[kcount++] = strdup(tmp);
            }
        }
        free(kmers_str);
        free(copy);

        if (count + 1 > capacity) {
            capacity = capacity ? capacity * 2 : 8;
            tasks = realloc(tasks, sizeof(RegionTask) * capacity);
        }
        tasks[count].start = start;
        tasks[count].end = end;
        tasks[count].kmer_count = kcount;
        tasks[count].kmers = kmers;
        count++;
    }
    free(line);
    fclose(f);
    *out_n = count;
    return tasks;
}

typedef struct {
    int task_index;
    RegionTask *tasks;
    char *sequence;
    const char *output_dir;
    const char *gene_name;
    int hd_percentage;
    double per_core_timeout_seconds;
    FILE *summary_file;
} Task;

typedef struct {
    Task *queue;
    int capacity;
    int size;
    int front;
    int rear;

    pthread_mutex_t lock;
    pthread_cond_t cond;
    int shutdown;
} TaskQueue;

void taskqueue_init(TaskQueue *q, int capacity) {
    q->queue = malloc(sizeof(Task) * capacity);
    q->capacity = capacity;
    q->size = 0;
    q->front = 0;
    q->rear = 0;
    pthread_mutex_init(&q->lock, NULL);
    pthread_cond_init(&q->cond, NULL);
    q->shutdown = 0;
}

void taskqueue_destroy(TaskQueue *q) {
    free(q->queue);
    pthread_mutex_destroy(&q->lock);
    pthread_cond_destroy(&q->cond);
}

void taskqueue_push(TaskQueue *q, Task task) {
    pthread_mutex_lock(&q->lock);
    while (q->size == q->capacity) {
        // optionally grow here
        pthread_cond_wait(&q->cond, &q->lock);
    }
    q->queue[q->rear] = task;
    q->rear = (q->rear + 1) % q->capacity;
    q->size++;
    pthread_cond_broadcast(&q->cond);
    pthread_mutex_unlock(&q->lock);
}

int taskqueue_pop(TaskQueue *q, Task *out) {
    pthread_mutex_lock(&q->lock);
    while (q->size == 0 && !q->shutdown) {
        pthread_cond_wait(&q->cond, &q->lock);
    }
    if (q->shutdown && q->size == 0) {
        pthread_mutex_unlock(&q->lock);
        return 0;
    }
    *out = q->queue[q->front];
    q->front = (q->front + 1) % q->capacity;
    q->size--;
    pthread_cond_broadcast(&q->cond);
    pthread_mutex_unlock(&q->lock);
    return 1;
}

void* worker_thread(void *arg) {
    TaskQueue *q = (TaskQueue*)arg;
    Task task;
    while (taskqueue_pop(q, &task)) {
        RegionTask *rt = &task.tasks[task.task_index];
        char outprefix[1024];
        snprintf(outprefix, sizeof(outprefix), "%s/region_%d", task.output_dir, task.task_index+1);
        printf("Size of queue: %d\n", q->size);
        process_row(rt->start, rt->end, task.task_index+1,
                    task.hd_percentage, task.sequence, outprefix,
                    task.gene_name, rt->kmers, rt->kmer_count, task.per_core_timeout_seconds, task.summary_file);
    }
    return NULL;
}

#define N_THREADS 64

int main(int argc, char *argv[]) {
    //if (argc < 4) {printf("Incorrect number of arguments, please give: input file, sequence file, output directory and max HD.\n");};
    if (argc != 5) {
        printf("\nWelcome to the Kmer Extension Tool!\n");
        printf("Usage: %s <input_file.tsv> <sequence.fasta> <output_dir> <max_HD>\n", argv[0]);
        printf("\nArguments:\n");
        printf("  input_file.tsv   - TSV file with region definitions\n");
        printf("  sequence.fasta   - FASTA file with the sequence\n");
        printf("  output_dir       - Directory where results will be saved\n");
        printf("  max_HD           - Maximum Hamming distance (integer percent)\n\n");
        return 1; // exit early
    }

    const char *input_file = argv[1];
    const char *fasta_file = argv[2];
    const char *output_dir = argv[3];
    int hd_percentage = (int)ceil(atof(argv[4]));

    // ensure output_dir exists
    if (mkdir(output_dir, 0755) != 0 && errno != EEXIST) {
        perror("mkdir");
        return 1;
    }
    
    char *sequence = read_fasta_single(fasta_file);

    const char *base = strrchr(fasta_file, '/');
    base = base ? base + 1 : fasta_file;
    char gene_name[32];
    if (strstr(base, "_seq.fasta")) {
        strncpy(gene_name, base, sizeof(gene_name)-1);
        gene_name[sizeof(gene_name)-1] = 0;
        char *p = strstr(gene_name, "_seq.fasta");
        if (p) *p = '\0';
    } else {
        strncpy(gene_name, base, sizeof(gene_name)-1);
        gene_name[sizeof(gene_name)-1] = 0;
        char *dot = strchr(gene_name, '.');
        if (dot) *dot = '\0';
    };
    //printf("This is the sequence: %s\n", sequence);
    //printf("Thisis the base %s \n", base);
    
    // load tasks
    int n_tasks = 0;
    RegionTask *tasks = load_tsv_tasks(input_file, &n_tasks);
    if (n_tasks == 0) { fprintf(stderr, "No regions found in input TSV\n"); return 1; }

    struct timespec tot_start, tot_end;
    clock_gettime(CLOCK_MONOTONIC, &tot_start);

    TaskQueue queue;
    taskqueue_init(&queue, n_tasks);
    double per_core_timeout = 300.0;
    
    char summary_path[1024];
    snprintf(summary_path, sizeof(summary_path), "%s/a_summary.tsv", output_dir);
    FILE *summary_file = fopen(summary_path, "w");
    if (!summary_file) {
        perror("fopen summary.tsv");
        return 1;
    }
    fprintf(summary_file, "gene_name\tstart\tend\tregion_id\tearly_stopping\n");

    pthread_t threads[N_THREADS];
    pthread_attr_t attr[N_THREADS];
    size_t stack_size = 8 * 1024 * 1024;

    for (int i = 0; i < N_THREADS; i++) {
        pthread_attr_init(&attr[i]);
        pthread_attr_setstacksize(&attr[i], stack_size);
        pthread_create(&threads[i], &attr[i], worker_thread, &queue);
    }

    // enqueue all tasks
    for (int i = 0; i < n_tasks; i++) {
        Task t = { .task_index = i, .tasks = tasks,
                   .sequence = sequence,
                   .output_dir = output_dir,
                   .gene_name = gene_name,
                   .hd_percentage = hd_percentage,
                   .per_core_timeout_seconds = per_core_timeout,
                   .summary_file = summary_file};
        taskqueue_push(&queue, t);
    }

    // shutdown
    pthread_mutex_lock(&queue.lock);
    queue.shutdown = 1;
    pthread_cond_broadcast(&queue.cond);
    pthread_mutex_unlock(&queue.lock);

    for (int i = 0; i < N_THREADS; i++) {
        pthread_join(threads[i], NULL);
        pthread_attr_destroy(&attr[i]);
    }

    taskqueue_destroy(&queue);
    fclose(summary_file);
    /*
    for (int i = 0; i < n_tasks; ++i) {
        char outprefix[1024];
        snprintf(outprefix, sizeof(outprefix), "%s/region_%d", output_dir, i+1);
        process_row(tasks[i].start, tasks[i].end, i+1, hd_percentage, sequence, outprefix, 
                    gene_name, tasks[i].kmers, tasks[i].kmer_count);
    }*/

    clock_gettime(CLOCK_MONOTONIC, &tot_end);
    double total_elapsed = (tot_end.tv_sec - tot_start.tv_sec) + (tot_end.tv_nsec - tot_start.tv_nsec)/1e9;
    fprintf(stderr, "\n✅ All regions processed successfully in %.2f minutes (%.2f seconds).\n",
            total_elapsed/60.0, total_elapsed);

    // free tasks
    for (int i = 0; i < n_tasks; ++i) {
        for (int j = 0; j < tasks[i].kmer_count; ++j) free(tasks[i].kmers[j]);
        free(tasks[i].kmers);
    }
    free(tasks);
    free(sequence);
    return 0;
};
