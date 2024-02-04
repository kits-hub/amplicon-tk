#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>

#include "kstring.h"
#include "khash.h"
#include "kseq.h"

typedef struct{
    int   parent;
    char *rank;
    char *name;
} node_t;


KHASH_MAP_INIT_INT(int64, int64_t *)
void kh_int64_destroy(khash_t(int64) *h);

KHASH_MAP_INIT_STR(reg, char *)
void kh_reg_destroy(khash_t(reg) *h);

KHASH_MAP_INIT_INT(int32, char *)
void kh_int32_destroy(khash_t(int32) *h);

KHASH_SET_INIT_STR(set)
void kh_set_destroy(khash_t(set) *h);

KHASH_MAP_INIT_INT(int, int)
KHASH_SET_INIT_INT(set32)

KHASH_MAP_INIT_INT(traversal, node_t *);

KSEQ_INIT(gzFile, gzread)
void print_kseq(kseq_t *s, FILE *file);
void rev_com(kstring_t  *ks);

void kh_taxon_init(khash_t(reg) *h);
void kh_traversal_destroy(khash_t(traversal) *h);