#include "utils.h"
#include <limits.h>

typedef struct{
    int  counts;
    int  rank;
    int  left;
    int  right;
    int  parent;
    kstring_t nodes;
} Node;

KHASH_MAP_INIT_INT(lca, Node *)

void kh_lca_destroy(khash_t(lca) *h);
void tree_label(khash_t(lca) *h, khash_t(int) *map);
int  lca(int max, int min, khash_t(lca) *h, khash_t(int) *map);
int  binning (kstring_t *kt, khash_t(lca) *h, khash_t(int) *map);

static int init  = 1;
static int value = 1;

int lca_main(int argc, char *argv[]){
    
    if ( optind == argc || argc != optind + 2) {

        fprintf(stderr, "\nUsage: amplicon-tk lca <taxon.map> <tsv>\n\n");
        return 1;
    }

    kstring_t  kt  = {0, 0, 0};

    int *fields, n, ret;

    khash_t(lca) *h;
    h = kh_init(lca);
    khint_t k, p;

    khash_t(int32) *db;
    db = kh_init(int32);
    kstring_t  kswab = {0, 0, 0};

    gzFile fp;
    fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r"); 
    kstream_t *ks;

    if (fp) {

       int node1, node2;
       ks  = ks_init(fp);
       while(ks_getuntil( ks, '\n', &kt, 0) >= 0 ){
            
            fields = ksplit(&kt, '\t', &n); 
            node1  = atoi( kt.s);
            node2  = atoi( kt.s + fields[1] );

            kswab.l = 0;
            kputs(kt.s + fields[2], &kswab);
            kputc('\t', &kswab);
            kputs(kt.s + fields[3], &kswab);

            k = kh_put(int32, db, node1, &ret);
            if(ret){
               kh_key(db, k) = node1;
               kh_val(db, k) = strdup(kswab.s);
            }

            k = kh_put(lca, h, node1, &ret);
            if(ret){
                Node  *node   = (Node *) calloc(1, sizeof(Node) );
                if( node != NULL ){

                    node->counts    =  0;
                    node->rank      =  0;
                    node->left      = -1;
                    node->right     = -1;
                    node->parent    = node2;
                    
                    kstring_t ktmp  = {0, 0, 0};
                    kputs(kt.s, &ktmp);
                    node->nodes     = ktmp;

                }else{
                    fprintf(stderr, "Node memory allocation failed!\n");
                    exit(-1);
                }

                kh_key(h, k)  = node1;
                kh_val(h, k)  = node;

            }

            if(node1 == 1) continue;

            p = kh_get(lca, h, node2);
            if(p == kh_end(h)){
                p = kh_put(lca, h, node2, &ret);       
                Node  *node   = (Node *) calloc(1, sizeof(Node) );
                if( node != NULL ){
                    node->counts    =  0;
                    node->rank      =  0;
                    node->left      = -1;
                    node->right     = -1;
                    node->parent    = -1;
                    
                    kstring_t ktmp  = {0, 0, 0};
                    kputs(kt.s + fields[1], &ktmp);
                    kputc('\t', &ktmp);
                    kputs(kt.s, &ktmp);
                    node->nodes     = ktmp;
                }else{
                    fprintf(stderr, "Node memory allocation failed!\n");
                    exit(-1);
                }

                kh_key(h, p)  = node2;
                kh_val(h, p)  = node;

            }else{
                kputc('\t', &kh_val(h, p)->nodes);
                kputs(kt.s, &kh_val(h, p)->nodes);
            }

            kh_val(h, p)->counts++;
            kh_val(h, k)->rank   = kh_val(h, p)->counts;
            kh_val(h, k)->parent = kh_key(h, p);

       }
       ks_destroy(ks);
       gzclose(fp);
    
    }else{
        fprintf(stderr, "[ERR]: can't open file %s\n", argv[optind]);
        exit(1);
    }

    khash_t(int) *map;
    map = kh_init(int);
    tree_label(h, map);

    fp = strcmp(argv[optind + 1], "-")? gzopen(argv[optind + 1], "r") : gzdopen(fileno(stdin), "r");

    if(fp){
        
       ks  = ks_init(fp);
       
       while( ks_getuntil( ks, '\n', &kt, 0) >=  0){       
            
            if(kt.l == 0) continue;
            fields  = ksplit(&kt, '\t', &n);
            
            int id = 0;

            if(atoi(kt.s + fields[2]) == 1) id = atoi(kt.s + fields[1]);
            else{
              kswab.l = 0;
              kputs(kt.s + fields[1], &kswab);
              id = binning(&kswab, h, map);
            }
            
            khint_t k = kh_get(int32, db, id);
            if(k == kh_end(db)){
               fprintf(stderr, "[ERR]: node %d can't map to a taxon id.\n", id);
               exit(-1);
            }

            printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", kt.s, id, kh_val(db, k), kt.s + fields[1], kt.s + fields[2], kt.s + fields[3], kt.s + fields[4], kt.s + fields[5], kt.s + fields[6], kt.s + fields[7]);
        
        }
    
    }else{
        fprintf(stderr, "[ERR]: can't open file %s\n", argv[optind + 1]);
        exit(1);
    }


    kh_destroy(int, map);
    kh_int32_destroy(db);
    kh_lca_destroy(h);
    return 0;
}

void kh_lca_destroy(khash_t(lca) *h){

    khint_t k;
    if (h == 0) return;
    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)){
            Node *node = kh_val(h, k);
            free(node->nodes.s);
            free(node);
        }
    }
    kh_destroy(lca, h);
}

void tree_label(khash_t(lca) *h, khash_t(int) *map){

     khint_t k, t;
     int ret;
     k = kh_get(lca, h, init);
     if( k != kh_end(h) ){
         
         kh_val(h, k)->left = value++ ;

         int *fields, n, i;
         fields = ksplit(&kh_val(h, k)->nodes, '\t', &n);
         for (i = 1; i < n; ++i){
              init  = atoi(kh_val(h, k)->nodes.s + fields[i]);
              tree_label(h, map);
         }

         kh_val(h, k)->right = value++;
         
         t = kh_put(int, map, kh_val(h, k)->right, &ret);
         if( ret ){
             kh_key(map, t) = kh_val(h, k)->right;
             kh_val(map, t) = kh_val(h, k)->parent;
         }

     }else{
        fprintf(stderr, "Node %d is obsolete! PASS...\n", init);
        exit(-1);
     }

}

int lca(int max, int min, khash_t(lca) *h, khash_t(int) *map){

    khint_t k, t;
    
    while(1){
        
        k = kh_get(int, map, min);
        if(k != kh_end(map)){

            t = kh_get(lca, h, kh_val(map, k) );
            if( t != kh_end(h) ){
                min = kh_val(h, t)->right;
                if( min >= max )
                    return kh_key(h, t);
            }else{
               fprintf(stderr, "Node %d map error!\n", kh_val(map, k) );
               return  -1;
            }
        
        }else{
            fprintf(stderr, "Index %d error!\n", min);
            return  -1;
        }

    }

}

int binning (kstring_t *kt, khash_t(lca) *h, khash_t(int) *map){

    int *fields, n, i;
    fields = ksplit(kt, ',', &n);

    if(n == 1) return atoi(kt->s);

    khint_t k;
    int max, min;
    max = -1;
    min = INT_MAX;

    for (i = 0; i < n; ++i){
        
        k = kh_get(lca, h, atoi(kt->s + fields[i]) );
        if( k != kh_end(h) ){
            if( max < kh_val(h, k)->right )  max = kh_val(h, k)->right;
            if( min > kh_val(h, k)->right )  min = kh_val(h, k)->right;
        }
    }

    if(max == min){
        return atoi(kt->s + fields[1]);
    }else{
        return lca(max, min, h, map);
    }

}