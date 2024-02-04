#include "utils.h"
#include "ksort.h"
#include "kvec.h"

static khash_t(int)     *maps;

#define khint_lt(a, b) ( kh_val(maps, (a) )  > kh_val(maps, (b)) )
KSORT_INIT(khint, khint_t, khint_lt)

int voting(kstring_t *kswab){
    
    int *fields, i, j, n, ret;
    for (i = 0; i < kswab->l; ++i){
        if(kswab->s[i] == ';') kswab->s[i] = ',';
    }

    fields   = ksplit(kswab, ',', &n);
    int num  = n/7; 

    khint_t k;
    for (i = 6; i >= 0; --i){
        
        khint_t bucket[num];
        memset(bucket, 0, num);
        int size = 0;

        for (j = 0; j < num; ++j){
            int node =  atoi(kswab->s + fields[i + j * 7]);
            if(node == 0) continue;

            k = kh_put(int, maps, node, &ret);
            kh_key(maps, k) = node;
            kh_val(maps, k) = 0;

        }

        for (j = 0; j < num; ++j){
            int node =  atoi(kswab->s + fields[i + j * 7]);
            if(node == 0) continue;
            k = kh_get(int, maps, node);
            if( k != kh_end(maps) ){
                if(kh_val(maps, k) == 0) bucket[size++] = k;
                kh_val(maps, k)++;
            }
        }

        if(size == 0) continue;    
        if(size == 1) return kh_key(maps, bucket[0]);
        else{
            ks_mergesort(khint, size, bucket, 0);
            if( kh_val(maps, bucket[0]) > kh_val(maps, bucket[1]) )
                return kh_key(maps, bucket[0]);
        }

    }
    
    return 1;

}

int voting_main(int argc, char *argv[]){
    
    if ( optind == argc || argc != optind + 2 ) {
        
        fprintf(stderr, "\nUsage: amplicon-tk voting <taxon.map> <tsv>\n\n");
        return 1;
    
    }

    kstring_t kt    = {0, 0, 0};
    kstring_t kswab = {0, 0, 0};

    int *fields, ret, n;
    khint_t k;

    khash_t(int32) *db;
    db = kh_init(int32);

    gzFile fp;
    fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r"); 
    kstream_t *ks;

    if (fp) {

       int node;
       ks  = ks_init(fp);
       while(ks_getuntil( ks, '\n', &kt, 0) >= 0 ){
            
          fields = ksplit(&kt, '\t', &n); 
          node   = atoi( kt.s);

          kswab.l = 0;
          kputs(kt.s + fields[2], &kswab);
          kputc('\t', &kswab);
          kputs(kt.s + fields[3], &kswab);
        
          k = kh_put(int32, db, node, &ret);
          if(ret){
            kh_key(db, k) = node;
            kh_val(db, k) = strdup(kswab.s);
          }


       }
       ks_destroy(ks);
       gzclose(fp);
    
    }else{
        fprintf(stderr, "[ERR]: can't open file %s\n", argv[optind]);
        exit(1);
    }

    
    maps = kh_init(int);

    fp = strcmp(argv[optind + 1], "-")? gzopen(argv[optind + 1], "r") : gzdopen(fileno(stdin), "r"); 
    if (fp) {
        
        ks  = ks_init(fp);       
        while( ks_getuntil( ks, '\n', &kt, 0) >=  0){
            
            if(kt.l == 0) continue;
            fields  = ksplit(&kt, '\t', &n);

            int id = 0;

            if(atoi(kt.s + fields[2]) == 1) id = atoi(kt.s + fields[4]);
            else{
              kswab.l = 0;
              kputs(kt.s + fields[1], &kswab);
              id = voting( &kswab );
            }

            k = kh_get(int32, db, id);
            if(k == kh_end(db)){
               fprintf(stderr, "[ERR]: node %d can't map to a taxon id.\n", id);
               exit(-1);
            }

            printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", kt.s, id, kh_val(db, k), kt.s + fields[1], kt.s + fields[2], kt.s + fields[3], kt.s + fields[4], kt.s + fields[5], kt.s + fields[6],  kt.s + fields[7]);
        
        }
        
        ks_destroy(ks);
        gzclose(fp);

    }else{
        fprintf(stderr, "[ERR]: can't open file %s\n", argv[optind + 1]);
        exit(1);
    }

    free(kt.s);
    free(kswab.s);
    kh_destroy(int, maps);

    return 0;
}