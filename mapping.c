#include "utils.h"

int mapping_main(int argc, char *argv[]){
    
    if ( optind == argc || argc != optind + 2 ) {
        
        fprintf(stderr, "\nUsage: amplicon-tk mapping <db> <tsv>\n\n");
        return 1;
    
    }

    kstring_t kt    = {0, 0, 0};
    int *fields, i, ret, n;

    khash_t(set) *maps;
    maps = kh_init(set);

    khash_t(reg) *h;
    h = kh_init(reg);
    khint_t k;


    kstream_t *ks;
    gzFile fp;
    fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r"); 
    int n_col = 0;

    if(fp){
        
        ks = ks_init(fp);
        while( ks_getuntil( ks, '\n', &kt, 0) >=  0){

            fields = ksplit(&kt, '\t', &n_col);
            k = kh_put(reg, h, kt.s, &ret);
            if(ret){
                kh_key(h, k)   = strdup(kt.s);
                kh_value(h, k) = strdup(kt.s + fields[1]);
            }
        
        }

        ks_destroy(ks);
        gzclose(fp);
    
    }else{
        fprintf(stderr, "[ERR]: can't open file %s\n", argv[optind]);
        exit(1);
    }
    
    fp = strcmp(argv[optind + 1], "-")? gzopen(argv[optind + 1], "r") : gzdopen(fileno(stdin), "r"); 
    if (fp) {
        
        ks  = ks_init(fp);
        kstring_t kswap = {0, 0, 0};
        kstring_t knode = {0, 0, 0};

        int *items, m;
       
        while( ks_getuntil( ks, '\n', &kt, 0) >=  0){
            
            fields = ksplit(&kt, '\t', &n);

            kswap.l= 0;
            kputs(kt.s + fields[1], &kswap);
            items  = ksplit(&kswap, ';', &m);

            knode.l = 0;
            
            for (i = 0; i < m; ++i){

                k = kh_get(reg, h, kswap.s + items[i]);
                if( k != kh_end(h) ){

                    khint_t s = kh_put(set, maps, kh_val(h, k), &ret);
                    if(ret) kh_key(maps, s) =  kh_val(h, k);

                }else{
                    fprintf(stderr, "[ERR]: item: %s has no taxon id.\n",  kswap.s + items[i]);
                    exit(-1);
                }

            }

            int num = 0;
            knode.l = 0;
            for (k = 0; k < kh_end(maps); ++k) {
                if (kh_exist(maps, k)) {
                    kputc(',' , &knode);
                    kputs(kh_key(maps, k) , &knode);
                    ++num;
                    kh_del(set, maps, k);
                }
            }

            k = kh_get(reg, h, kt.s + fields[2]);
            if( k == kh_end(h) ){
                fprintf(stderr, "[ERR]: item: %s has no taxon id.\n",  kt.s + fields[2]);
                exit(-1);
            }

            printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n", kt.s, knode.s + 1, num, kt.s + fields[2], kh_val(h, k), kt.s + fields[3], kt.s + fields[4], kt.s + fields[5]);
         
        }
        
        free(kswap.s);
        free(knode.s);
        free(kt.s);

        ks_destroy(ks);
        gzclose(fp);

    }else{
        fprintf(stderr, "[ERR]: can't open file %s\n", argv[optind + 1]);
        exit(1);
    }

    kh_reg_destroy( h );
    kh_destroy(set, maps);

    return 0;
}