#include "utils.h"

static unsigned int level_tab[] = {
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,  
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
      0,   0,   0,   6,   0,   0,   4,   3,   0,   0,  0,  8,  0,  0,  0,  5,
      7,   0,   0,   2,   1,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
};

int level_traversal(const khash_t(traversal) *h, const khash_t(reg) *map, int taxon);
int level_main(int argc, char *argv[]){

    if ( optind == argc || argc != optind + 2) {
        fprintf(stderr, "\nUsage: amplicon-tk level [options] <taxon.map> <tab>\n\n");
        return 1;
    }
 
    khash_t(traversal) *h;
    h = kh_init(traversal);

    khash_t(reg) *map;
    map = kh_init(reg);
    kh_taxon_init(map);

    khint_t    k;
    kstream_t *ks;

    int *fields, n;
    kstring_t kt    = {0, 0, 0};
    
    gzFile     fp;
    fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r"); 
    if (fp) {
       
       ks = ks_init(fp);  
       int ret, node1, node2;

       while( ks_getuntil( ks, '\n', &kt, 0) >=  0){
        
            fields = ksplit(&kt, '\t', &n);
            
            node1  = atoi( kt.s + fields[0] );
            node2  = atoi( kt.s + fields[1] );

            k = kh_put(traversal, h, node1, &ret);
            if(ret){

                kh_key(h, k) = node1;
                node_t *node   = (node_t *)malloc( sizeof(node_t) );
                
                if( node != NULL ){
                    node->parent  = node2;
                    
                    if( (kt.s + fields[2])[0] == 's'){
                        switch( (kt.s + fields[2])[2] ){
                            case 'b': 
                                        node->rank    = strdup( "typing" );
                                        break;
                            case 'r':    
                                        node->rank    = strdup( "typing" );
                                        break;
                            case 'e':    
                                        node->rank    = strdup( "species" );
                                        break;
                            case 'p': 
                                        node->rank    = strdup( "kingdom" );
                                        break;
                            default:
                                        node->rank    = strdup( kt.s + fields[2] );
                                        break;
                        }
                    }else
                        node->rank    = strdup( kt.s + fields[2] );
                    node->name    = strdup( kt.s + fields[3] );
                    kh_val(h, k) = node;
                }
            }
        }
        ks_destroy(ks);
        gzclose(fp);
    }else{
        fprintf(stderr, "[ERR]: %s :file streaming error\n", argv[optind]);
    }


    fp = strcmp(argv[optind + 1], "-")? gzopen(argv[optind + 1], "r") : gzdopen(fileno(stdin), "r");  
    if (fp) {
       
       ks = ks_init(fp);
       int i, taxon, node;

       while( ks_getuntil( ks, '\n', &kt, 0) >=  0){            
            
            if(kt.s[0] == '#'){
                puts(kt.s);
                continue;
            }

            fields = ksplit(&kt, '\t', &n);
            taxon = atoi( kt.s + fields[1] );
            node = level_traversal(h, map, taxon);
            
            khint_t k;
            k  = kh_get(traversal, h, node); 
            
            /*SRR12647582.16 29466   species Veillonella parvula 29466   1   CP019721.1:1196476-1199317  29466   98.0    strain*/

            if(k != kh_end(h) ) printf("%s\t%d\t%s\t%s", kt.s, kh_key(h, k) , kh_val(h, k)->rank, kh_val(h, k)->name);
            else printf("%s\t-\t-\t-", kt.s);

            for (i = 4; i < n; ++i) printf("\t%s",  kt.s + fields[i]);

            fputc('\n', stdout);
            
        
        }

        ks_destroy(ks);
        gzclose(fp);
    }else{        
        fprintf(stderr, "[ERR]: %s :file streaming error\n", argv[optind + 1]);       
    }

    free(kt.s);

    kh_reg_destroy( map );
    kh_traversal_destroy( h );

    return 0;
}

int level_traversal(const khash_t(traversal) *h, const khash_t(reg) *map, int taxon){
    
    khint_t k;
    k  = kh_get(traversal, h, taxon); 
    if(k == kh_end(h) ) return 0;

    while(1){
        

        if(level_tab[ (int)kh_val(h, k)->rank[0] ] != 0)
            return kh_key(h, k);
        
        if(kh_key(h, k) == 1) return 1; 
        
        k = kh_get(traversal, h, kh_val(h, k)->parent);
        if( k == kh_end( h ) ) return 0;

    }
}