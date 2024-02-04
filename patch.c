#include "utils.h"

int patch_main(int argc, char *argv[]){
    
    if ( optind == argc || argc != optind + 1 ) {
        
        fprintf(stderr, "\nUsage: amplicon-tk patch <tsv>\n\n");
        return 1;
    
    }

    kstring_t kt    = {0, 0, 0};
    
    int *fields,n, i;
    int *items, m;

    kstream_t *ks;
    gzFile fp;
    fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r"); 
    if (fp) {
        
        ks  = ks_init(fp);
        
        kstring_t kswap  = {0, 0, 0};
        kstring_t revise  = {0, 0, 0};

       
        while( ks_getuntil( ks, '\n', &kt, 0) >=  0){
            
            kswap.l = 0;
            kputs(kt.s, &kswap);
            fields = ksplit(&kt, '\t', &n);

            if(n == 6){
                if(strcmp(kt.s + fields[3], "genus") == 0  && strcmp(kt.s + fields[2], "species") == 0){
                    
                    revise.l = 0;
                    kputs(kt.s + fields[4], &revise);

                    for (i = revise.l - 1; i > 0; --i){                        
                        if(revise.s[i]  == ','){
                            revise.s[i] = '\0';
                            revise.l    = i + 1;
                            break;
                        }
                    }
                    printf("%s\t%s\n", kt.s, revise.s);
                }else printf("%s\t%s\n", kt.s, kt.s + fields[4]);

            }else{
                if(strcmp(kt.s + fields[6], "genus") == 0  && strcmp(kt.s + fields[1], "species") == 0){

                    revise.l = 0;
                    kputs(kt.s + fields[2], &revise);
                    items = ksplit(&revise, ' ', &m);
                    
                    fputs(kt.s, stdout);
                    if(strcmp(revise.s, "Candidatus") == 0) printf("\tgenus\t%s %s", revise.s, revise.s + items[1]);
                    else printf("\tgenus\t%s", revise.s);
                    printf("\t%s\t%s\t%s\t%s\t%s\n", kt.s + fields[3], kt.s + fields[4], kt.s + fields[5], kt.s + fields[6], kt.s + fields[7]);

                }else puts(kswap.s);
            }
        }
        
        free(kswap.s);
        free(kt.s);
        free(revise.s);

        ks_destroy(ks);
        gzclose(fp);

    }else{
        fprintf(stderr, "[ERR]: can't open file %s\n", argv[optind + 1]);
        exit(1);
    }

    return 0;
}