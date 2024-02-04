#include "utils.h"

int bin_main(int argc, char *argv[]){
    
    double sway     = 0.5;
    double coverage = 0.90;

    double strain   = 99.0;
    double species  = 97.0;
    double genus    = 95.0;
    int    num      = 5;

    double level_t[] = {
          0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
          0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,  
          0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
          0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
          0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
          0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
          0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
          0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
    };

   level_t['p'] = species;
   level_t['e'] = genus;
   level_t['t'] = strain;

   int c;
   while ( (c = getopt(argc, argv, "n:w:g:s:t:c:") ) >= 0) {
        if(c == 'n') num = atoi(optarg);
        else if(c == 'w') sway     = atof(optarg);
        else if(c == 'g') genus    = atof(optarg);
        else if(c == 's') species  = atof(optarg);
        else if(c == 't') strain   = atof(optarg);
        else if(c == 'c') coverage = atof(optarg);
    }

    if ( optind == argc || argc != optind + 1) {
        
        fprintf(stderr, "\nUsage: amplicon-tk bin [option] <tsv>\n\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -n INT   number of hits to print, default: [5]\n");
        fprintf(stderr, "  -w FOALT percent sway value:, default: [0.5]\n");
        fprintf(stderr, "           reference: github.com/Joseph7e/Assign-Taxonomy-with-BLAST\n");
        fprintf(stderr, "  -g FOALT genus level identity threshold, default: [95.0]\n");
        fprintf(stderr, "  -s FOALT species level identity threshold, default: [97.0]\n");
        fprintf(stderr, "  -t FOALT strain level identity threshold, default: [99.0]\n\n");

        fprintf(stderr, " Note:\n");
        fprintf(stderr, " usearch usearch_global userfields format:\n");
        fprintf(stderr, " \"query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+qcov\"\n\n");

        exit(1);
    
    }

    gzFile     fp;
    fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
    if (fp) {

        kstream_t *ks;
        ks = ks_init(fp);

        kstring_t kt   = {0, 0, 0};
        kstring_t bin  = {0, 0, 0};
        kstring_t kq   = {0, 0, 0};

        int *fields, n;
        kstring_t level = {0, 0, 0};
        kstring_t rep = {0, 0, 0};
        double prev = 0;
        int val  = 0;
        int ties = 0;

        while( ks_getuntil( ks, '\n', &kt, 0) >=  0){
            
            if(kt.l == 0 || kt.s[0] == '#') continue;
            
            fields = ksplit(&kt, '\t', &n);
            double identity = atof(kt.s + fields[2]);

            if( ( identity < genus) || (atof(kt.s + fields[12]) < coverage) ) continue;
            
            if( kq.s != NULL && strcmp(kt.s, kq.s) == 0){
                
                if( (prev - identity) <= sway && val < num && identity >= level_t[ (int)level.s[1] ]){
                   kputc(';',  &bin);
                   kputs(kt.s + fields[1],  &bin);
                   ++val;
                }

                if(prev == identity) ties++;
  
            }else{

                if(kq.s != NULL) printf("%s\t%s\t%s\t%.1lf\t%s\t%d\n", kq.s, bin.s, rep.s, prev, level.s, ties);
                
                level.l = 0;
                if(identity >= strain) kputs("strain", &level);
                else if(identity >= species) kputs("species", &level);
                else kputs("genus", &level);

                ties  = 1;
                rep.l = 0;
                kputs(kt.s + fields[1], &rep);

                val   = 1;
                kq.l  = 0;
                kputs(kt.s, &kq);
                bin.l = 0;
                kputs(kt.s + fields[1], &bin);
                prev = identity;

            }

        }

        /*last query*/
        if(kq.s != NULL) printf("%s\t%s\t%s\t%.1lf\t%s\t%d\n", kq.s, bin.s, rep.s, prev, level.s, ties);

        free(kt.s);
        free(kq.s);
        free(bin.s);
        free(level.s);
        free(rep.s);

        ks_destroy(ks);
        gzclose(fp);
    }else{
        fprintf(stderr, "[ERR]: can't open file %s .\n", argv[optind]);
        exit(1);
    }

    return 0;
}