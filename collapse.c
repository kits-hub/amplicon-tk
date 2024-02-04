#include "utils.h"

static unsigned int hierarchy_t[] = {
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,  
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
      0,   0,   0,   2,   0,   0,   4,   5,   0,   0,  0,  0,  0,  0,  0,  3,
      1,   0,   0,   6,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,
};

static khash_t(set) *maps;

static int nlevel = 6;

int levels (kstring_t *kswap, kstring_t *knode, kstring_t *flatten){

  int *fields, *items, i, j, n;
  fields   = ksplit(kswap, ',' , &n);
  if(n == 7) return 0;
  items    = ksplit(knode, ',', &n);

  int current, pad;
  current =  pad = 0;
  flatten->l = 0;

  for (i = 0; i < n; ++i){

     current = hierarchy_t[ (unsigned int)(kswap->s + fields[i])[0] ];

     if(current > pad){
        for (j = 0; j < current - pad; ++j) kputs(",0", flatten);
     }

     kputc(',', flatten);
     kputs(knode->s + items[i], flatten);
     pad = current + 1;

  }

  if(nlevel > current){
    for (j = 0; j < nlevel - current; ++j) kputs(",0", flatten);  
  }

  return 1;

}

int collapse_main(int argc, char *argv[]){

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
        fprintf(stderr, "  -n INT number of hits to print, default: [5]\n");
        fprintf(stderr, "  -w FOALT percent sway value:, default: [0.5]\n");
        fprintf(stderr, "           reference: github.com/Joseph7e/Assign-Taxonomy-with-BLAST\n");
        fprintf(stderr, "  -g FOALT genus level identity threshold, default: [95.0]\n");
        fprintf(stderr, "  -s FOALT species level identity threshold, default: [97.0]\n");
        fprintf(stderr, "  -t FOALT strain level identity threshold, default: [99.0]\n\n");

        fprintf(stderr, "Note:user defined fields format.\n\n");

        exit(1);
    
    }

    maps = kh_init(set);

    gzFile     fp;
    fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
    if (fp) {

        kstream_t *ks;
        ks = ks_init(fp);

        kstring_t kt      = {0, 0, 0};
        kstring_t bin     = {0, 0, 0};
        kstring_t kq      = {0, 0, 0};

        int *fields, n, f;
        kstring_t level   = {0, 0, 0};
        kstring_t rep     = {0, 0, 0};
        kstring_t knode   = {0, 0, 0};
        kstring_t flatten = {0, 0, 0};
        kstring_t kswap   = {0, 0, 0};

        int    taxon      = 0;
        int    val        = 0;
        double prev       = 0;
        int    ties       = 0;

        while( ks_getuntil( ks, '\n', &kt, 0) >=  0){
            
            if(kt.l == 0 || kt.s[0] == '#') continue;
            
            fields = ksplit(&kt, '\t', &n);
            double identity = atof(kt.s + fields[5]);

            if( ( identity < genus) || (atof(kt.s + fields[6]) < coverage) ) continue;
            
            kswap.l = knode.l = flatten.l = 0;
            kputs(kt.s + fields[3], &kswap);
            kputs(kt.s + fields[4], &knode);

            if( kq.s != NULL && strcmp(kt.s, kq.s) == 0){
                
                if( (prev - identity) <= sway && val < num && identity >= level_t[ (int)level.s[1] ]){
                  
                   f =  levels(&kswap, &knode, &flatten);
                   kputc(';',  &bin);
                   (f == 0) ? kputs(kt.s + fields[4],  &bin) : kputs(flatten.s + 1,  &bin);
                   ++val;
                }

                if(prev == identity) ties++;

            }else{

                if(kq.s != NULL) printf("%s\t%s\t%d\t%s\t%d\t%.1lf\t%s\t%d\n", kq.s, bin.s, val, rep.s, taxon, prev, level.s, ties);

                level.l = 0;
                if(identity >= strain) kputs("strain", &level);
                else if(identity >= species) kputs("species", &level);
                else kputs("genus", &level);

                rep.l = 0;
                kputs(kt.s + fields[1], &rep);

                val   = 1;
                kq.l  = 0;
                kputs(kt.s, &kq);
                
                ties  = 1;
                /*inital bin val*/
                bin.l = 0;
                f =  levels(&kswap, &knode, &flatten);
                (f == 0) ? kputs(kt.s + fields[4],  &bin) : kputs(flatten.s + 1,  &bin);
                
                prev = identity;
                taxon= atoi(kt.s + fields[2]);

            }

        }

        if(kq.s != NULL) printf("%s\t%s\t%d\t%s\t%d\t%.1lf\t%s\t%d\n", kq.s, bin.s, val, rep.s, taxon, prev, level.s, ties);

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

    kh_set_destroy(maps);
    return 0;
}