/*****************************************************************************
  amplicon-tk.c

  amplicon-tk command line interface.  

  (c) 2021-2024 - LEI ZHANG
  Logic Informatics Co.,Ltd.
  zhanglei@logicinformatics.com
  
  Licenced under The MIT License.
******************************************************************************/

#include <stdio.h>
#include <string.h>

int bin_main(int argc, char *argv[]);
int mapping_main(int argc, char *argv[]);
int lca_main(int argc, char *argv[]);

int collapse_main(int argc, char *argv[]);
int voting_main(int argc, char *argv[]);

int patch_main(int argc, char *argv[]);
int level_main(int argc, char *argv[]);

int uniques_main(int argc, char *argv[]);

static int usage(){
    
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   amplicon-tk <command> <arguments>\n");
    fprintf(stderr, "Version: 0.0.2\n\n");

    fprintf(stderr, "Command:\n");
    fprintf(stderr, "  -- LCA method.\n");
  
    fprintf(stderr, "     bin      bin same query sequence hits.\n");
    fprintf(stderr, "     mapping  mapping sequence id to taxon id.\n");
    fprintf(stderr, "     lca      taxon assignment using LCA strategy.\n\n");

    fprintf(stderr, "  -- Voting method.\n");
    fprintf(stderr, "     collapse collapse same query sequence hits.\n");
    fprintf(stderr, "     voting   taxon assignment using voting strategy.\n\n");

    fprintf(stderr, "  -- auxiliary utils.\n");
    fprintf(stderr, "     patch    patch the taxonomy level.\n");
    fprintf(stderr, "     level    bin node to latest major node.\n");
    fprintf(stderr, "     uniques  find unique reads from fasta/q.\n");
 
    fprintf(stderr, "\n");

    fprintf(stderr, "\nLicenced:\n");
    fprintf(stderr, "(c) 2021-2024 - LEI ZHANG\n");
    fprintf(stderr, "Logic Informatics Co.,Ltd.\n");
    fprintf(stderr, "zhanglei@logicinformatics.com\n");
    fprintf(stderr, "\n");

    return 1;

}

int main(int argc, char *argv[]){

    if (argc < 2) return usage();
    
    if (strcmp(argv[1], "bin") == 0) bin_main(argc - 1, argv + 1);
    else if (strcmp(argv[1], "mapping") == 0) mapping_main(argc - 1, argv + 1);
    else if (strcmp(argv[1], "lca") == 0) lca_main(argc - 1, argv + 1);
    else if (strcmp(argv[1], "collapse") == 0) collapse_main(argc - 1, argv + 1);
    else if (strcmp(argv[1], "voting") == 0) voting_main(argc - 1, argv + 1);
    else if (strcmp(argv[1], "patch") == 0) patch_main(argc - 1, argv + 1);
    else if (strcmp(argv[1], "level") == 0) level_main(argc - 1, argv + 1);
    else if (strcmp(argv[1], "uniques") == 0) uniques_main(argc - 1, argv + 1);
    else {
        fprintf(stderr, "[main] unrecognized command '%s'. Abort!\n", argv[1]);
        return 1;
    }
    return 0;

}