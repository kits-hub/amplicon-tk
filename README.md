###  amplicon-tk: tools for amplicon data analysis

<hr>

#### 1. install

```sh
git clone https://github.com/kits-hub/amplicon-tk
cd amplicon-tk
make
```

#### 2. interface

current: `version：0.0.2`


```text
Usage:   amplicon-tk <command> <arguments>
Version: 0.0.2

Command:
  -- LCA method.
     bin      bin same query sequence hits.
     mapping  mapping sequence id to taxon id.
     lca      taxon assignment using LCA strategy.

  -- Voting method.
     collapse collapse same query sequence hits.
     voting   taxon assignment using voting strategy.

  -- auxiliary utils.
     patch    patch the taxonomy level.
     level    bin node to latest major node.
     uniques  find unique reads from fasta/q.


Licenced:
(c) 2021-2024 - LEI ZHANG
Logic Informatics Co.,Ltd.
zhanglei@logicinformatics.com

```
