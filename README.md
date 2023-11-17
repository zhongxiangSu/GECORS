# GECORS
## Get diploid genome cricle map 

## positional arguments:
```
  g1                    Path of target genome_1 fasta file
  g2                    Path of target genome_2 fasta file
  i1                    Path of target genome_1 index file, you can use : $samtools faidx genome_2
  i2                    Path of target genome_2 index file, you can use : $samtools faidx genome_2
  l1                    Path of target genome_1 Chr_name_list file
  l2                    Path of target genome_2 Chr_name_list file
  gff1                  Path of target genome_1 gff file
  gff2                  Path of target genome_2 gff file
  map                   Path of two genome mapper file, You need to install mummer yourself.
```
## options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Path of output file

### Chr_name_list file
```
Chr1
Chr2
Chr3
Chr4
Chr5
Chr6
Chr7
Chr8
Chr9
Chr10
Chr11
```

### get two genome mapper file
```
#run the following commands: 
$nucmer --minmatch 2000 --maxgap=80 --mincluster=2000 -p G1_VS_G2 genome_1 genome_2;
$show-coords -rcl G1_VS_G2.delta > G1_VS_G2.coords
```
