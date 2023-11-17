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

### Chr_name_list_1 file
```
Chr1_1
Chr2_1
Chr3_1
Chr4_1
Chr5_1
Chr6_1
Chr7_1
Chr8_1
Chr9_1
Chr10_1
Chr11_1
```

### Chr_name_list_2 file
```
Chr1_2
Chr2_2
Chr3_2
Chr4_2
Chr5_2
Chr6_2
Chr7_2
Chr8_2
Chr9_2
Chr10_2
Chr11_2
```

### get two genome mapper file
```
#run the following commands: 
$nucmer --minmatch 2000 --maxgap=80 --mincluster=2000 -p G1_VS_G2 genome_1 genome_2;
$show-coords -rcl G1_VS_G2.delta > G1_VS_G2.coords
```
