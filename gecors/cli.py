import argparse
from gecors.pipelines import gecors_main

def main():
    # argument parse
    parser = argparse.ArgumentParser(
        prog='GECORS', description='Get diploid genome cricle map'
    )

    parser.add_argument("g1",
                        help="Path of target genome_1 fasta file", type=str)
    parser.add_argument("g2",
                        help="Path of target genome_2 fasta file", type=str)
    parser.add_argument("i1", 
                        help="Path of target genome_1 index file, you can use : $samtools faidx genome_2 ", type=str)
    parser.add_argument("i2", 
                        help="Path of target genome_2 index file, you can use : $samtools faidx genome_2", type=str)
    parser.add_argument("l1", 
                        help="Path of target genome_1 Chr_name_list file", type=str)
    parser.add_argument("l2",   
                        help="Path of target genome_2 Chr_name_list file", type=str)
    parser.add_argument("gff1", 
                        help="Path of target genome_1 gff file", type=str)
    parser.add_argument("gff2", 
                        help="Path of target genome_2 gff file", type=str)
    parser.add_argument("map", 
                        help="Path of two genome mapper file, You need to install mummer yourself and run the following commands: $nucmer --minmatch 2000 --maxgap=80 --mincluster=2000 -p G1_VS_G2  genome_1  genome_2; 2: $show-coords -rcl G1_VS_G2.delta > G1_VS_G2.coords ", type=str)
    parser.add_argument("-o", "--output",
                        help="Path of output file", type=str)
                        
    args = parser.parse_args()
    
    gecors_main(args)

if __name__ == '__main__':
    main()