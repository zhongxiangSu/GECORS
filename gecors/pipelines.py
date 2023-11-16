from gecors.src import *
def gecors_main(args):
    Chr_name_list_1 = get_Chr_list(args.Chr_name_list_1)
    Chr_name_list_2 = get_Chr_list(args.Chr_name_list_2)
    painter = GenomePainter(args.index_1_file, args.index_2_file, Chr_name_list_1, Chr_name_list_2, args.genome_file1, args.genome_file2, args.gff_file_1, args.gff_file_2, args.two_genome_mapper_file)
    painter.paint_genome_cricos()
    painter.paint_GC_ratio()
    painter.paint_CDS_denstiy()
    painter.paint_genome_mapping()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    pass