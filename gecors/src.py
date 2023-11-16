import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from matplotlib.path import Path
import matplotlib.patches as patches

class GenomePainter:
    def __init__(self, index_1_file, index_2_file, Chr_name_list_1, Chr_name_list_2,  genome_file_1, genome_file_2, gff_file_1, gff_file_2, genome_mapping_file, **kwargs):
        self.index_1_dict = self.read_index_to_dict(index_1_file)
        self.index_2_dict = self.read_index_to_dict(index_2_file)
        self.differece_len_multiple = sum(self.index_1_dict.values()) / sum(self.index_2_dict.values())
        self.index_2_dict = {k: v * self.differece_len_multiple for k, v in self.index_2_dict.items()}
        self.all_chr_index_dict = self.merge_two_dict(self.index_1_dict, self.index_2_dict)
        self.Chr_name_list_1 = Chr_name_list_1
        self.Chr_name_list_2 = Chr_name_list_2
        self.chr_list_name = Chr_name_list_1 + Chr_name_list_2
        self.chr_list_len = [self.all_chr_index_dict[n] for n in self.chr_list_name]
        self.genome_file_1 = genome_file_1
        self.genome_file_2 = genome_file_2
        self.chr_theta_positions = {}
        self.gff_file_1 = gff_file_1
        self.gff_file_2 = gff_file_2
        self.genome_mapping_file = genome_mapping_file 
        #default parameters
        self.Polar_coordinates_starting_direction = kwargs.get('Polar_coordinates_starting_direction', 'N')
        self.figure_size = kwargs.get('figure_size', (18, 18))
        self.color_line = kwargs.get('color_line', 'RdYlGn')
        self.R = kwargs.get('R', 1)
        self.gap_tuble = kwargs.get('gap_tuble', (0.75, 0.25))
        self.barh_edgecolors = kwargs.get('barh_edgecolors', 'black')
        self.barh_height = kwargs.get('barh_height', 0.02)
        self.label_radiu = kwargs.get('label_radiu', 0.1)
        self.scale_chr = kwargs.get('scale_chr', 1000000)
        #self.fig = plt.figure(figsize=(12, 12))
        #self.ax = self.fig.add_subplot(111, polar=True)

    
    @staticmethod
    def read_index_to_dict(index_file):
        index_dict = {} 
        with open(index_file, 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                index_dict[line[0]] = int(line[1])
        return index_dict
    
    @staticmethod
    def read_fasta_to_dict(fasta_file):
        output_dict = {}  

        for record in SeqIO.parse(fasta_file, "fasta"):
            output_dict[record.id] = str(record.seq)
        return output_dict

    
    @staticmethod
    def merge_two_dict(dict1, dict2):
 
        merged_dict = dict1.copy()  
        merged_dict.update(dict2)  
        return merged_dict

    @staticmethod
    def get_Chr_GC_ratio_dict(fasta_file, static_len=1000):
        Chr_dict = read_fasta_to_dict(fasta_file)
        Genome_ratio_dict = {}
        for chr in Chr_dict:
            Chr_GC_ratio_dict = {}
            seq =  Chr_dict[chr]
            for spilt_len in range(0, len(seq), static_len):
                spilt_seq = seq[spilt_len:spilt_len+static_len]
                spilt_sit = (spilt_len , spilt_len+static_len)
                GC_count = spilt_seq.count('G') + spilt_seq.count('C')
                GC_ratio = GC_count / static_len
                Chr_GC_ratio_dict[spilt_sit] = GC_ratio
            Genome_ratio_dict[chr] = Chr_GC_ratio_dict
        return Genome_ratio_dict
    
    @staticmethod
    def get_cds_sit(gff_file, genome_file1, static_len = 100000 ):
        in_handle = open(gff_file,'r')
        A37_cds_sit = {}
        CDS_density_dict = {}
        CDS_sort_dict = {}
        for rec in GFF.parse(in_handle):
            A37_cds_sit[rec.id] = {}
            for gene in rec.features:
                if gene.type == "gene":
                    for mRNA in gene.sub_features:
                        if mRNA.type == 'mRNA' and '.1' in mRNA.id:
                            for structure in mRNA.sub_features:
                                if structure.type == 'CDS':
                                    A37_cds_sit[rec.id][structure.location.start.position] = structure.location.end.position                   
        in_handle.close()

        for Chr in A37_cds_sit:
            CDS_sort_dict[Chr] = {}
            start_list = sorted(A37_cds_sit[Chr].keys())
            for start in start_list:
                CDS_sort_dict[Chr][start] = A37_cds_sit[Chr][start] - start   
        
        Chr_dict = read_fasta_to_dict(genome_file1)
        for chr in Chr_dict:
            Chr_GC_ratio_dict = {}
            seq =  Chr_dict[chr]
            for split_len in range(0, len(seq), static_len):
                split_start = split_len
                split_end = split_len+static_len
                split_cds = []
                for i in CDS_sort_dict[chr].keys():
                    if split_start < i < split_end:
                        split_cds.append(CDS_sort_dict[chr][i])
                Chr_GC_ratio_dict[(split_start, split_end)] = sum(split_cds) / static_len
            CDS_density_dict[chr] = Chr_GC_ratio_dict

        return CDS_density_dict 
    
    @staticmethod
    def get_two_genome_mappinginfo(genome_mapping_file, start=4):   
        two_genome_map_info_dict = {}        
        with open (genome_mapping_file, 'r') as f:
            lines = f.readlines()
            num = 1 
        for line in lines[start:]:
            line = line.strip()
            line = line.split('\t')
            two_genome_map_info_dict[num] = {'Genome1':[line[11],line[0], line[1]], 'Genome2':[line[12],line[2], line[3]]}
            num += 1
        return two_genome_map_info_dict


    def generate_category_colors(self):
        category_colors_origin = plt.get_cmap(self.color_line)(np.linspace(0.15, 0.85, len(self.Chr_name_list_1)))
        category_colors = {}
        Chr_name_list_2_color = self.Chr_name_list_2[::-1]
        for i, chr_name in enumerate(self.Chr_name_list_1):
            category_colors[chr_name] = category_colors_origin[i]
            category_colors[Chr_name_list_2_color[i]] = category_colors_origin[i]
        return category_colors
    

    
    def paint_genome_cricos(self):
        category_colors = self.generate_category_colors()        
        self.fig, self.ax = plt.subplots(figsize=self.figure_size)
        self.ax = plt.subplot(projection='polar')
        self.ax.set_theta_zero_location(self.Polar_coordinates_starting_direction)
        total_width = 2 * np.pi
        num_bars = len(self.all_chr_index_dict)
        cumulative_list = np.cumsum(self.chr_list_len).tolist()
        sit = 2 * np.pi * self.gap_tuble[1] / num_bars / 2
        for i, chr_name in enumerate(self.chr_list_name):           
            tick_bottom = self.R           
            if chr_name in self.Chr_name_list_1:
                length = self.all_chr_index_dict[chr_name]
                bar_width = length / sum(self.chr_list_len) * total_width * self.gap_tuble[0]
                cum_length = cumulative_list[i]
                starts = (cum_length / sum(self.chr_list_len) * total_width) - bar_width - sit
                colors = category_colors[chr_name]
                self.ax.barh(self.R, bar_width, left=starts, height=self.barh_height, color=colors, edgecolor=self.barh_edgecolors)
                label_angle = starts + (bar_width / 2)
                label_radius = self.R + self.label_radiu
                label_text = chr_name
                self.ax.text(label_angle, label_radius, label_text, ha='center', va='center', rotation=label_angle * 180 / np.pi)                
                tick_theta_positions_list = []
                x = 0
                chr_scale = int(length // self.scale_chr)
                self.chr_theta_positions[chr_name] = [starts, starts + bar_width]
                for a in range(chr_scale + 1):
                    tick_theta_positions_list.append((bar_width / length * self.scale_chr * a) + starts)
                
                for tick_theta in tick_theta_positions_list:
                    self.ax.plot([tick_theta, tick_theta], [tick_bottom + 0.01, tick_bottom + 0.015], color="k", linewidth=1)
                    self.ax.text(tick_theta, tick_bottom + 0.03, x, ha='center', va='center', rotation=label_angle * 180 / np.pi, size=6)
                    x += 1
            
            if chr_name in self.Chr_name_list_2:
                length_2 = self.all_chr_index_dict[chr_name]
                bar_width_2 = length_2 / sum(self.chr_list_len) * total_width * self.gap_tuble[0]
                cum_length = cumulative_list[i]
                starts = (cum_length / sum(self.chr_list_len) * total_width) - bar_width_2 - sit + bar_width_2
                colors = category_colors[chr_name]
                self.ax.barh(self.R, -bar_width_2, left=starts, height=self.barh_height, color=colors, edgecolor='black')
                label_angle = starts - (bar_width_2 / 2)
                label_radius = self.R + self.label_radiu
                label_text = chr_name
                self.ax.text(label_angle, label_radius, label_text, ha='center', va='center', rotation=label_angle * 180 / np.pi)
                
                tick_theta_positions_list = []
                x = 0
                chr_scale = int((length_2 / self.differece_len_multiple) // self.scale_chr)
                self.chr_theta_positions[chr_name] = [starts, starts - bar_width_2]
                for a in range(chr_scale + 1):
                    tick_theta_positions_list.append(starts - (bar_width_2 / length_2 * self.scale_chr * a) * self.differece_len_multiple)
                
                for tick_theta in tick_theta_positions_list:
                    self.ax.plot([tick_theta, tick_theta], [tick_bottom + 0.01, tick_bottom + 0.015], color="k", linewidth=1)
                    self.ax.text(tick_theta, tick_bottom + 0.03, x, ha='center', va='center', rotation=label_angle * 180 / np.pi, size=6)
                    self.ax.set_axis_off()
                    x += 1
        
        return self.fig, self.ax, self.chr_theta_positions
    
    def paint_GC_ratio(self):
        GC_cricos = self.R - 0.1
        genome_ratio_dict_1 = self.get_Chr_GC_ratio_dict(self.genome_file_1)
        genome_ratio_dict_2 = self.get_Chr_GC_ratio_dict(self.genome_file_2)
        genome_ratio_dict = merge_two_dict(genome_ratio_dict_1, genome_ratio_dict_2)

        for chr_name in self.chr_list_name:
            gc_length = []
            gc_ratio = []
            chr_ratio_dict = genome_ratio_dict[chr_name]
            if chr_name in self.Chr_name_list_1:               
                sit = (self.chr_theta_positions[chr_name][1] - self.chr_theta_positions[chr_name][0]) / len(chr_ratio_dict)
                gc_length.append(self.chr_theta_positions[chr_name][0])              
                for i in range(1, len(chr_ratio_dict)):
                    gc_length.append(gc_length[i-1] + sit)              
                for ratio in chr_ratio_dict.values():
                    gc_ratio.append((ratio/20) + GC_cricos)                
                x = np.linspace(min(gc_length), max(gc_length), 100)               
                y1 = np.full_like(x, GC_cricos)
                y2 = np.interp(x, gc_length, gc_ratio)
                self.ax.fill_between(x, y1, y2,  color='blue', alpha = 0.3)
            if chr_name in self.Chr_name_list_2:               
                sit = (self.chr_theta_positions[chr_name][0] - self.chr_theta_positions[chr_name][1]) / len(chr_ratio_dict)
                gc_length.append(self.chr_theta_positions[chr_name][1])               
                for i in range(1, len(chr_ratio_dict)):
                    gc_length.append(gc_length[i-1] + sit)
                for ratio in chr_ratio_dict.values():
                    gc_ratio.append((ratio/20) + GC_cricos)                
                x = np.linspace(min(gc_length), max(gc_length), 100)
                y1 = np.full_like(x, GC_cricos)
                y2 = np.interp(x, gc_length, gc_ratio)
                self.ax.fill_between(x, y1, y2,  color='blue', alpha = 0.3)

        self.ax.plot([0, 0], [1/20 + GC_cricos , GC_cricos], color="k", linewidth=1)

        tick_position = (1/20 + GC_cricos + GC_cricos) / 2
        self.ax.plot([ 0, 0.005], [1/20 + GC_cricos, 1/20 + GC_cricos], color="k", linewidth=1)
        self.ax.plot([ 0, 0.005], [tick_position, tick_position], color="k", linewidth=1)
        self.ax.plot([ 0, 0.005], [GC_cricos, GC_cricos], color="k", linewidth=1)
        self.ax.text(0.03, tick_position, '0.5', ha='left', va='center', fontsize=6)
  
    def paint_CDS_denstiy(self):

        CDS_cricos = self.R - 0.18
        genome_cds_dentisy_dict_1 = self.get_cds_sit(self.gff_file_1, self.genome_file_1)
        genome_cds_dentisy_dict_2 = self.get_cds_sit(self.gff_file_2, self.genome_file_2)
        genome_cds_dentisy_dict = merge_two_dict(genome_cds_dentisy_dict_1, genome_cds_dentisy_dict_2)

        for chr_name in self.chr_list_name:     
            cds_length = []
            cds_ratio = []
            chr_ratio_dict = genome_cds_dentisy_dict[chr_name]
            if chr_name in self.Chr_name_list_1:               
                sit = (self.chr_theta_positions[chr_name][1] - self.chr_theta_positions[chr_name][0]) / len(chr_ratio_dict)
                cds_length.append(self.chr_theta_positions[chr_name][0])               
                for i in range(1, len(chr_ratio_dict)):
                    cds_length.append(cds_length[i-1] + sit)            
                for ratio in chr_ratio_dict.values():
                    cds_ratio.append((ratio/10) + CDS_cricos)
                x = np.linspace(min(cds_length), max(cds_length), 100)
                y1 = np.full_like(x, CDS_cricos)
                y2 = np.interp(x, cds_length, cds_ratio)
                self.ax.fill_between(x, y1, y2,  color='gray', alpha = 0.5)
            if chr_name in self.Chr_name_list_2:     
                sit = (self.chr_theta_positions[chr_name][0] - self.chr_theta_positions[chr_name][1]) / len(chr_ratio_dict)
                cds_length.append(self.chr_theta_positions[chr_name][1])    
                for i in range(1, len(chr_ratio_dict)):
                    cds_length.append(cds_length[i-1] + sit)         
                for ratio in chr_ratio_dict.values():
                    cds_ratio.append((ratio/10) + CDS_cricos)
                x = np.linspace(min(cds_length), max(cds_length), 100)
                y1 = np.full_like(x, CDS_cricos)
                y2 = np.interp(x, cds_length, cds_ratio)
                self.ax.fill_between(x, y1, y2,  color='gray', alpha = 0.5)

    def paint_genome_mapping(self):
        line_wight = 1.5
        alphas = 0.8
        two_genome_mapper_dict = self.get_two_genome_mappinginfo(self.genome_mapping_file)
        category_colors = self.generate_category_colors()
        for region in two_genome_mapper_dict:
            chr1 = two_genome_mapper_dict[region]['Genome2'][0]
            start1 = two_genome_mapper_dict[region]['Genome2'][1]
            end1 = two_genome_mapper_dict[region]['Genome2'][2]
            chr2 = two_genome_mapper_dict[region]['Genome1'][0]
            start2 = two_genome_mapper_dict[region]['Genome1'][1]
            end2 = two_genome_mapper_dict[region]['Genome1'][2]
            angle1 = self.chr_theta_positions[chr1][0] + (self.chr_theta_positions[chr1][1] - self.chr_theta_positions[chr1][0]) * (int(start1) / self.all_chr_index_dict[chr1])
            angle2 = self.chr_theta_positions[chr2][1] - (self.chr_theta_positions[chr2][1] - self.chr_theta_positions[chr2][0]) * (int(start2) / self.all_chr_index_dict[chr2])  * self.differece_len_multiple
            r1 = self.R - 0.2
            r2 = self.R - 0.2
    

            if angle1 < np.pi:
                angle1_sit1 = (angle1, r1)
                angle2_sit2 = (angle2, r2)               
                sit1 =  np.pi/4 + np.pi/4 
                sit2 =  2 * np.pi * 6/8 
                ctrl1 = (sit1, 0.225)
                ctrl2 = (sit2, 0.225)
            else:
                angle1_sit1 = (angle1, r1)
                angle2_sit2 = (angle2, r2)
                sit1 =  np.pi * 3/4 
                sit2 =  2 * np.pi * 5/8
                ctrl1 = (sit1, 0.425)
                ctrl2 = (sit2, 0.425)
            verts = [angle1_sit1, ctrl1, ctrl2, angle2_sit2]
            codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
            path_line = Path(verts, codes)
            patch_line = patches.PathPatch(path_line, facecolor='none',ec = category_colors[chr1],lw= line_wight,alpha = alphas)
            self.ax.add_patch(patch_line)
            #self.ax.plot([angle1, angle2], [r1, r2], color='gray',alpha = 0.3)


def get_Chr_list():
    with open ('Chr_name_list.txt', 'r') as f:
        Chr_name_list = f.read().splitlines()
    return Chr_name_list