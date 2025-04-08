import argparse
import sys
import matplotlib.pyplot as plt
from statistics import mean
import matplotlib.patches as mpatches

def generate_contig_list(samtools_index_file):
    """Generate a list of contigs from samtools index file"""
    contig_list = []
    with open(samtools_index_file) as file:
        for line in file:
            fields = line.rstrip("\n").split("\t")
            contig_list.append(fields[0])
    return contig_list

def retrieve_busco_contigs(busco_full_table):
    """Returns a list of contigs, the number of contigs in the list indicates the number of BUSCO genes"""
    contigs_with_buscos = []
    with open(busco_full_table) as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                contigs_with_buscos.append(fields[2])
    return contigs_with_buscos

def generate_contig_dict(coverage_bed_file, contig_list, gc_content_file, samtools_index_file, busco_full_tables):
    """Generate a dict with keys as contigs and values as: mean cov, mean GC, contig size, and BUSCO counts"""
    contig_dict = {}

    # Mean coverage
    with open(coverage_bed_file) as file:
        contig = ""
        coverages = []
        for line in file:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in contig_list:
                if contig == fields[0]:
                    coverages.append(int(fields[3]))
                else:
                    if contig:
                        contig_dict[contig] = [mean(coverages)]
                    contig = fields[0]
                    coverages = [int(fields[3])]
        if contig and coverages:
            contig_dict[contig] = [mean(coverages)]

    # Add GC content
    with open(gc_content_file) as file:
        for line in file:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in contig_list and fields[0] in contig_dict:
                contig_dict[fields[0]].append(float(fields[1]))

    # Add contig size
    with open(samtools_index_file) as file:
        for line in file:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in contig_list and fields[0] in contig_dict:
                contig_dict[fields[0]].append(float(fields[1]))

    # Add BUSCO counts
    for table in busco_full_tables:
        busco_contig_list = retrieve_busco_contigs(table)
        for contig in contig_dict:
            contig_dict[contig].append(busco_contig_list.count(contig))

    return contig_dict

def main():
    print("""
        o               o               o
    o             o                               o
     ____        _     _     _      ____  _       _   
    | __ ) _   _| |__ | |__ | | ___|  _ \\| | ___ | |_  o
    |  _ \\| | | | '_ \\| '_ \\| |/ _ \\ |_) | |/ _ \\| __|
    | |_) | |_| | |_) | |_) | |  __/  __/| | (_) | |_ 
    |____/ \\__,_|_.__/|_.__/|_|\\___|_|   |_|\\___/ \\__| o
          o         o               o           o
                            o
    """)

    print("""
    Easily generate a GC content (%) vs Average Read Depth
    plot for your genome assembly.
          
    Please contact ak37@sanger.ac.uk for issues, questions 
    and requests.
    _______________________________________________________                
                    version 1.0.0
    """)

    parser = argparse.ArgumentParser()
    parser.add_argument("--samtools_index_file", metavar="FILE", required=True,
                        help="Samtools .fai file of fasta")
    parser.add_argument("--coverage_bed_file", metavar="FILE", required=True,
                        help="BED file with average read depth per contig")
    parser.add_argument("--gc_content_file", metavar="FILE", required=True,
                        help="GC content file generated using seqkit fx2tab -g -n")
    parser.add_argument("--busco_full_tables", metavar="LIST", required=True,
                        help="Comma-separated list of BUSCO full table files")
    parser.add_argument("--busco_labels", metavar="LIST", help="Comma-separated labels for BUSCO tables")
    parser.add_argument("--plot_title", metavar="STR", default=None)
    parser.add_argument("--output_prefix", metavar="STR", default="bubbleplot.asm")
    parser.add_argument("--min_gc", metavar="FLOAT", type=float, default=0)
    parser.add_argument("--max_gc", metavar="FLOAT", type=float, default=100)
    parser.add_argument("--min_cov", metavar="FLOAT", type=float, default=0)
    parser.add_argument("--max_cov", metavar="FLOAT", type=float, default=1000000)

    args = parser.parse_args()

    busco_full_tables = [table.strip() for table in args.busco_full_tables.split(",")]

    number_of_busco_full_tables = len(busco_full_tables)
    if number_of_busco_full_tables > 5:
        print("Error: Please provide no more than 5 BUSCO tables.")
        parser.print_help()
        sys.exit(1)

    taxon_labels = [taxon_label.strip() for taxon_label in args.busco_labels.split(",")]

    color_list = ["#66CC52", "#8F7DB3", "#CD527B", "#DE7B31", "#FFDE5A", "#629CF6"]

    contig_list = generate_contig_list(args.samtools_index_file)
    contig_dict = generate_contig_dict(
        coverage_bed_file=args.coverage_bed_file,
        contig_list=contig_list,
        gc_content_file=args.gc_content_file,
        samtools_index_file=args.samtools_index_file,
        busco_full_tables=busco_full_tables
    )

    output_file = "contig_list.tsv"
    
    plt.rcParams["figure.figsize"] = (10, 7)
    coverages = []

    for key, values in contig_dict.items():
        if len(values) < 3:
            continue
        cov, gc, size = values[:3]
        if args.min_cov <= cov <= args.max_cov and args.min_gc <= gc <= args.max_gc:
            coverages.append(cov)
            size_scaled = size / 2000
            if len(values) > 3 and max(values[3:]) > 0:
                busco_profile = max(values[3:])
                index = values.index(busco_profile) - 3
                color = color_list[index % len(color_list)]
                with open(output_file, "a+") as f:
                    output_line = key + "\t" + str(values[0]) + "\t" + str(values[1]) + "\t" + str(values[2]) + "\t" + str(taxon_labels[index])
                    f.write(output_line)
                    f.write("\n")
                f.close()
            elif max(values[3:]) == 0:
                color = "grey"
                with open(output_file, "a+") as f:
                    output_line = key + "\t" + str(values[0]) + "\t" + str(values[1]) + "\t" + str(values[2]) + "\t" + "Unidentified"
                    f.write(output_line)
                    f.write("\n")
                f.close()
            plt.scatter(gc, cov, s=size_scaled, alpha=0.5, color=color, edgecolors="black")


    plt.xlim(0, 100)
    if coverages:
        cov_range = max(coverages) - min(coverages)
        if cov_range < 100:
            plt.ylim(-10, 110)
        else:
            allowance = cov_range * 0.1
            plt.ylim(min(coverages) - allowance, max(coverages) + allowance)

    plt.xlabel("\nGC content (%)", size=15)
    plt.ylabel("Average read depth\n", size=15)

    if args.plot_title:
        plt.title(f"{args.plot_title}\n", size=20)

    legend_boxes = []
    for taxon in range(0, number_of_busco_full_tables):
        patch = mpatches.Patch(color=color_list[taxon], label=taxon_labels[taxon])
        legend_boxes.append(patch)

    patch = mpatches.Patch(color="grey", label="Unidentified")
    legend_boxes.append(patch)

    plt.legend(handles=legend_boxes)

    png_file = f"{args.output_prefix}.png"
    svg_file = f"{args.output_prefix}.svg"
    plt.savefig(png_file, dpi=500, format="png")
    plt.savefig(svg_file, dpi=500, format="svg")

if __name__ == "__main__":
    main()
