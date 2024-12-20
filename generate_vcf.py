
import argparse
import random

import pysam
import uuid

PREFIX = """
##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
"""

CHROMOSOMES = [
    ('chr1', 248956422),
    ('chr2', 242193529),
    ('chr3', 198295559),
    ('chr4', 190214555),
    ('chr5', 181538259),
    ('chr6', 170805979),
    ('chr7', 159345973),
    ('chr8', 145138636),
    ('chr9', 138394717),
    ('chr10', 133797422),
    ('chr11', 135086622),
    ('chr12', 133275309),
    ('chr13', 114364328),
    ('chr14', 107043718),
    ('chr15', 101991189),
    ('chr16', 90338345),
    ('chr17', 83257441),
    ('chr18', 80373285),
    ('chr19', 58617616),
    ('chr20', 64444167),
    ('chr21', 46709983),
    ('chr22', 50818468),
    ('chrX', 156040895),
    ('chrY', 57227415)
]

NUCLEOTIDES = ['A', 'T', 'C', 'G']

def generate_vcf_header():
    """ Generate VCF header.
    """
    header = pysam.VariantHeader()
    header.add_meta('fileformat', 'VCFv4.2')
    header.add_meta('source', 'SyntheticVCFGenerator')
    header.add_line('##contig=<ID=chr1,length=248956422>')
    header.add_line('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    
    sample_name = 'SAMPLE1'
    header.add_sample(sample_name)
    
    for chromosome, length in CHROMOSOMES:
        header.add_line(f'##contig=<ID={chromosome},length={length}>')

    return header

def generate_vcf_variants(num_variants):
    variants = []
    for _ in range(num_variants):
        chromosome, length = random.choice(CHROMOSOMES)
        # Generate a random position, reference allele, and alternative allele
        pos = random.randint(1, length)
        ref = random.choice(NUCLEOTIDES)
        alt = random.choice(NUCLEOTIDES)
        while alt == ref:  # Ensure ref and alt are not the same
            alt = random.choice(NUCLEOTIDES)

        # Add the variant information as a dictionary
        variant = {
            'contig': chromosome,
            'start': pos - 1,  # 0-based start position
            'stop': pos,       # 1-based end position
            'alleles': (ref, alt),
            'info': {'DP': random.randint(10, 100)},  # Example depth
            'samples': {'SAMPLE1': {'GT': (0, 1)}}  # Heterozygous genotype
        }
        variants.append(variant)

    chromosome_order = {chrom: i for i, (chrom, _) in enumerate(CHROMOSOMES)}
    variants.sort(key=lambda v: (chromosome_order[v['contig']], v['start']))
    return variants

def generate_vcf():
    """ Generate random VCF file.
    """
    random_id = str(uuid.uuid4())
    
    vcf_header = generate_vcf_header()
    vcf_variants = generate_vcf_variants(1000)
    
    with pysam.VariantFile(f'varsynth_{random_id}.vcf', 'w', header=vcf_header) as vcf_out:
        for variant in vcf_variants:
            record = vcf_out.new_record(
                contig=variant['contig'],
                start=variant['start'],
                stop=variant['stop'],
                alleles=variant['alleles'],
                info=variant['info']
            )
            record.samples['SAMPLE1']['GT'] = variant['samples']['SAMPLE1']['GT']
            vcf_out.write(record)
    

generate_vcf()

def main():
    parser = argparse.ArgumentParser(description="Generate a synthetic VCF file.")
    parser.add_argument(
        '-o', '--output', 
        type=str, 
        required=True, 
        help="Path to the output VCF file."
    )
    parser.add_argument(
        '-n', '--num_variants', 
        type=int, 
        default=100, 
        help="Number of variants to generate (default: 100)."
    )
    args = parser.parse_args()

    generate_vcf(args.output, args.num_variants)

if __name__ == "__main__":
    main()
