#!/usr/bin/env python3
import click
import pandas as pd
from collections import defaultdict

@click.command()
@click.option(
    '-g', 
    '--gene-value', 
    default='gene', 
    show_default=True, 
    help='Feature type in the GTF file to consider as a gene (default: gene).'
)
@click.option(
    '-l', 
    '--length', 
    default=100000, 
    show_default=True, 
    type=int, 
    help='Length of each interval for gene density calculation (default: 100,000).'
)
@click.option(
    '-i', 
    '--input-file', 
    required=True, 
    type=click.Path(exists=True), 
    help='Path to the input GTF file.'
)
@click.option(
    '-o', 
    '--output-file', 
    required=True, 
    type=click.Path(), 
    help='Path to the output file where gene density results will be saved.'
)
@click.version_option(version='0.2.0', prog_name='Gene Density Calculator')
def gene_density(gene_value, length, input_file, output_file):
    """
    Calculate gene density in intervals of specified length from a GTF file.
    """
    gene_intervals = defaultdict(list)

    # Read the GTF file and filter rows where the feature matches the specified gene value
    with open(input_file, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue
            feature_type = columns[2]
            if feature_type == gene_value:
                chrom = columns[0]
                start = int(columns[3])
                end = int(columns[4])
                gene_intervals[chrom].append((start, end))

    # Process intervals to calculate gene density
    density_results = []
    for chrom, intervals in gene_intervals.items():
        # Determine the maximum position in the chromosome
        max_position = max([end for _, end in intervals])
        interval_counts = [0] * ((max_position // length) + 1)

        # Count the number of genes in each interval
        for start, end in intervals:
            start_idx = start // length
            end_idx = end // length
            for idx in range(start_idx, end_idx + 1):
                interval_counts[idx] += 1

        # Collect results for each interval
        for idx, count in enumerate(interval_counts):
            start_pos = idx * length
            end_pos = min((idx + 1) * length, max_position)
            density_results.append([chrom, start_pos, end_pos, count])

    # Save results to the output file without a header
    df = pd.DataFrame(density_results, columns=['Chromosome', 'Start', 'End', 'Gene_Count'])
    df.to_csv(output_file, sep='\t', index=False, header=False)

    click.echo(f'Gene density results saved to {output_file}')

if __name__ == '__main__':
    gene_density()

