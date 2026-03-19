#!/usr/bin/env python3
import re
import click

@click.command()
@click.option('-i', '--input_file', required=True, type=click.Path(exists=True), help='Input FASTA file path.')
@click.option('-o', '--output_file', required=True, type=click.Path(), help='Output FASTA file path.')
def rename_fasta_headers(input_file, output_file):
    """
    Rename the headers of a FASTA file to a simpler format like chr1, chr2, etc.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        chrom_counter = 1
        for line in infile:
            if line.startswith('>'):
                # Extract chromosome number using regex and create new header
                new_header = f'>chr{chrom_counter}\n'
                chrom_counter += 1
                outfile.write(new_header)
            else:
                outfile.write(line)

    click.echo(f'Headers renamed and output saved to {output_file}')

if __name__ == '__main__':
    rename_fasta_headers()

