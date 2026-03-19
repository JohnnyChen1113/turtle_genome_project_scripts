#!/usr/bin/env python3
import click
import re

@click.command()
@click.option('-i', '--input_file', required=True, type=click.Path(exists=True), help='Input FASTA file path.')
@click.option('-o', '--output_file', required=True, type=click.Path(), help='Output GC content file path.')
@click.option('-l', '--length', default=10000, help='Length of each interval for calculating GC content.', type=int)
def calculate_gc_content(input_file, output_file, length):
    """
    Calculate GC content in intervals from the given FASTA file.
    """
    def gc_content(sequence):
        gc_count = sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
        return gc_count / len(sequence) if len(sequence) > 0 else 0

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_sequence = []
        current_chrom = None
        for line in infile:
            if line.startswith('>'):
                if current_chrom is not None:
                    # Process the previous chromosome
                    sequence = ''.join(current_sequence)
                    for i in range(0, len(sequence), length):
                        segment = sequence[i:i+length]
                        gc_percent = gc_content(segment) * 100
                        outfile.write(f'{current_chrom}\t{i+1}\t{i+len(segment)}\t{gc_percent:.2f}\n')
                # Start a new chromosome
                current_chrom = line[1:].strip().split()[0]  # Extract chromosome name
                current_sequence = []
            else:
                current_sequence.append(line.strip())

        # Process the last chromosome
        if current_chrom is not None:
            sequence = ''.join(current_sequence)
            for i in range(0, len(sequence), length):
                segment = sequence[i:i+length]
                gc_percent = gc_content(segment) * 100
                outfile.write(f'{current_chrom}\t{i+1}\t{i+len(segment)}\t{gc_percent:.2f}\n')

    click.echo(f'GC content results saved to {output_file}')

if __name__ == '__main__':
    calculate_gc_content()

