
## Amino Acid Mutation Identifier and Counter

This Python script is designed for the identification and counting of amino acid mutations from .fasta sequences. It utilizes the BioPython library to interpret sequences.

## Features

- Identification of nucleotide mutations between reference and query sequences.
- Translation of nucleotide changes into amino acid mutations.
- Counting the frequency of specific amino acid mutations across multiple sequences.
- Export of mutation data and summaries into CSV files for further analysis.

## Requirements

To run this script, you will need the following libraries:

- Pandas
- BioPython
- diff_match_patch
- Collections


## Usage

1. Place your reference and query FASTA files in the designated directory.
2. Ensure the `referencepath` variable in the script is set to the location of your reference files.
3. Run the script
4. The script will process each file, identify mutations, and output the results in CSV format in the specified output directory.

## Input File Format

The script expects FASTA files for both reference and query sequences. Ensure that there are no spaces in the FASTA file names, as they are processed to remove spaces for consistency.

## Output

The script generates several CSV files containing detailed mutation data, including:

- A summary of amino acid changes.
- A count of specific amino acid mutations.
- Mutation details for each sequence analyzed.

## Customization

You can customize the script to fit your specific needs, such as adjusting the maximum length deviation for sequences or modifying the output file paths.

## License

This project is open-source and available under the MIT License.
