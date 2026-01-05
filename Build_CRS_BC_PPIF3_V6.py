import pandas as pd
import gzip
from tqdm import tqdm
import ahocorasick  # Aho-Corasick for faster multiple substring searching

def process_fastq(input_file, variant_df, output_file, collapsed_output_file):
    """Process the FASTQ file and write both non-collapsed and collapsed results to output files."""

    # Create a dictionary for quick lookup of MappingSequence to VariantID
    mapping_dict = dict(zip(variant_df['MappingSequence'], variant_df['VariantID']))

    # Initialize Aho-Corasick automaton for fast multiple sequence matching
    A = ahocorasick.Automaton()
    for mapping_seq in mapping_dict:
        A.add_word(mapping_seq, mapping_seq)
    A.make_automaton()

    # Buffer for results
    results = []

    with gzip.open(input_file, 'rt') as infile:
        total_lines = sum(1 for _ in infile)  # Count total lines in file
        infile.seek(0)  # Rewind file to start reading

        for _ in tqdm(range(total_lines // 4), desc="Processing reads"):
            header = infile.readline().strip()
            if not header:
                break
            sequence = infile.readline().strip()
            infile.readline()  # plus line
            infile.readline()  # quality line

            # Skip if sequence length is outside the range
            if len(sequence) < 200 or len(sequence) > 1000:
                continue

            found_match = False
            for end_idx, mapping_seq in A.iter(sequence):
                mapping_length = len(mapping_seq)
                match_start = end_idx - mapping_length + 1
                match_end = end_idx + 1
                
                # Ensure match is full length and exact
                if match_end <= len(sequence):
                    barcode = sequence[match_end:match_end + 15]  # Extract barcode (15 bp immediately after mapping sequence)
                    
                    # Only process if a full-length barcode is found
                    if len(barcode) == 15:
                        variant_id = mapping_dict[mapping_seq]
                        results.append([header, mapping_seq, barcode, variant_id])
                        found_match = True
                        break  # Stop after the first match

    # Process the results DataFrame
    if results:
        output_df = pd.DataFrame(results, columns=['ReadID', 'MappingSequence', 'Barcode', 'VariantID'])
        
        # Debugging: Print header of the non-collapsed DataFrame
        print("Non-collapsed DataFrame header:")
        print(output_df.head())

        # Save non-collapsed DataFrame to file
        output_df.to_csv(output_file, sep='\t', index=False, mode='w', header=True)
        
        # Group by 'Barcode' and 'VariantID', and collapse rows
        collapsed_df = output_df.groupby(['Barcode', 'VariantID']).size().reset_index(name='count')

        # Debugging: Print header of the collapsed DataFrame
        print("Collapsed DataFrame header:")
        print(collapsed_df.head())

        # Save collapsed DataFrame to file
        collapsed_df.to_csv(collapsed_output_file, sep='\t', index=False, mode='w', header=True)
        
        # Debugging: Print file write confirmation
        print(f"Writing non-collapsed results to: {output_file}")
        print(f"Total non-collapsed rows written: {len(output_df)}")
        
        print(f"Writing collapsed results to: {collapsed_output_file}")
        print(f"Total collapsed rows written: {len(collapsed_df)}")

        print(f"Files '{output_file}' and '{collapsed_output_file}' created successfully.")

# Parameters
input_fastq = "fastq/MONTGOMERY_CUSTOM_508/MONTGOMERY_508_3.fastq.gz"
non_collapsed_output_file = "PPIF3_results/PPIF3_CRS-BC_Library_V6.txt"
collapsed_output_file = "PPIF3_results/PPIF3_CRS-BC_Library_V6_collapsed.txt"
variant_list_path = "design_files/V6_Promoter_Edits_OligoList_240920.tsv"

# Load Variant List DataFrame
try:
    variant_df = pd.read_csv(variant_list_path, delimiter='\t')
    # Process the FASTQ file and save both non-collapsed and collapsed versions
    process_fastq(input_fastq, variant_df, non_collapsed_output_file, collapsed_output_file)
except Exception as e:
    print(f"An error occurred: {e}")
