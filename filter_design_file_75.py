import pandas as pd
import argparse

# Main function to handle command-line inputs and process the design file
def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Process design file by filtering counts and handling duplicate barcodes.")
    parser.add_argument("-d", "--design_file", required=True, help="Path to the design file.")
    parser.add_argument("-o", "--output_file", required=True, help="Output file to save the results.")
    
    args = parser.parse_args()

    # Load the design file into a DataFrame
    df = pd.read_csv(args.design_file, sep='\t')

    # Drop rows where count is less than 25
    rows_before_count_filter = df.shape[0]
    df = df[df['count'] >= 5]
    rows_after_count_filter = df.shape[0]
    print(f"Initial number of rows before count filter: {rows_before_count_filter}")
    print(f"Initial number of rows after count filter: {rows_after_count_filter}")

    # Handle duplicate barcodes based on the 75% rule
    rows_before_dup_filter = df.shape[0]
    print(f"Initial number of rows before duplicate filter: {rows_before_dup_filter}")

    # Group by 'Barcode' and calculate the total count per barcode
    barcode_groups = df.groupby('Barcode')

    # Define a function to process each group of barcodes
    def filter_barcode_group(group):
        total_count = group['count'].sum()
        group['count_fraction'] = group['count'] / total_count  # Calculate the fraction of total counts for each VariantID

        # Keep only the rows where a VariantID has at least 75% of the total reads for that barcode
        return group[group['count_fraction'] >= 0.75]

    # Apply the filter function to each barcode group
    df_filtered = barcode_groups.apply(filter_barcode_group)

    # Reset index after groupby
    df_filtered = df_filtered.reset_index(drop=True)

    # Drop the temporary 'count_fraction' column
    df_filtered = df_filtered.drop(columns=['count_fraction'])

    # Print the final number of rows after filtering
    rows_after_dup_filter = df_filtered.shape[0]
    print(f"Number of rows after filtering: {rows_after_dup_filter}")

    # Export the updated design file
    df_filtered.to_csv(args.output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
