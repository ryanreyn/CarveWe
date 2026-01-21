#!/usr/bin/env python3
"""
This script consolidates individual per-genome sensitivity prediction files
generated using batch mode CarveWe into one composite long-format CSV.

Input format (per genome):
    - Wide format with arrays: each row is a metabolite category,
      each column is a genome, values are arrays of sensitivity scores

Output format (consolidated):
    - Long format: genome_id, model_number, and columns for each of the
      13 metabolite categories

Usage:
    python consolidate-batch-predictions.py -i output_dir/ -o final_matrix.csv
"""

import pandas as pd
import argparse
from pathlib import Path
import glob
import ast
import sys


def parse_sensitivity_file(filepath):
    """
    Parse a single sensitivity file and convert to long format.

    Transforms from wide format with arrays:
        Category         Genome_ID
        Amino Acids      [0.99, 0.98, 0.97, ...]
        Carbohydrates    [0.5, 0.6, 0.7, ...]

    To long format:
        genome_id  model_number  Amino Acids  Carbohydrates  ...
        Genome_ID  0            0.99         0.5            ...
        Genome_ID  1            0.98         0.6            ...

    Args:
        filepath: Path to individual genome sensitivity CSV

    Returns:
        DataFrame in long format (genome_id, model_number, category columns)

    Raises:
        ValueError: If file format is invalid or arrays have inconsistent lengths
    """
    try:
        # Read the file
        df = pd.read_csv(filepath, index_col=0)

        # Extract genome ID from column name (first and only data column after index)
        if len(df.columns) == 0:
            raise ValueError(f"No data columns found in {filepath}")

        genome_id = df.columns[0]

        # Initialize storage for transformed data
        rows = []

        # Determine number of models from first array
        first_category = df.index[0]
        first_value = df.loc[first_category, genome_id]

        # Parse the array (stored as string representation of list)
        if isinstance(first_value, str):
            array_vals = ast.literal_eval(first_value)
        else:
            # Already a list/array
            array_vals = first_value if hasattr(first_value, '__iter__') else [first_value]

        n_models = len(array_vals)

        # For each model number
        for model_num in range(n_models):
            row = {
                'genome_id': genome_id,
                'model_number': model_num
            }

            # For each metabolite category
            for category in df.index:
                value = df.loc[category, genome_id]

                # Parse the array
                if isinstance(value, str):
                    array_vals = ast.literal_eval(value)
                elif hasattr(value, '__iter__'):
                    array_vals = value
                else:
                    array_vals = [value]

                # Validate array length
                if len(array_vals) != n_models:
                    raise ValueError(
                        f"Inconsistent array length in {filepath}: "
                        f"category '{category}' has {len(array_vals)} values, "
                        f"expected {n_models}"
                    )

                # Extract value for this model
                row[category] = array_vals[model_num]

            rows.append(row)

        # Convert to DataFrame
        result_df = pd.DataFrame(rows)

        return result_df

    except Exception as e:
        raise ValueError(f"Error parsing {filepath}: {str(e)}")


def consolidate_sensitivity_files(input_dir, output_file, pattern='*_sensitivity.csv',
                                   batch_size=None, verbose=True):
    """
    Consolidate all sensitivity files into one long-format CSV.

    Args:
        input_dir: Directory containing input sensitivity files
        output_file: Path for consolidated output CSV
        pattern: Glob pattern for matching input files (default: *_sensitivity.csv)
        batch_size: If specified, process in batches to save memory (default: None)
        verbose: Print progress messages (default: True)

    Returns:
        int: Number of genomes successfully processed
    """
    input_dir = Path(input_dir)

    # Find all files matching pattern
    search_pattern = str(input_dir / pattern)
    files = sorted(glob.glob(search_pattern))

    if not files:
        print(f"ERROR: No files found matching pattern: {search_pattern}", file=sys.stderr)
        return 0

    if verbose:
        print(f"Found {len(files)} sensitivity files to consolidate")

    # Choose processing method based on batch_size
    if batch_size:
        return consolidate_in_batches(files, output_file, batch_size, verbose)
    else:
        return consolidate_all_at_once(files, output_file, verbose)


def consolidate_all_at_once(files, output_file, verbose=True):
    """
    Process all files at once (simpler, faster for moderate datasets).

    Args:
        files: List of file paths to process
        output_file: Output file path
        verbose: Print progress messages

    Returns:
        int: Number of genomes successfully processed
    """
    all_dfs = []
    errors = []

    for i, filepath in enumerate(files, 1):
        try:
            df = parse_sensitivity_file(filepath)
            all_dfs.append(df)

            if verbose and i % 50 == 0:
                print(f"  Processed {i}/{len(files)} files...")

        except Exception as e:
            errors.append((filepath, str(e)))
            if verbose:
                print(f"  WARNING: Skipping {Path(filepath).name}: {e}", file=sys.stderr)
            continue

    if not all_dfs:
        print("ERROR: No files were successfully processed", file=sys.stderr)
        return 0

    # Concatenate all genomes
    if verbose:
        print("Concatenating all data...")

    final_df = pd.concat(all_dfs, ignore_index=True)

    # Sort for cleaner output
    final_df = final_df.sort_values(['genome_id', 'model_number']).reset_index(drop=True)

    # Save
    final_df.to_csv(output_file, index=False)

    # Print summary
    if verbose:
        print(f"\n✓ Successfully consolidated {len(all_dfs)} genomes → {output_file}")
        print(f"  Total rows: {len(final_df):,}")
        print(f"  Columns: {len(final_df.columns)}")
        print(f"  Unique genomes: {final_df['genome_id'].nunique()}")

        # Calculate average models per genome
        models_per_genome = final_df.groupby('genome_id')['model_number'].max() + 1
        print(f"  Average models per genome: {models_per_genome.mean():.1f}")

        if errors:
            print(f"\n⚠ {len(errors)} files had errors and were skipped")

    return len(all_dfs)


def consolidate_in_batches(files, output_file, batch_size, verbose=True):
    """
    Process files in batches to conserve memory (for very large datasets).

    Args:
        files: List of file paths to process
        output_file: Output file path
        batch_size: Number of files to process per batch
        verbose: Print progress messages

    Returns:
        int: Number of genomes successfully processed
    """
    total_processed = 0
    errors = []

    for batch_num, i in enumerate(range(0, len(files), batch_size), 1):
        batch_files = files[i:i+batch_size]
        batch_dfs = []

        for filepath in batch_files:
            try:
                df = parse_sensitivity_file(filepath)
                batch_dfs.append(df)
                total_processed += 1
            except Exception as e:
                errors.append((filepath, str(e)))
                if verbose:
                    print(f"  WARNING: Skipping {Path(filepath).name}: {e}", file=sys.stderr)
                continue

        if not batch_dfs:
            if verbose:
                print(f"  Batch {batch_num}: No valid files")
            continue

        # Concatenate batch
        batch_result = pd.concat(batch_dfs, ignore_index=True)
        batch_result = batch_result.sort_values(['genome_id', 'model_number']).reset_index(drop=True)

        # Write to output (append after first batch)
        if i == 0:
            batch_result.to_csv(output_file, index=False, mode='w')
        else:
            batch_result.to_csv(output_file, index=False, mode='a', header=False)

        if verbose:
            print(f"  Batch {batch_num}: Processed {len(batch_dfs)} genomes "
                  f"({i+len(batch_files)}/{len(files)} total files)")

    if verbose:
        print(f"\n✓ Successfully consolidated {total_processed} genomes → {output_file}")
        if errors:
            print(f"⚠ {len(errors)} files had errors and were skipped")

    return total_processed


def main():
    """Main entry point for command-line execution."""
    parser = argparse.ArgumentParser(
        description='Consolidate per-genome sensitivity files into long-format CSV',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (current directory)
  python consolidate-batch-predictions.py -i ./ -o final_sensitivity.csv

  # Specify custom pattern
  python consolidate-batch-predictions.py -i results/ -p "*_met_depends.csv" -o output.csv

  # Use batched processing for large datasets
  python consolidate-batch-predictions.py -i results/ -o output.csv -b 100
        """
    )

    parser.add_argument(
        '-i', '--input-dir',
        dest='input_dir',
        required=True,
        help='Directory containing individual sensitivity files'
    )

    parser.add_argument(
        '-o', '--output',
        dest='output_file',
        default='model_sensitivity_by_met_category.csv',
        help='Output file path (default: model_sensitivity_by_met_category.csv)'
    )

    parser.add_argument(
        '-p', '--pattern',
        dest='pattern',
        default='*_sensitivity.csv',
        help='Glob pattern for input files (default: *_sensitivity.csv)'
    )

    parser.add_argument(
        '-b', '--batch-size',
        dest='batch_size',
        type=int,
        default=None,
        help='Process in batches of N files (for memory efficiency, default: process all at once)'
    )

    parser.add_argument(
        '-q', '--quiet',
        dest='quiet',
        action='store_true',
        help='Suppress progress messages'
    )

    args = parser.parse_args()

    # Validate input directory exists
    if not Path(args.input_dir).exists():
        print(f"ERROR: Input directory does not exist: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    # Run consolidation
    n_processed = consolidate_sensitivity_files(
        input_dir=args.input_dir,
        output_file=args.output_file,
        pattern=args.pattern,
        batch_size=args.batch_size,
        verbose=not args.quiet
    )

    # Exit with appropriate status
    if n_processed == 0:
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
