#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch P2RANK processing for multiple PDB structures

Usage:
    python batch_p2rank.py --input-dir <pdb_directory> --output-dir <output>
    python batch_p2rank.py --input-list <file_list.txt> --output-dir <output>

For 1000+ structures:
    python batch_p2rank.py --input-dir alphafold_models/ --output-dir p2rank_results/

Author: ASMC
Date: 2025-11-20
"""

import sys
import io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import os
import argparse
import subprocess
from pathlib import Path
import time
from datetime import datetime


class BatchP2RANK:
    """Batch processing for P2RANK"""

    def __init__(self, p2rank_dir, java_home=None):
        self.p2rank_dir = Path(p2rank_dir)
        self.prank_script = self.p2rank_dir / "prank"

        if not self.prank_script.exists():
            raise FileNotFoundError(f"P2RANK not found: {self.prank_script}")

        # Set Java environment
        if java_home:
            os.environ['JAVA_HOME'] = java_home
            os.environ['PATH'] = f"{java_home}/bin{os.pathsep}{os.environ.get('PATH', '')}"

    def create_dataset_file(self, pdb_files, dataset_file):
        """
        Create P2RANK dataset file

        Args:
            pdb_files: List of PDB file paths
            dataset_file: Output dataset file path
        """
        with open(dataset_file, 'w', encoding='utf-8') as f:
            f.write("# P2RANK batch dataset\n")
            f.write(f"# Generated: {datetime.now()}\n")
            f.write(f"# Total structures: {len(pdb_files)}\n")
            f.write("\n")

            for pdb_file in pdb_files:
                # Convert to absolute path
                abs_path = Path(pdb_file).absolute()
                f.write(f"{abs_path}\n")

        print(f"Created dataset file: {dataset_file}")
        print(f"  Contains {len(pdb_files)} structures")

    def run_p2rank(self, dataset_file, output_dir, threads=4):
        """
        Run P2RANK on dataset

        Args:
            dataset_file: Path to dataset file
            output_dir: Output directory
            threads: Number of threads to use

        Returns:
            dict: Results summary
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            str(self.prank_script),
            "predict",
            str(dataset_file),
            "-o", str(output_dir),
            "-threads", str(threads)
        ]

        print(f"\nRunning P2RANK...")
        print(f"Command: {' '.join(cmd)}")
        print(f"Output: {output_dir}")
        print()

        start_time = time.time()

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(self.p2rank_dir)
            )

            elapsed = time.time() - start_time

            if result.returncode == 0:
                print(f"P2RANK completed successfully in {elapsed:.1f}s")
                return {
                    'success': True,
                    'elapsed': elapsed,
                    'output_dir': output_dir
                }
            else:
                print(f"P2RANK failed with error code {result.returncode}")
                print(f"STDERR: {result.stderr}")
                return {
                    'success': False,
                    'error': result.stderr
                }

        except Exception as e:
            print(f"Error running P2RANK: {e}")
            return {
                'success': False,
                'error': str(e)
            }

    def collect_results(self, output_dir):
        """
        Collect all prediction results

        Returns:
            list: [(pdb_name, predictions_csv_path), ...]
        """
        output_dir = Path(output_dir)
        results = []

        for csv_file in output_dir.glob("*_predictions.csv"):
            pdb_name = csv_file.stem.replace("_predictions", "")
            results.append((pdb_name, csv_file))

        return sorted(results)


def find_pdb_files(input_dir, pattern="*.pdb"):
    """Find all PDB files in directory"""
    input_path = Path(input_dir)

    if input_path.is_file():
        return [input_path]

    pdb_files = list(input_path.glob(pattern))
    pdb_files.extend(input_path.glob("**/" + pattern))

    return sorted(set(pdb_files))


def main():
    parser = argparse.ArgumentParser(
        description="Batch P2RANK processing for multiple structures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all PDB files in directory
  python batch_p2rank.py --input-dir alphafold_models/ --output-dir p2rank_results/

  # Process specific files from list
  python batch_p2rank.py --input-list structures.txt --output-dir results/

  # Use custom P2RANK installation
  python batch_p2rank.py --input-dir models/ --output-dir results/ \\
      --p2rank-dir /path/to/p2rank_2.4.2

  # Specify Java home
  python batch_p2rank.py --input-dir models/ --output-dir results/ \\
      --java-home "/c/Program Files/Java/jdk-11"
        """
    )

    parser.add_argument('--input-dir', '-i',
                       help='Input directory containing PDB files')
    parser.add_argument('--input-list', '-l',
                       help='Text file with list of PDB file paths')
    parser.add_argument('--output-dir', '-o', required=True,
                       help='Output directory for P2RANK results')
    parser.add_argument('--p2rank-dir', '-p',
                       default='p2rank_2.4.2',
                       help='P2RANK installation directory (default: p2rank_2.4.2)')
    parser.add_argument('--java-home',
                       default='/c/Program Files/JetBrains/PyCharm Community Edition 2024.2.1/jbr',
                       help='JAVA_HOME path')
    parser.add_argument('--threads', '-t', type=int, default=4,
                       help='Number of threads (default: 4)')
    parser.add_argument('--pattern', default='*.pdb',
                       help='File pattern for input directory (default: *.pdb)')

    args = parser.parse_args()

    # Validate inputs
    if not args.input_dir and not args.input_list:
        parser.error("Either --input-dir or --input-list must be specified")

    print("="*80)
    print("Batch P2RANK Processing")
    print("="*80)
    print()

    # Find PDB files
    if args.input_list:
        print(f"Reading file list: {args.input_list}")
        with open(args.input_list, 'r') as f:
            pdb_files = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    else:
        print(f"Searching for PDB files in: {args.input_dir}")
        pdb_files = find_pdb_files(args.input_dir, args.pattern)

    if not pdb_files:
        print("Error: No PDB files found!")
        sys.exit(1)

    print(f"Found {len(pdb_files)} PDB files")
    print()

    # Initialize P2RANK
    try:
        p2rank = BatchP2RANK(args.p2rank_dir, java_home=args.java_home)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Create dataset file
    dataset_file = Path(args.output_dir) / "batch_dataset.ds"
    dataset_file.parent.mkdir(parents=True, exist_ok=True)

    p2rank.create_dataset_file(pdb_files, dataset_file)
    print()

    # Run P2RANK
    result = p2rank.run_p2rank(
        dataset_file,
        args.output_dir,
        threads=args.threads
    )

    if not result['success']:
        print("P2RANK processing failed!")
        sys.exit(1)

    # Collect results
    print()
    print("="*80)
    print("Collecting Results")
    print("="*80)

    results = p2rank.collect_results(args.output_dir)

    print(f"Generated {len(results)} prediction files")
    print()

    if results:
        print("Sample results:")
        for pdb_name, csv_file in results[:5]:
            print(f"  {pdb_name}: {csv_file.name}")

        if len(results) > 5:
            print(f"  ... and {len(results) - 5} more")

    print()
    print("="*80)
    print("Batch Processing Complete!")
    print("="*80)
    print()
    print(f"Results saved to: {args.output_dir}")
    print(f"Total structures: {len(pdb_files)}")
    print(f"Elapsed time: {result['elapsed']:.1f}s")
    print(f"Average: {result['elapsed']/len(pdb_files):.2f}s per structure")
    print()
    print("Next steps:")
    print(f"  1. Analyze results: python parse_p2rank_pockets.py {args.output_dir}/<file>_predictions.csv")
    print(f"  2. Compare with known sites")
    print(f"  3. Use for ASMC analysis")
    print()


if __name__ == "__main__":
    main()
