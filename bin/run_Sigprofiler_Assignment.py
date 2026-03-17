#!/usr/bin/env python
"""
Run SigProfilerAssignment to fit reference signatures to samples.

Two modes:
1. COSMIC-only: Use built-in COSMIC database (default, or with --include_sigs to subset)
2. Custom: Provide --signatures_file (optionally with --include_sigs to subset)
"""

import sys
import os
import argparse
import pandas as pd
from SigProfilerAssignment import Analyzer as Analyze


def main():
    parser = argparse.ArgumentParser(
        description="Fit reference signatures to samples using SigProfilerAssignment"
    )
    parser.add_argument("--samples", required=True,
                        help="Path to samples matrix (SigProfiler format)")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for assignment results")
    parser.add_argument("--reference_genome", required=True,
                        help="Reference genome (GRCh37 or GRCh38)")
    
    parser.add_argument("--signatures_file",
                        help="Path to custom signatures matrix file (TSV/CSV). If not provided, uses COSMIC.")
    parser.add_argument("--include_sigs",
                        help="Comma-separated list of signatures to fit (e.g., SBS1,SBS5,SBS9). Works with both COSMIC and custom files.")
    
    parser.add_argument("--cosmic_version", type=float, default=3.4,
                        help="COSMIC signature version (default: 3.4)")
    parser.add_argument("--no_plots", action="store_true", default=False,
                        help="Disable plot generation")
    parser.add_argument("--export_probabilities", action="store_true", default=False,
                        help="Export mutation-level signature probabilities")
    
    args = parser.parse_args()
    print(args)
    
    # Parse include_sigs if provided
    include_sigs = None
    if args.include_sigs:
        include_sigs = [s.strip() for s in args.include_sigs.split(',')]
        print(f"\nSubsetting to {len(include_sigs)} signatures: {include_sigs}")

    if args.signatures_file:
        # Mode 2: Use custom signatures file
        print(f"\nFitting custom signatures from: {args.signatures_file}")
        
        sig_file = args.signatures_file
        
        # Subset signatures if requested
        if include_sigs:
            # Read full signatures matrix
            sep = '\t' if args.signatures_file.endswith('.tsv') else ','
            sigs_df = pd.read_csv(args.signatures_file, sep=sep, index_col=0)
            
            # Check all requested signatures exist
            missing = [s for s in include_sigs if s not in sigs_df.columns]
            if missing:
                raise ValueError(f"Signatures not found in file: {missing}")
            
            # Subset to requested signatures
            sigs_subset = sigs_df[include_sigs]
            
            # Write to temp file
            os.makedirs(args.output_dir, exist_ok=True)
            sig_file = os.path.join(args.output_dir, "reference_signatures_subset.tsv")
            sigs_subset.to_csv(sig_file, sep='\t')
            print(f"Wrote subset signatures to: {sig_file}")
        
        Analyze.cosmic_fit(
            samples=args.samples,
            output=args.output_dir,
            genome_build=args.reference_genome,
            cosmic_version=args.cosmic_version,
            signature_database=sig_file,
            nnls_add_penalty=0.05,
            nnls_remove_penalty=0.01,
            initial_remove_penalty=0.05,
            make_plots=not args.no_plots,
            collapse_to_SBS96=True,
            connected_sigs=False,
            verbose=True,
            exome=False,
            input_type="matrix",
            export_probabilities=args.export_probabilities
        )
    else:
        # Mode 1: Use COSMIC signatures
        if include_sigs:
            print(f"\nFitting {len(include_sigs)} COSMIC signatures: {include_sigs}")
            
            # Load COSMIC database and subset it
            from SigProfilerMatrixGenerator import install as genInstall
            cosmic_path = os.path.join(
                os.path.dirname(genInstall.__file__),
                "references", "chromosomes", "tsb", args.reference_genome,
                f"COSMIC_v{args.cosmic_version}_SBS_{args.reference_genome}.txt"
            )
            
            # Try alternate path if not found
            if not os.path.exists(cosmic_path):
                cosmic_path = f"../../reference/cosmic/COSMIC_v{args.cosmic_version}_SBS_{args.reference_genome}.txt"
            
            if not os.path.exists(cosmic_path):
                raise FileNotFoundError(f"COSMIC signatures file not found. Tried: {cosmic_path}")
            
            # Read and subset COSMIC
            cosmic_df = pd.read_csv(cosmic_path, sep='\t', index_col=0)
            missing = [s for s in include_sigs if s not in cosmic_df.columns]
            if missing:
                raise ValueError(f"Signatures not found in COSMIC: {missing}")
            
            cosmic_subset = cosmic_df[include_sigs]
            
            # Write to temp file
            os.makedirs(args.output_dir, exist_ok=True)
            sig_file = os.path.join(args.output_dir, "cosmic_signatures_subset.tsv")
            cosmic_subset.to_csv(sig_file, sep='\t')
            print(f"Wrote subset COSMIC signatures to: {sig_file}")
            
            Analyze.cosmic_fit(
                samples=args.samples,
                output=args.output_dir,
                genome_build=args.reference_genome,
                cosmic_version=args.cosmic_version,
                signature_database=sig_file,
                nnls_add_penalty=0.05,
                nnls_remove_penalty=0.01,
                initial_remove_penalty=0.05,
                make_plots=not args.no_plots,
                collapse_to_SBS96=True,
                connected_sigs=False,
                verbose=True,
                exome=False,
                input_type="matrix",
                export_probabilities=args.export_probabilities
            )
        else:
            print("\nFitting all COSMIC signatures")
            
            Analyze.cosmic_fit(
                samples=args.samples,
                output=args.output_dir,
                genome_build=args.reference_genome,
                cosmic_version=args.cosmic_version,
                signature_database=None,
                nnls_add_penalty=0.05,
                nnls_remove_penalty=0.01,
                initial_remove_penalty=0.05,
                make_plots=not args.no_plots,
                collapse_to_SBS96=True,
                connected_sigs=True,
                verbose=True,
                exclude_signature_subgroups=None,
                exome=False,
                input_type="matrix",
                export_probabilities=args.export_probabilities
            )


if __name__ == "__main__":
    main()
