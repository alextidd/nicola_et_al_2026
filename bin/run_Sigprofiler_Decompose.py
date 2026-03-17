#!/usr/bin/env python
"""
Decompose de novo signatures to COSMIC reference signatures.

NOTE: This function only supports decomposition to COSMIC signatures.
For fitting samples to custom signatures, use run_Sigprofiler_Assignment.py instead.
"""
import sys, argparse
from SigProfilerAssignment import Analyzer as Analyze

def main():
    parser = argparse.ArgumentParser(
        description="Decompose de novo signatures to COSMIC reference signatures"
    )
    parser.add_argument("--signatures", required=True,
                        help="Path to de novo signatures file (e.g., SBS96_S8_Signatures.txt)")
    parser.add_argument("--samples", required=True,
                        help="Path to samples matrix or directory containing matrices")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for decomposition results")
    parser.add_argument("--reference_genome", required=True,
                        help="Reference genome (GRCh37 or GRCh38)")
    parser.add_argument("--cosmic_version", type=float, default=3.4,
                        help="COSMIC signature version (default: 3.4)")
    parser.add_argument("--no_plots", action="store_true", default=False,
                        help="Disable plot generation (useful if TMB plotting fails)")
    parser.add_argument("--collapse_to_sbs96", action="store_true", default=True,
                        help="Collapse to SBS96 context (default: True)")
    args = parser.parse_args()
    print(args)

    Analyze.decompose_fit(
        samples=args.samples,
        output=args.output_dir,
        signatures=args.signatures,
        genome_build=args.reference_genome,
        cosmic_version=args.cosmic_version,
        make_plots=not args.no_plots,
        collapse_to_SBS96=args.collapse_to_sbs96,
        verbose=True
    )

if __name__ == "__main__":
    main()
