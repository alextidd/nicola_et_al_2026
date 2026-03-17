#!/usr/bin/env python
import sys, argparse
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--output_dir", required=True, help="output_directory")
  parser.add_argument("--input_data", required=True, help="input nucleotide context")
  parser.add_argument("--project", required=True, help="project title")
  parser.add_argument("--reference_genome", required=True, help="reference genome")
  parser.add_argument("--cosmic_version", required=False, default=3.4, type=float,
                      help="COSMIC signature version to use for decomposition (default: 3.4)")
  args = parser.parse_args()
  print(args)

  sig.sigProfilerExtractor(
      input_type="matrix",
      output=args.output_dir,
      input_data=args.input_data,
      reference_genome=args.reference_genome,
      context_type="default",
      exome=False,
      minimum_signatures=1,
      maximum_signatures=20,
      nmf_replicates=100,
      resample=True,
      batch_size=1,
      cpu=1,
      gpu=False,
      nmf_init="random",
      precision="single",
      matrix_normalization="gmm",
      seeds="random",
      min_nmf_iterations=10000,
      max_nmf_iterations=1000000,
      nmf_test_conv=10000,
      nmf_tolerance=1e-15,
      get_all_signature_matrices=False,
      cosmic_version=args.cosmic_version
  )

if __name__ == "__main__":
  main()