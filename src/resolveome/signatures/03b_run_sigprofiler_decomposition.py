#!/usr/bin/env python3
# module purge ; module load python-3.12.0/perl-5.38.0
# source ~/envs/sigprofiler_v1.1.3/bin/activate

# packages
from SigProfilerAssignment import Analyzer as Analyze

# input matrix (mutation types x samples)
samples = "out/resolveome/signatures/matrices/trinuc_mut_mat_sigpro.txt"

# custom signatures
database = "out/resolveome/signatures/ref_sigs/ref_sigs.tsv"

# recommended solution (n=4)
signatures = "out/resolveome/signatures/sigprofiler/extractor/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt"

# default
output = "out/resolveome/signatures/sigprofiler/decompose/n=4,default/"
Analyze.decompose_fit(samples, output, signatures = signatures,
                      genome_build = "GRCh38", verbose = False,
                      cosmic_version = 3.4, make_plots = True,
                      collapse_to_SBS96 = True)

# optimised
output = "out/resolveome/signatures/sigprofiler/decompose/n=4,optimised/"
to_exclude = ['MMR_deficiency_signatures', 'POL_deficiency_signatures',
              'HR_deficiency_signatures', 'BER_deficiency_signatures',
              'Chemotherapy_signatures', 'Immunosuppressants_signatures',
              'Treatment_signatures', 'APOBEC_signatures', 'Tobacco_signatures',
              'UV_signatures', 'AA_signatures', 'Colibactin_signatures',
              'Artifact_signatures']
Analyze.decompose_fit(
	samples, output, signatures = signatures, genome_build = "GRCh38",
	verbose = False, cosmic_version = 3.4, make_plots = True,
	collapse_to_SBS96 = True,
	exclude_signature_subgroups = to_exclude)

# custom
output = "out/resolveome/signatures/sigprofiler/decompose/n=4,custom/"
Analyze.decompose_fit(
	samples, output, signatures = signatures, genome_build = "GRCh38",
	verbose = False, cosmic_version = 3.4, make_plots = True,
	collapse_to_SBS96 = True, signature_database = database)

# alternative solution (n=5)
signatures = "out/resolveome/signatures/sigprofiler/extractor/SBS96/All_Solutions/SBS96_5_Signatures/Signatures/SBS96_S5_Signatures.txt"

# default
output = "out/resolveome/signatures/sigprofiler/decompose/n=5,default/"
Analyze.decompose_fit(samples, output, signatures = signatures,
                      genome_build = "GRCh38", verbose = False,
                      cosmic_version = 3.4, make_plots = True,
                      collapse_to_SBS96 = True)

# optimised
output = "out/resolveome/signatures/sigprofiler/decompose/n=5,optimised/"
to_exclude = ['MMR_deficiency_signatures', 'POL_deficiency_signatures',
              'HR_deficiency_signatures', 'BER_deficiency_signatures',
              'Chemotherapy_signatures', 'Immunosuppressants_signatures',
              'Treatment_signatures', 'APOBEC_signatures', 'Tobacco_signatures',
              'UV_signatures', 'AA_signatures', 'Colibactin_signatures',
              'Artifact_signatures']
Analyze.decompose_fit(
	samples, output, signatures = signatures, genome_build = "GRCh38",
	verbose = False, cosmic_version = 3.4, make_plots = True,
	collapse_to_SBS96 = True,
	exclude_signature_subgroups = to_exclude)

# custom
output = "out/resolveome/signatures/sigprofiler/decompose/n=5,custom/"
Analyze.decompose_fit(
	samples, output, signatures = signatures, genome_build = "GRCh38",
	verbose = False, cosmic_version = 3.4, make_plots = True,
	collapse_to_SBS96 = True, signature_database = database)
