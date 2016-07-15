# pVAC-Seq
Cancer immunotherapy has gained significant momentum from recent clinical successes of checkpoint blockade inhibition. Massively parallel sequence analysis suggests a connection between mutational load and response to this class of therapy. Methods to identify which tumor-specific mutant peptides (neoantigens) can elicit anti-tumor T cell immunity are needed to improve predictions of checkpoint therapy response and to identify targets for vaccines and adoptive T cell therapies. Here, we provide a cancer immunotherapy pipeline for the identification of **p**ersonalized **V**ariant **A**ntigens by **C**ancer **Seq**uencing (pVAC-Seq) that integrates tumor mutation and expression data (DNA- and RNA-Seq).
http://www.genomemedicine.com/content/8/1/11

## New in version 2.0.1
<ul>
<li>Supports inframe indels and frameshifts.</li>
<li>Supports VCF as the input file format.</li>
</ul>

## Citation
Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. <a href="http://www.genomemedicine.com/content/8/1/11">pVAC-Seq: A genome-guided in silico approach to identifying tumor neoantigens</a>. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/26825632">26825632</a>.

## License
This project is licensed under <a href="http://opensource.org/licenses/NPOSL-3.0">NPOSL-3.0</a>.

## Installation Instructions
pVAC-Seq requires Python 3.5. You can check your installed Python version by running `python -V`. If you don't have Python 3.5 installed, we recommend using <a href="http://conda.pydata.org/docs/py2or3.html">Conda</a> to emulate a Python 3.5. environment. We've encountered problems with users that already have Python 2.* installed when they also try to install Python 3.5. The defaults will not be set correctly in that case. If you already have Python 2.* installed we <b>strongly</b> recommmend using Conda instead of installing Python 3.5 locally.

Once you have set up your Python 3.5 environment correctly you can use `pip` to install pVAC-Seq. Make sure you have `pip` installed.  `pip` is generally included in python distributions, but may need to be upgraded before use.  See the <a href="https://packaging.python.org/en/latest/installing/#install-pip-setuptools-and-wheel">instructions</a> for installing or upgrading pip.

After you have pip installed, type the following command on your Terminal(for Mac and Linux users) or the command prompt (for Windows users): `pip install pvacseq`. You can check that pvacseq has been installed under the default environment by running `pip list`.

pip will fetch and install pVAC-Seq and its dependencies for you.  After installing, you can run `pvacseq` directly from the Terminal/command prompt.

##Prerequisites
###<b>VEP</b>

The input to the pVAC-Seq pipeline is a VEP annotated VCF. In addition to the standard VEP annotations, pVAC-Seq also requires the annotations provided by the Downstream and Wildtype VEP plugins. To create a VCF for use with pVAC-Seq follow these steps:
- Download and install the VEP command line tool following the instructions found <a href="http://useast.ensembl.org/info/docs/tools/vep/script/index.html">here</a>.
- Download the VEP_plugins from their <a href="https://github.com/Ensembl/VEP_plugins">Github repository</a>.
- Copy the Wildtype plugin provided with the pVAC-Seq package to the folder with the other VEP_plugins by running `pvacseq install_vep_plugin`.
- Run VEP on the input vcf with at least the following options:<br>
`--format vcf`<br>
`--vcf`<br>
`--symbol`<br>
`--plugin Downstream`<br>
`--plugin Wildtype`<br>
`--terms SO`<br>
The `--dir_plugins <VEP_plugins directory>` option may need to be set depending on where the VEP_plugins were installed to. Additional VEP options that might be desired can be found <a href="http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html">here</a>.<br>
<b>Example VEP Command</b><br>
`perl variant_effect_predictor.pl --input_file <input VCF> --format vcf --output_file <output VCF> --vcf --symbol --terms SO --plugin Downstream --plugin Wildtype [--dir_plugins <VEP_plugins directory>]`

###<b>NetMHC 3.4</b>

pVAC-Seq uses NetMHC 3.4 to predict binding affinities. <a href="http://www.cbs.dtu.dk/cgi-bin/sw_request?netMHC+3.4">NetMHC 3.4 can be downloaded here</a>. Once NetMHC is properly installed and tested, pVAC-Seq expects the path to the installation directory.

## Pipeline Overview
![alt text][overview]
[overview]:
https://raw.githubusercontent.com/wiki/griffithlab/pVAC-Seq/images/pvacseq-code-python.png

## pvacseq commands
### run
`pvacseq run <input VCF> <sample name> <NetMHC installation path> <allele name> <epitope length> <ouput directory> [-l peptide sequence length] [-b binding threshold] [-c minimum fold change]`<br>
Use this command to run the full pVAC-Seq pipeline.  This will internally call the other commands, passing data between them to generate an output TSV file of neoepitope predictions. Multiple alleles and epiope length can be specified as comma-separated lists.

<b>Required inputs</b><br>
<ul>
<li><code>input VCF</code>: A VEP-annotated VCF containing transcript, Wildtype protein sequence, and Downstream protein sequence information. (Please see above for instructions)</li>
<li><code>sample name</code>: The name of the sample being processed. This will be used as a prefix for output files.</li>
<li><code>NetMHC installation path</code>: The path to the NetMHC installation directory (please see above for installation instructions)</li>
<li><code>allele name</code>: Name of the allele to use for epitope prediction. Mutliple alles can be specified using a comma-separated list.</li>
<li><code>epitope length</code>: This refers to the length of subpeptides (neoepitopes) to predict. The pipeline can handle multiple lengths, which can be specified using a comma-separated list. Typical epitope lengths vary between 8-11.</li>
<li><code>Output directory</code>: The directory for writing all result files.</li>
</ul>

<b>Optional inputs</b><br>
<ul>
<li><code>peptide sequence length</code>: Length of the peptide sequence to use when creating the FASTA. See "Additional Information" for details.
<li><code>binding threshold</code>: The user can choose to report only epitopes where the mutant allele has IC50 binding scores below this value. By default, pvacseq uses a cutoff of 500.
<li><code>minimum fold change</code>: This parameter sets the minimum fold change between mutant binding score and wild-type score to use for filtering. The default is 0, which filters no results. Using 1 will require that binding is better to the MT than WT.</li>
</ul>

### convert_vcf
`pvacseq convert_vcf <input VCF> <output TSV file>`<br>
Run this command to generate a TSV file with annotated variants from a VEP-annotated VCF.

### generate_fasta
`pvacseq generate_fasta <input TSV file> <peptide sequence length> <output FASTA file>`<br>
Run this command to generate a FASTA file for wildtype(WT) and mutant(MT) amino acid sequences for MHC Class I epitope prediction. The length of the amino acid sequences is determined by the peptide sequence length specified. The input file is the properly formatted TSV file of annotated variants.

### generate_fasta_key
`pvacseq generate_fasta_key <input FASTA file> <output key file>`<br>
NetMHC strips off the name of the FASTA header. This command generates a key file to lookup each NetMHC output entry to its original entry in the FASTA file.

### parse_output
`pvacseq parse_output <NetMHC output file> <input TSV file> <FASTA key file> <output parsed file>`<br>
After running NetMHC 3.4, this command parses the output for MHC Class I epitope prediction. It uses a special key file to link each NetMHC result entry to the original entry from the input TSV file. The parsed TSV output file contains predictions for the mutant as well as the wildtype version of the epitope, and compares binding affinities for the same. It also contains gene and transcript information from the input TSV file.

### binding_filter
`pvacseq binding_filter <input TSV file> <output file> [-b binding threshold] [-c minimum fold change]`<br>
Takes a comma-separated list of parsed NetMHC files for different allele-length combinations and outputs best candidates per gene based on binding affinities.

### coverage_filter
`pvacseq coverage_filters <input TSV file> <output file> [--normal-cov normal coverage cutoff] [--tdna-cov tumor DNA coverage cutoff] [--trna-cov tumor RNA coverage cutoff] [--normal-vaf normal vaf cutoff] [--tdna-vaf tumor DNA vaf cutoff] [--trna-vaf tumor RNA vaf cutoff] [--expn-val gene expression (fpkm) cutoff]`<br>
Depending on the type(s) of sequencing data available, a variety of coverage and expression based filters can be installed. The input file should contain the predicted epitopes along with read counts appended as additional columns. If specific type of sequencing data is not available, the columns can be left off. Column order is not important.

The input TSV file contains the following columns in tab-separated format:<br>
Chromosome<br>
Start<br>
Stop<br>
Reference<br>
Variant<br>
Transcript<br>
Ensembl Gene ID<br>
Variant Type<br>
Mutation<br>
Protein Position<br>
Gene Name<br>
HLA Allele<br>
Peptide Length<br>
Sub-peptide Position<br>
MT score<br>
WT score<br>
MT epitope seq<br>
WT epitope seq<br>
Fold Change<br>
Normal Ref Count<br>
Normal Var Count<br>
Tumor DNA Ref Count<br>
Tumor DNA Var Count<br>
Tumor RNA Ref Count<br>
Tumor RNA Var Count<br>
Gene Exp FPKM<br>

### download_example_data
`pvacseq download_example_data <destination directory>`
Downloads a set of example data files to the directory specififed.

### install_vep_plugin
`pvacseq install_vep_plugin <vep plugins path>`
Installs the Wildtype VEP plugin into the specified directory.


## Additional Information

Since the goal of the pVAC-Seq pipeline is to predict putative 'neo'antigens, we only consider a sub-section of the transcript sequence encompassing the mutated amino acid.
<ol type="a">
<li>In the following figure, the amino acid FASTA sequence is built using 10 flanking amino acids on each side of the mutated amino acid. The preceding or succeeding 20 amino acids are taken if the mutation lies near the end or beginning of the transcript, respectively.</li>
<li>All predicted candidate peptides from epitope prediction software based on selected k-mer window size.</li>
<li>Only localized peptides (those containing the mutant amino acid) are considered to compare to wild-type counterpart.</li>
<li>The ‘best candidate’ (lowest MT binding score) per mutation is chosen across all specified k-mers that were installed as input.</li>
</ol>
 ![alt text][logo]
 [logo]:
 https://github.com/jhundal/src/blob/master/bin/images/Fig1_fastav2.png
