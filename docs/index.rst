.. pVAC-Seq documentation master file, created by
   sphinx-quickstart on Wed Aug 31 08:31:01 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pVAC-Seq
====================================

pVAC-Seq is a cancer immunotherapy pipeline for the identification of **p**\ ersonalized **V**\ ariant **A**\ ntigens by **C**\ ancer **Seq**\ uencing (pVAC-Seq) that integrates tumor mutation and expression data (DNA- and RNA-Seq). It enables cancer immunotherapy research by using massively parallel sequence data to predicting tumor-specific mutant peptides (neoantigens) that can elicit anti-tumor T cell immunity. It is being used in studies of checkpoint therapy response and to identify targets for cancer vaccines and adoptive T cell therapies. For more general information, see the `manuscript published in Genome Medicine <http://www.genomemedicine.com/content/8/1/11>`_.

.. toctree::
   :maxdepth: 3

   features
   install
   prerequisites
   run
   filter_commands
   additional_commands
   optional_downstream_analysis_tools
   contact

New in version |version|
------------------------

We added an :ref:`optional downstream analysis tool<optional_downstream_analysis_tools_label>` to generate an annotated fasta file from a VCF with protein sequences of mutations and matching wildtypes. This tool can be run with the ``pvacseq generate_protein_fasta`` command.

This release fixes a couple of errors that were introduced in the previous version which would occur during the processing of certain inframe indels.

This version also fixes an error that would occur if the number of variants to
process was a multiple of the chosen ``--fasta-size``.

Citation
--------

Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi L. Griffith, Elaine R. Mardis, and Malachi Griffith. `pVAC-Seq: A genome-guided in silico approach to identifying tumor neoantigens <http://www.genomemedicine.com/content/8/1/11>`_. Genome Medicine. 2016, 8:11, DOI: 10.1186/s13073-016-0264-5. PMID: `26825632 <http://www.ncbi.nlm.nih.gov/pubmed/26825632>`_.

License
-------
This project is licensed under `NPOSL-3.0 <http://opensource.org/licenses/NPOSL-3.0>`_.
