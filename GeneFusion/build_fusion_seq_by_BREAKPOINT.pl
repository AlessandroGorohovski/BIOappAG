#!/usr/bin/perl -w

use strict;
use warnings;

$|++;
my $VERSION = '1.5'; my $AUTHOR="Alessandro Gorohovski"; my $YEAR="2020-2023";

=head1 NAME
=encoding utf8

build_fusion_seq_by_BREAKPOINT.pl - tool for preparing gene fusion sequences using breakpoint location.

=head1 VERSION

1.5

=head1 USAGE

  ./build_fusion_seq_by_BREAKPOINT.pl --assembly <hg19|hg38> --breakpoint <geneA:chrA:strandA:locationA:geneB:chrB:strandB:locationB> --disease <disease> --pubmed <pubmed> [cell_line]

=over 3

=item 1.
Example (+)(+) strands
  ./build_fusion_seq_by_BREAKPOINT.pl --assembly hg38 --breakpoint VTI1A:10:+:112,464,657:TCF7L2:10:+:113,141,184 --disease GBM --pubmed 32414213

=item 2.
Example (+)(-) strands
  ./build_fusion_seq_by_BREAKPOINT.pl --assembly hg38 --breakpoint TGFB1:19:+:41,351,202:TAL1:1:-:47,229,351 --disease GBM --pubmed 32414213

=item 3.
Example (-)(+) strands
  ./build_fusion_seq_by_BREAKPOINT.pl --assembly hg38 --breakpoint CD74:5:-:150,402,152:GID8:20:+:62,944,739 --disease GBM --pubmed 32414213

=item 4.
Example (-)(-) strands
  ./build_fusion_seq_by_BREAKPOINT.pl --assembly hg38 --breakpoint CDCA7L:7:-:21,945,781:MLLT3:9:-:20,354,879 --disease GBM --pubmed 32414213

=back

=head2 OUTPUT

stdout | OUT_fusion_sequence.csv file

=head2 OPTIONS

=over 3

=item * C<hg38|hg19> - Required

=item * C<geneA:chrA:strandA:locationA:geneB:chrB:strandB:locationB> - Required

=item * C<disease> - Required

=item * C<pubmed> - Required

=item * C<[cell_line]> - Optional

=back


=head2 REQUIREMENTS

C<extractseq> EMBOSS tool. See L<http://emboss.open-bio.org/>

Sequences of human chromosomes.
Prepare reference files:

    $ wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz
    $ gunzip GCF_000001405.25_GRCh37.p13_genomic.fna.gz
    $ mkdir GCF_000001405.25_GRCh37.p13.chroms
    $ cd GCF_000001405.25_GRCh37.p13.chroms
    $ ln -s ../GCF_000001405.25_GRCh37.p13_genomic.fna .
    $ seqretsplit GCF_000001405.25_GRCh37.p13_genomic.fna

and rename all fasta to chr{1,2,3,etc}.fa

=head1 AUTHOR

Alessandro N. Gorohovski (AG), E<lt>an.gorohovski@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2020-2023 by Alessandro N. Gorohovski (AG). All Rights Reserved.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.0 or,
at your option, any later version of Perl 5 you may have available.

=cut


use File::Spec;
use Getopt::Long;

###
# Default Options
my $OUTPUT;
my $START_ID;
my $BREAKPOINT;
my $ASSEMBLY;	# hg19 , hg38
my $DISEASE;
my $PUBMED;
my $CELL_LINE = '';
my $DELTA = 90; #51;	# длина каждого фрагмента HEAD & TAIL
my $HEADER;
my $BSDIR = '/home/aleksander/BioInf/Programs/B/BLAT/Chimeras_search/genomes/human';
my $AUTO;

###
# Parse input data
GetOptions(
	'output=s'     => \$OUTPUT,
	'start_id=i'   => \$START_ID,
	'breakpoint=s' => \$BREAKPOINT,
	'assembly=s'   => \$ASSEMBLY,
	'disease=s'    => \$DISEASE,
	'pubmed=s'     => \$PUBMED,
	'cell_line=s'  => \$CELL_LINE,
	'delta=i'      => \$DELTA,
	'header'       => \$HEADER,
	'bsdir=s'      => \$BSDIR,
	'auto'         => \$AUTO,
) or die &usage();

$DISEASE || die &usage();
$PUBMED  || die &usage();

if( defined $AUTO ){
#	$OUTPUT ||= 'OUT_fusion_sequence.csv';
	$START_ID ||= 1000;
	$ASSEMBLY ||= 'hg19';
}

$OUTPUT = undef if $OUTPUT && $OUTPUT =~/^stdout$/i;
my $fusID = $START_ID;

if( defined $OUTPUT ){
	if( open RF, $OUTPUT ){
		while( <RF>){
			s/^\s+|\s+$//g;
			next if /^$/ || /^#/;

			my( $id ) = split /\t/;
			$id =~s/\D+//g;
			$fusID = $id + 1 if $id+0 >= $fusID;
		}
		close RF;
	}
}
$fusID || die &usage();

my %chrD = (
		'hg19' => 'GCF_000001405.25_GRCh37.p13.chroms',
		'hg38' => 'GCF_000001405.39_GRCh38.p13.chroms',
	);

$ASSEMBLY || die &usage(); # hg38 , hg19
die &usage() unless exists $chrD{ $ASSEMBLY };

my $baseCHR = "$BSDIR/$chrD{ $ASSEMBLY }";
die "$baseCHR does not exist: $!" unless -d $baseCHR;

# geneA:chrA:strandA:locationA:geneB:chrB:strandB:locationB
$BREAKPOINT || die &usage(); # e.g. VTI1A:10:+:112,464,657:TCF7L2:10:+:113,141,184

my( $geneA,$chrA,$strandA,$locA, $geneB,$chrB,$strandB,$locB ) = split /:/, $BREAKPOINT;

# Проверяем и Чистим chromosomes
for( $chrA, $chrB ){
	s/^ch?r?//i;

	if( /^([1-9XYM]|1\d|2[012]|MT)$/ ){
		$_ = $1;
	}else{
		die"ERROR: wrong chromosome $_";
	}
}

# Проверяем и Чистим strands
for( $strandA, $strandB ){
	s/[^\+\-]//g;
	die"ERROR: wrong strand!" unless length;
}

# Проверяем и Чистим locations
for( $locA, $locB ){
	s/\D+//g;
	die"ERROR: wrong location!" unless length;
}

my( $seqA ) = &HT_location('HEAD', $strandA, $locA, $DELTA, "$baseCHR/chr$chrA.fa", "\x1b[32m");
my( $seqB ) = &HT_location('TAIL', $strandB, $locB, $DELTA, "$baseCHR/chr$chrB.fa", "\x1b[36m");

my $head_out = join("\t", qw(fusID ncbi_A gene_A ncbi_B gene_B sequence_A sequence_B exon_A exon_B), 'disease(,s)', 'pubmed(;s)', 'cell_line(,s)');
my $out = join("\t", "FUS$fusID", '.', $geneA, '.', $geneB, $seqA, $seqB, '.', '.', $DISEASE, $PUBMED, $CELL_LINE );

# Save results
if( $OUTPUT ){
	open OFILE, ">>$OUTPUT";
	print OFILE "# $head_out\n" if $HEADER;
	print OFILE "$out\n";
	close OFILE;

}else{
	print "$head_out\n" if $HEADER;
	print "$out\n";
}

exit;


sub HT_location {
	my( $HT, $strand, $loc, $delta, $chro_seq, $color ) = @_;

	my( $sloc, $eloc );
	my $revseq = '';

	if( $HT eq 'HEAD'){
		if( $strand eq '+'){
			$sloc = $loc - $delta;
			$eloc = $loc;

		}elsif( $strand eq '-'){
			$sloc = $loc;
			$eloc = $loc + $delta;

			$revseq = '| revseq -filter';
		}

	}elsif( $HT eq 'TAIL'){
		if( $strand eq '+'){
			$sloc = $loc;
			$eloc = $loc + $delta;

		}elsif( $strand eq '-'){
			$sloc = $loc - $delta;
			$eloc = $loc;

			$revseq = '| revseq -filter';
		}

	}else{
		die"ERROR: NO HEAD/TAIL";
	}

	my $run = qq{extractseq -auto -sequence $chro_seq -regions "$sloc:$eloc" -stdout $revseq};
	print "$run\n";

	my $seq;
	for(`$run`){
		next if /^>/;

		print "$color$_\x1b[0m\n";

		s/^\s+|\s+$//g;
	$seq .= $_;
	}

	return $seq;
}


sub usage {
	my( $msg ) = @_;
	$msg .= $msg ? "\n" : '';

	my $script = "\x1b[32m" . File::Spec->splitpath($0) . "\x1b[0m";
	return"$msg
$script version $VERSION; $YEAR; $AUTHOR

USAGE:
    $script [OPTIONS] --assembly <hg38|hg19> --breakpoint <geneA:chrA:strandA:locationA:geneB:chrB:strandB:locationB> --disease <disease> --pubmed <pubmed>

EXAMPLE:
    $script -assembly hg38 -breakpoint VTI1A:10:+:112464657:TCF7L2:10:+:113141184 -disease GBM -pubmed 32414213 -cell_line HeLa
or
    $script --start_id 1000 -assembly hg38 -breakpoint VTI1A:10:+:112464657:TCF7L2:10:+:113141184 -disease GBM -pubmed 32414213 -output OUT_fusion_sequence.csv

OPTIONS:
    --auto                      --  autocomplete options: --start_id, --assembly
    --output <file.tsv|stdout>  --  output table. By default, STDOUT
    --header                    --  output header of table
    --start_id                  --  Start Fusion ID (fusID). If --auto then 1000
    --assembly                  --  Required 'hg19' or 'hg38'. If --auto then 'hg19'
    --breakpoint                --  Required <geneA:chrA:strandA:locationA:geneB:chrB:strandB:locationB>
                                    chr = ([1-9XYM]|1[0-9]|2[012]|MT)
    --disease                   --  Required Disease Acronym(,s).
    --pubmed                    --  Required PubMed_ID(;s)
    --cell_line                 --  Optional
    --delta                     --  The length of output nucleotide sequence HEAD and TAIL. Optional, default is $DELTA bp
    --bsdir                     --  Defaul is '/home/aleksander/BioInf/Programs/B/BLAT/Chimeras_search/genomes/human'

OUTPUT TABLE FORMAT:
  1. fusID                      --  Fusion ID (inner). Required from INPUT --start_id option
  2. ncbi_A                     --  NCBI ID of gene_A (HEAD). Default is '.'
  3. gene_A                     --  gene symbol_A (HEAD). Required from INPUT --breakpoint option
  4. ncbi_B                     --  NCBI ID of gene_B (TAIL). By default, '.'
  5. gene_B                     --  gene symbol_B (TAIL). Required from INPUT --breakpoint option
  6. sequence_A                 --  nucleotide sequence of a gene_A at a fusion junction. Required
  7. sequence_B                 --  nucleotide sequence of a gene_B at a fusion junction. Required
  8. exon_A                     --  exon number of a gene_A at a fusion junction. Default is '.'
  9. exon_B                     --  exon number of a gene_B at a fusion junction. Default is '.'
  10.disease(,s)                --  Required from INPUT --disease option
  11.pubmed(;s)                 --  Required from INPUT --pubmed option
  12.cell_line(,s)              --  From INPUT --cell_line option. Default is empty ('')
";
}

__END__

=head1 SUPPLEMENTARY

=head2 Disease Acronyms

  +-------------+------------------------------------------------------------------+
  | Acronym     | Name                                                             |
  +-------------+------------------------------------------------------------------+
  | ACC         | Adrenocortical carcinoma                                         |
  | ALL         | Acute Lymphocytic Leukemia                                       |
  | ASD         | Autism spectrum disorder                                         |
  | BCC         | Basal Cell Carcinoma                                             |
  | BLCA        | Bladder Urothelial Carcinoma                                     |
  | BRCA        | Breast Invasive Carcinoma                                        |
  | CESC        | Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma |
  | CFC         | Cranial fasciitis of childhood                                   |
  | CHOL        | Cholangiocarcinoma                                               |
  | CLL         | Chronic Lymphocytic Leukemia                                     |
  | CML         | Chronic Myelogenous Leukemia                                     |
  | CNTL        | Controls                                                         |
  | COAD        | Colon adenocarcinoma                                             |
  | COADREAD    | Colon Adenocarcinoma and Rectum Adenocarcinoma                   |
  | CTCL        | Cutaneous T-Cell Lymphoma                                        |
  | DLBC        | Lymphoid Neoplasm Diffuse Large B-Cell Lymphoma                  |
  | EHE         | Epithelioid Hemangioendothelioma                                 |
  | ESCA        | Esophageal Carcinoma                                             |
  | ETT         | Epithelioid trophoblastic tumor                                  |
  | FPPP        | FFPE Pilot Phase II                                              |
  | GBC         | Gallbladder cancer                                               |
  | GBM         | Glioblastoma Multiforme                                          |
  | GBMLGG      | GBM + LGG Glioma                                                 |
  | GIST        | Gastrointestinal Stromal Tumor                                   |
  | HES         | Hypereosinophilic syndrome                                       |
  | HGG         | Brain High Grade Glioma                                          |
  | HL          | Hodgkin Lymphoma                                                 |
  | HNSC        | Head and Neck Squamous Cell Carcinoma                            |
  | IMT         | Inflammatory Myofibroblastic Tumor                               |
  | KICH        | Kidney Chromophobe                                               |
  | KIPAN       | Pan-kidney cohort (KICH+KIRC+KIRP)                               |
  | KIRC        | Kidney Renal Clear Cell Carcinoma                                |
  | KIRP        | Kidney Renal Papillary Cell Carcinoma                            |
  | LAML        | Acute Myeloid Leukemia                                           |
  | LCML        | Chronic Myelogenous Leukemia                                     |
  | LGG         | Brain Lower Grade Glioma                                         |
  | LIHC        | Liver Hepatocellular Carcinoma                                   |
  | LUAD        | Lung Adenocarcinoma                                              |
  | LUSC        | Lung Squamous Cell Carcinoma                                     |
  | MB          | Medulloblastoma                                                  |
  | MESO        | Mesothelioma                                                     |
  | MISC        | Miscellaneous                                                    |
  | MM          | Multiple Myeloma                                                 |
  | MPD         | Myeloproliferative Disorder                                      |
  | MPNST       | Malignant Peripheral Nerve Sheath Tumor                          |
  | NB          | Neuroblastoma                                                    |
  | NEC         | No Evidence of Cancer                                            |
  | NSCLC       | Non-Small Cell Lung Cancer                                       |
  | OV          | Ovarian Serous Cystadenocarcinoma                                |
  | PA          | Pilocytic astrocytoma                                            |
  | PAAD        | Pancreatic Adenocarcinoma                                        |
  | PCPG        | Pheochromocytoma and Paraganglioma                               |
  | PRAD        | Prostate Adenocarcinoma                                          |
  | READ        | Rectum adenocarcinoma                                            |
  | SARC        | Sarcoma                                                          |
  | SKCM        | Skin Cutaneous Melanoma                                          |
  | STAD        | Stomach Adenocarcinoma                                           |
  | STES        | Stomach + Esophageal carcinoma (ESCA)                            |
  | TGCT        | Testicular Germ Cell Tumors                                      |
  | THCA        | Thyroid Carcinoma                                                |
  | THYM        | Thymoma                                                          |
  | UCEC        | Uterine Corpus Endometrial Carcinoma                             |
  | UCS         | Uterine carcinosarcoma                                           |
  | UVM         | Uveal Melanoma                                                   |

=cut
