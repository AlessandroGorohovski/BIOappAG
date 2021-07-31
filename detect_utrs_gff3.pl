#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

++$|;

our $VERSION = '0.03';

###
# Default Options
our $UNIQ;
our $UTR;
my $HELP;
my $HOST_GFF;
my $EMSG = "\x1b[31mERROR\x1b[0m:";

###
# Parse input data
GetOptions(
	'help'       => \$HELP,
	'uniq'       => \$UNIQ,
	'utr'        => \$UTR,
	'host_gff'   => \$HOST_GFF,
) or &usage();

&usage() if $HELP;

# input gff3 file
my $ingff = $ARGV[0] || &usage("$EMSG GFF3 file not specified for processing!");

my @GFF_fields = qw(seq_id source seq_type start end score strand phase attributes);

my %gene;
my %utr_data;
our $gff_line;

open INGFF, $ingff or &usage("$EMSG can't open file $ingff!");
if( ($_=<INGFF>) =~/^##gff\-version\s+3\D?/i){	# 1st line
	print unless $UTR;
}else{
	&usage("$EMSG input file '$ingff' is not GFF3 or of unknown type!") 
}
while( $gff_line = <INGFF>){

	if( eof(INGFF) or $gff_line =~/^###+/ ){ # END of FILE or SEGMENT
		&collect_utrs( \%gene, \@GFF_fields, \%utr_data ) if %gene;

		undef %gene;

	}elsif( $gff_line =~/^#/ ){ # COMMENT or PRAGMA
		unless( defined $HOST_GFF ){
			$HOST_GFF = uc($1) if $gff_line =~/^#!(?:processor\s+(NCBI)\s+annotwriter|genome\-build\s+(Ensembl)\s+)/i;
		}

	}elsif( my $gff_record = &create_gff_record( \@GFF_fields ) ){
		&analyze_gff_record( $gff_record, \%gene, \@GFF_fields, \%utr_data );
	}

	print $gff_line unless $UTR;
}
close INGFF;

exit;


sub collect_utrs {
	my( $gene, $GFF_fields, $utr_data ) = @_;
	return undef if ! exists $gene->{'transcripts'};

	my @utrs;

	# Анализируем транскрипты, чтобы найти их UTRs
	for my $id ( sort keys %{ $gene->{'transcripts'} } ){

		# самые крайние границы CDS
		my $cds_Lstart = $gene->{'transcripts'}{$id}{'cds_Lstart'};	# left
		my $cds_Rend   = $gene->{'transcripts'}{$id}{'cds_Rend'};	# right

		# по цепи (strand), обозначить тип UTR
		my( $left_utr_type, $right_utr_type ) = $gene->{'transcripts'}{$id}{'strand'} eq '+' ?
																('five_prime_UTR',  'three_prime_UTR') :
																('three_prime_UTR', 'five_prime_UTR');

		# Перебираем экзоны
		for my $exon ( @{ $gene->{'transcripts'}{$id}{'exons'} } ){

			# ищем экзоны, которые составляют левый UTR
			if( $exon->{'start'} < $cds_Lstart ){
				my $utr_start = $exon->{'start'};
				my $utr_end = ($exon->{'end'} >= $cds_Lstart) ? $cds_Lstart - 1 : $exon->{'end'};

				&add_utr( $utr_start, $utr_end, $left_utr_type, $exon, $utr_data, \@utrs );
			}

			# ищем экзоны, которые составляют правый UTR
			if( $exon->{'end'} > $cds_Rend ){
				my $utr_start = ($exon->{'start'} <= $cds_Rend) ? $cds_Rend + 1 : $exon->{'start'};
				my $utr_end = $exon->{'end'};

				&add_utr( $utr_start, $utr_end, $right_utr_type, $exon, $utr_data, \@utrs );
			}
		}

	}

	&output_utrs( \@utrs, $GFF_fields );

	# Reset gene
	undef %{ $gene };
}


=comment
	Create GFFRecords of UTRs from: utr_start, utr_end, utr_type, exon
	Add it to the list: utrs (so it appends to the external list variable, does not return anything)
---
	Создать GFFRecords UTR из: utr_start, utr_end, utr_type, exon
	Добавляет её в список: utrs (чтобы он добавлялся во внешнюю переменную списка и ничего не возвращал)
=cut

sub add_utr {
	my( $utr_start, $utr_end, $utr_type, $exon, $utr_data, $utrs ) = @_;

	my $exon_attrs = &parse_gff_attributes( $exon->{'attributes'} );

	my $id_str = 'ID=utr';
	$id_str .= $1 if $exon_attrs->{'ID'} =~/(\d+)/;

	my $k = $exon->{'seq_id'} . lc($utr_type) . $utr_start . $exon->{'strand'} . $utr_end;
	$k .= lc($id_str) unless $UNIQ;

	return if exists $utr_data->{$k};
	$utr_data->{$k} = undef;

	my @utr_attrs = ( $id_str );
	push @utr_attrs, map { exists( $exon_attrs->{$_} ) ? "$_=$exon_attrs->{$_}" : ( ) }
									qw( Parent transcript_id Dbxref partial start_range end_range );

	my $utr_attrs_str = join ';', @utr_attrs;

	push @{ $utrs }, &prepare_gff_record(
								$exon->{'seq_id'},
								$exon->{'source'},
								$utr_type,
								$utr_start,
								$utr_end,
								'.',	# score
								$exon->{'strand'},
								'.',	# phase
								$utr_attrs_str
							);

}


# Создать и возвратить GFF record из строки, только, если она содержит нужные записи
sub create_gff_record {
	my( $GFF_fields ) = @_;
	return undef if $gff_line =~/^\s*$/;	# Empty line

# seq_id       	source	seq_type    	start	end	score	strand	phase	attributes
# 0            	1     	2           	3  	4  	5	6	7	8
# NW_011647988.1	Gnomon	C_gene_segment	49105	51791	.	-	.	ID=rna11;Parent=gene8;Dbxref=GeneID:105180724,Genbank:XM_011136973.1;Name=XM_011136973.1;gbkey=mRNA;gene=LOC105180724;product=uncharacterized LOC105180724%2C transcript variant X1;transcript_id=XM_011136973.1

	my @fields = split /\t/, $gff_line;
	return undef if ~~@fields != ~~@{ $GFF_fields };	# количество полей не соответствует ожидаемому кол-ву @GFF_FIELDS

	# Проверка seq_type на присутствие одного из допустимых
	return undef if $fields[2] !~/^(?:gene|mRNA|CDS|exon|[CVDJ]_gene_segment|(?:fiv|thre)e_prime_UTR)$/i;

	return &prepare_gff_record( @fields );
}


#< Create GFF record from the provided variables
#> Создать запись GFF из предоставленных переменных
sub prepare_gff_record {
	my( $seq_id, $source, $seq_type, $start, $end, $score, $strand, $phase, $attributes ) = @_;
	$attributes =~s/^\s+|\s+$//g;

	return {
		'seq_id'   => ($seq_id eq '.')   ? undef : $seq_id,
		'source'   => ($source eq '.')   ? undef : $source,
		'seq_type' => ($seq_type eq '.') ? undef : $seq_type,
		'start'    => ($start eq '.')    ? undef : int( $start ),
		'end'      => ($end eq '.')      ? undef : int( $end ),
		'score'    => ($score eq '.')    ? undef : $score,
		'strand'   => ($strand eq '.')   ? undef : $strand,
		'phase'    => ($phase eq '.')    ? undef : $phase,
		'attributes' => defined( $attributes ) && length( $attributes ) ? $attributes : undef,
	};

}


#< Analyze GFF record
#> Проанализировать запись GFF
sub analyze_gff_record {
	my( $gff_record, $cur_gene, $GFF_fields, $utr_data ) = @_;
	my $seq_type = $gff_record->{'seq_type'};

	if( $seq_type =~/^gene$/i ){

		# собрать и вывести UTRs для уже найденного гена, если есть
		&collect_utrs( $cur_gene, $GFF_fields, $utr_data ) if exists $cur_gene->{'transcripts'};

		# create new gene
		$cur_gene->{'transcripts'} = undef;

	}elsif( $seq_type =~/^(?:mRNA|[CVDJ]_gene_segment)$/i ){
		# Это mRNA или эквивалентный тип записи: создаем transcript и добавляем его к транскриптам текущего гена
		&GFFtranscript( $gff_record, $cur_gene );

	}elsif( $seq_type =~/^exon$/i ){
		# Это exon: add it to current transcript's exons
		&GFFexon( $gff_record, $cur_gene->{'transcripts'} );

	}elsif( $seq_type =~/^CDS$/i ){
		# Это CDS: create 'cdss' record with leftmost 'Lstart', rightmost 'Rend' of current gene
		&GFFcds( $gff_record, $cur_gene->{'transcripts'} );

	}elsif( $seq_type =~/^(?:fiv|thre)e_prime_UTR$/i ){
		&GFFutr( $gff_record, $utr_data );
	}

}


sub GFFtranscript {
	my( $gff_record, $cur_gene ) = @_;

	my $retval = &parse_gff_attributes( $gff_record->{'attributes'} );
	return undef unless $retval;

	my $id = $retval->{'ID'};
	$cur_gene->{'transcripts'}{ $id }{'strand'} = $gff_record->{'strand'};
}


# Создать структуру данных для хранения экзотов транскрипта
sub GFFexon {
	my( $gff_record, $transcript ) = @_;

	my $retval = &parse_gff_attributes( $gff_record->{'attributes'} );
	return undef unless $retval;

	my $id = $retval->{'Parent'};
	push @{ $transcript->{ $id }{'exons'} }, $gff_record;
}


# Создать структуру данных для хранения крайних границ CDS
sub GFFcds {
	my( $gff_record, $transcript ) = @_;

	my $retval = &parse_gff_attributes( $gff_record->{'attributes'} );
	return undef unless $retval;

	my $id = $retval->{'Parent'};

	# получить самую крайнюю левую границу CDS (левый UTR)
	$transcript->{ $id }{'cds_Lstart'} = $gff_record->{'start'} if ! exists($transcript->{ $id }{'cds_Lstart'}) or
																						$transcript->{ $id }{'cds_Lstart'} > $gff_record->{'start'};

	# получить самую крайнюю правую границу CDS (правый UTR)
	$transcript->{ $id }{'cds_Rend'} = $gff_record->{'end'} if ! exists($transcript->{ $id }{'cds_Rend'}) or
																					$transcript->{ $id }{'cds_Rend'} < $gff_record->{'end'};
}


sub GFFutr {
	my( $gff_record, $utr_data ) = @_;

	my $retval = &parse_gff_attributes( $gff_record->{'attributes'} );
	return undef unless $retval;

	my $id_str = 'id=' . $retval->{'ID'};

	my $k = $gff_record->{'seq_id'} . lc($gff_record->{'seq_type'}) . join '', @{$gff_record}{'start','strand','end'};
	$k .= lc($id_str) unless $UNIQ;

	return if exists $utr_data->{$k};
	$utr_data->{$k} = undef;

	print $gff_line if $UTR;
}


# Разбираем столбец атрибутов GFF3 и возвращаем словарь со всеми атрибутами.
sub parse_gff_attributes {
	my( $attr_str ) = @_;

	return undef if ! defined($attr_str) or length($attr_str) < 4;

	# hash of return values
	my %retval;

# ID=cds8;Parent=rna11;Dbxref=GeneID:105180724,Genbank:XP_011135275.1;Name=XP_011135275.1;gbkey=CDS;gene=LOC105180724;product=uncharacterized protein LOC105180724;protein_id=XP_011135275.1
	for( split /;/, $attr_str ){
		next if length($_) < 3;

		my( $key, $value ) = split /\s*=\s*/;
		next unless $value;

		$retval{ $key } = $value;
	}

	return \%retval;
}


sub output_utrs {
	my( $utrs, $GFF_fields ) = @_;
	return unless $utrs;

	for my $utr ( @{ $utrs } ){
		print join("\t", map{ $utr->{$_} || '.'} @{ $GFF_fields } ), "\n";
	}
}


sub usage {
	my( $msg ) = @_;
	$msg .= $msg ? "\n" : '';

	my $script = "\x1b[32m" . $0 . "\x1b[0m";

	die <<EOF;
$msg
$script version $VERSION
   Detection of 5' and 3' UTR features from GFF3 data.

USAGE:
   $script <file.gff> [OPTIONS]

EXAMPLE:
   $script  input.gff > input_with_UTRs.gff

   $script  input.gff -utr > output_UTRs_only.gff

HERE:
   <file.gff>  --  input NCBI's or ENSEMBL's GFF3 file

OPTIONS:
   --uniq      --  find unique UTR features only
   --utr       --  filtered output of UTR features
   --host_gff  --  GFF3 file from 'NCBI' or 'ENSEMBL'
   --help
EOF

}
