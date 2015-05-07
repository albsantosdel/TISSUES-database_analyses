#################################################################################################################################
#		This code generates all the files necessary to reproduce the figures shown in the manuscript			#
#			Comprehensive comparison of large-scale tissue expression datasets					#
#	Alberto Santos, Kalliopi Tsafou, Christian Stolte, Sune Pletscher-Frankild, Seán I. O’Donoghue and Lars Juhl Jensen	#
#################################################################################################################################
	
#!/usr/bin/ perl
use warnings;
use strict;

#Data directory
my $data_dir = "data/";
#Input files: this files need to be downloaded from FigShare (http://figshare.com/s/cb788d0ef4bd11e4b5ea06ec4b8d1f61)
# and stored in the data/datasets/ directory
my %options = ("unigene"=> $data_dir."datasets/unigene.tsv",
                "gnf"=> $data_dir."datasets/gnf.tsv",
                "exon"=> $data_dir."datasets/exon.tsv",
                "hpa"=> $data_dir."datasets/hpa_ihc.tsv",
                "hpa_rna"=> $data_dir."datasets/hpa_rna.tsv",
                "rna"=> $data_dir."datasets/rna_seq.tsv",
                "tm"=> $data_dir."datasets/text_mining.tsv",
                "hpm"=> $data_dir."datasets/hpm.tsv");

#dictionary files: These files are used for backtracking tissues to their respective parent tissues
# to identify which genes are expressed in the tissues of interest (more details in the README file)
my $labels_file = $data_dir."dictionary/labels.tsv";
my $bto_entities = $data_dir."dictionary/bto_entities.tsv";
my $bto_groups = $data_dir."dictionary/bto_groups.tsv";

#common tissues: List of common tissues to all datasets
my %common_tissues =("heart"=>1, "liver"=>1, "intestine"=>1, "nervous system"=>1, "kidney"=>1);

#Get the dataset name we are going to analyze   
my $dataset_file = "";

#Gold standard
my $goldstandard = "uniprot";
my $goldstandard_file = $data_dir."datasets/uniprot.tsv";

#sliding window
my $window_size = 100;

#mRNA datasets necessary to generate the mRNA reference set: Gene-tissue associations in at least half of the datasets 
# and a expression higher that these cutoffs will be included in the mRNA reference set 
my %mRNA = ("unigene"=> 20,
       	"gnf"=> 250,
	"exon"=> 250,
	"hpa_rna"=> 20,
	"rna"=> 5);

#Cutoffs used in the publication: [low, medium, high]
my %cutoffs = ("unigene"=>[0.9,10,20],
	"gnf"=> [50,100,250],
	"exon"=> [50,100,250],
	"hpa"=>[0,5.9,9.9],
	"hpa_rna"=>[1,10,20],
	"rna"=>[0.5,1,5],
	"uniprot"=>[0,0,0],
	"tm"=>[0,2.5,3.5],
	"hpm"=>[0,5,15]);

#output files
my $major_tissues_file = $data_dir."datasets_major_tissues.tsv"; #Specifies which datasets have which major tissues within the ones specified in the labels.tsv file
my $expression_breadth_file = "cutoff_expression_breadth.tsv"; #Expression breadth for each of the datasets in each cutoff
my $consistency_file = "consistency_analysis.tsv"; #Gene-tissue associations per dataset per cutoff
my $mRNA_reference_set = $data_dir."mRNA_reference_set.tsv"; #Contains the gene-tissue pairs that form the mRNA reference set used in the comparisons for confirmation of trends
my $hpa_single_antibody_prots = $data_dir."hpa_single_antibody_proteins.tsv"; #Proteins with a single antibody in HPA Immunohistochemistry experiments

##############
#  Execution #
##############
my ($labels) = &parse_labels_file(\%common_tissues);
my ($entities) = &parse_bto_entities_file();
my ($groups) = &parse_bto_groups_file();
my ($child_labels_hash) = &get_child_label($labels,$entities,$groups);

#dictionary with the number of proteins per tissue
my %major_tissues = ();

my ($goldstandard_data,$gold_btos) = &parse_goldstandard_file();
my ($goldstandard_labels) = &convert_btos_labels($goldstandard_data,$child_labels_hash, \%major_tissues, $goldstandard);


#dictionary with the number of tissues per protein per level (cutoff: low, ,medium, high)
my %express_breadth = ();
#dictionary with the protein-tissue pairs per level (cutoff: low, medium, high) needed for the consistency analysis
my %consistency = ();
&get_exp_breadth_and_consistency_analyses_files($goldstandard_labels,\%express_breadth, \%consistency, $goldstandard, \%cutoffs);

#filter common tissues in gold-standard
my ($goldstandard_filtered) = &filter_common_tissues($goldstandard_labels, \%common_tissues);

#dictionary with the mRNA reference set
my %reference_set = ();
#dictionary with all the data available
my %dataset_scored_pairs = ();

foreach my $dataset_name (keys %options){
	$dataset_file = $options{$dataset_name};
	
	my ($dataset, $single_antibodies) = &parse_studied_dataset_file($dataset_name);
	if($dataset_name eq "hpa"){
		&print_2dim($single_antibodies, $hpa_single_antibody_prots);
	}

	my $dataset_labels = undef;
	my $major_tissues = undef;
	if($dataset_name eq "unigene"){
		($dataset_labels) = &get_label_values($dataset,$labels, \%major_tissues, $dataset_name);
	}
	else{
		($dataset_labels) = &convert_btos_labels($dataset,$child_labels_hash, \%major_tissues, $dataset_name);
	}
	
	&format_score_pairs($dataset_labels, \%common_tissues, \%dataset_scored_pairs, $dataset_name);
	&get_exp_breadth_and_consistency_analyses_files($dataset_labels,\%express_breadth, \%consistency, $dataset_name, \%cutoffs);

	if (exists $mRNA{$dataset_name}){
		&build_mRNA_references_set($dataset_labels, \%reference_set, $mRNA{$dataset_name});
	}
}

my ($reference_set_filtered) = &filter_common_tissues(\%reference_set, \%common_tissues);
my %gold_standards = ("uniprot" => $goldstandard_filtered, "mRNA_reference_combined" => $reference_set_filtered);
#calculate fold enrichment for the different datasets
foreach my $goldstandard (keys %gold_standards){
	foreach my $dataset_name (keys %dataset_scored_pairs){
		my $output_file = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis.tsv";	
		&calculate_fold_enrichment(\%{$dataset_scored_pairs{$dataset_name}}, $gold_standards{$goldstandard}, $output_file, $window_size);
	}
}

#Print out the file with: Dataset Tissue NumberProteins
&print_2dim(\%major_tissues, $major_tissues_file);
#Print out the file with: Dataset Protein NumberTissues
&print_3dim(\%express_breadth, $expression_breadth_file, $data_dir);
#Print out consistency files
&print_3dim(\%consistency, $consistency_file, $data_dir);
#Print out mRNA reference set
my $min_datasets = int((keys %mRNA)/2) + 1;#In at least more than half of the datasets
&print_mRNA_reference_set(\%reference_set, $mRNA_reference_set, $min_datasets);

##############
#  Functions #
##############
sub print_2dim(){
	my ($data) = $_[0];
	my ($output_file) = $_[1];
	
	open(OUT,">$output_file") or die "Unable to open the Output file $output_file\n";
	foreach my $dataset (keys %{$data}){
		foreach my $tissue (keys %{${$data}{$dataset}}){
			print OUT $dataset."\t".$tissue."\t".${$data}{$dataset}{$tissue}."\n";
		}
	}
	close(OUT);
}

sub print_3dim(){
	my ($data) = $_[0];
	my ($output_file) = $_[1];
	my ($data_dir) = $_[2];
	
	foreach my $level (keys %{$data}){
		my $full_path = $data_dir.$level."_".$output_file;
		open(OUT,">$full_path") or die "Unable to open the Output file $full_path\n";
		foreach my $dataset (keys %{${$data}{$level}}){
			foreach my $id (keys %{${$data}{$level}{$dataset}}){
				print OUT $dataset."\t".${$data}{$level}{$dataset}{$id}."\t".$id."\n";
			}
		}
	}
	close(OUT);
}

sub print_mRNA_reference_set(){
	my ($data) = $_[0];
	my ($output_file) = $_[1];
	my $cutoff = $_[2];
	open(OUT,">$output_file") or die "Unable to open the Output file $output_file\n";
		
	foreach my $a (keys %{$data}){
		foreach my $b (keys %{${$data}{$a}}){
			if (${$data}{$a}{$b} >= $cutoff){
				print OUT $a."\t".$b."\n";
			}
		}
	}
	close(OUT);
}

sub parse_labels_file(){
    print "- Parsing the Labels file $labels_file with the tissues to be studied\n";
    
    open(LABELS,"$labels_file") or die "Unable to open $labels_file file\n";
    
    my %labels = ();
    my $tissues_type_code = "-25";
    
    while(<LABELS>){
        s/\r?\n//;
        my ($type,$name,$code) = split("\t",$_);
        if($type eq $tissues_type_code){
            $labels{$code} = $name;
        }
    }
    print "- The labels file has been correctly parsed\n";
    return \%labels;
}

sub parse_bto_entities_file(){
    print "- Parsing the BTO entities file $bto_entities\n";
    
    open(ENTITIES,"$bto_entities") or die "Unable to open $bto_entities file\n";
    
    my %entities = ();
    my $tissues_type_code = "-25";
    
    while(<ENTITIES>){
        s/\r?\n//;
        my ($id,$type,$bto) = split("\t",$_);
        if($type eq $tissues_type_code){
            $entities{$id} = $bto;
        }
    }
    print "- The BTO entities file has been correctly parsed\n";
    return \%entities;
}

sub parse_bto_groups_file(){
    print "- Parsing the BTO groups file $bto_groups\n";
    
    open(GROUPS,"$bto_groups") or die "Unable to open $bto_groups file\n";
    
    my %groups = ();
    
    while(<GROUPS>){
        s/\r?\n//;
        my ($id,$parent_id) = split("\t",$_);
        $groups{$id}{$parent_id} = 1;
    }
    print "- The BTO groups file has been correctly parsed\n";
    return \%groups;
}

sub get_child_label(){
    my ($labels) = $_[0];
    my($entities) = $_[1];
    my ($groups) = $_[2];
    
    my %child_labels = ();
    foreach my $id (keys %{$entities}){
        my $bto = ${$entities}{$id};
        foreach my $parent_id (keys %{${$groups}{$id}}){
            my $parent_bto = ${$entities}{$parent_id};
            if(exists ${$labels}{$parent_bto}){
                my $label = ${$labels}{$parent_bto};
                $child_labels{$bto}{$label} = 1;
            }
        }
        if(exists ${$labels}{$bto}){
            $child_labels{$bto}{${$labels}{$bto}} =1;
        }
    }
    
    return \%child_labels;
}

sub parse_studied_dataset_file(){
    	my $dataset_name = $_[0];
    	print "- Parsing the Dataset file $dataset_file\n";
    	open (DATASET,"$dataset_file") or die "Unable to open the Dataset file $dataset_file\n";
    	my %dataset = ();
    	my %single_antibodies = ();#only applicable to HPA-IHC
    	while(<DATASET>){
		s/\r?\n//;
		
		my (undef,$ensp,undef,$bto,$txt_aux,$score_scoretype,$stars,undef) = split("\t",$_);
		my $score = $stars;
		if($dataset_name ne "hpa" and $dataset_name ne "tm"){
	    		($score,undef) = split('\s',$score_scoretype);
		}
		elsif($dataset_name eq "tm"){
	    		$score = $txt_aux;
	    		$stars = $score_scoretype;
		}
		else{
			if(index($score_scoretype, '\\n') == -1 && index($score_scoretype, 'antibodies') == -1){
				$single_antibodies{$ensp}{$bto} = $score_scoretype;
			}
			$score= $stars * 5.5;       
	       	}
		next if $score == 0;
		$dataset{$ensp}{$bto} = "$score\t$stars";        
    	}
    	
    	print "- The Dataset file has been parsed\n";
    
	return \%dataset,\%single_antibodies;
}

sub parse_goldstandard_file(){
    print "- Parsing the goldstandard file $goldstandard_file\n";
    
    open (GOLD,"$goldstandard_file") or die "Unable to open the UniProt file $goldstandard_file\n";
    
    my %goldstandard_data = ();
    my %gold_btos = (); 
    
    while(<GOLD>){
        s/\r?\n//;
        my (undef,$ensp,undef,$bto,undef) = split("\t",$_);
        $goldstandard_data{$ensp}{$bto} = 1;
        $gold_btos{$bto} = 1;
    }
    print "- The gold standard file has been parsed\n";
    return \%goldstandard_data,\%gold_btos;
}

sub get_label_values(){
    	my ($dataset) = $_[0];
    	my ($labels) = $_[1];
	my ($major_tissues) = $_[2];
	my ($dataset_name) = $_[3];

    	print "- Getting the labeled btos of the UniGene data\n";
    	my %dataset_labels = ();
    	foreach my $ensp (keys %{$dataset}){
		foreach my $bto (keys %{${$dataset}{$ensp}}){
	    		if(exists ${$labels}{$bto}){
				my $label = ${$labels}{$bto};
				my ($score,$star) = split("\t",${$dataset}{$ensp}{$bto});
				if(exists $dataset_labels{$ensp}{$label}){
		    			my ($p_score,$p_stars) = split("\t",$dataset_labels{$ensp}{$label});
		    			$dataset_labels{$ensp}{$label} = ($dataset_labels{$ensp}{$label}, ${$dataset}{$ensp}{$bto})[$p_score<$score];
					
				}
				else{
					$dataset_labels{$ensp}{$label} = ${$dataset}{$ensp}{$bto};
				}
				${$major_tissues}{$dataset_name}{$label}++;
	    		}
		}
    	}
    	print "- The parent btos have been associated to each ENSP\n";
    
	return \%dataset_labels;
}

sub convert_btos_labels(){
    my ($data) = $_[0];
    my ($btos_labels) = $_[1];
    my ($major_tissues) = $_[2];
    my ($dataset_name) = $_[3];

    my %ensp_labels = ();
    my %labels = ();
    foreach my $ensp (keys %{$data}){
        foreach my $bto (keys %{${$data}{$ensp}}){
            foreach my $label (keys %{${$btos_labels}{$bto}}){
                my ($score,$star) = split("\t",${$data}{$ensp}{$bto});
                if(exists $ensp_labels{$ensp}{$label}){
                    my ($p_score,$p_stars) = split("\t",$ensp_labels{$ensp}{$label});
                    $ensp_labels{$ensp}{$label} = ($ensp_labels{$ensp}{$label},${$data}{$ensp}{$bto})[$p_score<$score];
                }
                else{
                    $ensp_labels{$ensp}{$label} = ${$data}{$ensp}{$bto};
                }
                $labels{$label}= 1;
		${$major_tissues}{$dataset_name}{$label}++;
            }
        }
    }
    
    return \%ensp_labels;
}

sub format_score_pairs(){
    	my ($dataset_labels) = $_[0];
	my ($common_tissues) = $_[1];
	my ($score_pairs) = $_[2];
	my ($dataset_name) = $_[3];

    	foreach my $ensp (keys %{$dataset_labels}){
		foreach my $label (keys %{${$dataset_labels}{$ensp}}){
			if (exists ${$common_tissues}{$label}){
				my ($score,$stars) = split("\t",${$dataset_labels}{$ensp}{$label});
				${$score_pairs}{$dataset_name}{$score}{$ensp}{$label} = $stars;
			}
		}
    	}
}

sub filter_common_tissues(){
	my ($data) = $_[0];
	my ($common_tissues) = $_[1];
	
	my %filtered = ();

	foreach my $ensp (keys %{$data}){
		foreach my $label (keys %{${$data}{$ensp}}){
			if(exists ${$common_tissues}{$label}){
				$filtered{$ensp}{$label} = 1;
			}
		}
	}

	return \%filtered;
}

sub get_exp_breadth_and_consistency_analyses_files(){
	my ($dataset_labels) = $_[0];
	my ($express_breadth) = $_[1];
	my ($consistency) = $_[2];
	my ($dataset) = $_[3];
	my ($cutoffs) = $_[4];

	my @cutoffs = @{${$cutoffs}{$dataset}};
	
	foreach my $ensp (keys %{$dataset_labels}){
		foreach my $label (keys %{${$dataset_labels}{$ensp}}){
			my ($score, $stars) = split("\t",${$dataset_labels}{$ensp}{$label});
			if ($score > $cutoffs[0]){
				${$express_breadth}{"low"}{$dataset}{$ensp}++;
				${$consistency}{"low"}{$dataset}{$label."\t".$ensp} = $score;
			}
			if ($score > $cutoffs[1]){
				${$express_breadth}{"medium"}{$dataset}{$ensp}++;
				${$consistency}{"medium"}{$dataset}{$label."\t".$ensp} = $score;
			}
			if ($score > $cutoffs[2]){
				${$express_breadth}{"high"}{$dataset}{$ensp}++;
				${$consistency}{"high"}{$dataset}{$label."\t".$ensp} = $score;
			}
		}
	}
}

#The mRNA references set is formed by all the high confident gene-tissue pairs that are found in at least half 
# ofthe transcriptomics datasets
sub build_mRNA_references_set(){
	my ($dataset_labels) = $_[0];
	my ($reference_set) = $_[1];
	my $cutoff = $_[2];

	foreach my $ensp (keys %{$dataset_labels}){
		foreach my $label (keys %{${$dataset_labels}{$ensp}}){
			my ($score, $stars) = split("\t",${$dataset_labels}{$ensp}{$label});
			if ($score > $cutoff){
				${$reference_set}{$ensp}{$label}++;
			}
		}
	}
}

#Calculate the fold enrichment for each dataset using a moving average window over the score
#This allowed making the datasets comparable
sub calculate_fold_enrichment(){
    	my ($dataset) = $_[0];
    	my ($gold) = $_[1];
    	my ($output_file) = $_[2];
    	my ($window_size) = $_[3];
    	
    	print "- Calculating precision at the level of score for the dataset\n";
    	my @pos_neg_array = ();
    	my @scores = ();
    	my @stars = ();
    	my $dataset_total_pairs = 0;
    	my %gold_valid_prots = ();
    	foreach my $score (sort {$a <=> $b} keys %{$dataset}){
		foreach my $ensp (keys %{${$dataset}{$score}}){
	    		if(exists ${$gold}{$ensp}){ 
		 		foreach my $label (keys %{${$dataset}{$score}{$ensp}}){
					$gold_valid_prots{$ensp} = 1;
					$dataset_total_pairs++;
					my $star = ${$dataset}{$score}{$ensp}{$label};
					if(exists ${$gold}{$ensp}{$label}){
						push @pos_neg_array,1;
					}
					else{
						push @pos_neg_array,0;
					}
					push @scores,$score;
					push @stars,$star;
				}
	    		}
		}
    	}
    	print scalar keys %gold_valid_prots;
    	print "- The dataset has $dataset_total_pairs pairs in total\n";
    	
	
    	#calculate gold_standard pairs and proteins
	my $gold_num_pairs = 0;
 	my %gold_tissues = ();
    	foreach my $ensp (keys %gold_valid_prots){
		foreach my $label (keys %{${$gold}{$ensp}}){
			$gold_tissues{$label} = 1;
	    		$gold_num_pairs++;	
		}
    	}
    	my $gold_num_prots = scalar keys %gold_valid_prots;
    	my $gold_num_tissues = scalar keys %gold_tissues;
    	my $gold_prob = $gold_num_pairs/($gold_num_tissues*$gold_num_prots);
     
	
    	print "- Gold standard number pairs: $gold_num_pairs\n";
    	print "- Gold standard number proteins: $gold_num_prots\n";
    	print "- Gold standard number tissues: $gold_num_tissues\n";
    	print "- Gold standard probability: $gold_prob\n";
    	
    	use List::Util qw(sum);
    	open(OUT,">$output_file") or die "Unable to open the Output file $output_file\n";
    	
    	my $i = 0;
    	my $j = $window_size -1;
    	my $last = 0;
    	
    	my $sum = sum(@pos_neg_array[$i..$j]);
    	my $fold_enrichment = ($sum/$window_size)/$gold_prob;
    	my $mean_stars = sum(@stars[$i..$j])/$window_size;
    	my $mean_score = sum(@scores[$i..$j])/$window_size;
    	
    	print OUT "$mean_stars\t$fold_enrichment\t$mean_score\n";
    	
    	$i = $j + 1;
    	if ($j + $window_size>(scalar @pos_neg_array-1)){
		$j = scalar @pos_neg_array-1;
		$last = 1;
    	}
    	else{
		$j = $j + $window_size;
    	}    
    	while($j<(scalar @pos_neg_array-1)){
		$sum = sum(@pos_neg_array[$i..$j]);
		$fold_enrichment = ($sum/($j-$i+1))/$gold_prob;
		$mean_stars = sum(@stars[$i..$j])/($j-$i+1);
		$mean_score = sum(@scores[$i..$j])/($j-$i+1);
		$i=$j+1;
		if ($j + $window_size>(scalar @pos_neg_array-1) and !$last){
	    		$j = scalar @pos_neg_array-1;
	    		$last = 1;
		}
		else{
	    		$j = $j + $window_size;
		}   
		print OUT "$mean_stars\t$fold_enrichment\t$mean_score\n";
    	}
    	close(OUT);
	
    	print "- The fold_enrichment values have been calculated and stored in the $output_file file\n";
}
