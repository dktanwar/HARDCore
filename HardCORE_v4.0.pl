#! /usr/bin/perl

use Cwd;
use List::Util qw(first);

sub ProkkaPrep{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	my @prokka_files = glob "*faa.prokka";
	foreach(@prokka_files){
		my $Genome = $_;
		$Genome =~ s{\.[^.]+$}{}; #strip .faa
		$Genome =~ s{\.[^.]+$}{}; #strip .prokka
		open(IN, "<$_") or die;
		open(OUT,">$Genome.mod.faa") or die;

		while(<IN>){
			chomp();
			my $line = $_;
			if($_ =~ m/>/){
				my @split = split(/\_/,$line);
				my @split2 = split(" ",$split[1]);
				my $prokka_id = $Genome . "_" . $split2[0];
				my $new = ">gi|$prokka_id|ref|1234| $prokka_id "."[$Genome]";
				print OUT $new . "\n";
			}#end if
			else{print OUT $line . "\n";}
		}#end while IN
		close(IN);
		close(OUT);
		system("rm $Genome.faa.prokka");
		system("mv $Genome.mod.faa $Genome.faa");
	}


}#end ProkkaPrep

sub ProkkaPrep_ffn{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	my @prokka_files = glob "*ffn.prokka";
	foreach(@prokka_files){
		my $Genome = $_;
		$Genome =~ s{\.[^.]+$}{}; #strip .faa
		$Genome =~ s{\.[^.]+$}{}; #strip .prokka
		open(IN, "<$_") or die;
		open(OUT,">$Genome.mod.ffn") or die;

		while(<IN>){
			chomp();
			my $line = $_;
			if($_ =~ m/>/){
				my @split = split(/\_/,$line);
				my @split2 = split(" ",$split[1]);
				my $prokka_id = $Genome . "_" . $split2[0];
				my $new = ">gi|$prokka_id|ref|1234| $prokka_id "."[$Genome]";
				print OUT $new . "\n";
			}#end if
			else{print OUT $line . "\n";}
		}#end while IN
		close(IN);
		close(OUT);
		system("rm $Genome.ffn.prokka");
		system("mv $Genome.mod.ffn $Genome.ffn");
	}


}#end ProkkaPrep_ffn

sub GeneMark{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	my @prokka_files = glob "*-GeneMark.faa";
	foreach(@prokka_files){
		my $Genome = $_;
		$Genome =~ s{\.[^.]+$}{}; #strip .faa
		$Genome =~ s{\-[^-]+$}{}; #strip .prokka
		open(IN, "<$_") or die;
		open(OUT,">$Genome.faa") or die;

		while(<IN>){
			chomp();
			my $line = $_;
			if($_ =~ m/>/){
				my @split = split(/\|/,$line);
				my @split2 = split(/\_/,$split[0]);
				my $Gene_mark_id = $Genome . "_" . $split2[1];
				my $new = ">gi|$Gene_mark_id|ref|1234| $Gene_mark_id "."[$Genome]";
				print OUT $new . "\n";
			}#end if
			else{print OUT $line . "\n";}
		}#end while IN
		close(IN);
		close(OUT);
		system("rm $Genome" . "-GeneMark.faa");
	}


}#end GeneMark

sub GeneMark_ffn{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	my @prokka_files = glob "*-GeneMark.ffn";
	foreach(@prokka_files){
		my $Genome = $_;
		$Genome =~ s{\.[^.]+$}{}; #strip .ffn
		$Genome =~ s{\-[^-]+$}{}; #strip .prokka
		open(IN, "<$_") or die;
		open(OUT,">$Genome.ffn") or die;

		while(<IN>){
			chomp();
			my $line = $_;
			if($_ =~ m/>/){
				my @split = split(/\|/,$line);
				my @split2 = split(/\_/,$split[0]);
				my $Gene_mark_id = $Genome . "_" . $split2[1];
				my $new = ">gi|$Gene_mark_id|ref|1234| $Gene_mark_id "."[$Genome]";
				print OUT $new . "\n";
			}#end if
			else{print OUT $line . "\n";}
		}#end while IN
		close(IN);
		close(OUT);
		system("rm $Genome" . "-GeneMark.ffn");
	}


}#end GeneMark_ffn

sub Glimmer{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	my @prokka_files = glob "*.glimmer.faa";
	foreach(@prokka_files){
		my $Genome = $_;
		$Genome =~ s{\.[^.]+$}{}; #strip .faa
		$Genome =~ s{\.[^.]+$}{}; #strip .prokka
		open(IN, "<$_") or die;
		open(OUT,">$Genome.faa") or die;

		while(<IN>){
			chomp();
			my $line = $_;
			if($_ =~ m/>/){
				my @split = split(/\>/,$line);
				my @split2 = split(" ",$split[1]);
				my $Gene_mark_id = $Genome . "_" . $split2[0];
				my $new = ">gi|$Gene_mark_id|ref|1234| $Gene_mark_id "."[$Genome]";
				print OUT $new . "\n";
			}#end if
			else{print OUT $line . "\n";}
		}#end while IN
		close(IN);
		close(OUT);
		system("rm $Genome" . "*.glimmer.faa");
	}


}#end Glimmer

sub Glimmer_ffn{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	my @prokka_files = glob "*.glimmer.ffn";
	foreach(@prokka_files){
		my $Genome = $_;
		$Genome =~ s{\.[^.]+$}{}; #strip .faa
		$Genome =~ s{\.[^.]+$}{}; #strip .prokka
		open(IN, "<$_") or die;
		open(OUT,">$Genome.ffn") or die;

		while(<IN>){
			chomp();
			my $line = $_;
			if($_ =~ m/>/){
				my @split = split(/\>/,$line);
				my @split2 = split(" ",$split[1]);
				my $Gene_mark_id = $Genome . "_" . $split2[0];
				my $new = ">gi|$Gene_mark_id|ref|1234| $Gene_mark_id "."[$Genome]";
				print OUT $new . "\n";
			}#end if
			else{print OUT $line . "\n";}
		}#end while IN
		close(IN);
		close(OUT);
		system("rm $Genome" . "*.glimmer.ffn");
	}


}#end Glimmer_ffn

sub AllvsAll
{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	#Get Protein List
	open(AllProtList,"<AllProteinList.all") or die;
	open(OUT,">All.ProteinList") or die;
	while(<AllProtList>){
		chomp();
		if($_ =~ m/>/){print OUT "$_\n";}
	}	
	close(AllProtList); close(OUT);
		open(ALL, "<AllProteinList.all") or die print "Could not open AllProteinList.all";
		system("blastp -db AllProteinList.all -query AllProteinList.all -outfmt '6 qseqid stitle qlen slen length nident' -num_threads 4 > AllProteinList.blastp");
		close(ALL);
}#end sub AllvsAll

sub AllvsAll_Nucleotide
{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	#Get Protein List
	open(AllProtList,"<AllProteinList.all") or die;
	open(OUT,">All.ProteinList") or die;
	while(<AllProtList>){
		chomp();
		if($_ =~ m/>/){print OUT "$_\n";}
	}	
	close(AllProtList); close(OUT);
		open(ALL, "<AllProteinList.all") or die print "Could not open AllProteinList.all";
		system("blastn -db AllProteinList.all -query AllProteinList.all -outfmt '6 qseqid stitle qlen slen length nident' -num_threads 4 > AllProteinList.blastp");
		close(ALL);
}#end sub AllvsAll

sub Filter
{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
		my $blast_file = "AllProteinList.blastp";
		open(IN,"<$blast_file") or die;
		open(OUT,">$blast_file.FILTER") or die;
		my $previous="empty";
		while(<IN>){
			chomp();
			my $line = $_;
			my @tab_split = split(/\t/,$line);
			my $SUBJECT = $tab_split[1];
			my @query_split = split(/\|/,$tab_split[0]);
			my $QUERY = $tab_split[0];
			my $query_gi = $query_split[1];
			my $qlen = $tab_split[2];
			my $slen = $tab_split[3];
			my $len = $tab_split[4];
			my $ident = $tab_split[5]/$len;
			my $Percent_len = $len/$slen;
			my $result = index($SUBJECT,$query_gi); #if this is -1 it is not there
			#ADJUST THE CUT-OFF HERE			
			if($Percent_len >= 0.8 && $ident >= 0.6){
					print OUT $QUERY . "\t" . $SUBJECT . "\n";
			}#end if
		}#end while(<IN>)
	close(IN);
	close(OUT);
}#End sub Filter

sub ModAllProteinDB{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	my $faa = $_[1];
	my $file = $faa;
	$file =~ s{\.[^.]+$}{};
	my $ScoreFilter = "$file.AllProteinList.blastp.FILTER";
	open (ALL,"<AllProteinList.all") or die;
	#Make a hash where the key is the SeqID and the value is the sequence
	my %AllProteins;
	my $key; #this will store the genome title
	my $value;
	my $count=0; #need a count to make sure I don't push nothing onto the hash
	while(<ALL>){
		my $line = $_;
		if($line =~ m/>/){
				if($count > 0){
					#print $key . $value;
					$AllProteins{$key}=$value;
				}#end count if
		$key = $line;
		chomp($key);
		$value=""; #re-initialize value
		}#end line if
		else{
			$value = $value . $line;
			$count=1;
		}#end else
	}#end While ALL
	close(ALL);
=c
	open(OUT,">AllProteinHashPrintOUT") or die;
	while (my ($key, $value) = each %AllProteins) {
				print OUT $key . $value;
			}#end hash print While
	close(OUT);
=cut
	#Now have AllProteins.all stored in %AllProteins hash.
	my @Hits;
	open (ScoreFilter,"<$ScoreFilter") or die;
	my $QUERY; #from the AllProteins list in the ScoreFilter file
	my $HIT; #from a reference genome in the ScoreFilter file
	while(<ScoreFilter>){
	chomp();
		my $line = $_;
		my @split=split(/\t/,$_);
		$QUERY=$split[0];
		my @QUERY_GI = split(/\|/,$QUERY);
		my $QUERY_GI = $QUERY_GI[1];
		$HIT=$split[1];
		my @HIT_GI = split(/\|/,$HIT);
		my $HIT_GI = $HIT_GI[1];
		push(@Hits,"$QUERY_GI\t$HIT_GI\t$HIT");
	}#end While(ScoreFilter)
	close(ScoreFilter);

	@Hits = sort { $a <=> $b } @Hits;
	foreach(@Hits){
		my $element = $_;
		my @split = split(/\t/,$_);
		my $QUERY_GI = $split[0];
		my $HIT_GI = $split[1];
		my $HIT = $split[2];
		my $found = "False";
		foreach(keys %AllProteins){
				if($_ =~ m/$HIT_GI/ && $_ !~ m/\$/){$found = "True";}  
		}#end foreach
		foreach(keys %AllProteins){
			my $line = $_;
			if($_ !~ m/$QUERY_GI/ && $_ =~ m/$HIT_GI/){$found = "True";}
			if($line =~ m/$QUERY_GI/ && $line !~ m/$HIT_GI/ && $found eq "True"){
					my $new_key = $line . "\$" . $HIT;
					$AllProteins{$new_key}=$AllProteins{$line};
					#print "$new_key\n";
					#print "$AllProteins{$new_key}\n";
					delete $AllProteins{$line};
					foreach(keys %AllProteins){
						my $nested = $_;
						if($nested =~ m/$HIT_GI/ && $nested !~ m/\$/){
							delete $AllProteins{$nested};
						}
					}#end 1st nested foreach
					foreach(keys %AllProteins){
						my $nested = $_;
						if($nested =~ m/$QUERY_GI/ && $nested !~ m/\$/){
							delete $AllProteins{$nested};
						}
					}#end 2nd nested foreach
			}#end if
		}#end foreach(AllProteins)
	}#end foreach keys

	system("rm AllProteinList.all");
	open(OUT,">AllProteinList.all") or die;
	while (my ($key, $value) = each %AllProteins){
		print OUT $key . "\n" . $value. "\n";
	}#end hash print While
	close(OUT);

}#end Sub ModAllProteinDB

sub CreateProteinList2
{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);
	#print getcwd . "\n";
	open(FAA,"<AllProteinList.all") or die;
	open(OUT,">AllProteinList.Curated.ProteinList") or die;
	while(<FAA>){
	chomp();
		if($_ =~ m/>/){
			print OUT "$_\n";						
		}
	}	
	close(FAA);
	close(OUT);
}#end CreateProteinList2

sub ParseBlast{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir($home);

	open (ALL,"<AllProteinList.all") or die;
	#Make a hash where the key is the GI SeqID and the value is the "NOT_YET"
	my %AllProteins;
	my $key; #this will store the genome title
	my $value;
	my $count=0; #need a count to make sure I don't push nothing onto the hash
	while(<ALL>){
		my $line = $_;
		if($line =~ m/>/){
				my @split = split(/\|/,$line);
				my $gi = $split[1];
				$AllProteins{$gi}="NOT_YET";
		}#end count if
	}#end While ALL
	close(ALL);

	open(OUT, ">HashTest") or die;
	while (my ($key, $value) = each %AllProteins) {
				print OUT "$key $value\n";
			}#end hash print While
	close(OUT);

	#Creates a hash of all the strain names
	open(STRAINS,"<$_[1]") or die;
	my %STRAINS;
	my $Hash_strain_size = keys %STRAINS;
	while(<STRAINS>){
		chomp();
		$STRAINS{$_}="NOT_DONE";	
	}
	close(STRAINS);


	system("mkdir Pan_Genome");
	open (ALL,"<AllProteinList.blastp.FILTER") or die;
	chdir("$home/Pan_Genome");
	my $previous;
	my $burned_GI;
	while(<ALL>){
		chomp();
		my $line = $_;
		my @Tab_split = split(/\t/,$line);
		my @Query_split = split(/\|/,$Tab_split[0]);
		my $Query_GI = $Query_split[1];
		my @Hit_split= split(/\|/,$Tab_split[1]);
		my $Hit_GI = $Hit_split[1];
		my $Hit_title = $Hit_split[4];
		my @Hit_strain = split(/\[/,$Hit_title);
		my $Hit_strain = $Hit_strain[1]; #This will store the strain of the hit
		chop($Hit_strain);
		#print $Hit_strain . "\n";
		my $Hit_out = ">$Hit_GI\$$Hit_title\n";
		#This first if will only create a new file for the hits if its the first instance of the GI
		if($AllProteins{$Query_GI} ne "DONE" && $previous ne $Query_GI ){ #&& $burned_GI ne $Query_GI
			system("touch $Query_GI");
			$previous = $Query_GI;
			if($Query_GI eq $Hit_GI){
				#$STRAINS{$Hit_strain}="DONE";
				open(OUT, ">>$Query_GI") or die;
				print OUT "$Hit_out";
				$AllProteins{$Query_GI}="DONE";
				close(OUT);
			}#end if that puts self hit into file
		}#end if to create new files
		if($AllProteins{$Hit_GI} ne "DONE" && $STRAINS{$Hit_strain} ne "DONE" ){ #&& $burned_GI ne $Query_GI
				#$STRAINS{$Hit_strain}="DONE";
				open(OUT, ">>$Query_GI") or die;
				print OUT "$Hit_out";
				$AllProteins{$Hit_GI}="DONE";
				close(OUT);
		}


	}#end while

=c		#Once all elements of the hash are found, we want to move onto the next GI		
	
		my $Hash_strains_count=0;
		foreach(keys %STRAINS){
			if($STRAINS{$_} eq "DONE"){
				$Hash_strains_count++;
			}
		}
		
		if($Hash_strains_count == $Hash_strain_size){
			$burned_GI = $Query_GI;
			foreach(keys %STRAINS){
				$STRAINS{$_} = "NOT_YET";
			}
		}


	}#end while
	
chdir($home);
open(OUT, ">HashTest.after") or die;
	while (my ($key, $value) = each %AllProteins) {
				print OUT "$key $value\n";
			}#end hash print While
	close(OUT);
=cut


}#end sub ParseBlast

#This is where the core is found
sub GetStrains{ 
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	chdir("$home");
	open(Strains,"<$_[1]") or die;
	my %STRAINS;
	while(<Strains>){
		chomp();
		$STRAINS{$_}="NOT_FOUND";
	}
	close(Strains);
	$Hash_size = keys %STRAINS;
	foreach(@STRAINS){print $_ . "\n";}
	chdir("$home");
	system("mkdir Core_Genome");
	chdir("$home/Pan_Genome");
	my @Files = glob "*";
	foreach my $file(@Files){
		open(IN,"<$file") or die;
		while (<IN>){
			my $line = $_;
			foreach my $strain(keys %STRAINS){
					if($line =~ m/$strain/){
						$STRAINS{$strain}="FOUND";
					}#end if
			}#end foreach(@STRAINS)
		}#end while(<IN>)		
		close(IN);
	my $count;
	foreach my $key(keys %STRAINS){
		if($STRAINS{$key} eq "FOUND"){
			$count++;
		}#end if
	}	

	if($count == $Hash_size){
		#print "$file\n";
		system("cp $file /$home/Core_Genome");
		
	}

	foreach(keys %STRAINS){
		$STRAINS{$_}="NOT_FOUND";
	}

	}#end foreach(@files)

}#end GetStrains


sub Uniques{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	system("mkdir Accessory");
	open(Strains,"<$_[1]") or die;
	my @STRAINS;
	while(<Strains>){
		chomp();
		push(@STRAINS,$_);
	}
	close(Strains);
	chdir("$home/Accessory");
	foreach(@STRAINS){
		system("mkdir $_");
	}#end foreach strains
	chdir("$home/Pan_Genome");
	my @Files = glob "*";
	foreach(@Files){
		my $file = $_;
			open(IN, "<$file") or die;
			my $line_count=0;
			$line_count++ while (<IN>);
			close(IN);
		if($line_count == 1){
			open(IN, "<$file") or die;
			my $genome;
			while(<IN>){
				chomp();
				$line = $_;
				my $gi;
				my @split = split(/\$/,$line);
				$gi = $split[0];
				$file2 = ">$file";
				foreach(@STRAINS){
					if (index($line,$_) != -1 && $file2 eq $gi){
						$genome = $_;
						system("cp $file $home/Accessory/$genome");
					}#end if
				}#end foreach STRAINS
			}#end while IN

			close(IN);
		}#end if line_count == 1

	}


}#end sub Uniques


sub Elim_Duplicates
{
	my $home = $_[0];
	chdir($home);
	$home = getcwd; #if the user enters './' as a dir this will solve issues
	system("mkdir Unique_Core");
	open(Strains,"<$_[1]") or die;
	my %STRAINS;
	while(<Strains>){
		chomp();
		$STRAINS{$_}="Not_Yet";
	}
	close(Strains);
	#while (my ($key, $value) = each %STRAINS) {
	#			print "$key $value\n";
	#}#end hash print While
	$Hash_size = keys %STRAINS;
	chdir("$home/Core_Genome");
	my @Core_Genes = glob "*";
	foreach(@Core_Genes){
		my $file = $_;

		open(IN, "<$file") or die;
			my $line_count=0;
			$line_count++ while (<IN>);
		close(IN);

		open(IN, "<$file") or die;
		while(<IN>){
		chomp();
		my @Split = split(/\[/,$_);
		my $strain = $Split[1];
		chop($strain);
		my $done_count=0;
		foreach(keys %STRAINS){
			if(index($strain,$_) != -1 && $STRAINS{$_} ne "DONE"){
				#print "$strain $_ \n";
				$STRAINS{$_}="DONE";
			}
			if($STRAINS{$_} eq "DONE"){
				$done_count++;
			}
		}#end foreach(keys %STRAINS)

		if($done_count eq $Hash_size && $line_count eq $Hash_size){
			system("cp $file $home/Unique_Core");
		}

	 }#end while(IN)
	close(IN);
	#Reset the hash		
	foreach(keys %STRAINS){
		$STRAINS{$_}="NOT_YET";
	}
	
	}#end foreach(@Core_Genes)

}


##############################################
################# MAIN #######################
##############################################

#Strain Names should not contain spaces

sub MAIN{

	my $usage="\nUSAGE: HardCORE_v4.0.pl /path/to/*.faa Strain_list_file \n";
	my $dir = $ARGV[0] or die ($usage);
	my $strains = $ARGV[1] or die ($usage);

	ProkkaPrep($dir); #File names must be 08-5578.faa.prokka
	#GeneMark($dir);   #File names must be 08-5578-GeneMark.faa
	#Glimmer($dir); #File names must be 08-5578.glimmer.faa


	#ProkkaPrep_ffn($dir); #File names must be ------.ffn.prokka
	#GeneMark_ffn($dir);   #File names must be 08-5578-GeneMark.ffn
	#Glimmer_ffn($dir); #File names must be 08-5578.glimmer.ffn

	system("cat *.faa > AllProteinList.all"); #Making a global list of proteins
	system("formatdb -i AllProteinList.all");
	print "AT AllvsAll STAGE\n";
	AllvsAll($dir);
	print "AT FILTER STAGE\n";
	Filter($dir);
	print "AT Parse STAGE\n";
	ParseBlast($dir,$strains);
	GetStrains($dir,$strains);
	Uniques($dir,$strains);
	Elim_Duplicates($dir, $strains);

}


MAIN();


=c
 _______   .__          __    
 \      \  |__|  ____  |  | __
 /   |   \ |  |_/ ___\ |  |/ /
/    |    \|  |\  \___ |    < 
\____|__  /|__| \___  >|__|_ \
        \/          \/      \/
__________          __                                   .__   .__           
\______   \  ____ _/  |_ _______   ____    ____    ____  |  |  |  |  _____   
 |     ___/_/ __ \\   __\\_  __ \ /  _ \  /    \ _/ __ \ |  |  |  |  \__  \  
 |    |    \  ___/ |  |   |  | \/(  <_> )|   |  \\  ___/ |  |__|  |__ / __ \_
 |____|     \___  >|__|   |__|    \____/ |___|  / \___  >|____/|____/(____  /
                \/                            \/      \/                  \/ 
=cut




