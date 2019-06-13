#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $version = 0.9;


(my $usage = <<OUT) =~ s/\t+//g;
Somatic variant calling pipeline 
Pipeline version: $version

Usage: perl $0  --rdir --log --t --q 


<rdir> = full path of the folder holding files for this sequence run (user must provide)
<log> = full path of the folder for saving log file; usually upper folder of rdir
<type> = somatic or germline
<q> which queue for submitting job; research-hpc, ding-lab, long (default)
 
OUT

my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m";

my $help = 0;
my $q_name="";
my $run_dir="";
my $log_dir="";
my $type="";

#__PARSE COMMAND LINE
my $status = &GetOptions (
      "rdir=s" => \$run_dir,
      "log=s"  => \$log_dir,
      "q=s" => \$q_name,
      "t=s" => \$type,
      "help" => \$help,
    );


if ($help || $run_dir eq "" || $log_dir eq "" || $type eq "") 
{
      print $usage;
      exit;
  }

print "run dir=",$run_dir,"\n";
print "log dir=",$log_dir,"\n";
print "type=",$type,"\n";
print "queue name=",$q_name,"\n";

my $working_name= (split(/\//,$run_dir))[-1];
my $HOME1=$log_dir;
#store job files here
if (! -d $HOME1)
{
`mkdir $HOME1`;
}
if (! -d $HOME1."/tmpcleaner") {
    `mkdir $HOME1"/tmpcleaner"`;
}
my $job_files_dir = $HOME1."/tmpcleaner";
#store SGE output and error files here
if (! -d $HOME1."/LSF_DIR_CLEANER") {
    `mkdir $HOME1"/LSF_DIR_CLEANER"`;
}

my $lsf_file_dir = $HOME1."/LSF_DIR_CLEANER";
my $run_script_path =`echo \$PWD`;
chomp $run_script_path;
my $script_dir=$run_script_path;
print $script_dir,"\n";

$run_script_path = "/gsc/bin/perl ".$run_script_path."/";

print $run_script_path,"\n";
my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";

opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
        $sample_name = $sample_dir_list[$i];
        if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
            $sample_full_path = $run_dir."/".$sample_name;
            if (-d $sample_full_path) { # is a full path directory containing a sample
                print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
            #   sleep 1;
                $current_job_file="";
			## gzip pindel big files ##
                   &bsub_clean();
			}
}

}

sub bsub_clean{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j1_clean_".$sample_name.".sh";
    
	my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    #`rm $current_job_file`;

    open(CLEAN, ">$job_files_dir/$current_job_file") or die $!;
    print CLEAN "#!/bin/bash\n";
    print CLEAN "RP=".$sample_full_path."/pindel/".$sample_name."_RP\n";
	print CLEAN "D=".$sample_full_path."/pindel/".$sample_name."_D\n";
	print CLEAN "TD=".$sample_full_path."/pindel/".$sample_name."_TD\n";
    print CLEAN "INV=".$sample_full_path."/pindel/".$sample_name."_INV\n";
    print CLEAN "SI=".$sample_full_path."/pindel/".$sample_name."_SI\n";	
	print CLEAN "RAW=".$sample_full_path."/pindel/pindel.out.raw\n";  	
	print CLEAN "FAIL=".$sample_full_path."/pindel/pindel.out.raw.CvgVafStrand_fail\n";
    print CLEAN "if [ -f \${RP} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${RP}\n";
    print CLEAN "fi\n";

    print CLEAN "if [ -f \${D} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${D}\n";
    print CLEAN "fi\n";
    
    print CLEAN "if [ -f \${TD} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${TD}\n";
    print CLEAN "fi\n";
    

    print CLEAN "if [ -f \${INV} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${INV}\n";
    print CLEAN "fi\n";
    

    print CLEAN "if [ -f \${SI} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${SI}\n";
    print CLEAN "fi\n";
    
    print CLEAN "if [ -f \${RAW} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${RAW}\n";
    print CLEAN "fi\n";

    print CLEAN "if [ -f \${FAIL} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${FAIL}\n";
    print CLEAN "fi\n";


    close CLEAN;

    my $sh_file=$job_files_dir."/".$current_job_file;

    if($q_name eq "research-hpc")
    {
    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -o $lsf_out -e $lsf_err bash $sh_file\n";     }
    else {
        $bsub_com = "bsub -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -o $lsf_out -e $lsf_err bash $sh_file\n";
    }

    print $bsub_com;
    system ($bsub_com);

}
