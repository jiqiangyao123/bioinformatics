#!/usr/bin/perl -w

use strict;
use Benchmark;
use Time::HiRes qw(gettimeofday tv_interval);
use Scalar::Util qw(looks_like_number);

###################################################################################################
my ($time0,$time1);start_time();
###################################################################################################

my ($line,$lncp,$big,$bigline,$newline,$size,$siz,$temp,$tmp,$sample,$samp,$tran,$flag);
my ($ave,$tot,$total,$fnm,$fnam,$fname,$prob,$gen,$gene,$len,$min,$max,$count,$counter);
my ($fr,$fw,$infile,$otfile,$inf,$otf,$dir,$otdir,$file);
my ($id,$type,$rscript,$cmd,$cmd1);
my ($chr,$chrm,$ch,$chrom,$head,$header,$pos,$cod,$cood,$coord,$cov,$cover,$coverage,$val);
my ($per,$perc,$percentage,$pval,$num,$star,$stop,$pre,$cur,$com,$dif,$diff,$stat,$name,$job);
my (@temp,@tmp,@samp,@sample,@samples,@gene,@genes,@id,@ids,@com,@diff,@dif,@lines,@files);
my (%norepeat,%temp,%tmp,%tot,%loc);
my ($c,$d,$e,$f,$g,$h,$i,$j,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
my ($A,$B,$C,$D,$E,$F,$G,$H,$I,$J,$K,$L,$M,$N,$O,$P,$Q,$R,$S,$T,$U,$V,$W,$X,$Y,$Z);
my (@a,@b,@c,@d,@e,@f,@g,@h,@i,@j,@k,@l,@m,@n,@o,@p,@q,@r,@s,@t,@u,@v,@w,@x,@y,@z);
my (@A,@B,@C,@D,@E,@F,@G,@H,@I,@J,@K,@L,@M,@N,@O,@P,@Q,@R,@S,@T,@U,@V,@W,@X,@Y,@Z);
my ($MIN,$MAX);
my $PRJDIR = "/share/dept_bbsr/Projects/Ahmed_Kamran/4213_19662_RNAseq_2023/2nd_analysis";
#### step1 ########################################################################################
  
     step1();

print_time();

###################################################################################################

sub step1
{
        my ($dir,$datdir,$wrkdir,$outdir,$inf,$otf,$big,@lines,@files,$size,@temp,$siz,%genum,@gen);
        my (%dup,%geninfo,@dup,%dup1,%avexp,$ensg,%done,@head,%rankensg,%rankgen,@group,$grp,$tmpf);
        $tmpf  = "temp1.csv";
        @group = (1,2,3,4);
        foreach $grp (@group)
        {
        $big      = "";
        %avexp    = %rankensg = %rankgen = %genum = %dup = %dup1 = %done = ();
        $inf      = "$PRJDIR/work/differential_expression/all_samples_RSEM_gene_count_Group$grp\_DE_allgenes.csv";
        %rankensg = get_rank_ensg($inf);ph('%rankensg',%rankensg);
	%avexp    = get_average_gene_expression($inf);ph1('%avexp',%avexp);
        @lines    = get_lines($inf);
        $size     = @lines;
        for($i=1;$i<$size;$i++)
	{
           $lines[$i] =~ s/\"//g;
	   @temp = split /\,/, $lines[$i];$siz  = @temp;
           $y    = 0;
           for($j=11;$j<$siz;$j++)
           {
               $x = $temp[$j];if($x==0){$y++;}
           }   if($y >= $siz/2){next;}
           $genum{$temp[3]}++;$geninfo{$temp[3]} .= "$temp[1],";
	}  @gen = sort keys %genum;
        foreach $g (@gen){if($genum{$g}>1){print "dup_gen: $g\n";$dup{$g}=$geninfo{$g};}}
        @dup = sort keys %dup;
        foreach $g (@dup)
        {
           print "$g\t$dup{$g}\n";
           $temp = $dup{$g};
           $temp =~ s/,$//;
           $ensg = fine_best_gene($temp,$inf,\%avexp);
           $dup1{$g} = $ensg;
        }  ph('%dup1',%dup1);

        my ($tot,$smp,$body,$symb,$head,@head1,$head1,$line1);
        $big  = "";
        $head = $lines[0];
        $head =~ s/\"//g;
        for($i=1;$i<$size;$i++)
	{
           $lines[$i] =~ s/\"//g;#print "line = $lines[$i]\n";
	   @temp = split /\,/, $lines[$i];$siz  = @temp;
           $y    = 0;
           for($j=11;$j<$siz;$j++)
           {
               $x = $temp[$j];if($x==0){$y++;}
           }   if($y >= $siz/2){next;}
           $ensg = $temp[1];
           $symb = $temp[10];$rankgen{$symb}= $rankensg{$ensg};
           if($dup1{$symb})
           {
              if($done{$symb}){next;}
              if($ensg eq $dup1{$symb})
              {
                $big .= "$lines[$i]\n";
                $done{$symb}=1;$rankgen{$symb}= $rankensg{$ensg};
              }
           }
           else{$big .= "$lines[$i]\n";$done{$symb}=1;}
	}  $big = "$big";
        $otf   = $tmpf;write_to_file($otf,$big);
        @lines = get_lines($otf);
        $size  = @lines;
        @head  = split /\,/, $head;
        $temp  = @head;
        $smp   = 0; foreach $h (@head){if($h =~ /S\d+/){$smp++;}}
        @head1 = ();foreach $h (@head){if($h =~ /S\d+/){push(@head1,$h);}}
        $head1 = join("\t",@head1);
        $line1 = "NAME\tDESCRIPTION\t$head1\n";

        $otf   = "Group$grp.gct";
        $head  = "#1.2\n$size\t$smp\n";
        $big   = "$head$line1";
        foreach $l (@lines)
        {
             $l =~ s/\"//g;
             for($j=0;$j<10;$j++){$l =~ s/^.*?,//;}
             $l =~ s/,/\t/g;
             $l =~ s/^(.*?)\t/$1\tNA\t/;
             $big .= "$l\n";
        }    write_to_file($otf,$big);

        
        $otf   = "Group$grp.rnk";
        @tmp   = sort keys %done;py('@tmp',@tmp);
        $big   = "";
        foreach $t (@tmp){$big .= "$t\t$rankgen{$t}\n";}
        write_to_file($otf,$big);
        }system("rm -rf $tmpf");
}

sub get_rank_ensg
{
     my $inf   = shift;#print "$inf\n\n";
     my (@lines,@tmp,$gen,$pva,$lgf,$val,%tmp);
        @lines = get_lines($inf);
	foreach my $l (@lines)
	{
           $l   =~ s/\"//g;if($l =~ /GencodeID/ || $l =~ /GeneName/){next;}
	   @tmp = split /\,/, $l;#print "l = $l\n";
           $gen = $tmp[1];
           $pva = $tmp[7];
           $lgf = $tmp[4];
           $val = (-1)*log10($pva)*$lgf;
           $tmp{$gen} = round($val,6);
	}  return %tmp;
}

sub fine_best_gene
{
    my ($genlst,$inf,$avexpref) = @_;
    my (@genes,$max,$find);
       @genes = split /,/, $genlst;
       $max   = 0;
       foreach $g (@genes)
       {
           if($avexpref->{$g} > $max){$max = $avexpref->{$g};$find=$g;}
       }   return $find;
}

sub get_average_gene_expression
{
    my $inf   = shift;print "inf = $inf\n";
    my (@lines,$size,%avexp,$i,$j,$l,@l,$x,$y,$s,@t);
       @lines = get_lines($inf);
       $size  = @lines;
       for($i=1;$i<$size;$i++)
       {
           $l = $lines[$i];
           $l =~ s/"//g;
           @l = split /,/, $l;
           $s = @l;
           @t = ();
           for($j=11;$j<$s;$j++){push(@t,$l[$j]);}
           $x = average(@t);
           $y = $l[1];
           $avexp{$y} = $x;
        }  return %avexp;
}
####################################################################################
########################  regular ##################################################
####################################################################################

sub start_time
{
   	$time0 = new Benchmark;    
}


sub print_time
{
	$time1 = new Benchmark;
        my $diff= timestr(timediff($time1,$time0));
        print "Time used: $diff\n";    
}

sub get_cur_date
{
    my ($hr,$min,$sec,$date,$mon,$yr,$time);
    my $curtime = localtime;
     $curtime =~ /^\w+\s+(\w+)\s+(\d+)\s+(.*?)\s+(\d+)/;
     $yr = $4; $mon = uc $1; $date = $2; $time = $3;
     
       if($mon eq "JAN"){$mon = "01";}
    elsif($mon eq "FEB"){$mon = "02";}
    elsif($mon eq "MAR"){$mon = "03";}
    elsif($mon eq "APR"){$mon = "04";}
    elsif($mon eq "MAY"){$mon = "05";}
    elsif($mon eq "JUN"){$mon = "06";}
    elsif($mon eq "JUL"){$mon = "07";}
    elsif($mon eq "AUG"){$mon = "08";}
    elsif($mon eq "SEP"){$mon = "09";}
    elsif($mon eq "OCT"){$mon = "10";}
    elsif($mon eq "NOV"){$mon = "11";}
    elsif($mon eq "DEC"){$mon = "12";}

    if($date < 10){$date = "0$date";}
    $time =~ /(\d+):(\d+):(\d+)/;
    $hr = $1; $min = $2; $sec = $3;
    $curtime = "$yr$mon$date";
    return $curtime;
}

sub txt2csv
{
        my ($inf,$otf) = @_;
        my $cmd = "cat $inf | sed -E 's/\t+\$//g' | tr '\t' ',' > $otf";runcmd($cmd,1);
}

sub openw
{
    my $file = shift;

    my $FW;
    open($FW, ">$file") || die "Couldn't open $file\n";
    return $FW;
}


sub openr
{
    my $file = shift;

    if(-e $file){}else{die "file not exist: $file\n";}

    my $FR;
    open($FR, "<$file") || die "Couldn't open $file\n";
    return $FR;
}

sub get_lines
{
        my $inf      = shift;
        my $fr       = openr($inf);
        my @temp     = ();
           while(!eof($fr))
           {
                 $line = <$fr>;chomp $line;if(!$line){next;}
                 push(@temp, $line);
           }close($fr);
           return @temp;
}

sub print_counter
{
    my ($nam,$c,$num) = @_;
  
       if(!$c){return;}
  
    my $y = 1;
    my $u = "M";

    if(!$num)
    {
       $y = 1000000;
    }
    else
    {
       $y = $num;
    }

    if($y == 1000000){$u = "M";}
    if($y == 10000  ){$u = "W";}
    if($y == 1000   ){$u = "K";}
    if($y == 100    ){$u = "H";}
    if($y == 10     ){$u = "T";}
    if($y == 1      ){$u = "";}

    if($c%$y == 0){my $v=$c/$y; print "$nam $v $u\n";}
}

sub create_directory
{
    my $dir = $_[0];

    if(-e $dir){return;}
    else
    {
       system("mkdir -p $dir");
       print "created: $dir\n";
    }
}

sub create_directory1
{
    my $dir = $_[0];

    if(-e $dir){system("rm -rf $dir");}
   
    system("mkdir -p $dir");
    print "created: $dir\n";

}

sub get_directory_name
{
     my ($dir,$pat)  = @_;if(!$pat){$pat = '\w';}
     my @temp = ();
       opendir(DIR,$dir) || die "Error in opening dir $dir\n";
     while(my $file = readdir(DIR))
     {
        if($file =~ /$pat/)
        {
          if(-d "$dir/$file"){push(@temp,$file);}
        }
     }close(DIR);
     @temp = sort @temp;
     return @temp;
}

sub read_directory
{
     my ($dir,$pat)  = @_;
     my @temp = ();
       opendir(DIR,$dir) || die "Error in opening dir $dir\n";
     while(my $file = readdir(DIR))
     {
        
        if($file =~ /$pat/i)
        {
          if(-d "$dir/$file"){} 
          else{push(@temp,$file);}
        }
     }close(DIR);
     @temp = sort @temp;
     return @temp;
}

sub extractfilename
{
    my $fullname = shift;
       $fullname =~ s/\/$//;
       $fullname =~ s/.+\///;
       return $fullname;
}

sub extractdirname
{
    my $dirname = shift;
       $dirname =~ s/(.+)\/.+/$1/;
       return $dirname;
}

sub commify 
{
    my $text = reverse $_[0];
       $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
       return scalar reverse $text;
}

sub print_scalar
{
    my ($fname,$value) = @_;
       print "$fname = $value\n";
}

sub print_scalar1
{
    my ($fname,$value) = @_;
       print_scalar($fname,$value);
       exit;
}

sub print_array
{
    my ($fname,@arr) = @_;
    my $size    = @arr;
    my $counter = 0;
   
       foreach my $a (@arr)
       {
          $counter++;
          print "$fname $counter $a\n";
       }
}

sub print_array1
{
    my ($fname,@arr) = @_;
       print_array($fname,@arr);
       exit;
}

sub print_hash
{
    my ($fname,%hash) = @_;
    my @arr  = sort keys %hash;

    my $counter = 0;
   
       foreach my $a (@arr)
       {
          $counter++;
          print "$fname $counter $a -> $hash{$a}\n";
       }
}

sub print_hash1
{
    my ($fname,%hash) = @_;
    
       print_hash($fname,%hash);
       exit;
}

sub print_hash_2D
{
    my ($fname,%ahash) = @_;

    my @arr1 = sort keys %ahash;#print_array('@arr1',@arr1);

    my $counter = 0;
   
       foreach my $a1 (@arr1)
       {
          my @arr2 = sort keys %{$ahash{$a1}};
       foreach my $a2 (@arr2)
       {
          $counter++;
          if($ahash{$a1}{$a2})
          {
          print "$fname $counter $a1 $a2 -> $ahash{$a1}{$a2}\n";
          }
          else
          {
          print "$fname $counter $a1 $a2 -> na\n";
          }
       }
       }
}

sub print_hash_2D1
{
    my ($fname,%ahash) = @_;

       print_hash_2D($fname,%ahash);
       exit;
}

sub print_hash_3D
{
    my ($fname,%ahash) = @_;

    my @arr1 = sort keys %ahash;#print_array('@arr1',@arr1);

    my $counter = 0;
   
       foreach my $a1 (@arr1)
       {
          my @arr2 = sort keys %{$ahash{$a1}};
       foreach my $a2 (@arr2)
       {
          my @arr3 = sort keys %{$ahash{$a1}{$a2}};
       foreach my $a3 (@arr3)
       {
          $counter++;
          print "$fname $counter $a1 $a2 $a3-> $ahash{$a1}{$a2}{$a3}\n";
       }
       }
       }
}

sub print_hash_3D1
{
    my ($fname,%ahash) = @_;

       print_hash_3D($fname,%ahash);
       exit;
}

sub ps
{
    my ($fname,$value) = @_;
        print_scalar($fname,$value);
}

sub ps1
{
    my ($fname,$value) = @_;
        print_scalar1($fname,$value);
}

sub py
{
    my ($fname,@arr) = @_;
        print_array($fname,@arr);
}

sub py1
{
    my ($fname,@arr) = @_;
        print_array1($fname,@arr);
}

sub ph
{
    my ($fname,%ahash) = @_;
        print_hash($fname,%ahash);
}

sub ph1
{
    my ($fname,%ahash) = @_;
        print_hash1($fname,%ahash);
}

sub ph2d
{
    my ($fname,%ahash) = @_;
        print_hash_2D($fname,%ahash);
}

sub ph2d1
{
    my ($fname,%ahash) = @_;
        print_hash_2D1($fname,%ahash);
}

sub ph3d
{
    my ($fname,%ahash) = @_;
        print_hash_3D($fname,%ahash);
}

sub ph3d1
{
    my ($fname,%ahash) = @_;
        print_hash_3D1($fname,%ahash);
}

sub write_array
{
    my ($fname,@arr) = @_;
    my $fw = openw($fname);
    my $c  = 0;
       foreach my $a (@arr){$c++;print $fw "$c\t$a\n";}
       close($fw);
}

sub write_array1
{
    my ($fname,@arr) = @_;
    my $fw = openw($fname);
    my $c  = 0;
       foreach my $a (@arr){$c++;print $fw "$c\t$a\n";}
       close($fw);

       exit;
}

sub write_hash
{
    my ($fname,%hash) = @_;
    my @arr  = sort keys %hash;

    my $counter = 0;
   
    my $fw = openw($fname);
       foreach my $a (@arr)
       {
          $counter++;
          print $fw "$fname $counter $a -> $hash{$a}\n";
       }  close($fw);
}

sub write_hash1
{
    my ($fname,%hash) = @_;
    my @arr  = sort keys %hash;

    my $counter = 0;
   
    my $fw = openw($fname);
       foreach my $a (@arr)
       {
          $counter++;
          print $fw "$fname $counter $a -> $hash{$a}\n";
       }  close($fw);

       exit;
}

sub write_hash_2D
{
    my ($fname,%ahash) = @_;

    my $fw   = openw($fname);
    my @arr1 = sort keys %ahash;#print_array('@arr1',@arr1);

    my $counter = 0;
   
       foreach my $a1 (@arr1)
       {
          my @arr2 = sort keys %{$ahash{$a1}};
          my $big  = "";
       foreach my $a2 (@arr2)
       {
           $counter++;
           $big .= "$fname $counter $a1 $a2 -> $ahash{$a1}{$a2}\n";
       }   print $fw $big;
       }close($fw);
}

sub write_hash_2D1
{
    my ($fname,%ahash) = @_;

    write_hash_2D($fname,%ahash);

    exit;
}

sub write_hash_3D
{
    my ($fname,%ahash) = @_;

    my $fw   = openw($fname);
    my @arr1 = sort keys %ahash;#print_array('@arr1',@arr1);

    my $counter = 0;
   
       foreach my $a1 (@arr1)
       {
          my @arr2 = sort keys %{$ahash{$a1}};
       foreach my $a2 (@arr2)
       {
          my @arr3 = sort keys %{$ahash{$a1}{$a2}};
          my $big  = "";
       foreach my $a3 (@arr3)
       {
          $counter++;
          $big .= "$fname $counter $a1 $a2 $a3-> $ahash{$a1}{$a2}{$a3}\n";
       }  print $fw $big;
       }
       }close($fw);
}

sub write_hash_3D1
{
    my ($fname,%ahash) = @_;

    write_hash_3D($fname,%ahash);

    exit;
}

sub wy
{
    my ($fname,@arr) = @_;
        write_array($fname,@arr);
}

sub wy1
{
    my ($fname,@arr) = @_;
        write_array1($fname,@arr);
}

sub wh
{
    my ($fname,%ahash) = @_;
        write_hash($fname,%ahash);
}

sub wh1
{
    my ($fname,%ahash) = @_;
        write_hash1($fname,%ahash);
}

sub wh2d
{
    my ($fname,%ahash) = @_;
        write_hash_2D($fname,%ahash);
}

sub wh2d1
{
    my ($fname,%ahash) = @_;
        write_hash_2D1($fname,%ahash);
}

sub wh3d
{
    my ($fname,%ahash) = @_;
        write_hash_3D($fname,%ahash);
}

sub wh3d1
{
    my ($fname,%ahash) = @_;
        write_hash_3D1($fname,%ahash);
}

sub round
{
    my ($v, $n) = @_;
        $n = "%.$n"."f";

    my $x = sprintf($n,$v);
    return $x;
}

sub add_to_hash
{
      my ($e,$n,$href) = @_;
         if($href->{$n}){$href->{$n} .= ",$e";}
         else{$href->{$n} = $e;}
}

sub printpercentage
{
    my ($v,$tot) = @_;

    my  $perc = ($v*100)/$tot;
        $perc = sprintf("%.2f", $perc);
        $perc = "($perc\%)";
        return $perc;
}

sub printpercentage1
{
    my ($v,$tot) = @_;

    my  $perc = ($v*100)/$tot;
        $perc = sprintf("%.2f", $perc);
        return $perc;
}

sub log2 
{
    my $n = shift;
    return log($n)/log(2);
}

sub log10 
{
    my $n = shift;
    return log($n)/log(10);
}

sub get_log_value
{
    my ($v, $b) = @_;

    if($v !~ /^[0-9\.]+$/){die "worng format 1 in get_log_value: $v\n";}
    if($b !~ /^\d++$/    ){die "worng format 2 in get_log_value: $b\n";}

    if(!$v){return 0;}

    return log($v)/log($b);
    
}

sub sci2float
{
    my $x = shift;if(!$x){print "error_sci2float: no value for \$x\n";exit;}

       if($x !~ /[Ee]/){return $x;}

       if($x =~ /^([Ee]\-)(\d+)/){$x = "1$x";}

       if($x =~ /^(\d+[Ee]\-)(\d+)/)
       {
          my $n = $2;

          if($n >100){$n = 100;$x = "$1"."100";}
       }

    my $y = sprintf("%.100f", $x);

       return $y
}

sub countfilelinenum
{
      my $inf     = shift;
      my $fr      = openr($inf);
      my $count   = 0;
         while(!eof($fr))
         {
            my $line = <$fr>;$count++; 
         }close($fr);

         return $count;
}

sub dec2bin 
{
    my $bits = shift;
    my $str = unpack("B8", pack("c", $bits));
    return $str;
}

sub bin2dec 
{
    return unpack("c", pack("B8", substr("0" x 8 . shift , -8)));
}


sub pdf2jpg
{
       my ($src, $des) = @_;

       my $ps = $des; $ps =~ s/\.jpg/\.ps/;

       my $cmd = "pdf2ps $src $ps";print "cmd =1= $cmd\n";system($cmd);
          $cmd = "gs -sDEVICE=jpeg -r300 -sPAPERSIZE=a4 -dBATCH -dNOPAUSE -sOutputFile=$des $ps";print "cmd =2= $cmd\n";system($cmd);
          $cmd = "rm -rf $ps";system($cmd);
}

sub unique
{
       my @data   = @_;
       my @unique = do { my %seen; grep { !$seen{$_}++ } @data };
          return @unique;
}


sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub average
{
       my @data = @_;
       my $size = @data;

       my $sum  = 0;
          foreach my $d (@data) 
          {
             $sum += $d;
          } 
       my $ave = $sum/$size;
          return $ave;
}

sub stdev
{
       my @data = @_;
       my $size = @data;
       my $mean = average(@data);

       my $sqsum = 0;
          foreach my $d (@data) 
          {
          $x     = $d - $mean;
          $sqsum += $x*$x;
          } 
          $sqsum /= ($size - 1);
       my $stdev = sqrt($sqsum);
          return $stdev;
}

sub zscore
{
       my ($x,@data) = @_;

       my $u   = average(@data);#print "u   = $u\n";
       my $sd  = stdev(@data);   #print "sd  = $sd\n";
       my $zsc = ($x - $u)/$sd; #print "zscore = $zsc\n\n";
        return $zsc;    
}

sub rmnewline #remove newline character
{
       my $line = shift;

          $line =~ s/[:a;N;$!ba;s]+/\n/g;

          return $line;
}

sub write_to_file
{
    my ($otf,$bigline) = @_;
    my $fw = openw($otf);
       if($bigline){print $fw $bigline;}
       close($fw);
}

sub clear_file
{
     my $file   = shift;
     my $cmd    = "rm -rf $file";

        if(-e $file){system($cmd);print "$file is cleared\n";}
}

sub excel_genename_correct
{
    my $gene = shift;
       $gene =~ s/(\d+)-(Mar)/MARCH$1/;
       $gene =~ s/(\d+)-(Sep)/SEPT$1/;
       $gene =~ s/(\d+)-(Dec)/DEC$1/;
       $gene =~ s/(\d+)-(Feb)/FEB$1/;
       $gene =~ s/1-Apr/MAGEH1/;
       $gene =~ s/2-Apr/FAM215A/;
       $gene =~ s/3-Apr/ATRAID/;
       $gene =~ s/1-Nov/C11orf40/;
       $gene =~ s/2-Nov/CTGF/;
       $gene =~ s/(\d+)-(Oct)/Oct-$1/;
       return $gene;
}

sub runcmd
{
     my ($cmd,$mode) = @_;print "\ncmd = $cmd\n";

     if($mode){system($cmd);}
}

sub addleadingzero
{
        my ($n,$m) = @_; 

        my $x = $n;
        my $y = $n;
        my $z = int($y/$m);
        while($z == 0)
        {
              $x  = "0$x";
              $y *= 10;
              $z  = int($y/$m);
        }
        return $x;
}

sub add_value_to_hash
{
    my ($key,$val,$hashref) = @_;
  
       if($hashref->{$key})
       {
          if($hashref->{$key} !~ /$val/)
          {
             $hashref->{$key} .= ",$val";
          }
       }
       else{$hashref->{$key} = $val;}
}

sub check_file_exist
{
     my $file  = shift;
     
     if(-e $file){}
     else
     {
         die "\nnot exist: $file\n\n";
     }
}

sub compare_list_size
{
     my ($listref1,$listref2) = @_;
 
     my @list1 = @{$listref1};
     my @list2 = @{$listref2};

     my (%common,%uniq1,%uniq2);

       foreach my $x (@list1)
       {
       foreach my $y (@list2)
       {
          if($x eq $y)
          {
             $common{$x} = 1;
          }
       }
       }

       foreach my $x (@list1)
       {
          if(!$common{$x})
          {
             $uniq1{$x} = 1;
          }
       }

       foreach my $y (@list2)
       {
          if(!$common{$y})
          {
             $uniq2{$y} = 1;
          }
       }


     my @common  = sort keys %common;
     my @unq1    = sort keys %uniq1;
     my @unq2    = sort keys %uniq2;
     my $comsize = @common;
     my $unqsiz1 = @unq1;
     my $unqsiz2 = @unq2;

        return ($comsize,$unqsiz1,$unqsiz2);
}


sub compare_list_file
{
     my ($listref1,$listref2,$otf) = @_;
 
     my @list1 = @{$listref1};
     my @list2 = @{$listref2};

     my (%common,%uniq1,%uniq2);

       foreach my $x (@list1)
       {
       foreach my $y (@list2)
       {
          if($x eq $y)
          {
             $common{$x} = 1;
          }
       }
       }

       foreach my $x (@list1)
       {
          if(!$common{$x})
          {
             $uniq1{$x} = 1;
          }
       }

       foreach my $y (@list2)
       {
          if(!$common{$y})
          {
             $uniq2{$y} = 1;
          }
       }


     my @common  = sort keys %common;
     my @uniq1   = sort keys %uniq1;
     my @uniq2   = sort keys %uniq2;
     my $comsize = @common;
     my $unqsiz1 = @uniq1;
     my $unqsiz2 = @uniq2;

     my $comon   = join("\n",@common);
     my $uniq1   = join("\n",@uniq1);
     my $uniq2   = join("\n",@uniq2);

     my $otfile  = "compare_result.txt";
     my $big     = "common\tuniq1\tuniq2\n";

     my @tmp     = ();push(@tmp,$comsize);push(@tmp,$unqsiz1);push(@tmp,$unqsiz2);
     my $max     = max($comsize,$unqsiz1,$unqsiz2);

        for(my $i=0;$i<$max;$i++)
        {
            $x = ".";if($common[$i]){$x = $common[$i];}
            $y = ".";if($uniq1[$i] ){$y = $uniq1[$i];}
            $z = ".";if($uniq2[$i] ){$z = $uniq2[$i];}
            $line = "$x\t$y\t$z\n";
            $big .= $line;
        }
        write_to_file($otf,$big);
}

sub read_fasta_file
{
       my $inf  = shift;
     
       my ($id,$seq,$line);
       my %seq  = ();
       my $fr   = openr($inf);
       my $stop = 0;
          while(!eof($fr))
          {
            if($stop==0){$line = <$fr>;chomp $line;}
            if($line =~ />(.+)/)
            {
               $stop = 0;
               $id   = $1;
               $line = <$fr>;chomp $line;
               $seq  = "";
               while(!eof($fr) && $line !~ /^>/)
               {
                     $seq  .= $line;
                     $line  = <$fr>;chomp $line;
               }
               if(eof($fr) && $line){$seq  .= $line;}
               $stop=1;
               $seq{$id} = $seq;
            }
          }close($fr);#ph1('%seq',%seq);
          return %seq;
}

sub get_seq_from_coordinate
{
    my $coordinate = shift;#ps('$coordinate',$coordinate);

       if($coordinate !~ /chr\w_\d+_\d+/ && $coordinate !~ /chr\w\s+\d+\s+\d+/){die "not good format: $coordinate\n";}
      
    my @cord = split /_/, $coordinate;
    my $size = @cord;
       if($size == 0){@cord = split /\s+/, $coordinate;$size = @cord;}

    my $chr  = $cord[0];
    my $sta  = $cord[$size - 2] - 1;
    my $stp  = $cord[$size - 1];
    my $len  = $stp - $sta + 2;#print "$coordinate  $len\n";

    my $otf  = "$coordinate.txt";
    my $cmd  = "twoBitToFa /home/yaoj/projects/public_data/grch38/hg38.2bit $otf -seq=$chr -start=$sta -end=$stp";
       runcmd($cmd,1);
    my @lines = get_lines($otf);
    my $tmp   = "";
       foreach my $l (@lines)
       {
            if($l =~ /^\>/){next;}
            $tmp .= $l;
       }
       system("rm -rf $otf");#print "$tmp\n";
       return $tmp;
}

sub search_files1
{
    my $path  = $_[0];
    my $patn  = $_[1];
    my @files = search_files($path);
    my @temp  = ();
       foreach my $f (@files)
       {
          if($f =~ /$patn/){push(@temp,$f);}
       }
       return @temp;
}

sub search_files
{
    my $path  = shift;
    my $files = search_dir($path);
    my @files = split /\n/, $files;
       return @files;
}

sub search_dir 
{
    my $path  = shift;
    my $files = shift;if(!$files){$files = "\n";}
    my @dir_entries = glob("$path/*");
    foreach my $entry (@dir_entries) 
    {

        if(-f $entry)
        {
           $files = "$files\n$entry";
        }
        $files = search_dir($entry,$files) if -d $entry;
    }
    return $files;
}

sub get_chromosomes
{
    my @chrs = ();
    for($i=1;$i<23;$i++)
    {
        push(@chrs,$i);
    }
    push(@chrs,'X');
    push(@chrs,'Y');
    return @chrs;
}

sub sort_chromosome
{
        my @chroms = unique(@_);

        my @temp = ();
        my @chrm = ();
 
           foreach my $ch (@chroms)
           {
               $ch =~ s/^chr//i;
               $ch =~ s/^ch//i;

               if($ch =~ /^\d+$/)
               {
                  push(@temp,$ch);
               }
           }
           foreach my $ch (@chroms)
           {
               $ch =~ s/^chr//i;
               $ch =~ s/^ch//i;

               if($ch =~ /^[XYxy]/)
               {
                  push(@temp,$ch);
               }
           }                           
           foreach my $ch (@temp)
           {
               $ch =~ s/^chr//i;
               $ch =~ s/^ch//i;

               if($ch =~ /^[Mm]/)
               {
                  push(@temp,$ch);
               }
           }
           @chrm = unique(@temp);

           @temp = ();
           for($i=1;$i<=22;$i++)
           {
           foreach $j (@chrm)
           {
              if($i eq $j){push(@temp,$j);}
           }}
           foreach $j (@chrm)
           {
              if($j =~ /[Xx]/){push(@temp,$j);}
           }
           foreach $j (@chrm)
           {
              if($j =~ /[Yy]/){push(@temp,$j);}
           }
           
           return @temp;
}

sub dna2aa
{
        my $dna = shift;
        my @dna = split //, $dna;
        my $siz = @dna;
        my $quo = $siz%3;if($quo > 0){die "wrong seq size\n";}

        my %getcod = get_genetic_code();
        my ($aa,$nuc);
           while($dna =~ /^(\w\w\w)(.+)/)
           {
                 $nuc  = $1;
                 $dna  = $2;
                 $aa  .= $getcod{$nuc};
           }
           return $aa;
}

sub get_genetic_code
{
	my %tmp = ("TGA" => "*",
                   "TAA" => "*",
                   "TAG" => "*",
                   "TTT" => "F",
                   "TTC" => "F",
                   "TTG" => "L",
                   "TTA" => "L",
                   "CTT" => "L",
                   "CTC" => "L",
                   "CTA" => "L",
                   "CTG" => "L",
                   "ATT" => "I",
                   "ATC" => "I",
                   "ATA" => "I",
                   "ATG" => "M",
                   "GTT" => "V",
                   "GTG" => "V",
                   "GTC" => "V",
                   "GTA" => "V",
                   "TCT" => "S",
                   "TCC" => "S",
                   "TCA" => "S",
                   "TCG" => "S",
                   "CCT" => "P",
                   "CCC" => "P",
                   "CCA" => "P",
                   "CCG" => "P",
                   "ACT" => "T",
                   "ACC" => "T",
                   "ACA" => "T",
                   "ACG" => "T",
                   "GCT" => "A",
                   "GCC" => "A",
                   "GCA" => "A",
                   "GCG" => "A",
                   "TAT" => "Y",
                   "TAC" => "Y",
                   "CAT" => "H",
                   "CAC" => "H",
                   "CAA" => "Q",
                   "CAG" => "Q",
                   "AAT" => "N",
                   "AAC" => "N",
                   "AAA" => "K",
                   "AAG" => "K",
                   "GAT" => "D",
                   "GAC" => "D",
                   "GAA" => "E",
                   "GAG" => "E",
                   "TGT" => "C",
                   "TGC" => "C",
                   "TGG" => "W",
                   "CGT" => "R",
                   "CGC" => "R",
                   "CGA" => "R",
                   "CGG" => "R",
                   "AGA" => "R",
                   "AGG" => "R",
                   "AGT" => "S",
                   "AGC" => "S",
                   "GGT" => "G",
                   "GGC" => "G",
                   "GGA" => "G",
                   "GGG" => "G");
           return %tmp;
}

sub check_overlap
{
    my ($id1,$id2) = @_;
    my @tmp1 = split /\_/, $id1;
    my @tmp2 = split /\_/, $id2;

    my $perc = 0;
    if($tmp1[0] ne $tmp2[0]){return 0;}
    
    my $sta1 = $tmp1[1];
    my $stp1 = $tmp1[2];
    my $sta2 = $tmp2[1];
    my $stp2 = $tmp2[2];

    my $flag = is_overlap_coordinates($sta1,$stp1,$sta2,$stp2);
    if($flag){return 1;}else{return 0;}
}

sub is_overlap_coordinates
{
    my ($x1,$y1,$x2,$y2) = @_;
   
    my $x    = "$x1 $y1";
    my $y    = "$x2 $y2";
    my $flag = is_overlap($x,$y);

    return $flag;
}

sub is_overlap
{
    my ($x,$y) = @_;
   
    my @x = split /\s+/, $x;
    my @y = split /\s+/, $y;

    if($x[0] <= $y[0] and $y[0] <= $x[1]){return 1;}
    if($y[0] <= $x[0] and $x[0] <= $y[1]){return 1;}

    return 0;
}

sub max
{
    my @vals = sort {$b <=> $a} @_;
       return $vals[0];
}

sub min
{
    my @vals = sort {$a <=> $b} @_;
       return $vals[0];
}

sub ttest
{
         my ($xref,$yref) = @_;
         my @x  = @{$xref};py('@x',@x);
         my @y  = @{$yref};py('@y',@y);

         my $x   = "c(";foreach my $e (@x){$x .= "$e,";} $x =~ s/\,$/\)/;
         my $y   = "c(";foreach my $e (@y){$y .= "$e,";} $y =~ s/\,$/\)/;
         my $r   = "t = t.test(x = $x, y = $y,paired = TRUE,na.action=na.omi)\nt\$p.value\n";
         my $otf = "ttest.r.txt";write_to_file($otf,$r);
         my $cmd = "R CMD BATCH $otf";print "$cmd\n";system($cmd);
         my $pval= "";
         my $fr  = openr("$otf.Rout");
            while(!eof($fr))
            {
                $line = <$fr>;
                if($line =~ /p.value/)
                {
                   $line = <$fr>;chomp $line;
                   $line =~ s/.+\s+(.+)/$1/;
                   $pval = $line;
                }
            }close($fr);

            $cmd = "rm -rf $otf";system($cmd);
            $cmd = "rm -rf $otf.Rout";system($cmd);

            return $pval;
}

sub t_test
{
         my ($xref,$yref) = @_;
             ttest($xref,$yref);
}

sub wilson_sign_test
{
     my ($xref,$yref) = @_;

     my (@x,@y,@z,$big,$otf,$x,$y,$p);

        @x = @{$xref};
        @y = @{$yref};

        $otf = "_tmp_.txt";
        $x   = join(",",@x);
        $y   = join(",",@y);
        $big = "$x\n$y\n";
        write_to_file($otf,$big);

       ####### write R script ######
       
        $big  = "raw = read.csv(\"_tmp_.txt\",header=F)\n";
        $big .= "x   = as.numeric(raw[1,])\n";
        $big .= "y   = as.numeric(raw[2,])\n";
        $big .= "z   = wilcox.test(x, y, paired = TRUE)\n";
        $big .= "write.table(z\$p.value,file=\"_out_.txt\",quote=F,col.names=F)\n";
        write_to_file("_run_.r",$big);
       
        $cmd = "Rscript _run_.r";runcmd($cmd,1);

        @z   = get_lines("_out_.txt");
        @m   = split /\s+/, $z[0];
        $p   = $m[1];

        #print "$big$p\n";
        $cmd = "rm -rf _tmp_.txt";runcmd($cmd,1);
        $cmd = "rm -rf _out_.txt";runcmd($cmd,1);
        $cmd = "rm -rf _run_.txt";runcmd($cmd,1);
        return $p;
}

sub is_number
{
    # use what Perl thinks is a number first
    # this is purely for speed, since the more complicated REGEX below should
    #  correctly handle all numeric cases
    if (looks_like_number($_[0]))
    {
        # Perl treats infinities as numbers, Excel does not
        if ($_[0] =~ /^[-+]*inf/)
        {
            return 0;
        }
        
        return 1;
    }

    # optional + or - sign at beginning
    # then require either:
    #  a number followed by optional comma stuff, then optional decimal stuff
    #  mandatory decimal, followed by optional digits
    # then optional exponent stuff
    #
    # Perl cannot handle American comma separators within long numbers.
    # Excel does, so we have to check for it.
    # Excel doesn't handle European dot separators, at least not when it is
    #  set to the US locale (my test environment).  I am going to leave this
    #  unsupported for now.
    #
    if ($_[0] =~ /^([-+]?)([0-9]+(,[0-9]{3,})*\.?[0-9]*|\.[0-9]*)([Ee]([-+]?[0-9]+))?$/)
    {
        # current REGEX can treat '.' as a number, check for that
        if ($_[0] eq '.')
        {
            return 0;
        }
        
        return 1;
    }
    
    return 0;
}

sub get_iupac_symbol
{
       my %temp = ();

          $temp{'A'} = "A";
          $temp{'C'} = "C";
          $temp{'G'} = "G";
          $temp{'T'} = "T";
          $temp{'U'} = "U";
          $temp{'W'} = "A,T";
          $temp{'S'} = "C,G";
          $temp{'M'} = "A,C";
          $temp{'K'} = "G,T";
          $temp{'R'} = "A,G";
          $temp{'Y'} = "C,T";
          $temp{'B'} = "C,G,T";
          $temp{'D'} = "A,G,T";
          $temp{'H'} = "A,C,T";
          $temp{'V'} = "A,C,G";
          $temp{'N'} = "A,C,G,T";
        return %temp;
}



sub run_commands_in_cluster
{
        my ($typ,$mode,$mem,$module,%cmds)    = @_;
        
        my (@cmds,@cmds1,$marker);

        @cmds       = sort keys %cmds;
        $size       = @cmds;
        if($size    > 0)
        {
        @cmds1      = get_cluster_commands($typ,$mem,$module,%cmds);if($mode eq 'print'){py1('@cmds1',@cmds1);}  
        if($mode eq 'none'){foreach $c (@cmds){runcmd($c,1);}}
        else{
        $marker     = get_marker($cmds1[0]);ps('$marker',$marker);
        foreach my $c (@cmds1)
        {
           runcmd($c,1);#control_job_number($marker);
        } #wait_for_job_done($marker);
        }
        }
}

sub get_marker
{
        my $cmd     = shift;#print "cmd = $cmd\n";
        my $marker  = 'na';
        if($cmd     =~ /\-j\s+(.*?)\_\d+\s+/)
        {
           $marker = $1;
        }
        if($marker eq 'na'){die "wrong marker value = $marker\n";}
        return $marker;
}

sub get_cluster_commands
{
        my ($typ,$mem,$module,%cmds)    = @_;
        
        my ($marker,$counter,%allcmds,$c,$log,$job,$cmd,$cmd1,@cmds,$mod);
           $counter = 0;
           $marker  = generate_batch_marker();
           @cmds    = sort keys %cmds;py('@cmds',@cmds);
           foreach my $cmd (@cmds)
           {
                   $counter++;

                   $job  = substr("$marker\_$counter",0,10);

                   $log  = $cmds{$cmd};

                   $mod  = "java-1.6,bwa/0.5.9,bowtie2,perl-5.22,python/2.7.6,R/3.1.1,fusioncatcher/20170816,samtools/0.1.18,star/2.5.4b,mhc_i/2.17,mhc_ii/2.16.3,htseq/0.6.1";

                   if($typ =~ /slm/)
                   {
                   $cmd1 = "submitslm -g $log -d $module -m $mem -j $job \'$cmd\'";$allcmds{$cmd1}=1;
                   }
                   else
                   {
                   $cmd1 = "submitjob -g $log -d $module -m $mem -j $job \'$cmd\'";$allcmds{$cmd1}=1;
                   }
           }

           @cmds = sort keys %allcmds;

           return @cmds;
}

sub generate_batch_marker
{
     my $marker = "";

        for(my $i=0;$i<4;$i++)
        {
            my $c       = get_random_character();
               $marker .= $c;
        }
        return $marker;
}

sub get_random_character
{
    my ($seconds, $microseconds) = gettimeofday;
    my $jobname     = $seconds;

    my @letters = get_character_set();
   
    my $random  = int(rand(61)); 
    
       return $letters[$random];       
}

sub get_character_set
{
     my @letters = ();
        for(my $i=0;$i<10;$i++){push(@letters,$i);}
        for(my $i=65;$i<91;$i++){push(@letters,chr($i));}
        for(my $i=97;$i<123;$i++){push(@letters,chr($i));}
        return @letters;
}

sub control_job_number
{
    my $marker        = shift;
    my $maxallowedjob = 50;

    my $runjob  = getrunjob($marker);
       while($runjob && $runjob >= $maxallowedjob)
       {
             print "\nrunning jobs = $runjob above the maximum number. Waiting for nodes......\n";
             sleep(180);$runjob  = getrunjob($marker);
       }
}
sub getrunjob
{
      my $marker = shift;

      my $user      = "";
      my $runjob    = 0;
      my $tmp       = ".__tmp__$marker.txt";
      my $cmd       = "/usr/local/bin/qstat > $tmp 2>&1";system($cmd);sleep(3);
      my $fr        = openr($tmp);
      my %alljobs   = ();
      my %watjobs   = ();
      my %runjobs   = ();
      my %comjobs   = ();

         while(!eof($fr))
         {
            my $line = <$fr>;

               if($line =~ /^(\d+)\..+\s+$marker.+\d+\s+[Q]\s+batch/)
               {
                  $watjobs{$1} = 1;
                  $alljobs{$1} = 1;
                  $runjob++;
                  print "$line";
               }  
               elsif($line =~ /^(\d+)\..+\s+$marker.+\d+\s+[R]\s+batch/)
               {
                  $runjobs{$1} = 1;
                  $alljobs{$1} = 1;
                  $runjob++;
                  print "$line";
               }
               elsif($line =~ /^(\d+)\..+\s+$marker.+\d+\s+[C]\s+batch/)
               {
                  $comjobs{$1} = 1;
                  $alljobs{$1} = 1;
                 #$runjob++;#this line will permit all jobs status disappear before submit new jobs
                  print "$line";
               }
         }close($fr);

      system("rm -rf $tmp");



      ################# check memory #####################################
      
      my ($memused,$memrequired,$flag,$jobid,@jobs,%memoryused,%memoryrequired,$line,$stop);

      $cmd    = "/usr/local/bin/qstat -f > $tmp 2>&1";system($cmd);
      $fr     = openr($tmp);

         while(!eof($fr))
         {
               if(!$stop){$line = <$fr>;}

               if($line =~ /^Job\s+Id:\s+(\d+)\./)
               {
                  $stop         = 0;
                  $jobid        = $1;
                  if($runjobs{$jobid} || $comjobs{$jobid})
                  {
                     $memused      = 0;
                     $memrequired  = 0;
                     $flag         = 0;

                     while(!eof($fr) && $line =~ /\w+/)
                     {
                          $line = <$fr>;

                          if($line =~ /resources_used.mem\s+=\s+(\d+)kb/)
                          {  
                             $memused = round($1/1000000,4);$flag=1;
                             $memoryused{$jobid} = $memused;
                          }
                          if($line =~ /Resource_List.mem\s+=\s+(\d+)gb/)
                          {  
                             $memrequired            = $1;
                             $memoryrequired{$jobid} = $memrequired;
                          }
                     }
                  }
               }  
         }close($fr);

      system("rm -rf $tmp");

      @jobs = sort keys %runjobs;

      foreach my $j (@jobs)
      {
            if(!$memoryused{$j}){next;}

            print "job=$j mem_used=$memoryused{$j} Gb mem_req=$memoryrequired{$j} Gb\n";
      }

      @jobs = sort keys %comjobs;

      foreach my $j (@jobs)
      {
            print "job=$j mem_used=$memoryused{$j} Gb mem_req=$memoryrequired{$j} Gb\n";
      }

      if($runjob == 0){print "no running jobs.\n";return;}

      @jobs = sort keys %runjobs;

      $flag = 0;
      foreach my $j (@jobs)
      {
          if($memused > $memrequired)
          {
             $flag = 1;
             print "memory used ($memused) >  memory required ($memrequired) marker = $marker\n";
             $cmd  = "/usr/local/bin/qdel $j";system($cmd);print "$j killed\n";
          }
      }

      if($flag){print "memory used ($memused) >  memory required ($memrequired) marker = $marker\n";exit;}


      return $runjob;
}

sub wait_for_job_done
{
     my $marker  = shift;

     my $runjob  = getrunjob($marker);

     my $rest    = 60*3;

     my $time    = localtime();

     while($runjob)
     {
          my $sep    = "";
          if($runjob < 1000){$sep .= " ";}
          if($runjob < 100 ){$sep .= " ";}
          if($runjob < 10  ){$sep .= " ";}

          print "There are $sep$runjob jobs running ($time).\n"; 
          sleep($rest);
          $runjob = getrunjob($marker);
          $time   = localtime();
     }
}

sub merge_peak
{
    my ($id1,$id2) = @_;
    my @tmp1 = split /\_/, $id1;
    my @tmp2 = split /\_/, $id2;
    if($tmp1[0] ne $tmp2[0]){die "chr are different: $id1 $id2\n";}
    
    my $sta1 = $tmp1[1];
    my $stp1 = $tmp1[2];
    my $sta2 = $tmp2[1];
    my $stp2 = $tmp2[2];

    my $min  = $sta1;if($sta2 < $min){$min = $sta2;}
    my $max  = $stp1;if($stp2 > $max){$max = $stp2;}

    my $tmp = "$tmp1[0]\_$min\_$max";
       return $tmp;
}

sub md5sum
{
    my $inf = shift;
    my $cmd = "md5sum $inf > _temp.txt";system($cmd);
    my $fr  = openr("_temp.txt");
    my $line=<$fr>;close($fr);
    my @temp=split /\s+/, $line;
    my $val = $temp[0];
       system("rm -rf _temp.txt"); 
       return $val;
}
####################################################################################
