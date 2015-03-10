package eForge::ePlot;

use 5.010;
use strict;
use warnings FATAL => 'all';
use Sort::Naturally;

=head1 NAME

eForge::ePlot - Plotting utilities for eForge

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

our (@ISA, @EXPORT, @EXPORT_OK);
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(Chart dChart table); # Symbols to export by default
@EXPORT_OK = qw(rChart);

=head1 SYNOPSIS

Provise plotting utilities for different plots to eForge

=head1 EXPORT

Chart
dChart
table

=head1 SUBROUTINES/METHODS

=head2 Chart

This is the original code using standard R plot to generate a static pdf.

=cut


sub Chart{
    print "Making static chart.\n";
    my ($filename, $lab, $resultsdir, $tissues, $cells, $label, $t1, $t2, $data) = @_;
    my $Rdir = $resultsdir;
    my $chart = "$lab.chart.pdf";
    my $rfile = "$Rdir/$lab.chart.R";
    #set some colors
    my ($sig, $msig, $ns, $abline, $tline) = qw(red palevioletred1 steelblue3 lightpink1 burlywood3); #alternate msig = pink2


    open my $rfh, ">", $rfile;
#results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t1, $t2, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE) # 99 and 95% CIs 1, 2, 3
$t1 = sprintf("%.2f", $t1);
$t2 = sprintf("%.2f", $t2);
    print $rfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\",header=TRUE,sep=\"\\t\")
results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t1, $t2, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)
pdf(\"$chart\", width=22.4, height=8)
ymin1= min(results\$Pvalue, na.rm=TRUE)*1.1
ymax1 = max(results\$Pvalue, na.rm=TRUE)*1.1
ymax = max(c(abs(ymin1),ymax1))
ymin = -ymin1
par(mar=c(11.5,4,3,1)+0.1)
tissue.cell.order <- unique(results[, c(\"Tissue\", \"Cell\")])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

palette(c(\"$sig\",\"$msig\",\"white\"))

plot(NA,ylab=\"-log10 binomial P\",xlab=\"\",main=\"MVPs in DNase1 sites (probably TF sites) in cell lines for $data $label\",ylim=c(ymin,ymax), las=2, pch=19, col = results\$Class2, xaxt='n', xlim=c(0,length(levels(results\$Cell))))

abline(h=par('yaxp')[1]:par('yaxp')[2],lty=1, lwd=0.1, col='#e0e0e0')

points(results\$Cell, results\$Pvalue, pch=19, col = results\$Class2, xaxt='n')

palette(c('black', '$msig', '$ns'))

points(results\$Cell, results\$Pvalue, pch=1, col = results\$Class2, xaxt='n')

palette(c(\"$ns\",\"$msig\",\"$sig\"))

axis(1, seq(1,length(levels(results\$Cell))),labels=levels(results\$Cell), las=2, cex.axis=0.7)

mtext(1,text=\"Cell\",line=7,cex=1.2)

legend('topleft', pch=19, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.005, 0.025), box.col='white', title='FDR q-value', text.col='white', bg='white')

palette(c('$ns', '$msig', 'black'))

legend('topleft', pch=1, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='darkgrey', title='FDR q-value')

tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))

abline(v=tissues[2:(length(tissues)-1)]+0.5,lty=6,col='#a0a0a0')

text((tissues[1:(length(tissues)-1)] + tissues[2:length(tissues)]) / 2 + 0.5, ymax, names(tissues[2:length(tissues)]), col=\"$tline\", adj=1, srt=90, cex=0.8) 

palette(\"default\")
dev.off()
";

#run the R code
    system("R", "--no-save", "--quiet", "--slave", "--file=$rfile");
  }

=head1 SUBROUTINES/METHODS

=head2 rChart

Makes a polycharts plot : note X axis labelling is problematical

=cut

sub rChart{

    print "Making rChart.\n";
    my ($filename, $lab, $resultsdir, $data, $label, $t1, $t2) = @_;
    my $chart = "$lab.rchart.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.rChart.R";
    open my $rcfh, ">", $rfile;
    print $rcfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\", header = TRUE, sep=\"\\t\")
results\$Colour<- 0 + (results\$Pvalue < $t2) + (results\$Pvalue < $t1)  # 99 and 95% CIs
require(rCharts)
r1 <- rPlot(Pvalue ~ Cell, data=results, color=\"bin(Colour, 0.25)\", type=\"point\", tooltip = \"function(item){ return (item.Pvalue + '\\\\n' + item.Cell + '\\\\n' + item.Tissue + '\\\\n' + item.File + '\\\\n' + item.Probe + '\\\\n' + item.Accession + '\\\\n')}\")
#r1\$guides(color=list(scale = list(type = \'gradient\', lower = \'\#CCC\', upper = \'\#000\'))) # optional code to make a grey scale
r1\$addParams(width = 2000, height=600, title=\"$label overlaps with $data DHS\")
ymin1 = min(results\$Pvalue, na.rm=TRUE)*1.2
ymax1 = max(results\$Pvalue, na.rm=TRUE)*1.2
ymax = max(c(abs(ymin1),ymax1))
ymin = -ymax
r1\$guides(x = list(numticks = length(unique(results\$Cell)), levels=results\$Cell), y = list(min = ymin, max = ymax))
r1\$save('$chart', cdn = F)
##r1\$show() #makes a temp file\n";

system("R", "--no-save", "--quiet", "--slave", "--file=$rfile");
  }

=head1 SUBROUTINES/METHODS

=head2 dChart

Make dimple interactive chart.

=cut

sub dChart{
    my ($filename, $lab, $resultsdir, $data, $label, $t1, $t2, $web) = @_;

    print "Making dChart.\n";
    my $chart = "$lab.dchart.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.dChart.R";
   open my $rcfh, ">", $rfile;
    print $rcfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\", header = TRUE, sep=\"\\t\")
results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t1, $t2, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE) # 99 and 95% CIs 1, 2, 3
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)
tissue.cell.order <- unique(results[, c(\"Tissue\", \"Cell\")])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))
require(rCharts)
d1 <- dPlot(
  y = \"Pvalue\",
  x = c(\"Cell\"),
  groups = c(\"Tissue\", \"Cell\", \"Probe\", \"Number\", \"Accession\", \"Pvalue\"),
  data = results,
  type = \"bubble\",
  width = 2000,
  height = 1500,
  bounds = list(x=90,y=50,height=600,width=1850),
  id = \"chart.$lab\"
)\n";
    if ($data =~ /erc/){
    print $rcfh "d1\$xAxis( type = \"addCategoryAxis\", grouporderRule = \"Cell\", orderRule = tissue.cell.order[,2])\n";
  }
else {
    print $rcfh "d1\$xAxis( type = \"addCategoryAxis\", grouporderRule = \"Tissue\", orderRule = \"Number\")\n";
  }

print $rcfh "d1\$yAxis( type = \"addMeasureAxis\" )
d1\$colorAxis(
   type = \"addColorAxis\",
   colorSeries = \"Class2\",
   palette = c('red', 'pink', 'lightblue'))

d1\$defaultColors(rep('lightgrey', 2))

d1\$layer(
  x = 'tissue',
  y = 'y',
  groups = c('y', 'tissue'),
  data = data.frame(tissue = rep(tissues[2:(length(tissues)-1)]/max(tissues)*100, each=2), y = rep(c(-10,20), length(tissues)-2), col=rep(1, 2*(length(tissues)-2)))
  ,type=\"line\"
  ,lineWeight=2
  ,xAxis = list(   type = 'addMeasureAxis', overrideMin=0, overrideMax=100)
  ,yAxis = list(   type = 'addMeasureAxis', overrideMin=0, overrideMax=10)
)

d1\$addParams(title=\"$label overlaps with $data DHS\")
d1\$save('$chart', cdn = F)\n";

    system("R", "--no-save", "--quiet", "--slave", "--file=$rfile");
    if ($web) {
        system("sed", "-e", "s/src='.*\\\/js/src='\\\/libraries\\\/dimple\\\/js/", "-i\"\"", "$resultsdir/$chart");
      }

  }

=head1 SUBROUTINES/METHODS

=head2 table

=cut

sub table{
    my ($filename, $lab, $resultsdir, $web) = @_;

    # Make Datatables table
    print "Making Table.\n";
    my $chart = "$lab.table.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.table.R";
    open my $rcfh, ">", $rfile;
    print $rcfh "setwd(\"$Rdir\")
    data<-read.table(\"$filename\", header = TRUE, sep=\"\\t\")
    results<-data.frame(data\$Cell, data\$Tissue, data\$Accession, data\$Pvalue, data\$Probe)
    names(results)<-c(\"Cell\", \"Tissue\", \"Accession\", \"Pvalue\", \"Probe\")
    require(rCharts)
    dt <- dTable(
      results,
      sScrollY= \"600\",
      bPaginate= F,
      sScrollX= \"100%\",
      width= \"680px\"
    )
    dt\$save('$chart', cdn = F)";
    system("R", "--no-save", "--quiet", "--slave", "--file=$rfile");
    if ($web) {
        system "sed -e \"s/href='.*\\\/css/href='\\\/libraries\\\/datatables\\\/css/; s/src='.*\\\/js/src='\\\/libraries\\\/datatables\\\/js/\" -i\"\" $resultsdir/$chart";
      }
  }

=head1 SUBROUTINES/METHODS

=head2 hChart

highcharts interface has problems with plotting - not used

=cut

#
#sub hChart{
#    print "Making hChart.\n";
#    my ($filename, $label) = @_;
#    my $chart = $lab . ".hchart.html";
#    my $Rdir = $cwd;
#    my $rfile = $lab . ".hChart.R";
#    open my $rcfh, ">", "$Rdir/$rfile";
#    print $rcfh "setwd(\"$Rdir\")
#results<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
#results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t1, $t2, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE) # 99 and 95% CIs 1, 2, 3
#require(rCharts)
#h1 <- hPlot(Pvalue ~ Number, data=results, type=\"scatter\", radius=5)
#h1\$addParams(width = 2000, height=600, title=list(\"$label overlaps with $data DHS\"))
#h1\$xAxis(title = list(text = \"Cell\"), labels = list(rotation=-90, align=\"right\"), categories = results\$Cell\)
#h1\$save('$chart', cdn = F)";
#
#system "R --no-save --quiet --slave < $rfile";
#}



=head1 AUTHOR

Charles Breeze, C<< <cbreeze at ebi.ac.uk> >>
Charles Breeze, C<< <c.breeze at ucl.ac.uk> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-eforge at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=eForge>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc eForge::ePlot


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Charles Breeze.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; version 2 dated June, 1991 or at your option
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available in the source tree;
if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA


=cut

1; # End of ePlot
