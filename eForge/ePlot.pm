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
    my ($sig, $msig, $ns, $abline, $tline, $tsline) = ('red', 'palevioletred1', 'steelblue3', 'lightpink1', 'burlywood3', '#e0e0e0'); #alternate msig = pink2


    open my $rfh, ">", $rfile;
#results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t1, $t2, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE) # 99 and 95% CIs 1, 2, 3
$t1 = sprintf("%.2f", $t1);
$t2 = sprintf("%.2f", $t2);
    print $rfh "setwd('$Rdir')
results<-read.table('$filename', header=TRUE,sep='\\t')

# Class splits the data into non-significant, marginally significant and significant according to $t1 and $t2 (in -log10 scale)
results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t1, $t2, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Plot an empty chart first
pdf('$chart', width=22.4, height=8)
ymax = max(results\$Pvalue, na.rm=TRUE)*1.1
ymin = -0.1
par(mar=c(11.5,4,3,1)+0.1)
plot(NA,ylab='-log10 binomial P', xlab='', main='MVPs in DNase1 sites (probably TF sites) in cell lines for $data $label',
    ylim=c(ymin,ymax), las=2, pch=19, col = results\$Class2, xaxt='n', xlim=c(0,length(levels(results\$Cell))))

# Add horizontal guide lines for the Y-axis
abline(h=par('yaxp')[1]:par('yaxp')[2],lty=1, lwd=0.1, col='#e0e0e0')

# Add points (internal color first)
palette(c('$sig', '$msig', 'white'))
points(results\$Cell, results\$Pvalue, pch=19, col = results\$Class2, xaxt='n')

# Add contour to the points
palette(c('black', '$msig', '$ns'))
points(results\$Cell, results\$Pvalue, pch=1, col = results\$Class2, xaxt='n')

# Add X-axis
axis(1, seq(1,length(levels(results\$Cell))), labels=levels(results\$Cell), las=2, cex.axis=0.7)
mtext(1, text='Cell', line=7, cex=1.2)

# Add legend (internal color first)
palette(c('white', '$msig', '$sig'))
legend('topleft', pch=19, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='white', title='FDR q-value', text.col='white', bg='white')

# Add contour to the points in the legend
palette(c('$ns', '$msig', 'black'))
legend('topleft', pch=1, legend=c('q < 0.01', 'q < 0.05', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='darkgrey', title='FDR q-value')

# Add vertical lines and labels to separate the tissues
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))
abline(v=tissues[2:(length(tissues)-1)]+0.5, lty=6, col='$tsline')
text((tissues[1:(length(tissues)-1)] + tissues[2:length(tissues)]) / 2 + 0.5, ymax, names(tissues[2:length(tissues)]), col='$tline', adj=1, srt=90, cex=0.8) 

palette('default')
dev.off()
";

#run the R code
    system("R", "--no-save", "--quiet", "--slave", "--file=$rfile");
  }


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

# Class splits the data into non-significant, marginally significant and significant according to $t1 and $t2 (in -log10 scale)
results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t1, $t2, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, 0.01, 0.05, 1), labels=FALSE, include.lowest=TRUE)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))

require(rCharts)

# Create a dimple plot, showing p-value vs cell, split data by tissue, cell, probe, etc to see individual points instead of aggregate avg
d1 <- dPlot(
  y = 'Pvalue',
  x = c('Cell'),
  groups = c('Tissue', 'Cell', 'Probe', 'Number', 'Accession', 'Pvalue'),
  data = results,
  type = 'bubble',
  width = 2000,
  height = 1500,
  bounds = list(x=90,y=50,height=600,width=1850),
  id = 'chart.$lab'
)

# Force the order on te X-axis
d1\$xAxis( type = 'addCategoryAxis', grouporderRule = 'Cell', orderRule = tissue.cell.order[,2])

d1\$yAxis( type = 'addMeasureAxis' )

# Color points according to the q-value
d1\$colorAxis(
   type = 'addColorAxis',
   colorSeries = 'Class2',
   palette = c('red', 'pink', 'lightblue'))

# Set default colors for the lines separating the tissue (we need an array here)
d1\$defaultColors(rep('#e0e0e0', 2))

# Include a layer with new axes to show the lines separating the tissues
d1\$layer(
  x = 'tissue',
  y = 'y',
  groups = c('y', 'tissue'),
  data = data.frame(tissue = rep(tissues[2:(length(tissues))]/max(tissues)*100, each=2), y = rep(c(-10,20), length(tissues)-1), col=rep(1, 2*(length(tissues)-1))),
  type = 'line',
  lineWeight = 1,
  xAxis = list(   type = 'addMeasureAxis', overrideMin=0, overrideMax=100),
  yAxis = list(   type = 'addMeasureAxis', overrideMin=0, overrideMax=10)
)

# Include title
d1\$addParams(title='$label overlaps with $data DHS')

d1\$save('$chart', cdn = F)\n";

    system("R", "--no-save", "--quiet", "--slave", "--file=$rfile");
    if ($web) {
        system("sed", "-e", "s/src='.*\\\/js/src='\\\/libraries\\\/dimple\\\/js/", "-i\"\"", "$resultsdir/$chart");
      }

  }


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
    print $rcfh "setwd('$Rdir')
results <- read.table('$filename', header = TRUE, sep='\\t')
results <- subset(results, T, select = c('Cell', 'Tissue', 'Accession', 'Pvalue', 'Qvalue', 'Probe'))
require(rCharts)
dt <- dTable(
    results,
    sScrollY = '600',
    bPaginate = F,
    sScrollX = '100%',
    width = '680px'
)
dt\$save('$chart', cdn = F)\n";
    system("R", "--no-save", "--quiet", "--slave", "--file=$rfile");
    if ($web) {
        system "sed -e \"s/href='.*\\\/css/href='\\\/libraries\\\/datatables\\\/css/; s/src='.*\\\/js/src='\\\/libraries\\\/datatables\\\/js/\" -i\"\" $resultsdir/$chart";
      }
  }


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
