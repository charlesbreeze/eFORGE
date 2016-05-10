package eForge::ePlot;

=head1 NAME

eForge::ePlot - Plotting utilities for eForge

=head1 VERSION

Version 0.01

=head1 LICENCE AND COPYRIGHT

Copyright (C) [2014-2015] EMBL - European Bioinformatics Institute and University College London

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

=head1 CONTACT

Charles Breeze, C<< <c.breeze at ucl.ac.uk> >>

Javier Herrero, C<< <javier.herrero at ucl.ac.uk> >>

=head1 ACKNOWLEDGEMENTS

This software is based on the FORGE tool developed by Ian Dunham at the EMBL-EBI

=cut

use 5.010;
use strict;
use warnings FATAL => 'all';
use Sort::Naturally;

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
    my ($filename, $lab, $resultsdir, $tissues, $cells, $label, $t_marginal, $t_strict, $data) = @_;
    my $Rdir = $resultsdir;
    my $chart = "$lab.chart.pdf";
    my $rfile = "$Rdir/$lab.chart.R";
    #set some colors
    my ($sig, $msig, $ns, $abline, $tline) = qw(red palevioletred1 steelblue3 lightpink1 brown); #alternate msig = pink2


    open my $rfh, ">", $rfile;
#results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t_marginal, $t_strict, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE) # 99 and 95% CIs 1, 2, 3
$t_marginal = sprintf("%.2f", $t_marginal);
$t_strict = sprintf("%.2f", $t_strict);
    print $rfh "setwd('$Rdir')
results<-read.table('$filename', header=TRUE,sep='\\t')

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
# Collapse into a single string (to support same cell type in different tissues)
tissue.cell.order2 <- apply(tissue.cell.order, 1, paste, collapse = ' -- ')
results\$TissueCell <- apply(results[, c('Tissue', 'Cell')], 1, paste, collapse = ' -- ')
results\$TissueCell <- factor(results\$TissueCell, levels=tissue.cell.order2)

# Plot an empty chart first
pdf('$chart', width=22.4, height=8)
ymax = max(-log10(results\$Pvalue), na.rm=TRUE)*1.1
ymin = -0.1
par(mar=c(15.5,4,3,1)+0.1)
plot(NA,ylab='', xlab='', main='MVPs in DNase1 sites (probably TF sites) in cell lines for $data $label',
    ylim=c(ymin,ymax), las=2, pch=19, col = results\$Class2, xaxt='n', xlim=c(0,length(levels(results\$TissueCell))), cex.main=2)

# Add horizontal guide lines for the Y-axis
abline(h=par('yaxp')[1]:par('yaxp')[2],lty=1, lwd=0.1, col='#e0e0e0')

# Add vertical lines and labels to separate the tissues
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))
abline(v=tissues[2:(length(tissues)-1)]+0.5, lty=6, col='$tline')
text((tissues[1:(length(tissues)-1)] + tissues[2:length(tissues)]) / 2 + 0.5, ymax, names(tissues[2:length(tissues)]), col='$tline', adj=1, srt=90, cex=1.4) 

# Add points (internal color first)
palette(c('$sig', '$msig', 'white'))
points(results\$TissueCell, -log10(results\$Pvalue), pch=19, col = results\$Class2, xaxt='n')

# Add contour to the points
palette(c('black', '$msig', '$ns'))
points(results\$TissueCell, -log10(results\$Pvalue), pch=1, col = results\$Class2, xaxt='n')

# Add X-axis (use cell name only and not TissueCell)
axis(1, seq(1,length(tissue.cell.order[,2])), labels=tissue.cell.order[,2], las=2, cex.axis=0.9)
mtext(1, text='Cell', line=14, cex=1.4)
mtext(2, text='-log10 binomial p-value', line=2, cex=1.4)

# Add legend (internal color first)
palette(c('white', '$msig', '$sig'))
legend('topleft', pch=19, legend=c('q < $t_strict', 'q < $t_marginal', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='white', title='FDR q-value', text.col='white', bg='white')

# Add contour to the points in the legend
palette(c('$ns', '$msig', 'black'))
legend('topleft', pch=1, legend=c('q < $t_strict', 'q < $t_marginal', 'non-sig'), col = 3:1, cex=0.8, inset=c(0.001, 0.005), box.col='darkgrey', title='FDR q-value')

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
    my ($filename, $lab, $resultsdir, $data, $label, $t_marginal, $t_strict, $web) = @_;

    print "Making dChart.\n";
    my $chart = "$lab.dchart.html";
    my $Rdir = $resultsdir;
    my $rfile = "$Rdir/$lab.dChart.R";
   open my $rcfh, ">", $rfile;
    print $rcfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\", header = TRUE, sep=\"\\t\")

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)
color.axis.palette = c();
if (length(which(results\$Class2 == 1)) > 0 ) {
    color.axis.palette = c('red');
}
if (length(which(results\$Class2 == 2)) > 0 ) {
    color.axis.palette = c(color.axis.palette, '#FF82ab');
}
color.axis.palette = c(color.axis.palette, 'lightblue');
if (length(color.axis.palette) < 2) {
    color.axis.palette = c(color.axis.palette, 'lightblue'); # Add it twice to force the color if only non-significant values
}

results\$log10pvalue <- -log10(results\$Pvalue)

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
# Collapse into a single string (to support same cell type in different tissues)
tissue.cell.order2 <- apply(tissue.cell.order, 1, paste, collapse = ' -- ')
results\$TissueCell <- apply(results[, c('Tissue', 'Cell')], 1, paste, collapse = ' -- ')
results\$TissueCell <- factor(results\$TissueCell, levels=tissue.cell.order2)

# Count number of cell types for each tissue (to be able to draw the vertical separation lines afterwards
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))

require(rCharts)

dplot.height=900
dplot.width=2000
bounds.x=60
bounds.y=50
bounds.height=dplot.height - 300
bounds.width=dplot.width - bounds.x - 20

# Create a dimple plot, showing p-value vs cell, split data by tissue, cell, probe, etc to see individual points instead of aggregate avg
d1 <- dPlot(
  y = 'log10pvalue',
  x = c('TissueCell'),
  groups = c('TissueCell', 'Accession', 'Pvalue', 'Qvalue', 'Datatype', 'Probe'),
  data = results,
  type = 'bubble',
  width = dplot.width,
  height = dplot.height,
  bounds = list(x=bounds.x, y=bounds.y, height=bounds.height, width=bounds.width),
  id = 'chart.$lab'
)

# Force the order on the X-axis
d1\$xAxis( type = 'addCategoryAxis', grouporderRule = 'Cell', orderRule = tissue.cell.order[,2])
d1\$xAxis( type = 'addCategoryAxis', grouporderRule = 'TissueCell', orderRule = as.factor(tissue.cell.order2))

d1\$yAxis( type = 'addMeasureAxis' )

# Color points according to the q-value
d1\$colorAxis(
   type = 'addColorAxis',
   colorSeries = 'Class2',
   palette = color.axis.palette)

# Builds a JS string to add labels for tissues
labels.string = paste(paste0(\"
    // Adds labels for tissues
    myChart.svg.insert('text', 'g')
      .attr('x', 0)
      .attr('y', 0)
      .attr('font-size', 20)
      .attr('font-family', 'Arial')
      .style('fill', 'brown')
      .attr('transform', 'translate(\", (bounds.x - 5 + bounds.width * (tissues[1:(length(tissues)-1)] +  tissues[2:length(tissues)]) / (2 * max(tissues))), \",60) rotate(90)')
      .attr('text-anchor', 'top')
      .text('\", names(tissues[2:length(tissues)]), \"')
\"), collapse='')

# Builds a JS string to add vertical lines to separate tissues
lines.string = paste(paste0(\"
    // Adds vertical lines between tissues
    myChart.svg.append('line')
      .attr('x1', \", (bounds.x + bounds.width * tissues[2:(length(tissues)-1)]/ max(tissues)), \")
      .attr('y1', 50)
      .attr('x2', \", (bounds.x + bounds.width * tissues[2:(length(tissues)-1)]/ max(tissues)), \")
      .attr('y2', \", (50 + bounds.height), \")
      .style('stroke', 'brown')
      .style('stroke-dasharray', '10,3,3,3')
\"), collapse='')

# Adds some JS to be run after building the plot to get the image we want
d1\$setTemplate(afterScript = paste0(\"
  <script>
    myChart.draw()

    // Substitutes TissueCell labels in X-axis by Cell labels
    myChart.axes[1].shapes
      .selectAll('text')
        .text(function (d) { 
           var i;
           for (i = 0; i < data.length; i += 1) {
               if (data[i].TissueCell === d) {
                   return data[i].Cell;
               }
           }
       });

    // Adds title for X-axis
    myChart.axes[1].titleShape
        .style('font-size', 20)

    // Adds title for Y-axis
    myChart.axes[2].titleShape
        .style('font-size', 20)
        .text('-log10 binomial p-value')

    // Adds main title
    myChart.svg.append('text')
      .attr('x', \", (dplot.width / 2), \")
      .attr('y', \", (bounds.y / 2), \")
      .attr('font-size', 24)
      .attr('font-family', 'Arial')
      .attr('font-weight', 'bold')
      .style('fill', 'black')
      .attr('text-anchor', 'middle')
      .text('MVPs in DNase1 sites (probably TF sites) in cell lines for $data $label')
    \", labels.string, \"
    \", lines.string, \"
    // Adds vertical line at the far right of the plot
    myChart.svg.append('line')
      .attr('x1', \", (bounds.x + bounds.width), \")
      .attr('y1', \", bounds.y, \")
      .attr('x2', \", (bounds.x + bounds.width), \")
      .attr('y2', \", (bounds.y + bounds.height), \")
      .style('stroke', 'rgb(0,0,0)')
  </script>
\"))


d1\$save('$chart', cdn = F)\n";

    system("R", "--no-save", "--quiet", "--slave", "--file=$rfile");

    if ($web) {
        $web =~ s/\/$//;
        open(FILE, "$resultsdir/$chart") or die;
        my @lines = <FILE>;
        close(FILE);
        open(FILE, ">", "$resultsdir/$chart") or die;        
        foreach my $line (@lines) {
            $line =~ s/src='.*\/js/src='$web\/libraries\/dimple\/js/;
            print FILE $line;
        }
        close(FILE);
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
results <- subset(results, T, select = c('Cell', 'Tissue', 'Datatype', 'Accession', 'Pvalue', 'Qvalue', 'Probe'))
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
        $web =~ s/\/$//;
        open(FILE, "$resultsdir/$chart") or die;
        my @lines = <FILE>;
        close(FILE);
        open(FILE, ">", "$resultsdir/$chart") or die;        
        foreach my $line (@lines) {
            $line =~ s/href='.*\/css/href='$web\/libraries\/datatables\/css/;
            $line =~ s/src='.*\/js/src='$web\/libraries\/datatables\/js/;
            print FILE $line;
        }
        close(FILE);
    }
}

1;
