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
@EXPORT = qw(Chart ChartAllH3Marks ChartAllChromatin15StateMarks dChart table); # Symbols to export by default

my $sig_Rcolor_default =    "red";
my $msig_Rcolor_default =   "palevioletred1";
my $ns_Rcolor_default =     "steelblue3";
my $abline_Rcolor_default = "lightpink1";
my $tline_Rcolor_default =  "grey";
my $tline_hex_default =     "'#A9A9A9'";

=head1 SYNOPSIS

Provise plotting utilities for different plots to eForge

=head1 EXPORT

Chart
ChartAllH3Marks
ChartAllChromatin15StateMarks
dChart
table

=head1 SUBROUTINES/METHODS

=head2 Chart

This is the original code using standard R plot to generate a static pdf.

=cut

sub Chart{
    print "Making generic static chart.\n";
    my ($filename, $lab, $resultsdir, $tissues, $cells, $label, $t_marginal, $t_strict, $data) = @_;
    my $Rdir = $resultsdir;
    my $chart = "$lab.chart.pdf";
    my $rfile = "$Rdir/$lab.chart.R";
    #set some colors
    my ($sig, $msig, $ns, $abline, $tline) = ($sig_Rcolor_default, $msig_Rcolor_default, $ns_Rcolor_default, $abline_Rcolor_default, $tline_Rcolor_default); #alternate msig = pink2


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
plot(NA,ylab='', xlab='', main='DMPs in DNase I sites (probably TF sites) in cell lines for $data $label',
    ylim=c(ymin,ymax), las=2, pch=19, col = results\$Class2, xaxt='n', xlim=c(0,length(levels(results\$TissueCell))), cex.main=2)

# Add horizontal guide lines for the Y-axis
abline(h=par('yaxp')[1]:par('yaxp')[2],lty=1, lwd=0.1, col='#e0e0e0')

# Add vertical lines and labels to separate the tissues
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))
abline(v=tissues[2:(length(tissues)-1)]+0.5, lty=6, col='$tline')
text((tissues[1:(length(tissues)-1)] + tissues[2:length(tissues)]) / 2 + 0.5, ymax, names(tissues[2:length(tissues)]), col='$tline', adj=1, srt=90, cex=1.2) 

# Add points (internal color first)
palette(c('$sig', '$msig', 'white'))
points(results\$TissueCell, -log10(results\$Pvalue), pch=19, col = results\$Class2, xaxt='n')

# Add contour to the points
palette(c('black', '$msig', '$ns'))
points(results\$TissueCell, -log10(results\$Pvalue), pch=1, col = results\$Class2, xaxt='n')

# Add X-axis (use cell name only and not TissueCell)
axis(1, seq(1,length(tissue.cell.order[,2])), labels=tissue.cell.order[,2], las=2, cex.axis=0.67)
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
  
  
=head2 ChartAllH3Marks

This is the original code using standard R plot to generate a static
pdf, which is modified to plot when the "All H3 Marks" dataset is 
selected from Analysis Options

=cut


sub ChartAllH3Marks{
    print "Making all-H3-marks static chart.\n";
    my ($filename, $lab, $resultsdir, $tissues, $cells, $label, $t_marginal, $t_strict, $data) = @_;
    my $Rdir = $resultsdir;
    my $chart = "$lab.chart.pdf";
    my $rfile = "$Rdir/$lab.chart.R";
    #set some colors
    my ($sig, $msig, $ns, $abline, $tline) = ($sig_Rcolor_default, $msig_Rcolor_default, $ns_Rcolor_default, $abline_Rcolor_default, $tline_Rcolor_default); #alternate msig = pink2


    open my $rfh, ">", $rfile;
#results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t_marginal, $t_strict, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE) # 99 and 95% CIs 1, 2, 3
$t_marginal = sprintf("%.2f", $t_marginal);
$t_strict = sprintf("%.2f", $t_strict);
    print $rfh "setwd('$Rdir')
 
#
# Initialize plot
#

all_results <- read.table('$filename', header=TRUE,sep='\\t')
    
# Plot an empty chart first
pdf('$chart', width=22.4, height=8)
ymax = max(-log10(all_results\$Pvalue), na.rm=TRUE)*1.1
ymin = -0.1
par(mar=c(15.5, 4, 3, 1)+0.1)
plot(NA,ylab='', xlab='', main='DMPs in DNase I sites (probably TF sites) in cell lines for $data $label',
    ylim=c(ymin,ymax), las=2, pch=19, col = all_results\$Class2, xaxt='n', xlim=c(0,length(levels(all_results\$Cell))), cex.main=2)

# Add horizontal guide lines for the Y-axis
abline(h=par('yaxp')[1]:par('yaxp')[2],lty=1, lwd=0.1, col='#e0e0e0')

# Add vertical lines and labels to separate the tissues
tissue.cell.order <- unique(all_results[, c('Tissue', 'Cell')])
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))
abline(v=tissues[2:(length(tissues)-1)]+0.5, lty=6, col='$tline')
text((tissues[1:(length(tissues)-1)] + tissues[2:length(tissues)]) / 2 + 0.5, ymax, names(tissues[2:length(tissues)]), col='$tline', adj=1, srt=90, cex=1.4)
    
#
# H3K4me1
#

results <- all_results[all_results\$Datatype == 'H3K4me1',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 19

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- '#e41a1c'

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2]) 

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col = col_by_qval, xaxt='n', lwd=0.5)

# Add contour to the points
# palette(c('black', 'steelblue3', 'steelblue3'))
# points(results\$Cell, -log10(results\$Pvalue), pch=1, col = results\$Class2, xaxt='n', lwd=0.5)

# Add X-axis
axis(1, seq(1,length(levels(results\$Cell))), labels=levels(results\$Cell), las=2, cex.axis=0.67)
mtext(1, text='Cell', line=14, cex=1.4)
mtext(2, text='-log10 binomial p-value', line=2, cex=1.4)

palette('default')

#
# H3K4me3
#

results <- all_results[all_results\$Datatype == 'H3K4me3',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 19

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- '#4daf4a'

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col = col_by_qval, xaxt='n', lwd=0.5)

# Add contour to the points
# palette(c('black', 'steelblue3', 'steelblue3'))
# points(results\$Cell, -log10(results\$Pvalue), pch=1, col = results\$Class2, xaxt='n', lwd=0.5)

palette('default')

#
# H3K27me3
#

results <- all_results[all_results\$Datatype == 'H3K27me3',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 19

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- '#e5e500'

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col = col_by_qval, xaxt='n', lwd=0.5)

# Add contour to the points
# palette(c('black', 'steelblue3', 'steelblue3'))
# points(results\$Cell, -log10(results\$Pvalue), pch=1, col = results\$Class2, xaxt='n', lwd=0.5)

palette('default')

#
# H3K36me3
#

results <- all_results[all_results\$Datatype == 'H3K36me3',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 19

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- '#984ea3'

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, xaxt='n', lwd=0.5)

# Add contour to the points
# palette(c('black', 'steelblue3', 'steelblue3'))
# points(results\$Cell, -log10(results\$Pvalue), pch=1, col = results\$Class2, xaxt='n', lwd=0.5)

palette('default')

#
# H3K9me3
#

results <- all_results[all_results\$Datatype == 'H3K9me3',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 19

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- '#ff7f00'

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, xaxt='n', lwd=0.5)

# Add contour to the points
# palette(c('black', 'steelblue3', 'steelblue3'))
# points(results\$Cell, -log10(results\$Pvalue), pch=1, col = results\$Class2, xaxt='n', lwd=0.5)

palette('default')

#
# Wrap-up
#

# strict
l1 <- legend(
  'topleft'
  , legend=c('H3K27me3'
    , 'H3K4me1'
    , 'H3K4me3'
    , 'H3K36me3'
    , 'H3K9me3')
  , title = expression(bold('FDR q < $t_strict')),
  , col=c('#e5e500', '#e41a1c', '#4daf4a', '#984ea3', '#ff7f00')
  , pch=19
  , inset=c(0.001, 0.005)
  , cex=0.8
  , box.col='white'
  , bty='n')

# marginally-significant
l2 <- legend(
  x = l1\$rect\$left, y = with(l1\$rect, top - h)
  , legend=c('H3K27me3'
    , 'H3K4me1'
    , 'H3K4me3'
    , 'H3K36me3'
    , 'H3K9me3')
  , title = expression(bold('FDR q < $t_marginal')),
  , col=c('#e5e500', '#e41a1c', '#4daf4a', '#984ea3', '#ff7f00')
  , pch=1
  , inset=c(0.001, 0.005)
  , cex=0.8
  , box.col='white'
  , bty='n')
  
# other
l3 <- legend(
  x = l2\$rect\$left, y = with(l2\$rect, top - h)
  , legend=c('Non-signif.')
  , title = expression(bold('Other')),
  , col='$ns'
  , pch=1
  , inset=c(0.001, 0.005)
  , cex=0.8
  , box.col='white'
  , bty='n')
  
# add the rectangle around the legend
rect(xleft = l1\$rect\$left, ybottom = with(l3\$rect, top - h), 
     xright = l1\$rect\$left + max(l1\$rect\$w, l2\$rect\$w, l3\$rect\$w), ytop = l1\$rect\$top,
     col='white')

# redraw atop filled rectangle that is correctly sized

# strict
l1 <- legend(
  'topleft'
  , legend=c('H3K27me3'
    , 'H3K4me1'
    , 'H3K4me3'
    , 'H3K36me3'
    , 'H3K9me3')
  , title = expression(bold('FDR q < $t_strict')),
  , col=c('#e5e500', '#e41a1c', '#4daf4a', '#984ea3', '#ff7f00')
  , pch=19
  , inset=c(0.001, 0.005)
  , cex=0.8
  , box.col='white'
  , bty='n')

# marginally-significant
l2 <- legend(
  x = l1\$rect\$left, y = with(l1\$rect, top - h)
  , legend=c('H3K27me3'
    , 'H3K4me1'
    , 'H3K4me3'
    , 'H3K36me3'
    , 'H3K9me3')
  , title = expression(bold('FDR q < $t_marginal')),
  , col=c('#e5e500', '#e41a1c', '#4daf4a', '#984ea3', '#ff7f00')
  , pch=1
  , inset=c(0.001, 0.005)
  , cex=0.8
  , box.col='white'
  , bty='n')
  
# other
l3 <- legend(
  x = l2\$rect\$left, y = with(l2\$rect, top - h)
  , legend=c('Non-signif.')
  , title = expression(bold('Other')),
  , col='$ns'
  , pch=1
  , inset=c(0.001, 0.005)
  , cex=0.8
  , box.col='white'
  , bty='n')
    
# add the rectangle around the legend
rect(xleft = l1\$rect\$left, ybottom = with(l3\$rect, top - h), 
     xright = l1\$rect\$left + max(l1\$rect\$w, l2\$rect\$w, l3\$rect\$w), ytop = l1\$rect\$top, 
     lwd=0.5)

dev.off()
";

#run the R code
    system("R", "--no-save", "--quiet", "--slave", "--file=$rfile");
  }
  
sub ChartAllChromatin15StateMarks{
    print "Making all-chromatin-15-state-marks static chart.\n";
    my ($filename, $lab, $resultsdir, $tissues, $cells, $label, $t_marginal, $t_strict, $data) = @_;
    my $Rdir = $resultsdir;
    my $chart = "$lab.chart.pdf";
    my $rfile = "$Rdir/$lab.chart.R";
    #set some colors
    my ($sig, $msig, $ns, $abline, $tline) = ($sig_Rcolor_default, $msig_Rcolor_default, $ns_Rcolor_default, $abline_Rcolor_default, $tline_Rcolor_default); #alternate msig = pink2


    open my $rfh, ">", $rfile;
#results\$Class<-cut(results\$Pvalue, breaks =c(min(results\$Pvalue), $t_marginal, $t_strict, max(results\$Pvalue)), labels=FALSE, include.lowest=TRUE) # 99 and 95% CIs 1, 2, 3
$t_marginal = sprintf("%.2f", $t_marginal);
$t_strict = sprintf("%.2f", $t_strict);
    print $rfh "setwd('$Rdir')
 
#
# Initialize plot
#

all_results <- read.table('$filename', header=TRUE,sep='\\t')

#
# Add functions to alter colors
#

darken <- function(color, factor=1.5){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


lighten <- function(color, factor=1.1){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  col
}
    
# Plot an empty chart first
pdf('$chart', width=40, height=8)
ymax = max(-log10(all_results\$Pvalue), na.rm=TRUE)*1.1
ymin = -0.1
par(mar=c(15.5, 4, 3, 1)+0.1)
plot(NA,ylab='', xlab='', main='DMPs in DNase I sites (probably TF sites) in cell lines for $data $label',
    ylim=c(ymin,ymax), las=2, pch=19, col = all_results\$Class2, xaxt='n', xlim=c(0,length(levels(all_results\$Cell))), cex.main=2)

# Add horizontal guide lines for the Y-axis
abline(h=par('yaxp')[1]:par('yaxp')[2],lty=1, lwd=0.1, col='#e0e0e0')

# Add vertical lines and labels to separate the tissues
tissue.cell.order <- unique(all_results[, c('Tissue', 'Cell')])
tissues <- c(0, cumsum(summary(tissue.cell.order[,'Tissue'])))
abline(v=tissues[2:(length(tissues)-1)]+0.5, lty=6, col='$tline')
text((tissues[1:(length(tissues)-1)] + tissues[2:length(tissues)]) / 2 + 0.5, ymax, names(tissues[2:length(tissues)]), col='$tline', adj=1, srt=90, cex=1.4)
    
#
# TssA
#

results <- all_results[all_results\$Datatype == 'TssA',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'Red'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('Red')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2]) 

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col = col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

# Add X-axis
axis(1, seq(1,length(levels(results\$Cell))), labels=levels(results\$Cell), las=2, cex.axis=0.67)
mtext(1, text='Cell', line=14, cex=1.4)
mtext(2, text='-log10 binomial p-value', line=2, cex=1.4)

palette('default')

#
# TssAFlnk
#

results <- all_results[all_results\$Datatype == 'TssAFlnk',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'OrangeRed'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('OrangeRed')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col = col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# TxFlnk
#

results <- all_results[all_results\$Datatype == 'TxFlnk',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'LimeGreen'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('LimeGreen')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col = col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# Tx
#

results <- all_results[all_results\$Datatype == 'Tx',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'Green'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('Green')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# TxWk
#

results <- all_results[all_results\$Datatype == 'TxWk',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'DarkGreen'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('DarkGreen')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# EnhG
#

results <- all_results[all_results\$Datatype == 'EnhG',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'GreenYellow'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('GreenYellow')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')


#
# Enh
#

results <- all_results[all_results\$Datatype == 'Enh',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- '#E5E500'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('#E5E500')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# ZNF-Rpts
#

results <- all_results[all_results\$Datatype == 'ZNF-Rpts',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- '#3D7B66'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('#3D7B66')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

#
# Het
#

results <- all_results[all_results\$Datatype == 'Het',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'PaleTurquoise'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('PaleTurquoise')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# TssBiv
#

results <- all_results[all_results\$Datatype == 'TssBiv',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'IndianRed'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('IndianRed')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# BivFlnk
#

results <- all_results[all_results\$Datatype == 'BivFlnk',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'DarkSalmon'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('DarkSalmon')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# EnhBiv
#

results <- all_results[all_results\$Datatype == 'EnhBiv',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'DarkKhaki'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('DarkKhaki')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# ReprPC
#

results <- all_results[all_results\$Datatype == 'ReprPC',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'gray50'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('gray50')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# ReprPCWk
#

results <- all_results[all_results\$Datatype == 'ReprPCWk',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'Gainsboro'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('Gainsboro')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# Quies
#

results <- all_results[all_results\$Datatype == 'Quies',]

# Class splits the data into non-significant, marginally significant and significant according to $t_marginal and $t_strict (in -log10 scale)
results\$Class <- cut(results\$Pvalue, breaks =c(0, $t_strict, $t_marginal, 1)/length(unique(results[,'Tissue'])), labels=FALSE, include.lowest=TRUE)

# Class splits the data into non-significant, marginally significant and significant according to q-value (B-Y FDR adjusted)
results\$Class2 <- cut(results\$Qvalue, breaks =c(0, $t_strict, $t_marginal, 1), labels=FALSE, include.lowest=TRUE)

pch_by_qval <- rep.int(1, length(results\$Class2))
pch_by_qval[which(results\$Qvalue < $t_strict)] <- 21

col_by_qval <- rep('steelblue3', length(results\$Class2))
col_by_qval[which(results\$Qvalue < $t_marginal)] <- 'whitesmoke'

bg_by_qval <- rep('white', length(results\$Class2))
bg_by_qval[which(results\$Qvalue < $t_marginal)] <- lighten('whitesmoke')

# Re-order the entries according to tissue first and then cell type/line
tissue.cell.order <- unique(results[, c('Tissue', 'Cell')])
tissue.cell.order <- tissue.cell.order[order(tissue.cell.order[,1], tissue.cell.order[,2]), ]
results\$Cell <- factor(results\$Cell, levels=tissue.cell.order[,2])

# Add points (internal color first)
points(results\$Cell, -log10(results\$Pvalue), pch=pch_by_qval, col=col_by_qval, bg = bg_by_qval, xaxt='n', lwd=0.5)

palette('default')

#
# Wrap-up
#

# strict
l1 <- legend(
  'topleft'
  , legend=c('Active TSS'
    , 'Flanking Active TSS'
    , 'Transcr. at gene 5p and 3p'
    , 'Strong transcription'
    , 'Weak transcription'
    , 'Genic enhancers'
    , 'Enhancers'
    , 'ZNF genes & repeats'
    , 'Heterochromatin'
    , 'Bivalent/Poised TSS'
    , 'Flanking Bivalent TSS/Enh'
    , 'Bivalent Enhancer'
    , 'Repressed PolyComb'
    , 'Weak Repressed PolyComb'
    , 'Quiescent/Low')
  , title = expression(bold('FDR q < $t_strict')),
  , col=sapply(c('Red', 'OrangeRed', 'LimeGreen', 'Green', 'DarkGreen', 'GreenYellow', '#E5E500', '#3D7B66', 'PaleTurquoise', 'IndianRed', 'DarkSalmon', 'DarkKhaki', 'gray50', 'Gainsboro', 'whitesmoke'), darken)
  , pt.bg=sapply(c('Red', 'OrangeRed', 'LimeGreen', 'Green', 'DarkGreen', 'GreenYellow', '#E5E500', '#3D7B66', 'PaleTurquoise', 'IndianRed', 'DarkSalmon', 'DarkKhaki', 'gray50', 'Gainsboro', 'whitesmoke'), lighten)
  , pch=21
  , inset=c(0.001, 0.005)
  , cex=0.55
  , box.col='white'
  , bty='n')

# marginally-significant
l2 <- legend(
  x = l1\$rect\$left, y = with(l1\$rect, top - h)
  , legend=c('Active TSS'
    , 'Flanking Active TSS'
    , 'Transcr. at gene 5p and 3p'
    , 'Strong transcription'
    , 'Weak transcription'
    , 'Genic enhancers'
    , 'Enhancers'
    , 'ZNF genes & repeats'
    , 'Heterochromatin'
    , 'Bivalent/Poised TSS'
    , 'Flanking Bivalent TSS/Enh'
    , 'Bivalent Enhancer'
    , 'Repressed PolyComb'
    , 'Weak Repressed PolyComb'
    , 'Quiescent/Low')
  , title = expression(bold('FDR q < $t_marginal')),
  , col=c('Red', 'OrangeRed', 'LimeGreen', 'Green', 'DarkGreen', 'GreenYellow', '#E5E500', '#3D7B66', 'PaleTurquoise', 'IndianRed', 'DarkSalmon', 'DarkKhaki', 'gray50', 'Gainsboro', 'whitesmoke')
  , pch=1
  , inset=c(0.001, 0.005)
  , cex=0.55
  , box.col='white'
  , bty='n')
  
# other
l3 <- legend(
  x = l2\$rect\$left, y = with(l2\$rect, top - h)
  , legend=c('Non-signif.')
  , title = expression(bold('Other')),
  , col='$ns'
  , pch=1
  , inset=c(0.001, 0.005)
  , cex=0.55
  , box.col='white'
  , bty='n')
  
# add the rectangle around the legend
rect(xleft = l1\$rect\$left, ybottom = with(l3\$rect, top - h), 
     xright = l1\$rect\$left + max(l1\$rect\$w, l2\$rect\$w, l3\$rect\$w), ytop = l1\$rect\$top,
     col='white')

# redraw atop filled rectangle that is correctly sized

# strict
l1 <- legend(
  'topleft'
  , legend=c('Active TSS'
    , 'Flanking Active TSS'
    , 'Transcr. at gene 5p and 3p'
    , 'Strong transcription'
    , 'Weak transcription'
    , 'Genic enhancers'
    , 'Enhancers'
    , 'ZNF genes & repeats'
    , 'Heterochromatin'
    , 'Bivalent/Poised TSS'
    , 'Flanking Bivalent TSS/Enh'
    , 'Bivalent Enhancer'
    , 'Repressed PolyComb'
    , 'Weak Repressed PolyComb'
    , 'Quiescent/Low')
  , title = expression(bold('FDR q < $t_strict')),
  , col=sapply(c('Red', 'OrangeRed', 'LimeGreen', 'Green', 'DarkGreen', 'GreenYellow', '#E5E500', '#3D7B66', 'PaleTurquoise', 'IndianRed', 'DarkSalmon', 'DarkKhaki', 'gray50', 'Gainsboro', 'whitesmoke'), darken)
  , pt.bg=sapply(c('Red', 'OrangeRed', 'LimeGreen', 'Green', 'DarkGreen', 'GreenYellow', '#E5E500', '#3D7B66', 'PaleTurquoise', 'IndianRed', 'DarkSalmon', 'DarkKhaki', 'gray50', 'Gainsboro', 'whitesmoke'), lighten)
  , pch=21
  , inset=c(0.001, 0.005)
  , cex=0.55
  , box.col='white'
  , bty='n')

# marginally-significant
l2 <- legend(
  x = l1\$rect\$left, y = with(l1\$rect, top - h)
  , legend=c('Active TSS'
    , 'Flanking Active TSS'
    , 'Transcr. at gene 5p and 3p'
    , 'Strong transcription'
    , 'Weak transcription'
    , 'Genic enhancers'
    , 'Enhancers'
    , 'ZNF genes & repeats'
    , 'Heterochromatin'
    , 'Bivalent/Poised TSS'
    , 'Flanking Bivalent TSS/Enh'
    , 'Bivalent Enhancer'
    , 'Repressed PolyComb'
    , 'Weak Repressed PolyComb'
    , 'Quiescent/Low')
  , title = expression(bold('FDR q < $t_marginal')),
  , col=c('Red', 'OrangeRed', 'LimeGreen', 'Green', 'DarkGreen', 'GreenYellow', '#E5E500', '#3D7B66', 'PaleTurquoise', 'IndianRed', 'DarkSalmon', 'DarkKhaki', 'gray50', 'Gainsboro', 'whitesmoke')
  , pch=1
  , inset=c(0.001, 0.005)
  , cex=0.55
  , box.col='white'
  , bty='n')
  
# other
l3 <- legend(
  x = l2\$rect\$left, y = with(l2\$rect, top - h)
  , legend=c('Non-signif.')
  , title = expression(bold('Other')),
  , col='$ns'
  , pch=1
  , inset=c(0.001, 0.005)
  , cex=0.55
  , box.col='white'
  , bty='n')
    
# add the rectangle around the legend
rect(xleft = l1\$rect\$left, ybottom = with(l3\$rect, top - h), 
     xright = l1\$rect\$left + max(l1\$rect\$w, l2\$rect\$w, l3\$rect\$w), ytop = l1\$rect\$top, 
     lwd=0.5)

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
    my $xAxisLabelShift = -0.7;
    if ($data eq "erc2-chromatin15state-all") {
      $xAxisLabelShift = -100.0;
    }
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
if (length(color.axis.palette) <= 2) {
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
  groups = c('TissueCell', 'Accession', 'Pvalue', 'Qvalue', 'Datatype', 'Probe', 'Class2'),
  data = results,
  type = 'bubble',
  width = dplot.width,
  height = dplot.height,
  bounds = list(x=bounds.x, y=bounds.y, height=bounds.height, width=bounds.width),
  id = 'chart.$lab'
)

# Force the order on the X-axis
d1\$xAxis( type = 'addCategoryAxis', grouporderRule = 'Cell',       orderRule = tissue.cell.order[,2])
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
      .attr('font-size', 16)
      .attr('font-family', 'Arial')
      .style('fill', $tline_hex_default)
      .attr('transform', 'translate(\", (bounds.x + 5 + bounds.width * (tissues[1:(length(tissues)-1)] +  tissues[2:length(tissues)]) / (2 * max(tissues))), \", 60) rotate(-90)')
      .attr('text-anchor', 'end')
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
      .style('stroke', $tline_hex_default)
      .style('stroke-dasharray', '10,3,3,3')
\"), collapse='')

# Adds some JS to be run after building the plot to get the image we want
d1\$setTemplate(afterScript = paste0(\"
  <script>
    myChart.draw()
    
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return (d['aggField'][3] < $t_strict); })
      .style('fill', function(d) {
        var datatype = d['aggField'][4];
        // http://pinetools.com/darken-color to darken hex values for chromatin states
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#FF0000';
          case 'TssAFlnk':
            return '#FF4500';
          case 'TxFlnk':
            return '#32CD32';
          case 'Tx':
            return '#008000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#C2E105';
          case 'Enh':
            return '#E5E500';
          case 'ZNF-Rpts':
            return '#3D7B66';
          case 'Het':
            return '#8A91D0';
          case 'TssBiv':
            return '#CD5C5C';
          case 'BivFlnk':
            return '#E9967A';
          case 'EnhBiv':
            return '#BDB76B';
          case 'ReprPC':
            return '#808080';
          case 'ReprPCWk':
            return '#C0C0C0';
          case 'Quies':
            return '#F8F8F8';
          default:
            return '#ff0000';
        }
      })
      .style('stroke', function(d) {
        var datatype = d['aggField'][4];
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#7f0000';
          case 'TssAFlnk':
            return '#7f2200';
          case 'TxFlnk':
            return '#186618';
          case 'Tx':
            return '#004000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#003200';
          case 'Enh':
            return '#727200';
          case 'ZNF-Rpts':
            return '#1e3d32';
          case 'Het':
            return '#31387b';
          case 'TssBiv':
            return '#712222';
          case 'BivFlnk':
            return '#983919';
          case 'EnhBiv':
            return '#66622d';
          case 'ReprPC':
            return '#404040';
          case 'ReprPCWk':
            return '#606060';
          case 'Quies':
            return '#7c7c7c';
          default:
            return '#ff0000';
        }
      });
      
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return ((d['aggField'][3] >= $t_strict) && (d['aggField'][3] < $t_marginal)); })
      .style('stroke', function(d) {
        var datatype = d['aggField'][4];
        switch (datatype) {
          case 'H3K27me3':
            return '#e5e500';
          case 'H3K4me1':
            return '#e41a1c';
          case 'H3K4me3':
            return '#4daf4a';
          case 'H3K36me3':
            return '#984ea3';
          case 'H3K9me3':
            return '#ff7f00';
          case 'TssA':
            return '#FF0000';
          case 'TssAFlnk':
            return '#FF4500';
          case 'TxFlnk':
            return '#32CD32';
          case 'Tx':
            return '#008000';
          case 'TxWk':
            return '#006400';
          case 'EnhG':
            return '#C2E105';
          case 'Enh':
            return '#E5E500';
          case 'ZNF-Rpts':
            return '#3D7B66';
          case 'Het':
            return '#8A91D0';
          case 'TssBiv':
            return '#CD5C5C';
          case 'BivFlnk':
            return '#E9967A';
          case 'EnhBiv':
            return '#BDB76B';
          case 'ReprPC':
            return '#808080';
          case 'ReprPCWk':
            return '#C0C0C0';
          case 'Quies':
            return '#F8F8F8';
          default:
            return '#FF82ab';
        }
      })
      .style('stroke-width', '2')
      .style('fill', 'white');
      
    myChart.svg
      .selectAll('.dimple-bubble')
      .filter(function(d) { return (d['aggField'][3] >= $t_marginal); })
      .style('stroke', 'rgb(173, 216, 230)')
      .style('fill', 'white');

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
      })
      .style('text-anchor', 'start');

//       .attr('transform', function() {
//         var t = d3.transform(d3.select(this).attr('transform'));
//         return 'rotate(-90) translate(' + -1*parseFloat(t.translate[0]) + ', ' + $xAxisLabelShift*parseFloat(t.translate[1]) + ')';
//       });

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
      .text('DMPs in DNase I sites (probably TF sites) in cell lines for $data $label')
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
