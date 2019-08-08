package eForge::eStats;

use Data::Dumper;

=head1 NAME

eForge::eStats - Stats for use in eForge

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

our $VERSION = '0.01';

our (@ISA, @EXPORT);
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(mean variance std zscore log10 binomial factorial fdr); # Symbols to export by default



=head1 SYNOPSIS

Provide various stats for eForge to do its stuff


=head1 EXPORT

mean
variance
std
log10
binomial
factorial
fdr

=head1 SUBROUTINES/METHODS

=head2 mean

Calculates the biased mean of an array

pass it a float array and it will return the mean
reused from Ben Brown

=cut


sub mean {
    if (($#_+1) == 0) {
      print Dumper "divide by zero!";
      print Dumper \$_;
    }
    my $sum = 0;
    foreach (@_){
        $sum+= $_;
      }
    return $sum/($#_+1);
  }

=head2 variance

Calculates the biased variance of an array

Pass it a float array and it will return the variance

Reused from Ben Brown

=cut

sub variance {
    my $ev = mean(@_);
    my $sum = 0;
    foreach (@_) { $sum += ($_ - $ev)**2 };

    return $sum/($#_+1);
  }

=head2 std

Calulates the standard deviation of an array: this is just the sqrt of the var

=cut

sub std { sqrt(variance(@_)) }

=head2 log10

log 10 since perl doesn't have

=cut

sub log10 {
    my $n = shift;
    return log($n)/log(10);
  }

=head2 zscore

Calculates the z-score for a given result and an array of values

=cut

sub zscore {
    my ($teststat, $values) = @_;
    my $zscore;
    
    if (! @$values) {
      print Dumper "undefined values!";
    }
    
    my $mean = mean(@$values);
    my $sd = std(@$values);
    if ($sd == 0) {
        $zscore = "NA";
    } else {
        $zscore = sprintf("%.3f", ($teststat-$mean)/$sd);
    }

    return $zscore;
}

=head2 binomial

Exact solution of binomial probability for k picks out of n, for n or greater need to sum for each k up to n

=cut

sub binomial {

    my ($k, $n, $p) = @_;

    my $prob = exp($k*log($p) + ($n-$k)*log(1 - $p) + log_factorial($n) - log_factorial($k) - log_factorial($n - $k));

    return $prob;
  }

=head2 log_factorial

Calculate log(N!). Required for binomial.

Uses a cache to speed up the calculation.

=cut

my %_log_factorial_cache;

sub log_factorial{
    my ($n) = shift;
    return 0 if($n <=1 ); # log(1) = 0
    return $_log_factorial_cache{$n} if (exists($_log_factorial_cache{$n}));
    my $result = 0; # log(1) = 0;
    for (my $i = $n; $i > 1; $i--) {
        $result += log($i);
    }
    $_log_factorial_cache{$n} = $result;
    return $_log_factorial_cache{$n};
  }

=head2 fdr

Empirical false discovery rate = FP/TP+FP.

Need to modify this now have switched to binomial p values.

=cut


sub fdr{
    my ($tp, $dmps, $cells) = @_;
    if ($tp == 0){
        return "NA";
      }
    else{
        my $fpr = 0.0085 * exp(-0.04201 * $dmps) + 0.00187; # from simulations of random data (Forge, not eForge->) 0.0085*exp(-0.04201. SNPs) + 0.00187
        my $fdr = ($cells * $fpr) / $tp;
        return $fdr;
      }
  }

1;
