package eForge::eStats;

use 5.010;
use strict;
use warnings FATAL => 'all';
use Math::BigInt;
use Math::BigFloat;

=head1 NAME

eStats - Stats for use in eForge

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

our (@ISA, @EXPORT);
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(mean variance std log10 binomial factorial fdr); # Symbols to export by default



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

=head2 binomial

Exact solution of binomial probability for k picks out of n, for n or greater need to sum for each k up to n

=cut

sub binomial {

    my ($k, $n, $p) = @_;
    my $prob = Math::BigFloat->new(($p**$k) * ((1 - $p)**($n - $k))) * factorial($n) / (factorial($k) * factorial($n - $k));
    $prob = sprintf("%.30f", $prob);
    return $prob;
  }

=head2 factorial

Calculate N!. Required for binomial.

Uses a cache to speed up the calculation.

=cut

my %_factorial_cache;

sub factorial{
    my ($n) = shift;
    return 1 if($n <=1 );
    return $_factorial_cache{$n} if (exists($_factorial_cache{$n}));
    Math::BigInt->new($n);
    $_factorial_cache{$n} = Math::BigInt->bfac($n);
    return $_factorial_cache{$n};
  }

=head2 fdr

Empirical false discovery rate = FP/TP+FP.

Need to modify this now have switched to binomial p values.

=cut


sub fdr{
    my ($tp, $mvps, $cells) = @_;
    if ($tp == 0){
        return "NA";
      }
    else{
        my $fpr = 0.0085 * exp(-0.04201 * $mvps) + 0.00187; # from simulations of random data (Forge, not eForge->) 0.0085*exp(-0.04201. SNPs) + 0.00187
        my $fdr = ($cells * $fpr) / $tp;
        return $fdr;
      }
  }





=head1 AUTHOR

Charles Breeze, C<< <cbreeze at ebi.ac.uk> >>
Charles Breeze, C<< <c.breeze at ucl.ac.uk> >>

=head1 BUGS
#not yet implemented below, it refers to eForge but it is not yet implemented (for cpan) even for Forge
Please report any bugs or feature requests to C<bug-eforge at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=eForge>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Forge::eStats



=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

eForge::eStats.pm

Copyright (C) 2014  EMBL - European Bioinformatics Institute

This program is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version. This program is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. Neither
the institution name nor the name eforge.pl can be used to endorse or promote products derived from
this software without prior written permission. For written permission, please contact
cbreeze@ebi.ac.uk or c.breeze@ucl.ac.uk . Products derived from this software may not be called eforge.pl nor may eforge.pl
appear in their names without prior written permission of the developers. You should have received
a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses/.

=cut

1; # End of eStats
