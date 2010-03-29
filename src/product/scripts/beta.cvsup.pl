#!/usr/bin/env perl

#	$Id: beta.cvsup.pl,v 1.2 2010/03/19 16:12:01 hiram Exp $
#
# beta.cvsup.pl - do a cvs update on -rbeta branch, report results concisely
#

use warnings;
use strict;

sub usage() {
    printf STDERR "beta.cvsup.pl - run the cvs command:\n";
    printf STDERR "    cvs up -rbeta -d -P <source tree directory>\n";
    printf STDERR "on the specified <source tree directory>\n";
    printf STDERR "usage: beta.cvsup.pl <source tree directory>\n";
    printf STDERR "where <source tree directory> is usually simply .\n";
    printf STDERR "to update the current directory.  The output to stdout\n";
    printf STDERR "is a convenient summary of what changed in the source\n";
    printf STDERR "tree during the update.\n";
    exit 255;
}

my $argCount = scalar(@ARGV);
if ($argCount < 1) {
    usage;
}

# verify we have a valid CVSROOT environment
my $cvsRoot = $ENV{'CVSROOT'};
if ( length($cvsRoot) < 3) {
    printf STDERR "beta.cvsup.pl ERROR: can not find CVSROOT in shell ";
    printf STDERR "environment.\nPlease set a CVSROOT in your environment.\n";
    usage;
}

my $cvsdir = $ARGV[0] || ".";
my $cvscmd = "cvs up -rbeta -d -P $cvsdir";

#
# isInteresting: return true if it looks like source that should be checked in.
#
my @interestingFiles = ('makefile', 'Makefile');
my @interestingSuffixes = ('c', 'h', 'as', 'sql', 'html', 'ra');
sub isInteresting {
    my $fname = shift(@_);
    foreach my $i (@interestingFiles) {
	if ($fname eq $i) {
	    return 1;
	}
    }
    if ($fname =~ m/\.(.+)$/) {
	my $suffix = $1;
	foreach my $i (@interestingSuffixes) {
	    if ($suffix eq $i) {
		return 1;
	    }
	}
    }
    return 0;
}

open(CVSUP, "$cvscmd 2>&1 |") ||
  die "Couldn't open pipe from \"$cvscmd\": $!\n";

my $date = `date`;
chop($date);
printf "in directory: %s", `pwd`;
print "Updating $cvsdir at $date\n\n";

# Pick out the interesting lines of output
my @additions;
my @updates;
my @modifieds;
my @removals;
my @patches;
my @conflicts;
my @removed;
my @noncvsInteresting;
my @noncvs;
my @unrecognized;
my $aborted;
while (<CVSUP>) {
  if      (/^A (\S+)/) {
    push @additions, $1;
  } elsif (/^U (\S+)/) {
    push @updates, $1;
  } elsif (/^M (\S+)/) {
    push @modifieds, $1;
  } elsif (/^R (\S+)/) {
    push @removals, $1;
  } elsif (/^P (\S+)/) {
    push @patches, $1;
  } elsif (/^C (\S+)/) {
    push @conflicts, $1;
  } elsif (/^\? (\S+)/) {
    my $f = $1;
    if (isInteresting($f)) {
      push @noncvsInteresting, $f;
    } else {
      push @noncvs, $f;
    }
  } elsif (/^(\S \S+)/) {
    push @unrecognized, $1;
  } elsif (/^cvs update: (\S+) is no longer in the repository/) {
    push @removed, $1;
  } elsif (/^cvs (update|server): Updating/ ||
	   /^RCS file: / ||
	   /^retrieving revision / ||
	   /^rcsmerge: warning: conflicts during merge/ ||
	   /^cvs update: conflicts found in/ ||
	   /^Merging differences between/) {
    next;
  } else {
    print $_;
    if (/^cvs \[update aborted\]:/) {
      $aborted = $_;
    }
  }
}

# Report them.
if (defined $aborted) {
  print "\nABORTED:\n$aborted\n";
}
if (scalar(@unrecognized) > 0) {
  print "Unrecognized update types:\n";
  foreach my $f (@unrecognized) {
    print "  $f\n";
  }
  print "\n";
}
if (scalar(@conflicts) > 0) {
  print "CONFLICTS -- manually edit to resolve:\n";
  foreach my $f (@conflicts) {
    print "  $f\n";
  }
  print "\n";
}
if (scalar(@additions) > 0) {
  print "Added (not checked in) files:\n";
  foreach my $f (@additions) {
    print "  $f\n";
  }
  print "\n";
}
if (scalar(@removals) > 0) {
  print "Removed (not checked in) files:\n";
  foreach my $f (@removals) {
    print "  $f\n";
  }
  print "\n";
}
if (scalar(@modifieds) > 0) {
  print "Modified (possibly merged) files:\n";
  foreach my $f (@modifieds) {
    print "  $f\n";
  }
  print "\n";
}
if (scalar(@patches) > 0) {
  print "Patched files:\n";
  foreach my $f (@patches) {
    print "  $f\n";
  }
  print "\n";
}
if (scalar(@updates) > 0) {
  print "Updated files:\n";
  foreach my $f (@updates) {
    print "  $f\n";
  }
  print "\n";
}
if (scalar(@removed) > 0) {
  print "Files removed by others:\n";
  foreach my $f (@removed) {
    print "  $f\n";
  }
  print "\n";
}
if (scalar(@noncvsInteresting) > 0) {
  print "Files/directories not checked in to CVS that look like source:\n";
  foreach my $f (@noncvsInteresting) {
    print "  $f\n";
  }
  print "\n";
}
if (scalar(@noncvs) > 0) {
  print "Files/directories not checked in to CVS that don't look like source:\n";
  foreach my $f (@noncvs) {
    print "  $f\n";
  }
  print "\n";
}
if (scalar(@conflicts) > 0) {
  print "CONFLICTS -- manually edit to resolve:\n";
  foreach my $f (@conflicts) {
    print "  $f\n";
  }
  print "\n";
}
if (defined $aborted) {
  print "ABORTED:\n$aborted\n";
}

close(CVSUP);
# Unfortunately $? actually doesn't pass cvs up's exit status... maybe because 
# of sh's handling of exit statuses in pipes??  Oh well, hopefully error 
# messages will be clear enough.  
exit $?;