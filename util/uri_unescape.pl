#!/usr/bin/env perl

use strict;
use warnings;
use URI::Escape;

while (<>) {
	print uri_unescape($_);
}

exit(0);

