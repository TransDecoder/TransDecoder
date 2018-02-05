package PWM;

use strict;
use warnings;
use Carp;

my $PSEUDOCOUNT = 0.1;

srand(1); # reproducibility

###
sub new {
    my $packagename = shift;
    
    my $self = { pos_freqs => [],
                 pos_probs => [],
                 _pwm_built_flag => 0,
    };

    bless($self, $packagename);

    return($self);
}


sub is_pwm_built {
    my ($self) = @_;
    
    return($self->{_pwm_built_flag});

}

sub has_pos_freqs {
    my ($self) = @_;

    if (@{$self->{pos_freqs}}) {
        return(1);
    }
    else {
        return(0);
    }
}



sub add_feature_seq_to_pwm {
    my ($self, $feature_seq) = @_;

    $feature_seq = uc $feature_seq;
    
    my $pwm_len = $self->get_pwm_length();
    if ($pwm_len && length($feature_seq) != $pwm_len) {
        confess "Error, pwm_len: $pwm_len and feature_seq len: " . length($feature_seq) . " are unmatched.";
    }
    if ($feature_seq =~ /[^GATC]/) {
        print STDERR "Error, feature_seq: $feature_seq contains non-GATC chars... skipping\n";
        return;
    }
    
    my @chars = split(//, $feature_seq);
    for (my $i = 0; $i <= $#chars; $i++) {
        my $char = $chars[$i];

        $self->{pos_freqs}->[$i]->{$char}++;
    }

    return;
}


sub remove_feature_seq_from_pwm {
    my ($self, $feature_seq) = @_;
    
    if (! $self->has_pos_freqs()) {
        confess "Error, pwm obj doesn't have position frequencies set";
    }
    
    $feature_seq = uc $feature_seq;
    
    my $pwm_len = $self->get_pwm_length();
    if ($pwm_len && length($feature_seq) != $pwm_len) {
        confess "Error, pwm_len: $pwm_len and feature_seq len: " . length($feature_seq) . " are unmatched.";
    }
    if ($feature_seq =~ /[^GATC]/) {
        print STDERR "Warning, feature_seq: $feature_seq contains non-GATC chars... skipping.\n";
        return;
    }
    
    my @chars = split(//, $feature_seq);
    for (my $i = 0; $i <= $#chars; $i++) {
        my $char = $chars[$i];
        $self->{pos_freqs}->[$i]->{$char}--;
    }
    
    return;
}


    
sub get_pwm_length {
    my ($self) = @_;

    my $pwm_length;
    
    if ($self->is_pwm_built()) {
        $pwm_length = $#{$self->{pos_probs}} + 1;
    }
    else {
        $pwm_length = $#{$self->{pos_freqs}} + 1;
    }
    
    return($pwm_length);
}


sub simulate_feature {
    my ($self) = @_;

    if (! $self->is_pwm_built()) {
        confess "Error, pwm not built yet";
    }

    my $probs_aref = $self->{pos_probs};

    my @chars = qw(G A T C);
    my $feature_seq = "";

    my $pwm_length = $self->get_pwm_length();
    for (my $i = 0; $i < $pwm_length; $i++) {
        my $sum_prob = 0;
        my $selected_flag = 0;
        #print "rand val: $rand_val\n";
        my $tot_prob = 0;
        for (my $j = 0; $j <= $#chars; $j++) {
            my $char = $chars[$j];
            my $p = $probs_aref->[$i]->{$char};
            $tot_prob += $p;
        }
        my $rand_val = rand($tot_prob); # tot_prob should be 1 or very close...
        
        
        for (my $j = 0; $j <= $#chars; $j++) {
            my $char = $chars[$j];
            my $p = $probs_aref->[$i]->{$char};
            #print "char: $char, p: $p\n";
            $sum_prob += $p;
            if ($j == $#chars) {
                $sum_prob = 1;
            }
            
            if ($rand_val <= $sum_prob) {
                # choose char
                $feature_seq .= $char;
                $selected_flag = 1;
                last;
            }
        }
        if (! $selected_flag) {
            confess "Error, didn't select a random char for feature seq";
        }
    }

    return($feature_seq);

}



####
sub write_pwm_file {
    my ($self, $filename) = @_;

    unless ($self->is_pwm_built()) {
        confess "Error, pwm needs to be built first before writing to file";
    }
    
    open (my $ofh, ">$filename") or confess "Error, cannot write to file: $filename";

    my $pos_probs_aref = $self->{pos_probs};
    my $pwm_length = $self->get_pwm_length();
    
    my @chars = qw(G A T C);
    print $ofh join("\t", "pos", @chars) . "\n";
    for (my $i = 0; $i < $pwm_length; $i++) {
        
        my @vals = ($i);
        foreach my $char (@chars) {
            my $prob = $pos_probs_aref->[$i]->{$char} || 0;
            push (@vals, $prob);
        }
                
        print $ofh join("\t", @vals) . "\n";
    }
    
    close $ofh;
    
    return;
}



####
sub build_pwm {
    my ($self) = @_;

    my $pos_freqs_aref = $self->{pos_freqs};

    for (my $i = 0; $i <= $#{$pos_freqs_aref}; $i++) {
        my @chars = keys %{$pos_freqs_aref->[$i]};

        my $sum = 0;
        foreach my $char (@chars) {
            my $val = $pos_freqs_aref->[$i]->{$char};
            $sum += $val;
        }
        
        # now convert to relative freqs
        @chars = qw(G A T C); # the complete set of chars we care about.
        foreach my $char (@chars) {
            my $val = $pos_freqs_aref->[$i]->{$char} || 0;
            my $prob = sprintf("%.6f", ($val + $PSEUDOCOUNT) / ($sum + 4 * $PSEUDOCOUNT) );
            $self->{pos_probs}->[$i]->{$char} = $prob;
        }
        
    }

    $self->{_pwm_built_flag} = 1;

    return;
}


sub score_pwm_using_base_freqs {
    my ($self, $target_sequence, $base_freqs_href, %options ) = @_;

    ###  Options can include:
    ###
    ###   mask => [ coordA, coordB, coordC ],  # ignored pwm positions in scoring
    ###   pwm_range => [pwm_idx_start, pwm_idx_end], # start and end are inclusive
    ###
    
    for my $key (keys %options) {
        if (! grep {$_ eq $key} ('mask', 'pwm_range')) {
            confess "Error, option $key is not recognized";
        }
    }
        
    #print STDERR "target_seq: [$target_sequence]\n";
    
    $target_sequence = uc $target_sequence;
    unless ($target_sequence =~ /^[GATC]+$/) {
        # can only score GATC-containing sequences.
        return("NA");
    }
    
    if (! $self->is_pwm_built()) {
        confess("pwm not built yet!");
    }
    
    my $pwm_length = $self->get_pwm_length();

    if (length($target_sequence) != $pwm_length) {
        confess "Error, len(target_sequence)=" . length($target_sequence) . " and pwm length = $pwm_length";
    }

    my %mask;
    if (my $mask_positions_aref = $options{'mask'}) {
        %mask = map { + $_ => 1 } @$mask_positions_aref;
    }
    
    my $motif_score = 0;

    my @seqarray = split(//, $target_sequence);

    my ($pwm_start, $pwm_end) = (0, $pwm_length-1);
    if (my $pwm_range_aref = $options{'pwm_range'}) {
        ($pwm_start, $pwm_end) = @$pwm_range_aref;
    }
        
    for (my $i = $pwm_start; $i <= $pwm_end; $i++) {
        
        if ($mask{$i}) {
            next;
        }
        
        my $char = $seqarray[$i];
        my $prob = $self->{pos_probs}->[$i]->{$char};

        unless ($prob) {
            return("NA");
        }
        
        my $prob_rand = $base_freqs_href->{$char};
        unless ($prob_rand) {
            die "Error, no non-zero probability value specified for char [$char] ";
        }
        
        my $loglikelihood = log($prob/$prob_rand);
        $motif_score += $loglikelihood;

    }
    
    return($motif_score);

}


sub score_plus_minus_pwm {
    my ($self, $target_sequence, $pwm_minus, %options) = @_;
    
    ###  Options can include:
    ###
    ###   mask => [ coordA, coordB, coordC ],  # ignored pwm positions in scoring
    ###   pwm_range => [pwm_idx_start, pwm_idx_end], # start and end are inclusive
    ###
    
    for my $key (keys %options) {
        if (! grep {$_ eq $key} ('mask', 'pwm_range')) {
            confess "Error, option $key is not recognized";
        }
    }
        
    #print STDERR "target_seq: [$target_sequence]\n";
    
    $target_sequence = uc $target_sequence;
    unless ($target_sequence =~ /^[GATC]+$/) {
        # can only score GATC-containing sequences.
        return("NA");
    }
    
    if (! $self->is_pwm_built()) {
        confess("pwm not built yet!");
    }
    
    my $pwm_length = $self->get_pwm_length();

    if (length($target_sequence) != $pwm_length) {
        confess "Error, len(target_sequence)=" . length($target_sequence) . " and pwm length = $pwm_length";
    }

    my %mask;
    if (my $mask_positions_aref = $options{'mask'}) {
        %mask = map { + $_ => 1 } @$mask_positions_aref;
    }
    
    my $motif_score = 0;

    my @seqarray = split(//, $target_sequence);

    my ($pwm_start, $pwm_end) = (0, $pwm_length-1);
    if (my $pwm_range_aref = $options{'pwm_range'}) {
        ($pwm_start, $pwm_end) = @$pwm_range_aref;
    }
        
    for (my $i = $pwm_start; $i <= $pwm_end; $i++) {
        
        if ($mask{$i}) {
            #print STDERR "masking $i\n";
            next;
        }
        
        my $char = $seqarray[$i];
        my $prob = $self->{pos_probs}->[$i]->{$char};

        unless ($prob) {
            return("NA");
        }
        
        my $prob_rand = $pwm_minus->{pos_probs}->[$i]->{$char};
        unless ($prob_rand) {
            die "Error, no non-zero probability value specified for char [$char] ";
        }
        
        my $loglikelihood = log($prob/$prob_rand);
        $motif_score += $loglikelihood;
        
    }
    
    return($motif_score);

}

####
sub load_pwm_from_file {
    my ($self, $pwm_file) = @_;

    open(my $fh, $pwm_file) or confess "Error, cannot open file: $pwm_file";

    my $header = <$fh>;
    chomp $header;
    my @chars = split(/\t/, $header);
    shift @chars;
    
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $idx_pos = shift @x;
        for (my $i = 0; $i <= $#chars; $i++) {
            my $char = $chars[$i];
            my $prob = $x[$i];

            $self->{pos_probs}->[$idx_pos]->{$char} = $prob;
        }
    }
    close $fh;

    $self->{_pwm_built_flag} = 1;

    return;
}


1; #EOM
