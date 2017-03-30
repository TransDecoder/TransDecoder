package PWM;

use strict;
use warnings;
use Carp;


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


sub add_feature_seq_to_pwm {
    my ($self, $feature_seq) = @_;

    my @chars = split(//, $feature_seq);
    for (my $i = 0; $i <= $#chars; $i++) {
        my $char = $chars[$i];

        $self->{pos_freqs}->[$i]->{$char}++;
    }

    return;
}


####
sub write_pwm_file {
    my ($self, $filename) = @_;

    unless ($self->is_pwm_built()) {
        croak "Error, pwm needs to be built first before writing to file";
    }
    
    open (my $ofh, ">$filename") or croak "Error, cannot write to file: $filename";

    my $pos_probs_aref = $self->{pos_probs};
    my $pwm_length = $#{$pos_probs_aref} + 1;
    
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
        foreach my $char (@chars) {
            my $val = $pos_freqs_aref->[$i]->{$char};
            my $prob = sprintf("%.3f", $val / $sum);
            $self->{pos_probs}->[$i]->{$char} = $prob;
        }
        
    }

    $self->{_pwm_built_flag} = 1;

    return;
}


1; #EOM
