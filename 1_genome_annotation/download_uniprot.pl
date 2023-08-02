use strict;
use warnings;
use LWP::UserAgent;

my $list = $ARGV[0]; # File containg list of UniProt identifiers.

my $base = 'https://www.uniprot.org';
my $tool = 'uploadlists';

my $contact = 'jgoodheart@ucsd.edu'; # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
push @{$agent->requests_redirectable}, 'POST';

my $response = $agent->post("$base/$tool/",
                            [ 'file' => [$list],
                              'format' => 'tab',
                              'from' => 'ACC+ID',
                              'to' => 'ACC',
                              # If you want tab-separated output instead of full entries, specify 'tab' instead of 'txt' above
                              # and include the following optional line to include your preferred set of columns
                              # instead of the default from-to output:
                              'columns' => 'id,entry name,reviewed,protein names,genes,organism,length,database(RefSeq),go',

                            ],
                            'Content_Type' => 'form-data');

while (my $wait = $response->header('Retry-After')) {
  print STDERR "Waiting ($wait)...\n";
  sleep $wait;
  $response = $agent->get($response->base);
}

$response->is_success ?
  print $response->content :
  die 'Failed, got ' . $response->status_line .
    ' for ' . $response->request->uri . "\n";
