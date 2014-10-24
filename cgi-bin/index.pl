#! /usr/bin/env perl
use warnings;
use strict;

use CGI ':form', ':cgi', ':html';                             # load CGI routines
use lib("modules");
use IPC::Cmd qw[run];
use LWP::Simple qw(!head);
use Time::HiRes qw(gettimeofday tv_interval);
use DBI;
use Data::UUID;

use Template;

my $debug = 0;


##################################################################################################
##################################################################################################
##
##   W E B   S E R V E R   C O N F I G U R A T I O N
##
##################################################################################################
##################################################################################################
##

my $WEB_ROOT_OUTDIR = "/files";
my $LOG_FILE = "../log/server_log.db";
my $BIN_DIR = "../bin";
my $INPUT_DATAFILE = "input.txt";
my $STDOUT_FILE = "output.txt";

my $colour="bright-blue";
my $plot_colour="#29A6C9";

my $title = "eFORGE";

my $breadcrumbs = [
    {"UCL Cancer Institute" => "http://www.ucl.ac.uk/cancer"},
    {"Cancer Biology" => "http://www.ucl.ac.uk/cancer/rescancerbiol/cancerbiology"},
    {"Medical Genomics" => "http://www.ucl.ac.uk/cancer/medical-genomics/medgenhome"},
];

my $left_menu = [
    {"__logo__" => "/logo.jpg"},
    {"__title__" => "eFORGE"},
    {"Start" => "/"},
    {"Help" => "/?help"},
    {"Download" => "/?download"},
    {"About" => "/?about"},
    {"__title__" => "UCL Cancer Institute"},
    {"Home" => "http://www.ucl.ac.uk/cancer/"},
    {"Medical Genomics" => "http://www.ucl.ac.uk/cancer/medical-genomics/medgenhome"},
    {"Bill Lyons Informatics Centre" => "http://www.ucl.ac.uk/cancer/blic/"},
];

my $right_column = undef;

my $REFRESH_TIME = 10; # in seconds

##
##################################################################################################
##################################################################################################

my $_OUT_DIR; # Private global. Do not use outside of get_web_outdir()!!!

my $q = CGI->new;

if (param("keywords") and param("keywords") eq "help") {
    print_help_page();

} elsif (param("keywords") and param("keywords") eq "download") {
    print_download_page();

} elsif (param("keywords") and param("keywords") eq "about") {
    print_about_page();

} elsif (param("action") and param("action") eq "Run") {
    print_run_page();

} else {
    print_main_page()

}

exit(0);


=head2 get_web_outdir

 Arg[1]         : -none-
 Example        : my $web_outdir = get_web_outdir();
 Description    : Gets the location of the output directory w.r.t. web server document root, i.e.
                  the WEB_OUTDIR portion in http://mytool.cs.ucl.ac.uk/WEB_OUTDIR.
                  If the directory does not exist yet, it will create one using a UUID.
 Returns        : string $web_outdir
 Exceptions     : dies if it fails to create the new directory when needed.

=cut

sub get_web_outdir {
    if (!$_OUT_DIR) {
        my $ug = new Data::UUID;
        $_OUT_DIR = $WEB_ROOT_OUTDIR."/".$ug->to_hexstring($ug->create());
        mkdir($ENV{'DOCUMENT_ROOT'}.$_OUT_DIR) or die "Cannot create output directory for the run";
    }
    return $_OUT_DIR;
}


=head2 get_absolute_outdir

 Arg[1]         : -none-
 Example        : my $absolute_outdir = get_absolute_outdir();
 Description    : Gets the location of the output directory in the file system. For instance, it
                  will be something like /var/www/htdocs/files/0xG4214242AD2EEC1CC56354/
 Returns        : string $absolute_outdir
 Exceptions     : None

=cut

sub get_absolute_outdir {
    my $web_outdir = get_web_outdir();

    return $ENV{'DOCUMENT_ROOT'}.$web_outdir;
}

=head2 get_absolute_root_outdir

 Arg[1]         : -none-
 Example        : my $absolute_root_outdir = get_absolute_root_outdir();
 Description    : Gets the location of the root output directory in the file system. For instance,
                  it will be something like /var/www/htdocs/files
 Returns        : string $absolute_root_outdir
 Exceptions     : None

=cut

sub get_absolute_root_outdir {
    my $web_outdir = get_web_outdir();

    return $ENV{'DOCUMENT_ROOT'}.$WEB_ROOT_OUTDIR;
}

=head2 print_form

 Arg[1]         : -none-
 Example        : print_form();
 Description    : Prints the form on the main page
 Returns        : 
 Exceptions     : None

=cut

sub print_form {
    print start_multipart_form(-id=>'eforge');
    print table({-width=>"100%", -border=>"0", -cellspacing=>"0", -cellpadding=>"0"},
        Tr({-valign=>"TOP"}, [

            ## ================================================================
            ## Input data section
            ## ================================================================
            th(["<h6>Data</h6>", ""]),
            td(["Paste data:",
                textarea('data_text', '# Example with a filtered set of monocyte tDMPs from Jaffe AE and Irizarry RA, Genome Biol 2014, 15:R31.
cg13430807
cg10480329
cg06297318
cg19301114
cg23244761
cg26872907
cg18066690
cg04468741
cg16636767
cg10624395
cg20918393', 10, 60)]),
            td(["",
                "<a href=\"javascript: void(0);\" onClick=\"document.getElementById('eforge').".
                    "data_text.value ='';\">Clear box</a>"]),
            td(["Upload file:",
                filefield('data_file','starting value', 58)]),
            td(["Provide file URL:",
                textfield('data_url', '', 58)]),
            td(["Input file format:",
                popup_menu('data_format', ['probeid', 'bed'], 'probeid', {'probeid'=>'Probe list',
                    'bed'=>'BED file'})]),

            th({-colspan=>2}, ["<hr>"]),

            ## ================================================================
            ## Option section
            ## ================================================================
            th(["Options", ""]),
            td(["Name for this data (optional):",
                textfield('label', '', 58)]),
            td(["Analysis data from:",
                radio_group('ref_data', ['erc'], 'erc', 'true', {'erc'=>' Epigenome Roadmap',
                    'encode'=>' Encode'} )]),
            td(["",
                radio_group('ref_data', ['encode'], 'erc', 'true', {'erc'=>' Epigenome Roadmap',
                    'encode'=>' Encode'} )]),
            td(["Depletion:",
                checkbox('depletion', '', 'on', '')]),
            td(["Proximity:",
                popup_menu('proxy', ['No filter', '1kb'], '1kb', {'1kb'=>'1 kb window',
                    'none'=>'No proxy'})]),
            td(["Background repetitions (100-1000):",
                textfield('reps', '1000', 10)]),
            td(["Significance threshold:",
                ""]),
            td(["&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Min:",
                textfield('thresh1', '0.01', 10)]),
            td(["&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Max:",
                textfield('thresh2', '0.05', 10)]),

            th({-colspan=>2}, ["<hr>"]),

            ## ================================================================
            ## Submit
            ## ================================================================
            td(["",
                submit('action','Run')]),
        ])
        );
    print end_form();
}


=head2 validate_form

 Arg[1]         : -none-
 Example        : validate_form();
 Description    : Validates the input from the form to ensure the options are valid.
                  Fills in the $input_data and @validate_args global variables.
 Returns        : (arrayref $validated_args, int $input_data_size)
 Exceptions     : prints full web page with header, error message, new form and footer on error.

=cut

sub validate_form {
    my $data;
    my $validated_args;
    my @error_messages;
    
    my $input_data;
    if (param("data_file")) {
        my $data_fh = upload('data_file');
        $input_data = join("", <$data_fh>);
        close($data_fh);

    } elsif (param("data_url")) {
        $input_data=get(param("data_url"));

    } elsif (param("data_text")) {
        $input_data = param("data_text");

    } else {
        push(@error_messages, "Please upload a file, provide a URL or paste your data in the".
            " text box.");
    }
    my @lines = grep {/\w/} split(/[\r\n]+/, $input_data);
    my $input_data_size = grep {!/^#/} @lines;
    if ($input_data_size > 1000) {
        push(@error_messages, "Maximum number of probes is 1000.");
    }

    my $data_format = param("data_format");
    if (!$data_format or ($data_format ne "probeid" and $data_format ne "bed")) {
        $data_format =~ s/\W/_/g;
        push(@error_messages, "Data format &quot;$data_format&quot; is not valid. Please".
            " specify a valid one.");
    }
    push(@$validated_args, "--format", param("data_format"));

    my $label = param("label");
    if ($label and $label =~ /\W/) {
        $label=~s/\W/_/g;
    }
    push(@$validated_args, "--label", $label) if ($label);

    my $ref_data = param("ref_data");
    if (!$ref_data or ($ref_data ne "erc" and $ref_data ne "encode")) {
        $ref_data =~ s/\W/_/g;
        push(@error_messages, "Unknown set of analysis data &quot;$ref_data&quot;. Please".
            " specify a valid one.");
    }
    push(@$validated_args, "--data", param("ref_data"));

    my $depletion = param("depletion");
    if ($depletion and $depletion eq 'on') {
        push(@$validated_args, "--depletion");
    }

    my $proxy = param("proxy");
    if (!$proxy or ($proxy ne "1kb" and $proxy ne "none")) {
        $proxy =~ s/\W/_/g;
        push(@error_messages, "Unknown proximity option &quot;$proxy&quot; is not valid.".
            " Please specify a valid one.");
    }
    if ($proxy eq "none") {
        push(@$validated_args, "--noproxy");
    } else {
        push(@$validated_args, "--proxy", "1kb");
    }

    my $reps = param("reps");
    if (!$reps or $reps !~ /^\d+$/) {
        push(@error_messages, "Background repetitions must be a positive number.");
    } elsif ($reps < 100) {
        push(@error_messages, "Background repetitions cannot be fewer than 100.");
    } elsif ($reps > 1000) {
        push(@error_messages, "Background repetitions cannot be more than 1000.");
    }
    push(@$validated_args, "--reps", $reps);

    my $thresh1 = param("thresh1");
    if (!$thresh1 or $thresh1 !~ /^\d+(\.\d*)?$/) {
        push(@error_messages, "Min. significance threshold must be a positive number.");
    } elsif ($thresh1 >= 1) {
        push(@error_messages, "Min. significance threshold must be less than 1.");
    }

    my $thresh2 = param("thresh2");
    if (!$thresh2 or $thresh2 !~ /^\d+(\.\d*)?$/) {
        push(@error_messages, "Max. significance threshold must be a positive number.");
    } elsif ($thresh2 >= 1) {
        push(@error_messages, "Max. significance threshold must be less than 1.");
    }

    if ($thresh1 > $thresh2) {
        print $q->header;
        push(@error_messages, "Min. significance threshold must be less than max. significance".
            " threshold.");
        return 0;
    }
    push(@$validated_args, "--thresh", "$thresh1,$thresh2");

    if (@error_messages) {
        print $q->header;
        print Template::start($title, $breadcrumbs, $left_menu, $colour, $right_column);
        print Template::header("Error");
        foreach my $this_error_message (@error_messages) {
            print Template::error_box($this_error_message);
        }
        print_form();
        print Template::content_box_1("Params", map {$_.": ".$q->param($_)} $q->param) if ($debug);
        print_credits_box();
        print Template::end;
        exit(0);
    }

    ## It seems like all the options are valid, so we can now store the input data in the output
    ## directory
    my $absolute_outdir = get_absolute_outdir();
    open(INPUT, ">$absolute_outdir/$INPUT_DATAFILE") or
        die "Cannot open $absolute_outdir/$INPUT_DATAFILE";
    foreach my $this_line (@lines) {
        print INPUT $this_line, "\n";
    }
#     print INPUT "# This file contains $input_data_size probes (before checking for duplicates)\n";
    close(INPUT);
    
    return ($validated_args, $input_data_size);
}


=head2 run_tool

 Arg[1]         : arrayref string $validated_args
 Arg[2]         : int $input_size (for logging purposes)
 Example        : run_tool($validated_args, 3212);
 Description    : Forks the process, closes STDIN, STDOUT and STDERR on the child so the parent can
                  finish and tell Apache that the page is complete and runs the tool in the
                  background.
                  When the process finishes, it overwrites the existing index.html waiting page with
                  the final result. It also forwards the STDOUT and STDIN to to the $STDOUT_FILE
                  in the data directory. That file is read by the AJAX magic to give some feedback
                  about the progress.
                  The "result" section of the page is also written to an output.html page. This is
                  used as a flag to tell the AJAX magic that the process has come to an end. The
                  waiting page will load that content and display it to the user.
                  Note that clients with JavaScript disabled still get some functionality through
                  automatic refreshes of the waiting page (until this is overwritten on success or
                  failure).
 Returns        : 
 Exceptions     : Prints an error message in an error_box if the fork fails.

=cut

sub run_tool {
    my ($validated_args, $input_size) = @_;

    if (my $pid = fork) { 
        # Parent, just continue...
        exit();
    } elsif (defined $pid) { 
        # Child. Run the tool
        open STDIN, "</dev/null";
        open STDOUT, ">/dev/null";
        open STDERR, ">/dev/null";

        my $absolute_outdir = get_absolute_outdir();

        ## ====================================================================
        ## Run the tool.
        ## Note that I am leaving the full input file name and the output
        ## directory out of the @$validated_args array to mask them from the
        ## command line shown to the user.
        ## ====================================================================        
        my $t0 = [gettimeofday];
        chdir($BIN_DIR);
        my ($ok, $err, $full_buff, $stdout_buff, $stderr_buff) = run(
            command => ["perl", "eforge.pl", "-f", "$absolute_outdir/$INPUT_DATAFILE", "-out_dir",
                $absolute_outdir, @$validated_args, ">", "$absolute_outdir/$STDOUT_FILE", "2>&1"]);
        my $t1 = [gettimeofday];
        ## ====================================================================        
        
        open(INDEX, ">$absolute_outdir/index.html") or
            die "Cannot open $absolute_outdir/index.html";
        open(OUTPUT, ">$absolute_outdir/output.html") or
            die "Cannot open $absolute_outdir/output.html";
        print INDEX Template::start($title, $breadcrumbs, $left_menu, $colour, $right_column);

        if (!$ok) {
            print INDEX Template::header("Error");
            print INDEX Template::error_box("Error while running eForge.pl: $err.",
                map({"ERR: $_"} @$stderr_buff),
                map({"OUT: $_"} @$stdout_buff),
                );

            print OUTPUT Template::header("Error");
            print OUTPUT Template::error_box("Error while running eForge.pl: $err.",
                map({"ERR: $_"} @$stderr_buff),
                map({"OUT: $_"} @$stdout_buff),
                );

        } else {
            my $web_outdir = get_web_outdir();
            my $url = $q->url(-base=>1)."$web_outdir/index.html";
            my @output = qx"cat $absolute_outdir/$STDOUT_FILE";
            print INDEX Template::header("eFORGE");
            print INDEX Template::content_box("Done.",
                "<strong>".join(" ", "perl", "eforge.pl", "-f", "input.txt", @$validated_args).
                    "</strong>",
                textarea('output', join("", @output), 10, 90),
                "URL: <a href=\"$url\">$url</a>");
            print_result(*INDEX);
            print_result(*OUTPUT);
        }

        print Template::content_box_1("Params", map {$_.": ".$q->param($_)} $q->param) if ($debug);

        print INDEX Template::end;
        close(INDEX);

        log_usage(tv_interval($t0, $t1), $input_size, join(" ", @$validated_args), $err);

    } else { 
        # Error while forking
        print Template::error_box("Error while attempting to fork the process.");
    }
    
}


=head2 print_result

 Arg[1]         : ref *filehandle
 Example        : print_result(*STDIN);
 Description    : Writes out the two result boxes, with the outptu files and the dynamic table
                  respectively.
 Returns        : 
 Exceptions     : Prints an error message in an error_box if it cannot find the table file.

=cut

sub print_result {
    my ($fh) = @_;
    
    my $absolute_outdir = get_absolute_outdir();
    my $web_outdir = get_web_outdir();

    opendir(DIR, $absolute_outdir);
    my @files = grep {/(.pdf|.html|.tsv|.R)$/} readdir(DIR);
    closedir(DIR);
    my $table_file = (grep {/.table.html$/} @files)[0];
    my $table_R = (grep {/.table.R$/i} @files)[0];
    my $dchart_file = (grep {/.dchart.html$/} @files)[0];
    my $dchart_R = (grep {/.dchart.R$/i} @files)[0];
    my $tsv_file = (grep {/.chart.tsv$/} @files)[0];
    my $pdf_file = (grep {/.chart.pdf$/} @files)[0];
    my $pdf_R = (grep {/.chart.R$/i} @files)[0];

    print $fh Template::content_box_1("Results",
        "<a href=\"$web_outdir/$INPUT_DATAFILE\">Input data (txt)</a>",
        "<a href=\"$web_outdir/$tsv_file\">Raw data (tsv)</a>",
        "<a href=\"$web_outdir/$pdf_file\">Static chart (PDF)</a>",
        "<a href=\"$web_outdir/$dchart_file\">Interactive chart (HTML)</a>",
        "<a href=\"$web_outdir/$table_file\">Interactive table (HTML)</a>",
        "R code for re-generating these files: [<a href=\"$web_outdir/$pdf_R\">PDF</a>]
        [<a href=\"$web_outdir/$dchart_R\">chart</a>] [<a href=\"$web_outdir/$table_R\">table</a>]",
        );
    if (-e "$absolute_outdir/$table_file") {
        my @table_file_content = qx"cat $absolute_outdir/$table_file";
#         open(JS, ">$absolute_outdir/table.js") or die;
#         my $print_mode = 0;
#         foreach my $this_line (@table_file_content) {
#             if ($this_line =~ /<script/ and $this_line !~ /<\/script/) {
#                 $print_mode = 1;
#             } elsif ($this_line =~ /<\/script/) {
#                 $print_mode = 0;
#             } elsif ($print_mode) {
#                 print JS $this_line;
#             }
#         }
#         close(JS);
        print $fh Template::content_box("Interactive table",
            "<noscript>This table is only visible if your browser supports JavaScript.</noscript>".
            join("", @table_file_content));
    } else {
        print $fh Template::error_box("Cannot find the table file...");
    }

}

=head2 print_run_page

 Arg[1]         : arrayref string $validated_args
 Arg[2]         : int $input_size (for logging purposes)
 Example        : print_waiting_page(*STDOUT, 10);
 Description    : Prints the waiting page. This is used to create the output of the CGI when running
                  the tool and to create a static waiting page as well. Both versions refresh the
                  page to the static page, i.e. the output of the CGI redirects to the static page
                  the first time and that one refreshes itself until it is rewritten by the child
                  process running the tool.
 Returns        : 
 Exceptions     : None

=cut

sub print_run_page {
    my ($validated_args, $input_size) = validate_form();

    my $absolute_outdir = get_absolute_outdir();

    open(INDEX, ">$absolute_outdir/index.html") or die;
    print_waiting_page(*INDEX, $validated_args);
    close(INDEX);
    
    print $q->header;
    print_waiting_page(*STDOUT, $validated_args);

    run_tool($validated_args, $input_size);
}


=head2 print_waiting_page

 Arg[1]         : $filehandle for printing
 Arg[2]         : arrayref string $validated_args
 Example        : print_waiting_page(*STDOUT, $validated_args);
 Description    : Prints the waiting page. This is used to create the output of the CGI when running
                  the tool and to create a static waiting page as well. If JavaScript is disabled,
                  both versions refresh the page to the static page, i.e. the output of the CGI
                  redirects to the static page the first time and that one refreshes itself until it
                  is rewritten by the child process running the tool (see run_tool()).
                  If JavaScript is enabled, the page checks for the progress on the server every
                  second and shows the output of the tool in a textbox. When the tool has finished
                  (i.e. when output.html exists in the directory), AJAX loads that content on the
                  page.
 Returns        : 
 Exceptions     : None

=cut

sub print_waiting_page {
    my ($fh, $validated_args) = @_;

    my $web_outdir = get_web_outdir();
    my $url = $q->url(-base=>1)."$web_outdir/index.html";

    my $start = Template::start($title, $breadcrumbs, $left_menu, $colour, $right_column);
    $start =~ s/(\s+<title>)/<noscript><META HTTP-EQUIV="refresh" CONTENT="$REFRESH_TIME;URL=$url" id="reload"><\/noscript>$1/i;
    print $fh $start;
    print $fh Template::header("eFORGE");
    print $fh Template::content_box("<div id=\"status_title\">Running...</div>",
        "<strong>".join(" ", "perl", "eforge.pl", "-f", "input.txt", @$validated_args).
            "</strong>",
        textarea(-name=>'output', -default=>'', -rows=>10, -columns=>90, -id=>'output'),
        "URL: <a href=\"$url\">$url</a>",
        "<noscript>Your browser does not support JavaScript. This page refreshes automatically".
        " every 10 sec until the result is ready.</noscript>");
    print $fh "<div id=\"result\"></div>";
    print $fh <<EOF;
<script type="text/javascript">

// From http://stackoverflow.com/questions/511088/use-javascript-to-place-cursor-at-end-of-text-in-text-input-element
(function(\$)
{
    jQuery.fn.putCursorAtEnd = function()
    {
    return this.each(function()
    {
        \$(this).focus()

        // If this function exists...
        if (this.setSelectionRange)
        {
        // ... then use it
        // (Doesn't work in IE)

        // Double the length because Opera is inconsistent about whether a carriage return is one character or two. Sigh.
        var len = \$(this).val().length * 2;
        this.setSelectionRange(len, len);
        }
        else
        {
        // ... otherwise replace the contents with itself
        // (Doesn't work in Google Chrome)
        \$(this).val(\$(this).val());
        }

        // Scroll to the bottom, in case we're in a tall textarea
        // (Necessary for Firefox and Google Chrome)
        this.scrollTop = 999999;
    });
    };
})(jQuery);

function loadXMLDoc()
{

var xmlhttp;
var xmlhttp2;
if (window.XMLHttpRequest)
  {// code for IE7+, Firefox, Chrome, Opera, Safari
  xmlhttp=new XMLHttpRequest();
  xmlhttp2=new XMLHttpRequest();
  }
else
  {// code for IE6, IE5
  xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
  xmlhttp2=new ActiveXObject("Microsoft.XMLHTTP");
  }
xmlhttp.onreadystatechange=function()
  {
  if (xmlhttp.readyState==4 && xmlhttp.status==200)
    {
    document.getElementById("output").innerHTML=xmlhttp.responseText;
    \$("#output").putCursorAtEnd();
    }
  }
xmlhttp2.onreadystatechange=function()
  {
  if (xmlhttp2.readyState==4 && xmlhttp2.status==200)
    {
    // document.getElementById("result").innerHTML=xmlhttp2.responseText;
    document.getElementById("status_title").innerHTML="Done.";
    \$("#result").load("$web_outdir/output.html")
    // \$.getScript("$web_outdir/table.js")
    clearInterval(myVar)
    location.replace("$web_outdir/index.html")
    }
  }
xmlhttp.open("GET","$web_outdir/output.txt",true);
xmlhttp.send();
xmlhttp2.open("GET","$web_outdir/output.html",true);
xmlhttp2.send();
}

loadXMLDoc();

var myVar = setInterval(function(){loadXMLDoc(); }, 1000);

//document.getElementById('reload').parentNode.removeChild(document.getElementById('reload'));

</script>

EOF

    print $fh Template::content_box_1("Params", map {$_.": ".$q->param($_)} $q->param) if ($debug);
    print $fh Template::end;

}


=head2 print_intro_box

 Arg[1]         : -none-
 Example        : print_intro_box();
 Description    : Prints the Intro box used in the main and in the help pages
 Returns        : 
 Exceptions     : None

=cut

sub print_intro_box {
    print Template::content_box_1("Description",
        "eFORGE is the epigenetic equivalent of
            <a href=\"http://www.1000genomes.org/forge-analysis\">FORGE</a>, using EWAS rather than
             GWAS data.
           <br \>",
        "eFORGE identifies tissue or cell type-specific signal by analysing a minimum set of 5
            differentially methylated positions (DMPs) for overlap with DNase 1 hypersensitive sites
            (DHSs) compared to matched background DMPs and provides both graphical and tabulated
            outputs.
           <br \>",
    );
}


=head2 print_credits_box

 Arg[1]         : -none-
 Example        : print_credits_box();
 Description    : Prints the Credits box
 Returns        : 
 Exceptions     : None

=cut

sub print_credits_box {
    print Template::content_box("Credits",
        "Breeze CE, Paul DS, Butcher LM, Herrero J, Birney E, Dunham I, Beck S. eFORGE: A tool for
           identifying tissue-specific signal in epigenomic data. <i>Manuscript in preparation.</i>
           <br \>",
    );
}


=head2 print_main_page

 Arg[1]         : -none-
 Example        : print_main_page();
 Description    : Prints the main page.
 Returns        : 
 Exceptions     : exit(0) always

=cut

sub print_main_page {
    print $q->header;
    print Template::start($title, $breadcrumbs, $left_menu, $colour, $right_column);
    print Template::header($title);

    print_intro_box();
    
    print_form();

    print_credits_box();

    print Template::content_box_1("Params", map {$_.": ".$q->param($_)} $q->param) if ($debug);

    print Template::end;
    exit(0);
}


=head2 print_help_page

 Arg[1]         : -none-
 Example        : print_help_page();
 Description    : Prints the Help page.
 Returns        : 
 Exceptions     : exit(0) always

=cut

sub print_help_page {
    print $q->header;
    print Template::start($title, $breadcrumbs, $left_menu, $colour, $right_column);
    print Template::header("eFORGE &gt; Help");

    print_intro_box();

    print Template::content_box_1("Input data",
       "<strong>Source</strong><br \><br \>
           You can upload a file, provide a URL for an existing file or simply
           copy &quot; paste your data on the text box provided.
           <br \>",
       "<strong>Input file format</strong><br \><br \>
           If f is specified, specify the file format as follow:
           <br \><br \>
           <strong>Probe list</strong>: list of mvps as probeids each on a separate line.
           Optionally can add other fields after the probeid which are
           ignored, unless the pvalue filter is specified, in which case
           eForge assumes that the second field is the minus log10 pvalue
           <br \><br \>
           <strong>BED file</strong>: File given is a bed file of locations (chr\tbeg\tend).  bed
           format should be 0 based and the chromosome should be given as
           chrN.  However we will also accept chomosomes as just N (ensembl)
           and 1-based format where beg and end are the same*.
           <br \>",
       "<strong>Maximum number of probes</strong><br \><br \>
           The web version is limited to 1000 probes. While it is possible to input a larger number
           of probes with the standalone version, however you must consider how the background
           selection occurs.
           <br \>",
    );

    print Template::content_box_1("Options",
       "<strong>Name for this data</strong><br \><br \>
           Supply a label that you want to use for the plotting titles, and
           filenames.
           <br \>",

       "<strong>Analysis data from</strong><br \><br \>
           Dataset to analyse. Either ENCODE data ('encode') or Roadmap
           Epigenome data ('erc'). erc by default.
           <br \>",

       "<strong>Depletion</strong><br \><br \>
           Analyse for DHS depletion pattern instead of the default DHS
           enrichment analysis. Use when dealing with datasets suspected not
           to overlap with DHS. Specifying depletion will be indicated on the
           label (the text \"Depletion Analysis\" will be added to the file
           label).
           <br \>",

       "<strong>Proximity</strong><br \><br \>
           Apply filter for DMPs in proximity (within 1 kb of another test
           DMP). With proximity filter specified, eFORGE will report DMPs
           removed due to proximity with another DMP in the list and will
           randomly pick one of the probes among the set of probes that are in
           proximity (within 1 kb of each other).
           <br \>",

       "<strong>Background repetitions</strong><br \><br \>
           The number of background matching sets to pick and analyse. (1-1000)
           <br \>",

       "<strong>Significance threshold</strong><br \><br \>
           Alter the default binomial p value thresholds. (0 < Min < Max < 1)
           <br \>",
    );

    print Template::content_box_1("Output",
       "<strong>Raw data</strong><br \><br \>
           A tab-separated values (tsv) file with columns for the Zscore, Pvalue, Cell, Tissue,
           File, Probe, Number and  Accession.
           <br \>",
       "<strong>Static chart</strong><br \><br \>
           A PDF chart of the data.
           <br \>",
       "<strong>Interactive chart</strong><br \><br \>
           A dimple (http://dimplejs.org) d3 interactive graphic using rCharts.
           <br \>",
       "<strong>Interactive table</strong><br \><br \>
           A table using the Datatables
           (<a href=\"https://datatables.net\">https://datatables.net</a>) plug-in for the jQuery
           Javascript library, again accessed through rChartsA dimple
           (<a href =\"http://dimplejs.org\">http://dimplejs.org</a>) d3 interactive graphic using
           rCharts.
           <br \>",
       "<strong>R code</strong><br \><br \>
           The source files for generating the charts and the table. You can download the raw data
           and use R (<a href=\"http://www.r-project.org\">http://www.r-project.org</a>)
           with these file to re-generate or modify the output.
           <br \>",
    );

    print Template::end;
    exit(0);
}


=head2 print_download_page

 Arg[1]         : -none-
 Example        : print_download_page();
 Description    : Prints the Download page.
 Returns        : 
 Exceptions     : exit(0) always

=cut

sub print_download_page {
    print $q->header;
    print Template::start($title, $breadcrumbs, $left_menu, $colour, $right_column);
    print Template::header("eFORGE &gt; Download");
    print Template::content_box("Download",
    "The code is available on GitHub:
    <a href=\"https://github.com/charlesbreeze/eFORGE\">https://github.com/charlesbreeze/eFORGE</a>");
    print Template::content_box("License",
    "<strong>eforge.pl</strong> Functional analysis of EWAS DMPs
       <br \><br \>
       Copyright (C) 2014  EMBL - European Bioinformatics Institute
       <br \><br \>
       This program is free software: you can redistribute it and/or modify it
       under the terms of the GNU General Public License as published by the
       Free Software Foundation, either version 3 of the License, or (at your
       option) any later version. This program is distributed in the hope that
       it will be useful, but WITHOUT ANY WARRANTY; without even the implied
       warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
       the GNU General Public License for more details. Neither the
       institution name nor the name eforge.pl can be used to endorse or
       promote products derived from this software without prior written
       permission. For written permission, please contact
       <strong>c.breeze(at)ucl.ac.uk</strong>.
       Products derived from this software may not be called eforge.pl nor may
       eforge.pl appear in their names without prior written permission of the
       developers. You should have received a copy of the GNU General Public
       License along with this program.  If not, see
       <a href=\"http://www.gnu.org/licenses/\">http://www.gnu.org/licenses/</a>");
    print Template::end;
    exit(0);
}


=head2 print_about_page

 Arg[1]         : -none-
 Example        : print_about_page();
 Description    : Prints the About page. If the request is internal, also prints the usage stats,
                  after refreshing them if necessary.
 Returns        : 
 Exceptions     : exit(0) always

=cut

sub print_about_page {
    print $q->header;
    print Template::start($title, $breadcrumbs, $left_menu, $colour, $right_column);
    print Template::header("eFORGE &gt; About");
    print Template::content_box("eFORGE",
    "eFORGE was developed by <a href=\"http://www.ucl.ac.uk/cancer/medical-genomics/mg_staff\">Charles
    Breeze</a> while on secondment at the
    <a href=\"http://www.ebi.ac.uk\">European Bioinformatics Institute</a> as part of the
    <a href=\"http://www.epitrain.eu\">EpiTrain Integrated Training
    Network</a>.",
    "The system is inspired by the <a href=\"http://www.1000genomes.org/forge-analysis\">FORGE
    tool</a> developed by <a href=\"http://www.ebi.ac.uk/about/people/ian-dunham\">Ian Dunham</a>.",
    "The web interface was developed by
    <a href=\"https://iris.ucl.ac.uk/iris/browse/profile?upi=JHERR07\">Javier Herrero</a> from the
    <a href=\"http://www.ucl.ac.uk/cancer/blic\">Bill Lyons Informatics Centre</a> team.",
    "The example data set corresponds to a filtered set of monocyte tDMPs from
    <a href=\"http://genomebiology.com/2014/15/2/R31\">Jaffe AE and Irizarry RA, Genome Biol 2014, 15:R31</a>.",
    "<p align=\"right\">(last updated on ".scalar(localtime((stat("$0"))[9])).")</p>",
    );
    print Template::content_box("System",
    "eForge is currently running on a virtual machine provided by <a href=\"".
    "http://www.cs.ucl.ac.uk/\">UCL Computer Sciences</a> department and maintained by the <a href=\"".
    "http://www.ucl.ac.uk/cancer/blic/\">Bill Lyons Informatics Centre</a>.");
    if ($q->remote_addr() =~ /^128\.40\.233\./ or $q->remote_addr() =~ /^127\.0\.0\.1$/) {
        my $absolute_root_outdir = get_absolute_root_outdir();
        refresh_usage_stats($absolute_root_outdir);
        print Template::content_box("Usage",
            "<img src=\"$WEB_ROOT_OUTDIR/last_week.png\">",
            "<img src=\"$WEB_ROOT_OUTDIR/last_month.png\">",
            "<img src=\"$WEB_ROOT_OUTDIR/scalability.png\">",
            "<strong>Disk usage:</strong> ".qx"cat $absolute_root_outdir/du.txt"." in ".
                qx"cat $absolute_root_outdir/num.txt"." folders",
            "You are accessing this server from ".$q->remote_addr(),
            "<p align=\"right\">(last updated on ".
                scalar(localtime((stat("$absolute_root_outdir/last_week.png"))[9])).")</p>");
    }
    print Template::content_box("Contact",
    "Email <a href=\"http://www.ucl.ac.uk/cancer/medical-genomics/mg_staff\">Charles Breeze</a>:
    c.breeze(at)ucl.ac.uk and
    <a href=\"https://iris.ucl.ac.uk/iris/browse/profile?upi=BECKX39\">Stephan Beck</a>:
    s.beck(at)ucl.ac.uk.",
    );
    print_credits_box();
    print Template::end;
    exit(0);
}


=head2 log_usage

 Arg[1]         : float $runtime_sec, the amount of seconds the program ran for.
 Arg[2]         : int $input_size, the size of the input file, typically determined with -s $input
 Arg[3]         : string $cmd_line, the command line to be stored
 Arg[4]         : string $err, the error message to be stored (if any)
 Example        : log_usage(21.34, 3424, "tool input.txt", undef);
 Description    : Stores the stats on a SQLite db for monitoring usage
 Returns        : 
 Exceptions     : Creates a text file $absolute_outdir/error.txt if an error occurs.

=cut

sub log_usage {
    my ($runtime_sec, $input_size, $cmd_line, $err) = @_;

    ## Note that by the time we run this bit of code, the STDERR and STDOUT are closed.
    ## For debugging purposes, this is run before the results page is closed.

    my $error_msg;
    
    my $dsn = "DBI:SQLite:dbname=$LOG_FILE";

    if (!-e $LOG_FILE) {
        my $dbh = DBI->connect($dsn, "", "") or die $DBI::errstr;

        $dbh->do("CREATE TABLE server_log (
            id INTEGER PRIMARY KEY,
            runtime_sec REAL NOT NULL,
            input_size INTEGER NOT NULL,
            cmd_line TEXT NOT NULL,
            error TEXT DEFAULT NULL,
            timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
            )");
        $dbh->disconnect;
    }

    my $dbh = DBI->connect($dsn, "", "") or
            $error_msg .= "Fail to connect to the server log DB ($DBI::errstr)\n";
    my $sth = $dbh->prepare("INSERT INTO server_log (runtime_sec, input_size, cmd_line, error)".
            " VALUES (?, ?, ?, ?)") or
            $error_msg .= "Fail to prepare query ($dbh->errstr)\n";
    $sth->execute($runtime_sec, $input_size, $cmd_line, $err) or
            $error_msg .= "Fail to write to the server log DB ($DBI::errstr)\n";;
    $dbh->disconnect;
    
    if ($error_msg) {
        my $absolute_outdir = get_absolute_outdir();
        open(ERROR, ">$absolute_outdir/error.txt") or die "Cannot open $absolute_outdir/error.txt";
        print ERROR $error_msg;
        close(ERROR);
    }
}


=head2 refresh_usage_stats

 Arg[1]         : -none-
 Example        : refresh_usage_stats();
 Description    : If the usage images are more than 10 min old, it will re-build them and refresh
                  the disk usage stats.
 Returns        :
 Exceptions     : prints an error webpage on error.

=cut

sub refresh_usage_stats {
    my ($absolute_root_outdir) = @_;
    my $img_file = "$absolute_root_outdir/last_week.png";

    return if (-e $img_file and (-M $img_file < 0.007)); # 0.007 of a day ~ 10 minutes

    my ($ok, $err, $full_buff, $stdout_buff, $stderr_buff) = run(
        command => ["Rscript", "$BIN_DIR/log.R", $LOG_FILE, "--outdir", $absolute_root_outdir,
            "--colour", $plot_colour, "--input_size", "lines"]);
    
    if (!$ok) {
        print Template::error_box("Error while running R to generate the usage plots: $err.",
            map({"ERR: $_"} @$stderr_buff),
            map({"OUT: $_"} @$stdout_buff),
            );
    }
    
    run(command => ["du", "-hs", $absolute_root_outdir, "|", "awk", "-e", "{print \$1}", ">",
        "$absolute_root_outdir/du.txt"]);
    run(command => ["ls", "-d", "$absolute_root_outdir/0x*", "|", "wc", "-l", ">",
        "$absolute_root_outdir/num.txt"]);
}


