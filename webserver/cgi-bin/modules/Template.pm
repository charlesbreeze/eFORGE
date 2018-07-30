package Template;

use strict;
use warnings;

=pod

=head1 Name

Template

=head1 Synopsis

    use CGI;
    use Template;

    my $q = CGI->new;
    print $q->header;

    my $title = "My Tool";

    my $breadcrumbs = [
        {"UCL Cancer Institute" => "http://www.ucl.ac.uk/cancer"},
        {"Cancer Biology" => "http://www.ucl.ac.uk/cancer/rescancerbiol/cancerbiology"},
    ];

    my $left_menu = [
        {"My Tools" => "/"},
        {"Help" => "/help"},
        {"Download" => "/download"},
        {"About" => "/about"},
    ];

    Template::start($title, $breadcrumbs, $left_menu, "pink");

    Template::header($title);

    Template::content_box_1(
        "About",
        "Some text goes here",
        "And some more");

    Template::end;

=head1 Methods

=head2 start

  Arg[1]        : string $title
  Arg[2]        : arrayref of hashes $breadcrumbs. Each element of the array is a hash with one
                  key representing the label and one value with URL. These are used to create the
                  breadcrumbs on the horizontal bar next to the search box
  Arg[3]        : arrayref of hashes $left_menu. Same as above, but for the left-hand menu.
  Arg[4]        : string $colour. This defines the colour scheme for the page. UCL provides several
                  ones like light-blue, bright-blue, orange, rich-red, red, mid-green, yellow or
                  pink.
  Description   : This method prints the header with all the information, including the JavaScript,
                  CSS, etc for this web page. It also includes the head banner with the search box
                  and the breadcrumbs, and the left-hand menu.

=head2 end

  Arg[1]        : N/A
  Description   : This method prints the footer of the page, with all the usual bits.

=head2 header


=head2 content_box

  Arg[1]        : string $title
  Arg[...]      : string $content

 +------------------------------------------------------------------------+
 | Title                                                                  |
 |                                                                        |
 | content1: Lorem ipsum dolor sit amet, consectetuer adipiscing elit,    |
 | sed diam nonummy nibh euismod                                          |
 |                                                                        |
 | content2: tincidunt ut laoreet dolore magna                            | 
 +------------------------------------------------------------------------+


=head2 content_box_1

  Arg[1]        : string $title
  Arg[...]      : string $content

 +------------------------------------------------------------------------+
 | Title       > content1: Lorem ipsum dolor sit amet, consectetuer       |
 |               adipiscing elit, sed diam nonummy nibh euismod           |
 |                                                                        |
 |             > content2: tincidunt ut laoreet dolore magna              |
 +------------------------------------------------------------------------+

=cut

sub error_box {
    my (@content) = @_;
    my $string;

$string .= qq{<div class="content-box">
<br \>
<div class="hiliteButton">
<H3 class="heading">ERROR</H3>
</div>
};
foreach my $this_content (@content) {

$string .= qq{
<span class="p"><strong>$this_content</strong></span>
<br />
<br />
};

}

$string .= qq{
</div>

};

    return $string;

}


sub content_box {
    my ($title, @content) = @_;
    my $string;

$string .= qq{<div class="content-box">
<table width="100%" border="0" cellspacing="0" cellpadding="0">
<tr><td valign="top"><H3 class="heading">$title</H3></td></tr>
};
foreach my $this_content (@content) {

$string .= qq{
<tr><td><span class="p">$this_content</span></td></tr>
};

}

$string .= qq{
</table></div>

};

    return $string;

}


sub content_box_1 {
    my ($title, @content) = @_;
    my $string;
    
    $string .= qq{<div class="content-box">
<table width="100%" border="0" cellspacing="0" cellpadding="0">
<tr>
<td valign="top" width="15%"><H3 class="heading">$title</H3></td>
<td width="5">&nbsp;</td>
<td>
};

foreach my $this_content (@content) {

    $string .= qq{
<span class="p"><strong>&gt;&nbsp;</strong>$this_content</span>
<br />
<br />
};

}

$string .= qq{</td>

<td width="5">&nbsp;</td>
<td valign="top">&nbsp;</td>
</tr>
</table>
</div>

};

    return $string;

}

sub content_box_2 {
    my (@content) = @_;
    my $string;

$string .= qq{<div class="content-box">};
foreach my $this_content (@content) {

$string .= qq{
<span class="p">$this_content</span>
<br />
<br />
};

}

$string .= qq{
</div>

};

    return $string;

}


sub header {
    my ($title) = @_;

    my $string = qq{
<h2 class="heading">
<!--  START-ENTRY group name -->
$title
<!--  END-ENTRY group name -->
</h2>

};

    return $string;
}


sub start {
    my ($title, $breadcrumbs, $left_menu, $colour, $right_column) = @_;

my $layout = "ucl_default_2col_layout.css";
if ($right_column) {
    $layout = "ucl_default_layout.css";
}

    my $string;
    $string .= qq{
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en"
      xml:lang="en">
            <head>

            <link rel="canonical"
                  href="http://www.ucl.ac.uk/cancer" />
            <link rel="dns-prefetch" href="//static.ucl.ac.uk" />
            <title>$title</title>
            <meta http-equiv="X-UA-Compatible" content="IE=9" />
            <link rel="shortcut icon" href="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/images/favicon2.ico" />
            


            <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />


<meta name="contact.name" content="Web Services" />
<meta name="contact.email" content="silva-support\@ucl.ac.uk" />


    






            <!--Start General Layout Styles-->
            
                
                
                    <!--Start UCL default layout 2 column styles-->
                    <link href="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/css/$layout" rel="stylesheet" type="text/css" />
                    <!--[if lte IE 8]><link rel="stylesheet" type="text/css" href="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/css/patches/patch_$layout" /><![endif]-->
                    <!--End UCL default layout 2 column styles-->
                
            
            
                <!--Corporate identity styling below-->
                <link href="http://static.ucl.ac.uk/silva/UCLDefaultLayoutV3/css/screen/corp-identity-$colour.css"
                      rel="stylesheet" type="text/css" />
                <!--Corporate identity styling above-->
            

            <!--End General Layout Styles-->

            


            

            <!--Start Print CSS-->
            <link href="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/css/print/print.css" media="print" rel="stylesheet" type="text/css" />
            <!--Set IE to be print smaller width, so that it doesn't crop content-->
            <!--[if lte IE 8]>
                <style type="text/css">
                    \@media print {body, div.page, #main{width:600px;}}
                </style>
            <![endif]-->
            <!--End Print CSS-->


            <!--Start IE CSS hacks-->
            
            <!--End IE CSS hacks-->

            <!-- Start Javascripts-->

    <!--Start Gallery javascripts-->
        
    <!--End Gallery javascripts-->

    <!--Start Ticker Javascript-->
        <script language="JavaScript" src="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/scripts/ticker.js" type="text/JavaScript"></script>
    <!--End Ticker Javascript-->

    <!--Start JQuery Javascript-->
        <script src="//ajax.googleapis.com/ajax/libs/jquery/1.3.2/jquery.min.js"></script>
        <script type="text/javascript">
            window.jQuery || document.write('<script src="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/scripts/jquery-1.3.2.min.js"><\\/script>')
        </script>
    <!--End JQuery Javascript-->
        
    

    <!--Start UCL Javascript-->
        <!--<script type="text/javascript" src="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/scripts/ucl.min.js"></script>-->
    <!--End UCL Javascript-->

    <!--Start JQuery Core User Interface Javascript-->
        <script type="text/javascript" src="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/scripts/ui.core.min.js"></script>
    <!--End JQuery Core User Interface Javascript-->

    <!--Start JQuery Font Resizer Javascript-->
        <script type="text/javascript" src="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/scripts/fontsizer.jquery.js"></script>
    <!--End JQuery Font Resizer Javascript-->

    <!--Start JQuery cookie Javascript-->
        <script type="text/javascript" src="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/scripts/jquery.cookie.min.js"></script>
    <!--End JQuery cookie Javascript-->

    <!--Start UCL alerts Javascript-->
        <script type="text/javascript" src="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/scripts/ucl.alert.js"></script>
    <!--End UCL alerts Javascript-->



            

    <script type="text/javascript" src="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/scripts/jquery.doc.ready.js"></script>


            <!--start Silva Style Sheets-->
            

            <!-- Stylesheet for high contrast style switching -->
            <link rel="stylesheet" type="text/css" id="contrastCSS" href="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/css/screen/ucl_normal_view.css" />
            <!--end Silva Style Sheets-->

            <!-- _______ Silva version 2.1.7 _______ -->
            </head>
            <body>
                <!-- alert message block -->
                

                <!-- alert message block: END -->
                <div class="page_margins">
                    <!-- Graphic Border - Begin Part 1 -->
                    <div id="border-top">
                        <div id="edge-tl"> </div>
                        <div id="edge-tr"> </div>
                    </div>
                    <!-- Graphic Border - End Part 1 -->
                    <div class="page">
                        <div id="printLogo"><img src="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/images/printlogo.gif" alt="UCL logo" /></div>
                        
                                <div id="header" class="tb-dark_grey">
      <div id="topnav">
        <!-- start: skip link navigation -->
        <a class="skip" href="#navigation" title="skip link">skip to navigation</a><span class="hideme">.</span>
        <a class="skip" href="#content" title="skip link">skip to content</a><span class="hideme">.</span>
        <!-- end: skip link navigation -->
      </div>
      <h1 class="section_header_white">UCL Cancer Institute</h1>
      <span class="section_subheader_white">School of Life and Medical Sciences</span>
      <div id="logo_holder"> 
            <a href="http://www.ucl.ac.uk/">
                <img src="//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/images/corp-identity-$colour.gif"
                     alt="UCL Home" />
            </a>
      </div> 
    </div>


                            <!-- Begin bradcrumbs row -->
    <div id="nav" class="noprint">
      <div class="hlist" id="corp-identity-$colour">
        <!--Start Search box-->
        <div id="search"> 
            <form action="http://search2.ucl.ac.uk/search/search.cgi" method="get" id="googlesearch">
            <input placeholder="Search UCL" name="query" id="query" type="text" value="" accesskey="q" value="Search UCL">
            <input type="hidden" name="collection" value="ucl-public-meta">
            <input type="submit" value="GO" class="submit" name="Submit">
            </form>
        </div>
        <!--End search Box-->
        

<ul>
};





    foreach my $this_breadcrumb ({"UCL Home" => "http://www.ucl.ac.uk"}, @$breadcrumbs) {
        my ($label, $url) = %$this_breadcrumb;
        $string .= qq{  <li>
    <a href="$url">$label</a>
    <!--[if lte IE 7]>&raquo;<![endif]-->
  </li>
};

    }





    $string .= qq{
</ul>

      </div>
    </div>
<!--End Breadcrumbs-->

                            
                            
                        
                        

                        <!-- begin: main content area #main -->
                        <div id="main">
                            <div class="cancer">
                                <!-- begin: #col1 - first float column -->
                                <div id="col1">
                                    <div id="col1_content" class="clearfix">
                                    <div id="col_top_left" class="col_top left_coltop">
    
</div>

                                    <a id="navigation" name="navigation"></a>
                                            
          <!--Start Generate Menu-->
              
          <!--End Generate Menu-->
        <!--Start left include-->
        
            
                <div id="left-silva-content"
                     class="leftcontainer">
<br />
};





    my $name = (keys(%{$left_menu->[0]}))[0];

    my $inside_list = 0;
    foreach my $this_menu_item (@$left_menu) {
        my ($label, $url) = %$this_menu_item;
        if ($label =~ /^__logo__$/) {
            if ($inside_list) {
                $string .= qq{</ul>\n};
                $inside_list = 0;
            }
            $string .= qq{<img src="$url" alt="logo" width="180">\n<br />\n};
            next;
        }
        if ($label =~ /^__title__$/) {
            if ($inside_list) {
                $string .= qq{</ul>\n};
                $inside_list = 0;
            }
            $string .= qq{<h6 class="vlist">$url</h6>\n<br />\n};
            next;
        }
        if (!$inside_list) {
            $string .= qq{<ul class="vlist">\n};
            $inside_list = 1;
        }
        $string .= qq{<li><a href="$url" target="_self" class="menuitem">$label</a></li>\n};

    }
    if ($inside_list) {
        $string .= qq{</ul>\n};
        $inside_list = 0;
    }

#     $string .= qq{<hr><li><a href="http://www.ucl.ac.uk/cancer/blic" target="_self" class="menuitem">Bill Lyons Informatics Centre</a></li>\n};




    $string .= qq{



<br />
<br />
<br />
<br />




</div>
            
        
        <!--End left include-->
      

                                    </div>
                                 </div>
                                <!-- end: #col1 -->
};





if ($right_column) {
    $string .= qq{
                                    <!-- begin: #col2 second float column -->
                                    <div id="col2">
                                        <div id="col2_content" class="clearfix">
                                            <div id="col_top_right" class="col_top right_coltop">
    
</div>

                                            <div class="textsize"> 
              <div class="fontResizer">
                <a href="#" class="smallFont">A</a>
                <a href="#" class="medFont">A</a>
                <a href="#" class="largeFont">A</a>
              </div>
</div>

                                            <!--Start rightcontent-->
<div id="rightcontent" class="container"> 
    
      
        
            <!--insert index_right-->
            <div id="right-column">
$right_column
</div>
        
      
      
    
    <!--Start image randomiser-->
    
    <!--End image randomiser-->
</div>
<!--End rightcontent-->

                                        </div>
                                    </div>
                                    <!-- end: #col2 -->
};

}





    $string .= qq{

                                <!-- begin: #col3 static column -->
                                <div id="col3">
                                    <div id="col3_content" class="clearfix">
                                        <a id="content" name="content"></a>
                                        <!-- skiplink anchor: Content -->
                                        
                                            
                                        
                                        <div id="toptabs-container" class="top-tabs">
    
</div>

                                        <!--Start center content area-->		
<!-- start main slot where silva content gets rendered--> 
};

    return $string;

}


sub end {
    my $string;
    $string .= qq{
<!-- end main slot where silva content gets rendered--> 
<!--End center content area-->

                                        <!--Start last modified info-->

    
	
	
    
    

<!--End last modified info-->




                                        <!--Start contacts include-->

<!--End contacts include-->

                                    </div>
                                    <div id="ie_clearing">&nbsp;</div>
                                <!-- End: IE Column Clearing -->
                                </div>
                                <!-- end: #col3 -->
                            </div>
                        </div>
                        <!-- end: #main -->
                        <!-- begin: #footer -->
                        <!-- =============== FOOTER =============== -->


 <div id="footer"> 
  <ul>
   <li><a href="http://www.ucl.ac.uk/disclaimer">Disclaimer</a></li>
   <li><a href="http://www.ucl.ac.uk/foi">Freedom of Information</a></li>
   <li><a href="http://www.ucl.ac.uk/accessibility">Accessibility</a></li>
   <li><a href="http://www.ucl.ac.uk/privacy">Privacy</a></li>
   <li><a href="http://www.ucl.ac.uk/cookies">Cookies</a></li>
   <li><a href="http://www.ucl.ac.uk/advanced-search">Advanced Search</a></li>
   <li><a href="http://www.ucl.ac.uk/contact-list/">Contact Us</a></li>
  </ul>
  
  <div>
   <address class="vcard">
    <span class="adr">
     University College London - Gower Street - London - WC1E 6BT
    </span>
    <span class="tel">
     <span class="type">Tel</span>:
     <span class="value">+44 (0)20 7679 2000</span>
    </span>
   </address>
   <p>&#169; UCL 1999&#8211;2014</p>
  </div>

 </div>

 
 <!-- end #footer -->

                        
                    </div>
                    <!-- Graphic Border - Begin Part 2 -->
                    <div id="border-bottom">
                        <div id="edge-bl"> </div>
                        <div id="edge-br"> </div>
                    </div>
                    <!-- Graphic Border - End Part 2 -->
                </div>
            </body>
</html>
};

    return $string;

}

1;