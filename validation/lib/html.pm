#!/usr/bin/env perl
=head1 NAME

    html.pm

=head1 SYNOPSIS

	use html.pm;
	html::indexHeader($index_html, $s_exp);
	html::indexBody($index_html, \%counts, \%sensitivity);
	html::indexFooter($index_html);

=head1 AUTHOR     Mike Halagan <mhalagan@nmdp.org>
    
    Associate Bioinformatics Scientist
    3001 Broadway Stree NE
    Minneapolis, MN 55413
    ext. 8225

=head1 DESCRIPTION

    This module creates the html files for the ngs-validation-report script.©

 =head1 LICENSE

    pipeline  Consensus assembly and allele interpretation pipeline.
    Copyright (c) 2015 National Marrow Donor Program (NMDP)

    This library is free software; you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation; either version 3 of the License, or (at
    your option) any later version.

    This library is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.
 
    You should have received a copy of the GNU Lesser General Public License
    along with this library;  if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.

    > http://www.gnu.org/licenses/lgpl.html    

=cut
################################################################################################################
package html;

use Data::Dumper;
use warnings;
use strict;


################################################################################################################
=head2 subjectsHeader

  Title: subjectsHeader
  Usage:
  Function:
  Returns:
  Args:

=cut
sub subjectsHeader{

  my ( $html, $s_exp, $s_ID, $rh_qc_verdict ) = @_;

  my $s_experiment = $s_exp;

  my $header = qq{
    <!DOCTYPE html>
  <html lang="en">
    <head>
      <meta charset="utf-8">
      <meta http-equiv="X-UA-Compatible" content="IE=edge">
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <meta name="description" content="">
      <meta name="author" content="">
      
      <title>Validation Report</title>

      <!-- Bootstrap core CSS -->
      <link href="http://netdna.bootstrapcdn.com/bootstrap/3.0.0/css/bootstrap.min.css" rel="stylesheet">
      
      <script src="../../js/Chart.js"></script>
      <link href="../../css/dashboard.css" rel="stylesheet">


        <style type='text/css'>
    .my-legend .legend-scale ul {
      margin: 0;
      padding: 0;
      float: left;
      list-style: none;
      }
    .my-legend .legend-scale ul li {
      display: block;
      float: left;
      width: 50px;
      margin-bottom: 6px;
      text-align: center;
      font-size: 80%;
      list-style: none;
      }
    .my-legend ul.legend-labels li span {
      display: block;
      float: left;
      height: 15px;
      width: 50px;
      }
    .my-legend .legend-source {
      font-size: 70%;
      color: #999;
      clear: both;
      }
    .my-legend a {
      color: #777;
      }

      #expanding{
        display:none;
      }


      </style>
    </head>

    <body>

      <nav class="navbar navbar-inverse navbar-fixed-top" role="navigation">
        <div class="container-fluid">
          <div class="navbar-header">
            <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
              <span class="sr-only">Toggle navigation</span>
              <span class="icon-bar"></span>
              <span class="icon-bar"></span>
              <span class="icon-bar"></span>
            </button>
            <a class="navbar-brand" href="../../experiment.html">Validation Report - $s_experiment </a>
          </div>
          <div id="navbar" class="navbar-collapse collapse">
            <ul class="nav navbar-nav navbar-right">
              <li><a href="../../experiment.html">Experiments</a></li>
              <li><a href="../../log.html">Log</a></li>
              <li><a href="../../help.html">Help</a></li>
            </ul>
          </div>
        </div>
      </nav>

      <div class="container-fluid">
        <div class="row">
          <div class="col-sm-3 col-md-2 sidebar">
            <ul class="nav nav-sidebar">
              <li><a href="../index.html">Overview</a></li>
              <li class="active"><a href="../results.html">Results<span class="sr-only">(current)</span></a></li>
              <li><a href="../fails.html">Failed</a></li>
              <li><a href="../errors.html">Dropout</a></li>
              <li><a href="../drbx.html">DRBX</a></li>
      };
      print $html $header;
      print $html "<li><a href=\"../qc.html\">QC</a></li>\n" if defined $$rh_qc_verdict{$s_exp};

      my $header2 = qq{
            </ul>
          </div>
          <div class="col-sm-9 col-sm-offset-3 col-md-10 col-md-offset-2 main">
            <h1 class="page-header">$s_ID Validation Results</h1>
        <div class="table-responsive">
              <table class="table table-striped">
                <thead>
                  <tr>
                    <th>Locus</th>
                    <th>Locus Verdict</th>
                    <th>Expected</th>
                    <th>Allele Call 1</th>
                    <th>Allele Call 1</th>
                  </tr>
                </thead>
              <tbody> 

          };
          print $html $header2;

}

################################################################################################################
=head2 htmlFooter

  Title: htmlFooter
  Usage:
  Function:
  Returns:
  Args:

=cut
sub htmlFooter{

  my ( $html ) = @_;

  my $footer = qq{

    <script type="text/javascript">
      function downloadcsv() {
      \$("table").toCSV();
      }
      </script>
      <script type="text/javascript">
      jQuery.fn.toCSV = function() {

        var i = 9;
        var checks     = \$('body').find('input');

        var data    = \$(this).first(); //Only one table
        var csvData = [];
        var tmpArr  = [];
        var tmpStr  = '';

        data.find("tr").each(function() {

              if(\$(this).find("th").length) {
                  \$(this).find("th").each(function() {
                    tmpStr = \$(this).text().replace(/"/g, '""');
                    tmpArr.push(tmpStr);
                  });
                  csvData.push(tmpArr);
              } else {

                if(typeof(checks[i]) !== 'undefined' && checks[i].checked == true){
                    tmpArr = [];
                       \$(this).find("td").each(function() {
                            if(\$(this).text().match(/^\\n/) && \$(this).text().match(/    /)) {
                                
                            } else {
                                tmpStr = \$(this).text().replace(/(\\r\\n|\\n|\\r)/gm,"");
                                tmpArr.push(tmpStr);
                            }
                       });
                    csvData.push(tmpArr.join(','));
                }
                i++;
              }
        });

        var output = csvData.join('\\n');
        var uri = 'data:application/csv;charset=UTF-8,' + encodeURIComponent(output);
        window.open(uri);
      }
      </script>
      
  <script>
       // If the user clicks in the window, set the background color of  to yellow
        function expand(obj) {
          if(obj.getElementsByTagName("span")[0].style.display == "table"){
            obj.getElementsByTagName("span")[0].style.display = "none";
          }else{
             obj.getElementsByTagName("span")[0].style.display = "table";
          }
        }

    
      </script>
              <script>
       // If the user clicks in the window, set the background color of  to yellow
        function Addfilter() {

          if(document.getElementsByClassName("tablesorter-filter-row")[0].style.display == "none"){
            \$('.tablesorter-filter-row').css('display', 'table-row');
          }else{
            \$('.tablesorter-filter-row').css('display', 'none');
          }
          
        }

    
      </script>

            <script>
       // If the user clicks in the window, set the background color of <body> to yellow
        function expand(obj) {
          if(obj.getElementsByTagName("span")[0].style.display == "table"){
            obj.getElementsByTagName("span")[0].style.display = "none";
          }else{
             obj.getElementsByTagName("span")[0].style.display = "table";
          }
        }
      </script>
            <footer>
        <div class="footer">
          <img src="../img/bethematch.jpeg" alt="img02">
          <p><a href="http://bioinformatics.nmdp.org/Copyright_Information.aspx">Copyright ©</a> 2006 - 2014 National Marrow Donor Program. All Rights Reserved.</p>
        </div>
      </footer> 
      <style type="text/css">
          .bargraph {
            list-style: none;
            padding-top: 20px;
            width:560px;
          }
          a.id_links {
        text-decoration: none;
        font:bold;
        color:black;
          }
          a.id_links:hover {
        color:#A8A8A8;
          }
          footer {
          color: #888;
          clear:both;
          position: relative; 
          bottom: 10px;
          left: 0; width: 100%; 
          text-align: center;
          font-size:9px;
          
          background: rgba(255, 255, 255, 0.6);
      }         

        .hovering{
        color:black;
        text-decoration: none;    
          }
          .hovering:hover{
        color:#AAA;
        text-decoration: none;
          }  
      </style>
      
      <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.0.0/js/bootstrap.min.js"></script>
    
    </body>
  </html>
  };
  print $html $footer;

}
################################################################################################################
=head2 htmlHeader

  Title:     htmlHeader
  Usage:
  Function:  Prints the header to the failed.html file
  Returns:
  Args:

=cut
sub htmlHeader{

  my ( $html, $s_exp, $rh_qc_verdict, $n_max_observed, $s_type, $b_drbx, $b_error, $b_failed  ) = @_;


  my $s_experiment = $s_exp;

  my $header = qq{
    <!DOCTYPE html>
  <html lang="en">
    <head>
      <meta charset="utf-8">
      <meta http-equiv="X-UA-Compatible" content="IE=edge">
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <meta name="description" content="">
      <meta name="author" content="">
      
      <title>Validation Report</title>

      <!-- Bootstrap core CSS -->
      <link href="http://netdna.bootstrapcdn.com/bootstrap/3.0.0/css/bootstrap.min.css" rel="stylesheet">

      <!-- Custom styles for this template -->
      <link   href="../css/dashboard.css" rel="stylesheet">
      <script type="text/javascript" src="../js/jquery-latest.min.js"></script>
      <script type="text/javascript" src="../js/jquery-ui.min.js"></script>
      <script type="text/javascript" src="../js/docs.min.js"></script>
      <script type="text/javascript" src="../js/Chart.js"></script>
 
      <!-- Tablesorter: theme -->
      <link class="theme default" rel="stylesheet" href="../css/theme.bootstrap.css">

      <!-- Tablesorter script: required -->
      <script src="http://mottie.github.io/tablesorter/js/jquery.tablesorter.js"></script>
      <script src="http://mottie.github.io/tablesorter/js/jquery.tablesorter.widgets.js"></script>
      <script src="http://mottie.github.io/tablesorter/addons/pager/jquery.tablesorter.pager.js"></script>

  <script id="js">\$(function(){

    \$(\'#table1\').tablesorter({
      widthFixed : true,
      showProcessing: true,
      headerTemplate : '{content} {icon}', // Add icon for various themes

      widgets: [ 'stickyHeaders', 'filter' ],

      widgetOptions: {

        // extra class name added to the sticky header row
        stickyHeaders : '',
        // number or jquery selector targeting the position:fixed element
        stickyHeaders_offset : 50,
        // added to table ID, if it exists
        stickyHeaders_cloneId : '-sticky',
        // trigger "resize" event on headers
        stickyHeaders_addResizeEvent : false,
        // if false and a caption exist, it won't be included in the sticky header
        stickyHeaders_includeCaption : false,
        // The zIndex of the stickyHeaders, allows the user to adjust this to their needs
        stickyHeaders_zIndex : 2,
        // jQuery selector or object to attach sticky header to
        stickyHeaders_attachTo : null,
        // jQuery selector or object to monitor horizontal scroll position (defaults: xScroll > attachTo > window)
        stickyHeaders_xScroll : null,
        // jQuery selector or object to monitor vertical scroll position (defaults: yScroll > attachTo > window)
        stickyHeaders_yScroll : null,
        // scroll table top into view after filtering
        stickyHeaders_filteredToTop: true,

      }
    });
    \$('input[type=\"checkbox\"][value=\"change\"]').change(function() {
         if(this.checked) {
            var checkboxes = \$('body').find('input');
            for(var i =0;i<checkboxes.length;i++){
              if(checkboxes[i].value !== "change"){
                checkboxes[i].checked = true;
              }
            }
         }
         if(!this.checked) {
            var checkboxes = \$('body').find('input');
            for(var i =0;i<checkboxes.length;i++){
              if(checkboxes[i].value !== "change"){
                checkboxes[i].checked = false;
              }
            }
         }

     });
    
});</script>
        <style type='text/css'>
    .my-legend .legend-scale ul {
      margin: 0;
      padding: 0;
      float: left;
      list-style: none;
      }
    .my-legend .legend-scale ul li {
      display: block;
      float: left;
      width: 50px;
      margin-bottom: 6px;
      text-align: center;
      font-size: 80%;
      list-style: none;
      }
    .my-legend ul.legend-labels li span {
      display: block;
      float: left;
      height: 15px;
      width: 50px;
      }
    .my-legend .legend-source {
      font-size: 70%;
      color: #999;
      clear: both;
      }
    .my-legend a {
      color: #777;
      }

      #expanding{
        display:none;
      }


      </style>
};
print $html $header;

my $header2 = qq{
    </head>

    <body>

      <nav class="navbar navbar-inverse navbar-fixed-top" role="navigation">
        <div class="container-fluid">
          <div class="navbar-header">
            <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
              <span class="sr-only">Toggle navigation</span>
              <span class="icon-bar"></span>
              <span class="icon-bar"></span>
              <span class="icon-bar"></span>
            </button>
            <a class="navbar-brand" href="../experiment.html">Validation Report - $s_experiment </a>
          </div>
          <div id="navbar" class="navbar-collapse collapse">
            <ul class="nav navbar-nav navbar-right">
              <li><a href="../experiment.html">Experiments</a></li>
              <li><a href="../log.html">Log</a></li>
              <li><a href="../help.html">Help</a></li>
            </ul>
          </div>
        </div>
      </nav>

      <div class="container-fluid">
        <div class="row">
          <div class="col-sm-3 col-md-2 sidebar">
            <ul class="nav nav-sidebar">
    };

 print $html $header2;

  my $index_results   = qq{
    <li ><a href="index.html">Overview</a></li>
  };
  print $html $index_results;
  my $result_inactive = qq{
    <li ><a href="results.html">Results</a></li>
  };
  print $html $result_inactive,"\n" if($s_type ne "results");    
  my $results_active   = qq{
    <li class="active"><a href="results.html">Results<span class="sr-only">(current)</span></a></li>
  };
  print $html $results_active,"\n" if($s_type eq "results");    

  my $failed_inactive = qq{
    <li ><a href="fails.html">Failed</a></li>
  };
  print $html $failed_inactive,"\n" if($s_type ne "failed" && $b_failed);
  my $fails_active   = qq{
    <li class="active"><a href="fails.html">Failed<span class="sr-only">(current)</span></a></li>
  };
  print $html $fails_active,"\n" if( $s_type eq "failed");

  my $error_inactive  = qq{
    <li ><a href="errors.html">Dropout</a></li>
  };
  print $html $error_inactive if($s_type ne "error" && $b_error);
  my $error_active   = qq{
    <li class="active"><a href="errors.html">Dropout<span class="sr-only">(current)</span></a></li>
  };
  print $html $error_active,"\n" if( $s_type eq "error");

  my $drbx_inactive   = qq{
    <li ><a href="drbx.html">DRBX</a></li>
  };
  print $html $drbx_inactive if($s_type ne "drbx" && $b_drbx);
  my $drbx_active    = qq{
    <li class="active"><a href="drbx.html">DRBX<span class="sr-only">(current)</span></a></li>
  };
  print $html $drbx_active,"\n" if( $s_type eq "drbx");
      
  print $html "<li><a href=\"qc.html\">QC</a></li>\n" if defined $$rh_qc_verdict{$s_exp};
  my $header3 = qq{
            </ul>
          </div>
          <div class="col-sm-9 col-sm-offset-3 col-md-10 col-md-offset-2 main">
          <div style="float:right;margin-top:20px;">
                <a class="hovering" >
              <span  onclick="Addfilter()" style="font-weight:bold;">Filter Rows</span>         
            </a>
            <span  style="font-weight:bold;"> | </span>
                <a class="hovering" >
                  <span  onclick="downloadcsv()" style="font-weight:bold;">Download to CSV</span>
                </a>
            </div>
            <h1 class="page-header">Validation Results</h1>
            <div class="table-responsive">
              <table class="table table-striped tablesorter-bootstrap" id="table1">
                <thead>
                  <tr>
                    <td style="border-color:#dddddd;border-style:solid;border-left-width:0px;border-right-width:0px;border-bottom-width:2px;" class="filter-false sorter-false">
                      <div style="margin-top:17px;"><input type = "checkbox" value="change"/></div>
                    </td>  
                    <th>Sample ID</th>
                    <th>Locus</th>
                    <th>Locus Verdict</th>
                    <th>Allele Verdict</th>
                    <th>Observed Seqs</th>
                    <th>Expected Allele</th>
          };
          print $html $header3;
          for(1..$n_max_observed){
            print $html "<th>Observed Allele Call ",$_,"</th>\n";
          }
          my $header4 = qq{
                  </tr>
                </thead>
                <tbody>
          };
          print $html $header4;

}

################################################################################################################
=head2 logHeader

	Title:     logHeader
	Usage:
	Function:  Prints the header to the log.html file
	Returns:
	Args:

=cut
sub logHeader{

	my $html = shift;

	my $header = 
		qq{
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="description" content="">
    <meta name="author" content="">
    
    <title>Validation Report</title>

    <!-- Bootstrap core CSS -->
    <link href="http://netdna.bootstrapcdn.com/bootstrap/3.0.0/css/bootstrap.min.css" rel="stylesheet">
    <script src="js/Chart.js"></script>
    <link href="css/dashboard.css" rel="stylesheet">


      <style type='text/css'>
  .my-legend .legend-scale ul {
    margin: 0;
    padding: 0;
    float: left;
    list-style: none;
    }
  .my-legend .legend-scale ul li {
    display: block;
    float: left;
    width: 50px;
    margin-bottom: 6px;
    text-align: center;
    font-size: 80%;
    list-style: none;
    }
  .my-legend ul.legend-labels li span {
    display: block;
    float: left;
    height: 15px;
    width: 50px;
    }
  .my-legend .legend-source {
    font-size: 70%;
    color: #999;
    clear: both;
    }
  .my-legend a {
    color: #777;
    }

    #expanding{
      display:none;
    }

    </style>
  </head>

  <body>

    <nav class="navbar navbar-inverse navbar-fixed-top" role="navigation">
      <div class="container-fluid">
        <div class="navbar-header">
          <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
            <span class="sr-only">Toggle navigation</span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a class="navbar-brand" href="experiment.html">Validation Report</a>
        </div>
        <div id="navbar" class="navbar-collapse collapse">
          <ul class="nav navbar-nav navbar-right">
            <li><a href="experiment.html">Experiments</a></li>
            <li><a href="log.html">Log</a></li>
            <li><a href="help.html">Help</a></li>
          </ul>
        </div>
      </div>
    </nav>

	<div class="container-fluid">
		<div class="row">
			<div class="col-sm-3 col-md-2 sidebar">

			</div>
        </div>
		<div class="col-sm-9 col-sm-offset-3 col-md-10 col-md-offset-2 main">
      		<h1 class="page-header">Log</h1>

        

		};
	print $html $header;

}
################################################################################################################
=head2 experimentsHtml

	Title:     experimentsHtml
	Usage:
	Function:  Prints the header to the log.html file
	Returns:
	Args:

=cut
sub experimentsHtml{

	my ($html, $rh_counts, $rh_experiments, $b_verbose) = @_;

	my %h_counts = %$rh_counts;
	my %h_experiments = %$rh_experiments;

	my $header = 
		qq{
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="description" content="">
    <meta name="author" content="">
    
    <title>Validation Report</title>

    <!-- Bootstrap core CSS -->
    <link href="http://netdna.bootstrapcdn.com/bootstrap/3.0.0/css/bootstrap.min.css" rel="stylesheet">
    <script src="js/Chart.js"></script>
    <link href="css/dashboard.css" rel="stylesheet">


      <style type='text/css'>
  .my-legend .legend-scale ul {
    margin: 0;
    padding: 0;
    float: left;
    list-style: none;
    }
  .my-legend .legend-scale ul li {
    display: block;
    float: left;
    width: 60px;
    margin-bottom: 6px;
    text-align: center;
    font-size: 80%;
    list-style: none;
    }
  .my-legend ul.legend-labels li span {
    display: block;
    float: left;
    height: 15px;
    width: 60px;
    }
  .my-legend .legend-source {
    font-size: 70%;
    color: #999;
    clear: both;
    }
  .my-legend a {
    color: #777;
    }

    #expanding{
      display:none;
    }


    </style>
};
print $html $header;

print $html "<script>\n";
print $html "\$(document).ready(function(){ \n";
print $html "\$(\"#myTable\").tablesorter(); \n";
print $html "}\n";
print $html ");\n";
print $html "</script>\n";

my $header2 = qq{
  </head>

  <body>
    <nav class="navbar navbar-inverse navbar-fixed-top" role="navigation">
      <div class="container-fluid">
        <div class="navbar-header">
          <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
            <span class="sr-only">Toggle navigation</span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a class="navbar-brand" href="experiment.html">Validation Report</a>
        </div>
        <div id="navbar" class="navbar-collapse collapse">
          <ul class="nav navbar-nav navbar-right">
            <li><a href="experiment.html">Experiments</a></li>
            <li><a href="log.html">Log</a></li>
            <li><a href="help.html">Help</a></li>
          </ul>
        </div>
      </div>
    </nav>
    };

    print $html $header2;

    my $sidebar = qq{
    <div class="container-fluid">
      <div class="row">
        <div class="col-sm-3 col-md-2 sidebar">
          <ul class="nav nav-sidebar">
    };
    print $html $sidebar;
    foreach my $s_exp (sort keys %h_experiments){
    	print $html "<li><a href=\"".$s_exp."/index.html\">".$s_exp."</a></li>\n";
    }
    print $html "</ul>\n</div>\n</div>\n";
	

      my $table = qq{
	        <div class="col-sm-9 col-sm-offset-3 col-md-10 col-md-offset-2 main">
              <div class='my-legend' style="float:right;">
                <div class='legend-scale'>
                  <ul class='legend-labels'>
                <li><span style='background:#46BFBD;'></span><b>PASS</b></li>
                <li><span style='background:#F7464A;'></span><b>FAIL</b></li>
                <li><span style='background:#D3D5D6;'></span><b>DROPOUT</b></li>
                <li><span style='background:#0d8fdb;'></span><b>DRBX</b></li>
                  </ul>
                </div>
              </div>
	          <h1 class="page-header">Experiments</h1>  
	 };
   print $html $table;


		 my $end_table = qq{
		 		              
		              </tbody>
		              </table>
		            </div>

		 };


########## TABLE 1 ##########
my $table2 = qq{
    <div class="table-responsive">
              <table class="table table-striped tablesorter" id="myTable">
                <thead>
                  <tr>
                    <th>Experiment</th>
                    <th># of subjects</th>
                    <th>Passed Subjects</th>
                    <th>Failed Subjects</th>
                    <th>Dropout Subjects</th>
                    <th>DRBX Subjects</th>
                  </tr>
                </thead>
                <tbody>
      
   };
   print $html $table2;

   foreach my $s_exp (sort keys %h_experiments){
    print $html "\t<tr>\n";
    print $html "\t\t<td>$s_exp</td>\n";
    for my $s_type (qw(SUBJECT LOCUS ALLELE)){
      for my $s_report_cde (qw(DRBX ERROR FAIL PASS)){
          $h_counts{$s_exp}{$s_type}{TOTAL}{$s_report_cde} = !defined $h_counts{$s_exp}{$s_type}{TOTAL}{$s_report_cde} ? 
            0 : $h_counts{$s_exp}{$s_type}{TOTAL}{$s_report_cde};
      }
    }

    my $total_subjects = $h_counts{$s_exp}{SUBJECT}{TOTAL}{PASS} +
      $h_counts{$s_exp}{SUBJECT}{TOTAL}{FAIL} + $h_counts{$s_exp}{SUBJECT}{TOTAL}{ERROR} + $h_counts{$s_exp}{SUBJECT}{TOTAL}{DRBX};;

    my $total_loci = $h_counts{$s_exp}{LOCUS}{TOTAL}{PASS} +
      $h_counts{$s_exp}{LOCUS}{TOTAL}{FAIL} + $h_counts{$s_exp}{LOCUS}{TOTAL}{ERROR} + $h_counts{$s_exp}{LOCUS}{TOTAL}{DRBX};

    my $total_allele = $h_counts{$s_exp}{ALLELE}{TOTAL}{PASS} +
      $h_counts{$s_exp}{ALLELE}{TOTAL}{FAIL} + $h_counts{$s_exp}{ALLELE}{TOTAL}{ERROR} + $h_counts{$s_exp}{ALLELE}{TOTAL}{DRBX};  

    print $html "\t\t<td>$total_subjects</td>\n";

    if($total_subjects == 0){ 
      print $s_exp," WFT!\n";
      print STDERR Dumper(%h_counts),"\n";
      exit;
    }
    #####  Subjects #####
    my $percent_passed = "(".sprintf("%2.1f",(($h_counts{$s_exp}{SUBJECT}{TOTAL}{PASS} / $total_subjects) * 100))."%)";
    print STDERR "------\n" if $b_verbose;
    print STDERR "PASSED: ",$s_exp,"\t",$h_counts{$s_exp}{SUBJECT}{TOTAL}{PASS}." ".$percent_passed,"\n"if $b_verbose;
    $percent_passed = length($h_counts{$s_exp}{SUBJECT}{TOTAL}{PASS}) == 1 ? "&nbsp;&nbsp;".$percent_passed : $percent_passed;
    print $html "\t\t<td>".$h_counts{$s_exp}{SUBJECT}{TOTAL}{PASS}." ".$percent_passed."</td>\n";

    my $percent_failed = "(".sprintf("%2.1f",(($h_counts{$s_exp}{SUBJECT}{TOTAL}{FAIL} / $total_subjects) * 100))."%)";
    print STDERR "FAILED: ",$s_exp,"\t",$h_counts{$s_exp}{SUBJECT}{TOTAL}{FAIL}." ".$percent_failed,"\n" if $b_verbose;
    $percent_failed = length($h_counts{$s_exp}{SUBJECT}{TOTAL}{FAIL}) == 1 ? "&nbsp;&nbsp;".$percent_failed : $percent_failed;
    print $html "\t\t<td>".$h_counts{$s_exp}{SUBJECT}{TOTAL}{FAIL}." ".$percent_failed."</td>\n";

    my $percent_error = "(".sprintf("%2.1f",(($h_counts{$s_exp}{SUBJECT}{TOTAL}{ERROR} / $total_subjects) * 100))."%)";
    print STDERR "ERROR: ",$s_exp,"\t",$h_counts{$s_exp}{SUBJECT}{TOTAL}{ERROR}."&nbsp;&nbsp;".$percent_error,"\n" if $b_verbose;
    $percent_error = length($h_counts{$s_exp}{SUBJECT}{TOTAL}{ERROR}) == 1 ? " ".$percent_error : $percent_error;
    print $html "\t\t<td>".$h_counts{$s_exp}{SUBJECT}{TOTAL}{ERROR}." ".$percent_error."</td>\n";
    print STDERR "------\n" if $b_verbose;

    my $percent_drbx = "(".sprintf("%2.1f",(($h_counts{$s_exp}{SUBJECT}{TOTAL}{DRBX} / $total_subjects) * 100))."%)";
    print STDERR "DRBX: ",$s_exp,"\t",$h_counts{$s_exp}{SUBJECT}{TOTAL}{DRBX}."&nbsp;&nbsp;".$percent_drbx,"\n" if $b_verbose;
    $percent_drbx = length($h_counts{$s_exp}{SUBJECT}{TOTAL}{DRBX}) == 1 ? " ".$percent_drbx : $percent_drbx;
    print $html "\t\t<td>".$h_counts{$s_exp}{SUBJECT}{TOTAL}{DRBX}." ".$percent_drbx."</td>\n";
    print STDERR "------\n" if $b_verbose;
  

   }


  print $html $end_table;

  
  ########## TABLE 1 ##########
my $table3 = qq{
    <div class="table-responsive">
              <table class="table table-striped tablesorter" id="myTable">
                <thead>
                  <tr>
                    <th>Experiment</th>
                    <th># of Loci</th>
                    <th>Passed Locus</th>
                    <th>Failed Locus</th>
                    <th>Dropout Locus</th>
                    <th>DRBX Locus</th>
                  </tr>
                </thead>
                <tbody>
      
   };
   print $html $table3;

   foreach my $s_exp (sort keys %h_experiments){
    print $html "\t<tr>\n";
    print $html "\t\t<td>$s_exp</td>\n";
    my $total_subjects = $h_counts{$s_exp}{SUBJECT}{TOTAL}{PASS} +
      $h_counts{$s_exp}{SUBJECT}{TOTAL}{FAIL} + $h_counts{$s_exp}{SUBJECT}{TOTAL}{ERROR} + $h_counts{$s_exp}{SUBJECT}{TOTAL}{DRBX};;

    my $total_loci = $h_counts{$s_exp}{LOCUS}{TOTAL}{PASS} +
      $h_counts{$s_exp}{LOCUS}{TOTAL}{FAIL} + $h_counts{$s_exp}{LOCUS}{TOTAL}{ERROR} + $h_counts{$s_exp}{LOCUS}{TOTAL}{DRBX};

    my $total_allele = $h_counts{$s_exp}{ALLELE}{TOTAL}{PASS} +
      $h_counts{$s_exp}{ALLELE}{TOTAL}{FAIL} + $h_counts{$s_exp}{ALLELE}{TOTAL}{ERROR} + $h_counts{$s_exp}{ALLELE}{TOTAL}{DRBX};  

    print $html "\t\t<td>$total_loci</td>\n";

    #####  Locus #####
    my $loci_passed = "(".sprintf("%2.1f",(($h_counts{$s_exp}{LOCUS}{TOTAL}{PASS} / $total_loci) * 100))."%)";
    print STDERR "PASSED: ",$s_exp,"\t",$h_counts{$s_exp}{LOCUS}{TOTAL}{PASS}." ".$loci_passed,"\n" if $b_verbose;
    $loci_passed = length($h_counts{$s_exp}{LOCUS}{TOTAL}{PASS}) == 1 ? "&nbsp;&nbsp;".$loci_passed : $loci_passed;
    print $html "\t\t<td>".$h_counts{$s_exp}{LOCUS}{TOTAL}{PASS}." ".$loci_passed."</td>\n";

    my $loci_failed = "(".sprintf("%2.1f",(($h_counts{$s_exp}{LOCUS}{TOTAL}{FAIL} / $total_loci) * 100))."%)";
    print STDERR "FAILED: ",$s_exp,"\t",$h_counts{$s_exp}{LOCUS}{TOTAL}{FAIL}." ".$loci_failed,"\n" if $b_verbose;
    $loci_failed = length($h_counts{$s_exp}{LOCUS}{TOTAL}{FAIL}) == 1 ? "&nbsp;&nbsp;".$loci_failed : $loci_failed;
    print $html "\t\t<td>".$h_counts{$s_exp}{LOCUS}{TOTAL}{FAIL}." ".$loci_failed."</td>\n";

    my $loci_error = "(".sprintf("%2.1f",(($h_counts{$s_exp}{LOCUS}{TOTAL}{ERROR} / $total_loci) * 100))."%)";
    print STDERR "ERROR: ",$s_exp,"\t",$h_counts{$s_exp}{LOCUS}{TOTAL}{ERROR}." ".$loci_error,"\n" if $b_verbose;
    $loci_error = length($h_counts{$s_exp}{LOCUS}{TOTAL}{ERROR}) == 1 ? "&nbsp;&nbsp;".$loci_error : $loci_error;
    print $html "\t\t<td>".$h_counts{$s_exp}{LOCUS}{TOTAL}{ERROR}." ".$loci_error."</td>\n";
    print STDERR "------\n" if $b_verbose;

    my $loci_drbx = "(".sprintf("%2.1f",(($h_counts{$s_exp}{LOCUS}{TOTAL}{DRBX} / $total_loci) * 100))."%)";
    print STDERR "DRBX: ",$s_exp,"\t",$h_counts{$s_exp}{LOCUS}{TOTAL}{DRBX}." ".$loci_drbx,"\n" if $b_verbose;
    $loci_drbx = length($h_counts{$s_exp}{LOCUS}{TOTAL}{DRBX}) == 1 ? "&nbsp;&nbsp;".$loci_drbx : $loci_drbx;
    print $html "\t\t<td>".$h_counts{$s_exp}{LOCUS}{TOTAL}{DRBX}." ".$loci_drbx."</td>\n";
    print STDERR "------\n" if $b_verbose;

   }


  print $html $end_table;
    
  

 ########## TABLE 1 ##########
my $table4 = qq{
    <div class="table-responsive">
              <table class="table table-striped tablesorter" id="myTable">
                <thead>
                  <tr>
                    <th>Experiment</th>
                    <th># of Alleles</th>
                    <th>Passed Allele</th>
                    <th>Failed Allele</th>
                    <th>Dropout Allele</th>
                    <th>DRBX Allele</th>
                  </tr>
                </thead>
                <tbody>
      
   };
   print $html $table4;

   foreach my $s_exp (sort keys %h_experiments){
    print $html "\t<tr>\n";
    print $html "\t\t<td>$s_exp</td>\n";

    $h_counts{$s_exp}{SUBJECT}{TOTAL}{PASS} = defined $h_counts{$s_exp}{SUBJECT}{TOTAL}{PASS} ? 0 : $h_counts{$s_exp}{SUBJECT}{TOTAL}{PASS};

    my $total_subjects = $h_counts{$s_exp}{SUBJECT}{TOTAL}{PASS} +
      $h_counts{$s_exp}{SUBJECT}{TOTAL}{FAIL} + $h_counts{$s_exp}{SUBJECT}{TOTAL}{ERROR} + $h_counts{$s_exp}{SUBJECT}{TOTAL}{DRBX};;

    my $total_loci = $h_counts{$s_exp}{LOCUS}{TOTAL}{PASS} +
      $h_counts{$s_exp}{LOCUS}{TOTAL}{FAIL} + $h_counts{$s_exp}{LOCUS}{TOTAL}{ERROR} + $h_counts{$s_exp}{LOCUS}{TOTAL}{DRBX};

    my $total_allele = $h_counts{$s_exp}{ALLELE}{TOTAL}{PASS} +
      $h_counts{$s_exp}{ALLELE}{TOTAL}{FAIL} + $h_counts{$s_exp}{ALLELE}{TOTAL}{ERROR} + $h_counts{$s_exp}{ALLELE}{TOTAL}{DRBX};  

    print $html "\t\t<td>$total_allele</td>\n";

    #####  Locus #####
    my $loci_passed = "(".sprintf("%2.1f",(($h_counts{$s_exp}{ALLELE}{TOTAL}{PASS} / $total_allele) * 100))."%)";
    print STDERR "PASSED: ",$s_exp,"\t",$h_counts{$s_exp}{ALLELE}{TOTAL}{PASS}." ".$loci_passed,"\n" if $b_verbose;
    $loci_passed = length($h_counts{$s_exp}{ALLELE}{TOTAL}{PASS}) == 1 ? "&nbsp;&nbsp;".$loci_passed : $loci_passed;
    print $html "\t\t<td>".$h_counts{$s_exp}{ALLELE}{TOTAL}{PASS}." ".$loci_passed."</td>\n";

    my $loci_failed = "(".sprintf("%2.1f",(($h_counts{$s_exp}{ALLELE}{TOTAL}{FAIL} / $total_allele) * 100))."%)";
    print STDERR "FAILED: ",$s_exp,"\t",$h_counts{$s_exp}{ALLELE}{TOTAL}{FAIL}." ".$loci_failed,"\n" if $b_verbose;
    $loci_failed = length($h_counts{$s_exp}{ALLELE}{TOTAL}{FAIL}) == 1 ? "&nbsp;&nbsp;".$loci_failed : $loci_failed;
    print $html "\t\t<td>".$h_counts{$s_exp}{ALLELE}{TOTAL}{FAIL}." ".$loci_failed."</td>\n";

    my $loci_error = "(".sprintf("%2.1f",(($h_counts{$s_exp}{ALLELE}{TOTAL}{ERROR} / $total_allele) * 100))."%)";
    print STDERR "ERROR: ",$s_exp,"\t",$h_counts{$s_exp}{ALLELE}{TOTAL}{ERROR}." ".$loci_error,"\n" if $b_verbose;
    $loci_error = length($h_counts{$s_exp}{ALLELE}{TOTAL}{ERROR}) == 1 ? "&nbsp;&nbsp;".$loci_error : $loci_error;
    print $html "\t\t<td>".$h_counts{$s_exp}{ALLELE}{TOTAL}{ERROR}." ".$loci_error."</td>\n";
    print STDERR "------\n" if $b_verbose;

    my $loci_drbx = "(".sprintf("%2.1f",(($h_counts{$s_exp}{ALLELE}{TOTAL}{DRBX} / $total_allele) * 100))."%)";
    print STDERR "DRBX: ",$s_exp,"\t",$h_counts{$s_exp}{ALLELE}{TOTAL}{DRBX}." ".$loci_drbx,"\n" if $b_verbose;
    $loci_drbx = length($h_counts{$s_exp}{ALLELE}{TOTAL}{DRBX}) == 1 ? "&nbsp;&nbsp;".$loci_drbx : $loci_drbx;
    print $html "\t\t<td>".$h_counts{$s_exp}{ALLELE}{TOTAL}{DRBX}." ".$loci_drbx."</td>\n";
    print STDERR "------\n" if $b_verbose;

   }


  print $html $end_table;
    
  
	my $charts = qq{

	          <div class="row placeholders">  
	              <div class="col-xs-6 col-sm-3 placeholder" style="margin-left: 20px;margin-top:20px;">
	              <div id="canvas-holder">
	                <canvas id="canvas" width="300" height="300"/>
	              </div>
	               <h4>Subjects</h4>
	               <span class="text-muted">Number of subjects pass and fails for each experiemnt</span>
	            </div>
	            <div class="col-xs-6 col-sm-3 placeholder" style="margin-left: 20px;margin-top:20px;">
	              <div id="canvas-holder">
	                <canvas id="canvas2" width="300" height="300" />
	              </div>
	               <h4>Locus</h4>
	               <span class="text-muted">Number of locus pass and fails for each experiemnt</span>
	            </div>    
	            <div class="col-xs-6 col-sm-3 placeholder" style="margin-left: 20px;margin-top:20px;">
	              <div id="canvas-holder">
	                <canvas id="canvas3" width="300" height="300"/>
	              </div>
	               <h4>Allele</h4>
	               <span class="text-muted">Number of allele pass and fails for each experiemnt</span>
	            </div>        

	          </div>
	         </div>

	     };


		print $html $charts;

		print $html "<SCRIPT>\n";

	    ######## 
	    # Subjects JS
	    ########
		print $html "var barChartData = {\n";
		print $html "\t\tlabels : [";
		my $label;
		foreach my $s_exp (sort keys %h_experiments){
			$label .= "\"$s_exp\",";
		}
		$label =~ s/,$//;
		print $html $label."],\n";
	    print $html "\t\tdatasets : [\n";
	   

	    #Pass columns
	    print $html "\t\t\t{\n";
		my $chart_pass = qq{
	            fillColor : "#46BFBD",
	            strokeColor : "#46BFBD",
	            highlightFill :  "#5AD3D1",
	            highlightStroke :  "#5AD3D1",
	        };
	    print $html $chart_pass;
	    print $html "\t\t\tdata : [";
	    my $subject_chart;
	    foreach my $s_exp (sort keys %h_experiments){
	    	$subject_chart .= $h_counts{$s_exp}{SUBJECT}{TOTAL}{PASS}.",";
	    }
	    $subject_chart =~ s/,$//;
	    print $html $subject_chart."],\n";
	    print $html "\t\t\t},\n";

	    #Fail columns
	    print $html "\t\t\t{\n";
	    my $chart_fail = qq{
	            fillColor : "#F7464A",
	            strokeColor : "#F7464A",
	            highlightFill: "#FF5A5E",
	            highlightStroke: "#FF5A5E",
	        };
	    print $html $chart_fail;    
	    print $html "\t\t\tdata : [";
	    my $fail_chart;
	    foreach my $s_exp (sort keys %h_experiments){
	    	$fail_chart .= $h_counts{$s_exp}{SUBJECT}{TOTAL}{FAIL}.",";
	    }
	    $fail_chart =~ s/,$//;
	    print $html $fail_chart."],\n";
	    print $html "\t\t\t},\n";

      #DRBX columns
      print $html "\t\t\t{\n";
      my $chart_drbx = qq{
              fillColor : "#0d8fdb",
              strokeColor : "#0d8fdb",
              highlightFill: "#3DA5E2",
              highlightStroke: "#3DA5E2",
          };
      print $html $chart_drbx;    
      print $html "\t\t\tdata : [";
      my $drbx_chart;
      foreach my $s_exp (sort keys %h_experiments){
        $drbx_chart .= $h_counts{$s_exp}{SUBJECT}{TOTAL}{DRBX}.",";
      }
      $drbx_chart =~ s/,$//;
      print $html $drbx_chart."],\n";
      print $html "\t\t\t},\n";

	    #Error charts
	    print $html "\t\t\t{\n";   
	    my $chart_error = qq{   
		        fillColor : "#D3D5D6",
		        strokeColor : "#D3D5D6",
		        highlightFill :  "#DCDDDE",
		        highlightStroke :  "#DCDDDE",
		    };
		print $html $chart_error;        
	    print $html "\t\t\tdata : [";
	    my $error_chart;
	    foreach my $s_exp (sort keys %h_experiments){
	    	$error_chart .= $h_counts{$s_exp}{SUBJECT}{TOTAL}{ERROR}.",";
	    }
	    $error_chart =~ s/,$//;
	    print $html $error_chart."],\n";
	    print $html "\t\t\t}\n";
	    print $html "\t\t]\n";
	    print $html "\t}\n";


	    ######## 
	    # Locus JS
	    ########
	   	print $html "var barChartData2 = {\n";
		  print $html "\t\tlabels : [";
		  print $html $label."],\n";
	    print $html "\t\tdatasets : [\n";
	   
	    #Pass columns
	    print $html "\t\t\t{\n";
	    print $html $chart_pass;
	    print $html "\t\t\tdata : [";
	    my $locus_chart;
	    foreach my $s_exp (sort keys %h_experiments){
	    	$locus_chart .= $h_counts{$s_exp}{LOCUS}{TOTAL}{PASS}.",";
	    }
	    $locus_chart =~ s/,$//;
	    print $html $locus_chart."],\n";
	    print $html "\t\t\t},\n";

	    #Fail columns
	    print $html "\t\t\t{\n";
	    print $html $chart_fail;    
	    print $html "\t\t\tdata : [";
	    my $fail_locus;
	    foreach my $s_exp (sort keys %h_experiments){
	    	$fail_locus .= $h_counts{$s_exp}{LOCUS}{TOTAL}{FAIL}.",";
	    }
	    $fail_locus =~ s/,$//;
	    print $html $fail_locus."],\n";
	    print $html "\t\t\t},\n";


      #DRBX columns
      print $html "\t\t\t{\n";
      print $html $chart_drbx;    
      print $html "\t\t\tdata : [";
      my $drbx_locus;
      foreach my $s_exp (sort keys %h_experiments){
        $drbx_locus .= $h_counts{$s_exp}{LOCUS}{TOTAL}{DRBX}.",";
      }
      $drbx_locus =~ s/,$//;
      print $html $drbx_locus."],\n";
      print $html "\t\t\t},\n";

	    #Error charts
	    print $html "\t\t\t{\n";   
		print $html $chart_error;        
	    print $html "\t\t\tdata : [";
	    my $error_locus;
	    foreach my $s_exp (sort keys %h_experiments){
	    	$error_locus .= $h_counts{$s_exp}{LOCUS}{TOTAL}{ERROR}.",";
	    }
	    $error_locus =~ s/,$//;
	    print $html $error_locus."],\n";
	    print $html "\t\t\t}\n";
	    print $html "\t\t]\n";
	    print $html "\t}\n";

	    ######## 
	    # Allele JS
	    ########
	   	print $html "var barChartData3 = {\n";
		  print $html "\t\tlabels : [";
		  print $html $label."],\n";
	    print $html "\t\tdatasets : [\n";
	   
	    #Pass columns
	    print $html "\t\t\t{\n";
	    print $html $chart_pass;
	    print $html "\t\t\tdata : [";
	    my $allele_chart;
	    foreach my $s_exp (sort keys %h_experiments){
	    	$allele_chart .= $h_counts{$s_exp}{ALLELE}{TOTAL}{PASS}.",";
	    }
	    $allele_chart =~ s/,$//;
	    print $html $allele_chart."],\n";
	    print $html "\t\t\t},\n";

	    #Fail columns
	    print $html "\t\t\t{\n";
	    print $html $chart_fail;    
	    print $html "\t\t\tdata : [";
	    my $fail_allele;
	    foreach my $s_exp (sort keys %h_experiments){
	    	$fail_allele .= $h_counts{$s_exp}{ALLELE}{TOTAL}{FAIL}.",";
	    }
	    $fail_allele =~ s/,$//;
	    print $html $fail_allele."],\n";
	    print $html "\t\t\t},\n";

      #DRBX columns
      print $html "\t\t\t{\n";
      print $html $chart_drbx;    
      print $html "\t\t\tdata : [";
      my $drbx_allele;
      foreach my $s_exp (sort keys %h_experiments){
        $drbx_allele .= $h_counts{$s_exp}{ALLELE}{TOTAL}{DRBX}.",";
      }
      $drbx_allele =~ s/,$//;
      print $html $drbx_allele."],\n";
      print $html "\t\t\t},\n";

	    #Error charts
	    print $html "\t\t\t{\n";   
		print $html $chart_error;        
	    print $html "\t\t\tdata : [";
	    my $error_allele;
	    foreach my $s_exp (sort keys %h_experiments){
	    	$error_allele .= $h_counts{$s_exp}{ALLELE}{TOTAL}{ERROR}.",";
	    }
	    $error_allele =~ s/,$//;
	    print $html $error_allele."],\n";
	    print $html "\t\t\t}\n";
	    print $html "\t\t]\n";
	    print $html "\t}\n";


		my $js = qq{
			window.onload = function(){

	        var ctx3 = document.getElementById("canvas").getContext("2d");
	        window.myBar = new Chart(ctx3).Bar(barChartData, {
	          responsive : true
	        });

	        var ctx5 = document.getElementById("canvas2").getContext("2d");
	        window.myBar = new Chart(ctx5).Bar(barChartData2, {
	          responsive : true
	        });

	        var ctx6 = document.getElementById("canvas3").getContext("2d");
	        window.myBar = new Chart(ctx6).Bar(barChartData3, {
	          responsive : true
	        });


	      };
	       </script>
		};
		print $html $js;
		
	

	my $footer = qq{
	<footer>
		<div class="footer">
			<img src="img/bethematch.jpeg" alt="img02">
      <p><a href="https://github.com/nmdp-bioinformatics/pipeline/blob/master/LICENSE">Copyright ©</a> 2014 - 2015 National Marrow Donor Program. All Rights Reserved.</p>		</div>
	</footer>	
	<style type="text/css">
	    .bargraph {
	      list-style: none;
	      padding-top: 20px;
	      width:560px;
	    }
	    footer {
				color: #888;
				clear:both;
				position: relative; 
				bottom: 10px;
				left: 0; width: 100%; 
				text-align: center;
				font-size:9px;
				margin-top: 300px;
				background: rgba(255, 255, 255, 0.6);
		}
			</style>

		    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.0.0/js/bootstrap.min.js"></script>
		    <script src="js/docs.min.js"></script>

  		</body>	
		};

	print $html $footer;

}
################################################################################################################
=head2 helpHtml

	Title:     helpHtml
	Usage:
	Function:  Prints the header to the log.html file
	Returns:
	Args:

=cut
sub helpHtml{

	my $html = shift;

	my $header = 
		qq{

<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="description" content="">
    <meta name="author" content="">
    
    <title>Validation Report</title>

    <!-- Bootstrap core CSS -->
    <link href="http://netdna.bootstrapcdn.com/bootstrap/3.0.0/css/bootstrap.min.css" rel="stylesheet">
    <script src="js/Chart.js"></script>


    <!-- Custom styles for this template -->
   <link href="css/dashboard.css" rel="stylesheet">
   <link rel="stylesheet" href="css/style.css">


    <style type='text/css'>

  .my-legend .legend-scale ul {
    margin: 0;
    padding: 0;
    float: left;
    list-style: none;
    }
  .my-legend .legend-scale ul li {
    display: block;
    float: left;
    width: 50px;
    margin-bottom: 6px;
    text-align: center;
    font-size: 80%;
    list-style: none;
    }
  .my-legend ul.legend-labels li span {
    display: block;
    float: left;
    height: 15px;
    width: 50px;
    }
  .my-legend .legend-source {
    font-size: 70%;
    color: #999;
    clear: both;
    }
  .my-legend a {
    color: #777;
    }

    #expanding{
      display:none;
    }


    </style>
  </head>

  <body>

    <nav class="navbar navbar-inverse navbar-fixed-top" role="navigation">
      <div class="container-fluid">
        <div class="navbar-header">
          <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
            <span class="sr-only">Toggle navigation</span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a class="navbar-brand" href="experiment.html">Validation Report</a>
        </div>
        <div id="navbar" class="navbar-collapse collapse">
          <ul class="nav navbar-nav navbar-right">
            <li><a href="experiment.html">Experiments</a></li>
            <li><a href="log.html">Log</a></li>
            <li><a href="help.html">Help</a></li>
          </ul>
        </div>
      </div>
    </nav>


    <div class="container-fluid">
    <div class="row">
      <div class="col-sm-3 col-md-2 sidebar">

      </div>

        </div>
         <div class="col-sm-9 col-sm-offset-3 col-md-10 col-md-offset-2 main">
          <h1 class="page-header">FAQs</h1>
          <div>
<ul id="mobile" class="cd-faq-group">
      <li>
        <a class="cd-faq-trigger  faq-format" href="#0">How do I run ngs-validation-report?</a>
        <div class="cd-faq-content">
          <p>Lorem ipsum dolor sit amet, consectetur adipisicing elit. Impedit quidem delectus rerum eligendi mollitia, repudiandae quae beatae. Et repellat quam atque corrupti iusto architecto impedit explicabo repudiandae qui similique aut iure ipsum quis inventore nulla error aliquid alias quia dolorem dolore, odio excepturi veniam odit veritatis. Quo iure magnam, et cum. Laudantium, eaque non? Tempore nihil corporis cumque dolor ipsum accusamus sapiente aliquid quis ea assumenda deserunt praesentium voluptatibus, accusantium a mollitia necessitatibus nostrum voluptatem numquam modi ab, sint rem.</p>
        </div> <!-- cd-faq-content -->
      </li>
      <li>
        <a class="cd-faq-trigger faq-format" href="#0">What does it mean for the typing to be designated as a fail?</a>
        <div class="cd-faq-content">
          <p>Lorem ipsum dolor sit amet, consectetur adipisicing elit. Blanditiis provident officiis, reprehenderit numquam. Praesentium veritatis eos tenetur magni debitis inventore fugit, magnam, reiciendis, saepe obcaecati ex vero quaerat distinctio velit.</p>
        </div> <!-- cd-faq-content -->
      </li>
      <li>
        <a class="cd-faq-trigger faq-format" href="#0">What does it mean for the typing to be designated as a dropout?</a>
        <div class="cd-faq-content">
          <p>Lorem ipsum dolor sit amet, consectetur adipisicing elit. Blanditiis provident officiis, reprehenderit numquam. Praesentium veritatis eos tenetur magni debitis inventore fugit, magnam, reiciendis, saepe obcaecati ex vero quaerat distinctio velit.</p>
        </div> <!-- cd-faq-content -->
      </li>
      <li>
        <a class="cd-faq-trigger faq-format" href="#0">What does it mean for the typing to be designated as a DRBX?</a>
        <div class="cd-faq-content">
          <p>Lorem ipsum dolor sit amet, consectetur adipisicing elit. Blanditiis provident officiis, reprehenderit numquam. Praesentium veritatis eos tenetur magni debitis inventore fugit, magnam, reiciendis, saepe obcaecati ex vero quaerat distinctio velit.</p>
        </div> <!-- cd-faq-content -->
      </li>
      <li>
        <a class="cd-faq-trigger faq-format" href="#0">Can a subject's DRB1 typing still pass if they have DRBX typing being called?</a>
        <div class="cd-faq-content">
          <p>Lorem ipsum dolor sit amet, consectetur adipisicing elit. Blanditiis provident officiis, reprehenderit numquam. Praesentium veritatis eos tenetur magni debitis inventore fugit, magnam, reiciendis, saepe obcaecati ex vero quaerat distinctio velit.</p>
        </div> <!-- cd-faq-content -->
      </li>
      <li>
        <a class="cd-faq-trigger faq-format" href="#0">How do I include my QC data?</a>
        <div class="cd-faq-content">
          <p>Lorem ipsum dolor sit amet, consectetur adipisicing elit. Nisi tempore, placeat quisquam rerum! Eligendi fugit dolorum tenetur modi fuga nisi rerum, autem officiis quaerat quos. Magni quia, quo quibusdam odio. Error magni aperiam amet architecto adipisci aspernatur! Officia, quaerat magni architecto nostrum magnam fuga nihil, ipsum laboriosam similique voluptatibus facilis nobis? Eius non asperiores, nesciunt suscipit veniam blanditiis veritatis provident possimus iusto voluptas, eveniet architecto quidem quos molestias, aperiam eum reprehenderit dolores ad deserunt eos amet. Vero molestiae commodi unde dolor dicta maxime alias, velit, nesciunt cum dolorem, ipsam soluta sint suscipit maiores mollitia assumenda ducimus aperiam neque enim! Quas culpa dolorum ipsam? Ipsum voluptatibus numquam natus? Eligendi explicabo eos, perferendis voluptatibus hic sed ipsam rerum maiores officia! Beatae, molestias!</p>
        </div> <!-- cd-faq-content -->
      </li>
      <li>
        <a class="cd-faq-trigger faq-format" href="#0">Why is there no log output?</a>
        <div class="cd-faq-content">
          <p>Lorem ipsum dolor sit amet, consectetur adipisicing elit. Blanditiis provident officiis, reprehenderit numquam. Praesentium veritatis eos tenetur magni debitis inventore fugit, magnam, reiciendis, saepe obcaecati ex vero quaerat distinctio velit.</p>
        </div> <!-- cd-faq-content -->
      </li>
      <li>
        <a class="cd-faq-trigger faq-format" href="#0">Why do the number of loci and alleles defer between experiments?</a>
        <div class="cd-faq-content">
          <p>Lorem ipsum dolor sit amet, consectetur adipisicing elit. Blanditiis provident officiis, reprehenderit numquam. Praesentium veritatis eos tenetur magni debitis inventore fugit, magnam, reiciendis, saepe obcaecati ex vero quaerat distinctio velit.</p>
        </div> <!-- cd-faq-content -->
      </li>

    </ul> <!-- cd-faq-group -->
            
          </div>
    </div>
    </div>
  
  <footer>
    <div class="footer">
      <img src="img/bethematch.jpeg" alt="img02">
      <p><a href="https://github.com/nmdp-bioinformatics/pipeline/blob/master/LICENSE">Copyright ©</a> 2014 - 2015 National Marrow Donor Program. All Rights Reserved.</p>    </div>
    </div>
  </footer> 
  <script src="js/main.js"></script> <!-- Resource jQuery -->
  <style type="text/css">
      .bargraph {
        list-style: none;
        padding-top: 20px;
        width:560px;
      }
      footer {
        color: #888;
        clear:both;
        position: relative; 
        bottom: 10px;
        left: 0; width: 100%; 
        text-align: center;
        font-size:9px;
        margin-top: 300px;
        background: rgba(255, 255, 255, 0.6);
    }
    .faq-format{
       text-decoration: none;
       color:black;

    }
    .cd-faq-group {
  /* hide group not selected */
  display: none;
  list-style-type: none;

}
.cd-faq-group.selected {
  display: block;
}
.cd-faq-group .cd-faq-title {
  background: transparent;
  box-shadow: none;

}


.no-js .cd-faq-group {
  display: block;

}
\@media only screen and (min-width: 768px) {
  .cd-faq-group {
    /* all groups visible */
    display: block;
     margin-left: 0px;
     width: 75%;
  }
  .cd-faq-group > li {
    background: #ffffff;
    margin-bottom: 6px;
    margin-left: -38px;
    box-shadow: 0 1px 2px rgba(0, 0, 0, 0.08);
    -webkit-transition: box-shadow 0.2s;
    -moz-transition: box-shadow 0.2s;
    transition: box-shadow 0.2s;

  }

  .no-touch .cd-faq-group > li:hover {
    box-shadow: 0 1px 10px rgba(108, 125, 142, 0.3);
  }
  .cd-faq-group .cd-faq-title {
    margin: 2em 0 1em;
  }
  .cd-faq-group:first-child .cd-faq-title {
    margin-top: 0;
  }
}

.cd-faq-trigger {
  position: relative;
  display: block;
  line-height: 1.2;
  text-decoration: none;
}
\@media only screen and (min-width: 768px) {
  
  .cd-faq-trigger {
    font-size: 24px;
    font-size: 1.5rem;

    font-weight: 300;
    padding: 24px 72px 24px 24px;
  }

  a.cd-faq-trigger {
    text-decoration: none;
    font-color:black;
  }
  .cd-faq-trigger::before, .cd-faq-trigger::after {
    /* arrow icon on the right */
    position: absolute;
    top: 50%;
    height: 2px;
    width: 13px;
    background: #cfdca0;
font-color:black;
    -webkit-backface-visibility: hidden;
    backface-visibility: hidden;
    -webkit-transition-property: -webkit-transform;
    -moz-transition-property: -moz-transform;
    transition-property: transform;
    -webkit-transition-duration: 0.2s;
    -moz-transition-duration: 0.2s;
    transition-duration: 0.2s;
  }
  .cd-faq-trigger::before {
    -webkit-transform: rotate(45deg);
    -moz-transform: rotate(45deg);
    -ms-transform: rotate(45deg);
    -o-transform: rotate(45deg);
    transform: rotate(45deg);
    
  }
  .cd-faq-trigger::after {
    -webkit-transform: rotate(-45deg);
    -moz-transform: rotate(-45deg);
    -ms-transform: rotate(-45deg);
    -o-transform: rotate(-45deg);
    transform: rotate(-45deg);
  }
  .content-visible .cd-faq-trigger::before {
    -webkit-transform: rotate(-45deg);
    -moz-transform: rotate(-45deg);
    -ms-transform: rotate(-45deg);
    -o-transform: rotate(-45deg);
    transform: rotate(-45deg);
  }
  .content-visible .cd-faq-trigger::after {
    -webkit-transform: rotate(45deg);
    -moz-transform: rotate(45deg);
    -ms-transform: rotate(45deg);
    -o-transform: rotate(45deg);
    transform: rotate(45deg);
  }
}

.cd-faq-content p {
  font-size: 14px;
  font-size: 1.2rem;
  line-height: 1.4;
  color: #6c7d8e;
}
\@media only screen and (min-width: 768px) {
  .cd-faq-content {
    display: none;
    padding: 0 24px 30px;
  }
  .cd-faq-content p {
    line-height: 1.6;
  }
  .no-js .cd-faq-content {
    display: block;
  }
}

      </style>
       <!-- <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script> -->
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.0.0/js/bootstrap.min.js"></script>
        <script src="js/docs.min.js"></script>

        <!-- IE10 viewport hack for Surface/desktop Windows 8 bug -->
        <script src="js/ie10-viewport-bug-workaround.js"></script>
      </body> 
    
		};

	print $html $header;

}
################################################################################################################
=head2 indexHeader

	Title:     indexHeader
	Usage:
	Function:  Prints the header to the index.html file
	Returns:
	Args:

=cut
sub indexHeader{

	my ($html, $s_exp, $rh_qc_verdict, $b_drbx, $b_error, $b_failed) = @_;

	my $s_experiment = $s_exp;
	my $header = 
		qq{
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="description" content="">
    <meta name="author" content="">
    
    <title>Validation Report</title>

    <!-- Bootstrap core CSS -->
    <link href="http://netdna.bootstrapcdn.com/bootstrap/3.0.0/css/bootstrap.min.css" rel="stylesheet">
    <script src="../js/Chart.js"></script>
    <!-- Custom styles for this template -->
    <link href="../css/dashboard.css" rel="stylesheet">

    <!--[if lt IE 9]>
      <script src="https://oss.maxcdn.com/html5shiv/3.7.2/html5shiv.min.js"></script>
      <script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
    <![endif]-->


      <style type='text/css'>
  .my-legend .legend-scale ul {
    margin: 0;
    padding: 0;
    float: left;
    list-style: none;
    }
  .my-legend .legend-scale ul li {
    display: block;
    float: left;
    width: 60px;
    margin-bottom: 6px;
    text-align: center;
    font-size: 80%;
    list-style: none;
    }
  .my-legend ul.legend-labels li span {
    display: block;
    float: left;
    height: 15px;
    width: 60px;
    }
  .my-legend .legend-source {
    font-size: 70%;
    color: #999;
    clear: both;
    }
  .my-legend a {
    color: #777;
    }

    #expanding{
      display:none;
    }


    </style>
  </head>

  <body>

    <nav class="navbar navbar-inverse navbar-fixed-top" role="navigation">
      <div class="container-fluid">
        <div class="navbar-header">
          <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
            <span class="sr-only">Toggle navigation</span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a class="navbar-brand" href="index.html">Validation Report - $s_experiment</a>
        </div>
        <div id="navbar" class="navbar-collapse collapse">
          <ul class="nav navbar-nav navbar-right">
            <li><a href="../experiment.html">Experiments</a></li>
            <li><a href="../log.html">Log</a></li>
            <li><a href="../help.html">Help</a></li>
          </ul>
        </div>
      </div>
    </nav>

    <div class="container-fluid">
      <div class="row">
        <div class="col-sm-3 col-md-2 sidebar">
          <ul class="nav nav-sidebar">
            <li class="active"><a href="#">Overview <span class="sr-only">(current)</span></a></li>
    };

  print $html $header;

  my $results_inactive = qq{
    <li ><a href="results.html">Results</a></li>
  };
  print $html $results_inactive;
  my $failed_inactive = qq{
    <li ><a href="fails.html">Failed</a></li>
  };
  print $html $failed_inactive if($b_failed);
  my $error_inactive  = qq{
    <li ><a href="errors.html">Dropout</a></li>
  };
  print $html $error_inactive if($b_error);
  my $drbx_inactive   = qq{
    <li ><a href="drbx.html">DRBX</a></li>
  };
  print $html $drbx_inactive if($b_drbx);
  print $html "<li><a href=\"qc.html\">QC</a></li>\n" if defined $$rh_qc_verdict{$s_exp};
    
    my $header2 = qq{      </ul>
        </div>
        <div class="col-sm-9 col-sm-offset-3 col-md-10 col-md-offset-2 main">
            <div class='my-legend' style="float:right;margin-right:20%;">
            <div class='legend-scale'>
              <ul class='legend-labels'>
                <li><span style='background:#46BFBD;'></span><b>PASS</b></li>
                <li><span style='background:#F7464A;'></span><b>FAIL</b></li>
                <li><span style='background:#D3D5D6;'></span><b>DROPOUT</b></li>
                <li><span style='background:#0d8fdb;'></span><b>DRBX</b></li>
              </ul>
            </div>
          </div>
          <h1 class="page-header" >Overview</h1>
          <div class="row placeholders">  

            <div class="col-xs-6 col-sm-3 placeholder">
              <div id="canvas-holder">
                <canvas id="chart-area3" width="300" height="300"/>
              </div>
               <h4>Subjects</h4>
              <span class="text-muted">Number of subjects that pass or fail</span>
            </div>
            <div class="col-xs-6 col-sm-3 placeholder" style="margin-left: 20px;">
              <div id="canvas-holder">
                <canvas id="chart-area" width="300" height="300"/>
              </div>
               <h4>Locus</h4>
              <span class="text-muted">Number of locus level pass and fails</span>
            </div>          
            <div class="col-xs-6 col-sm-3 placeholder">
              <div id="canvas-holder">
                <canvas id="chart-area2" width="300" height="300"/>
              </div>
               <h4>Allele</h4>
              <span class="text-muted">Number of allele level pass and fails</span>
            </div>            

          </div>

          <div class="row placeholders">  
                        <div class="col-xs-6 col-sm-3 placeholder" style="margin-left: 20px;">
              <div id="canvas-holder">
                <canvas id="canvas3" width="300" height="300"/>
              </div>
               <h4>Subjects</h4>
              <span class="text-muted">Number of subjects that pass or fail</span>
            </div>
            <div class="col-xs-6 col-sm-3 placeholder" style="margin-left: 20px;">
              <div id="canvas-holder">
                <canvas id="canvas" width="300" height="300" />
              </div>
               <h4>Locus</h4>
              <span class="text-muted">Number of allele level pass and fails by locus</span>
            </div>    
            <div class="col-xs-6 col-sm-3 placeholder" style="margin-left: 20px;">
              <div id="canvas-holder">
                <canvas id="canvas2" width="300" height="300"/>
              </div>
               <h4>Allele</h4>
              <span class="text-muted">Number of allele level pass and fails by Allele</span>
            </div>        

          </div>
			

      </div>
    </div>

		};
	print $html $header2;

}
################################################################################################################
=head2 indexFooter

	Title:     indexFooter
	Usage:
	Function:  Prints the footer to the index.html file
	Returns:
	Args:

=cut
sub indexFooter{

	my $html = shift;

	my $footer = 
		qq{
			<footer>
				<div class="footer">
					<img src="../img/bethematch.jpeg" alt="img02">
          <p><a href="https://github.com/nmdp-bioinformatics/pipeline/blob/master/LICENSE">Copyright ©</a> 2014 - 2015 National Marrow Donor Program. All Rights Reserved.</p>    </div>
				</div>
			</footer>	
	<style type="text/css">
    .bargraph {
      list-style: none;
      padding-top: 20px;
      width:560px;
    }
    footer {
			color: #888;
			clear:both;
			position: relative; 
			bottom: 10px;
			left: 0; width: 100%; 
			text-align: center;
			font-size:9px;
			
			background: rgba(255, 255, 255, 0.6);
	}

	</style>
	    <!-- <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script> -->
	    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.0.0/js/bootstrap.min.js"></script>
	    <script src="../js/docs.min.js"></script>

  </body>		
		};
	print $html $footer;
	close $html;

}
################################################################################################################
=head2 logFooter

	Title:     logFooter
	Usage:
	Function:  Prints the footer to the log.html file
	Returns:
	Args:

=cut
sub logFooter{

	my $html = shift;


	my $footer = 
		qq{
			</div>
			</div>
			<footer>
				<div class="footer">
					<img src="img/bethematch.jpeg" alt="img02">
          <p><a href="https://github.com/nmdp-bioinformatics/pipeline/blob/master/LICENSE">Copyright ©</a> 2014 - 2015 National Marrow Donor Program. All Rights Reserved.</p>    </div>
				</div>
			</footer>	
	<style type="text/css">
    .bargraph {
      list-style: none;
      padding-top: 20px;
      width:560px;
    }
    footer {
			color: #888;
			clear:both;
			position: relative; 
			bottom: 10px;
			left: 0; width: 100%; 
			text-align: center;
			font-size:9px;
			margin-top: 300px;
			background: rgba(255, 255, 255, 0.6);
	}

	</style>
	    <!-- <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script> -->
	    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.0.0/js/bootstrap.min.js"></script>
	    <script src="js/docs.min.js"></script>
	    <!-- IE10 viewport hack for Surface/desktop Windows 8 bug -->
	    <script src="js/ie10-viewport-bug-workaround.js"></script>
  
  </body>		
		};
	print $html $footer;
	close $html;

}
################################################################################################################
=head2 indexBody

	Title:     indexBody
	Usage:
	Function:  Prints the body to the index.html file
	Returns:
	Args:

=cut
sub indexBody{

	my($html, $ra_counts, $ra_sensitivity) = @_;

	my $pie_data = qq{
	<script>
		var pieData = [
        {
          value: $$ra_counts{LOCUS}{TOTAL}{PASS},
          color: "#46BFBD",
          highlight: "#5AD3D1",
          label: "Pass"
        },
        {
          value: $$ra_counts{LOCUS}{TOTAL}{FAIL},
          color:"#F7464A",
          highlight: "#FF5A5E",
          label: "Fail"
        },
        {
          value: $$ra_counts{LOCUS}{TOTAL}{ERROR},
          color: "#D3D5D6",
          highlight: "#DCDDDE",
          label: "Dropout"
        },
        {
          value: $$ra_counts{LOCUS}{TOTAL}{DRBX},
          color: "#0d8fdb",
          highlight: "#3DA5E2",
          label: "DRBX"
        }


      ];
		var pieData2 = [
       
        {
          value: $$ra_counts{ALLELE}{TOTAL}{PASS},
          color: "#46BFBD",
          highlight: "#5AD3D1",
          label: "Pass"
        },
        {
          value: $$ra_counts{ALLELE}{TOTAL}{FAIL},
          color:"#F7464A",
          highlight: "#FF5A5E",
          label: "Fail"
        },
        {
          value: $$ra_counts{ALLELE}{TOTAL}{ERROR},
          color: "#D3D5D6",
          highlight: "#DCDDDE",
          label: "Dropout"
        },
        {
          value: $$ra_counts{ALLELE}{TOTAL}{DRBX},
          color: "#0d8fdb",
          highlight: "#3DA5E2",
          label: "DRBX"
        }

      ];

		var pieData3 = [
        {
          value: $$ra_counts{SUBJECT}{TOTAL}{PASS},
          color: "#46BFBD",
          highlight: "#5AD3D1",
          label: "Pass"
        },
        {
          value: $$ra_counts{SUBJECT}{TOTAL}{FAIL},
          color:"#F7464A",
          highlight: "#FF5A5E",
          label: "Fail"
        },
        {
          value: $$ra_counts{SUBJECT}{TOTAL}{ERROR},
          color: "#D3D5D6",
          highlight: "#DCDDDE",
          label: "Dropout"
        },
        {
          value: $$ra_counts{SUBJECT}{TOTAL}{DRBX},
          color: "#0d8fdb",
          highlight: "#3DA5E2",
          label: "DRBX"
        }

      ];

	};
	print $html $pie_data;

	my $bar_chart_data = qq{
		      var barChartData = {
        labels : ["A","B","C","DRB1","DQB1"],
        datasets : [
          {
            fillColor : "#46BFBD",
            strokeColor : "#46BFBD",
            highlightFill :  "#5AD3D1",
            highlightStroke :  "#5AD3D1",
            data : [$$ra_counts{LOCUS}{"A"}{PASS},$$ra_counts{LOCUS}{"B"}{PASS},$$ra_counts{LOCUS}{"C"}{PASS},$$ra_counts{LOCUS}{"DRB1"}{PASS},$$ra_counts{LOCUS}{"DQB1"}{PASS}]
          },
            {
            fillColor : "#F7464A",
            strokeColor : "#F7464A",
            highlightFill: "#FF5A5E",
            highlightStroke: "#FF5A5E",
            data : [$$ra_counts{LOCUS}{"A"}{FAIL},$$ra_counts{LOCUS}{"B"}{FAIL},$$ra_counts{LOCUS}{"C"}{FAIL},$$ra_counts{LOCUS}{"DRB1"}{FAIL},$$ra_counts{LOCUS}{"DQB1"}{FAIL}]
          },
          {
            fillColor : "#D3D5D6",
            strokeColor : "#D3D5D6",
            highlightFill :  "#DCDDDE",
            highlightStroke :  "#DCDDDE",
            data : [$$ra_counts{LOCUS}{"A"}{ERROR},$$ra_counts{LOCUS}{"B"}{ERROR},$$ra_counts{LOCUS}{"C"}{ERROR},$$ra_counts{LOCUS}{"DRB1"}{ERROR},$$ra_counts{LOCUS}{"DQB1"}{ERROR}]
          },
          {
            fillColor : "#0d8fdb",
            strokeColor : "#0d8fdb",
            highlightFill :  "#3DA5E2",
            highlightStroke :  "#3DA5E2",
            data : [0,0,0,$$ra_counts{LOCUS}{"DRB1"}{DRBX},0]
          }      
        ]
      }

      var barChartData2 = {
        labels : ["A","B","C","DRB1","DQB1"],
        datasets : [
          
          {
            fillColor : "#46BFBD",
            strokeColor : "#46BFBD",
            highlightFill :  "#5AD3D1",
            highlightStroke :  "#5AD3D1",
            data : [$$ra_counts{ALLELE}{"A"}{PASS},$$ra_counts{ALLELE}{"B"}{PASS},$$ra_counts{ALLELE}{"C"}{PASS},$$ra_counts{ALLELE}{"DRB1"}{PASS},$$ra_counts{ALLELE}{"DQB1"}{PASS}]
          },
          {
            fillColor : "#F7464A",
            strokeColor : "#F7464A",
            highlightFill: "#FF5A5E",
            highlightStroke: "#FF5A5E",
            data : [$$ra_counts{ALLELE}{"A"}{FAIL},$$ra_counts{ALLELE}{"B"}{FAIL},$$ra_counts{ALLELE}{"C"}{FAIL},$$ra_counts{ALLELE}{"DRB1"}{FAIL},$$ra_counts{ALLELE}{"DQB1"}{FAIL}]
          },
          {
            fillColor : "#D3D5D6",
            strokeColor : "#D3D5D6",
            highlightFill: "#DCDDDE",
            highlightStroke: "#DCDDDE",
            data : [$$ra_counts{ALLELE}{"A"}{ERROR},$$ra_counts{ALLELE}{"B"}{ERROR},$$ra_counts{ALLELE}{"C"}{ERROR},$$ra_counts{ALLELE}{"DRB1"}{ERROR},$$ra_counts{ALLELE}{"DQB1"}{ERROR}]
          },
          {
            fillColor : "#0D8FDB",
            strokeColor : "#0D8FDB",
            highlightFill: "#3DA5E2",
            highlightStroke: "#3DA5E2",
            data : [0,0,0,$$ra_counts{ALLELE}{"DRB1"}{DRBX},0]
          }

        ]
      }

      var barChartData3 = {
        labels : ["Subject"],
        datasets : [
          {
            fillColor : "#46BFBD",
            strokeColor : "#46BFBD",
            highlightFill :  "#5AD3D1",
            highlightStroke :  "#5AD3D1",
            data : [$$ra_counts{SUBJECT}{TOTAL}{PASS}]
          },
          {
            fillColor : "#F7464A",
            strokeColor : "#F7464A",
            highlightFill: "#FF5A5E",
            highlightStroke: "#FF5A5E",
            data : [$$ra_counts{SUBJECT}{TOTAL}{FAIL}]
          },
          {
            fillColor : "#D3D5D6",
            strokeColor : "#D3D5D6",
            highlightFill: "#DCDDDE",
            highlightStroke: "#DCDDDE",
            data : [$$ra_counts{SUBJECT}{TOTAL}{ERROR}]
          },
          {
            fillColor : "#0d8fdb",
            strokeColor : "#0d8fdb",
            highlightFill: "#3DA5E2",
            highlightStroke: "#3DA5E2",
            data : [$$ra_counts{SUBJECT}{TOTAL}{DRBX}]
          }
        ]
      }




	};
	print $html $bar_chart_data;

	my $js = qq{
		window.onload = function(){
        var ctx = document.getElementById("chart-area").getContext("2d");
        window.myPie = new Chart(ctx).Pie(pieData, {
          responsive : true
        });
        
        var ctx2 = document.getElementById("chart-area2").getContext("2d");
        window.myPie = new Chart(ctx2).Pie(pieData2, {
          responsive : true
        });

        var ctx4 = document.getElementById("chart-area3").getContext("2d");
        window.myPie = new Chart(ctx4).Pie(pieData3, {
          responsive : true
        });

        var ctx3 = document.getElementById("canvas").getContext("2d");
        window.myBar = new Chart(ctx3).Bar(barChartData, {
          responsive : true
        });

        var ctx5 = document.getElementById("canvas2").getContext("2d");
        window.myBar = new Chart(ctx5).Bar(barChartData2, {
          responsive : true
        });

        var ctx6 = document.getElementById("canvas3").getContext("2d");
        window.myBar = new Chart(ctx6).Bar(barChartData3, {
          responsive : true
        });

        var ctx7 = document.getElementById("canvas4").getContext("2d");
        window.myBar = new Chart(ctx7).Bar(barChartData4, {
          responsive : true
        });

        var ctx8 = document.getElementById("canvas5").getContext("2d");
        window.myBar = new Chart(ctx8).Bar(barChartData5, {
          responsive : true
        });

      };
       </script>
	};
	print $html $js;


}
################################################################################################################
=head2 resultsBody

	Title:     resultsBody
	Usage:
	Function:  Prints the body to the results.html file
	Returns:
	Args:

=cut
sub resultsBody{

	my ($html, $s_exp, $rh_locus_verdict, $rh_qc_verdict, $rh_qc_expected, $rh_verified, $rh_observed, $ra_subject_pages, $rh_subject_errors, $rh_consensus ) = @_;

	my $s_experiment = $s_exp;

	my %h_subject_errors;
	foreach my $s_ID (keys %$rh_verified){
		foreach my $s_loc (keys %{$$rh_verified{$s_ID}}){

			foreach my $s_typing (keys %{$$rh_verified{$s_ID}{$s_loc}}){

				print $html "\t\t<tr>\n";

				my ($n_cutoff1,$n_cutoff2) = $s_loc =~ /^[A|B|C]/ ? (5,16) : (4,13);
				
				if(!defined $$ra_subject_pages{$s_ID}){
					print STDERR "Subject has no valid loci defined!! $s_ID\n" if !defined $rh_subject_errors || !defined $h_subject_errors{$s_ID};
					$h_subject_errors{$s_ID}++;
					next;
				}
                              
        my $s_checkbox = qq{
                  <td >
                    <div style="margin-left:4px;margin-top:10px;"><input type = "checkbox"
                       id = "myCheck"
                       value = "ham" />
                      </div>
                   </td>
                  };
        print $html $s_checkbox;

				my $id_link = "subjects/subject".$$ra_subject_pages{$s_ID}.".html";
				print $html  "\t\t\t<td><a href=\"$id_link\" class=\"id_links\">$s_ID</a></td>\n";
				
				print $html "\t\t\t<td>$s_loc</td>\n";
				if($$rh_locus_verdict{$s_ID}{$s_loc} eq "FAIL"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>FAIL</b></td>\n";
				}elsif($$rh_locus_verdict{$s_ID}{$s_loc} eq "ERROR"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>DROPOUT</b></td>\n";
				}elsif($$rh_locus_verdict{$s_ID}{$s_loc} eq "DRBX"){
          print $html "\t\t\t<td style=\"color:#F7464A\"><b>DRBX</b></td>\n";
        }else{
					print $html  "\t\t\t<td>PASS</td>\n";
				}

				
				if($$rh_verified{$s_ID}{$s_loc}{$s_typing} eq "FAIL" && $$rh_locus_verdict{$s_ID}{$s_loc} eq "ERROR"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>DROPOUT</b></td>\n";
				}elsif($$rh_verified{$s_ID}{$s_loc}{$s_typing} eq "FAIL"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>FAIL</b></td>\n";
				}elsif($$rh_verified{$s_ID}{$s_loc}{$s_typing} eq "DRBX"){
          print $html "\t\t\t<td style=\"color:#F7464A\"><b>DRBX</b></td>\n";
        }else{
					print $html "\t\t\t<td>$$rh_verified{$s_ID}{$s_loc}{$s_typing}</td>\n";
				}

        if(defined $$rh_consensus{$s_ID}{$s_loc}){
				  print $html "\t\t\t<td>$$rh_consensus{$s_ID}{$s_loc}</td>\n";
        }else{
          print $html "\t\t\t<td></td>\n";
        }

				print $html  "\t\t\t<td>$s_typing</td>\n";

				if(defined @{$$rh_observed{$s_ID}{$s_loc}}[0]){
          for(0..$#{$$rh_observed{$s_ID}{$s_loc}}){
            my $observed_gl = @{$$rh_observed{$s_ID}{$s_loc}}[$_];
            printObservedAlleles($html,$observed_gl,$s_typing,$n_cutoff2,$n_cutoff1);
          }

				}else{
					print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
					print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
				}

				print $html "</tr>\n";


			}
		}
	}

		my $end_body = qq{</tbody>
	            </table>
	          </div>
	        </div>
	      </div>
	    </div>
		};
		print $html $end_body,"\n";


}

################################################################################################################
=head2 printObservedAlleles

  Title:     printObservedAlleles
  Usage:
  Function:  Prints the body to the failed.html file
  Returns:
  Args:

=cut
sub printObservedAlleles{

    my( $html, $observed_gl, $s_typing, $n_cutoff2, $n_cutoff1) = @_;


    my %h_alleles = map{ $_ => 1 } split(/\//,$observed_gl);
    my $num_alleles = keys %h_alleles;

    if($num_alleles > 9){

      my @a_alleles = keys %h_alleles;
      for(my $i=1;$i<=$#a_alleles+1;$i++){
        print $html "<td><span  onclick=\"expand(this)\">\n" if $i == 1;
        my $s_end = $i == $#a_alleles+1 ? "" : "/"; 

        if($i == $n_cutoff2){

          if(g2p($a_alleles[$i-1]) eq g2p($s_typing)){
            print $html "...<br><span id=\"expanding\"><b>$a_alleles[$i-1]".$s_end."</b>";
          }else{
            print $html "...<br><span id=\"expanding\">$a_alleles[$i-1]".$s_end;
          }
        }elsif($i % $n_cutoff1 == 0) {
          if(g2p($a_alleles[$i-1]) eq g2p($s_typing)){
            print $html "<b>$a_alleles[$i-1]".$s_end."</b><br>";
          }else{
            print $html "$a_alleles[$i-1]".$s_end."<br>\n";
          }
        }else{
          if($a_alleles[$i-1] eq $s_typing){
            print $html "<b>$a_alleles[$i-1]".$s_end."</b>";
          }else{
            print $html "$a_alleles[$i-1]".$s_end;
          }
        }

      }

      print $html "</span></span></td>\n";
    }else{
      my @a_alleles = keys %h_alleles;
      print $html "<td>";
      for(my $i=1;$i<=$#a_alleles+1;$i++){
        my $s_end = $i == $#a_alleles+1 ? "" : "/"; 

        if($i % $n_cutoff1 == 0) {
          if(g2p($a_alleles[$i-1]) eq g2p($s_typing)){
            print $html "<b>$a_alleles[$i-1]".$s_end."</b><br>\n";
          }else{
            print $html "$a_alleles[$i-1]".$s_end."<br>\n";
          }
        }else{
          if(g2p($a_alleles[$i-1]) eq g2p($s_typing)){
            print $html "<b>$a_alleles[$i-1]".$s_end."</b>";
          }else{
            print $html "$a_alleles[$i-1]".$s_end;
          }
        }

      }
      print $html "</td>\n";
    }


}

################################################################################################################
=head2 failedBody

	Title:     failedBody
	Usage:
	Function:  Prints the body to the failed.html file
	Returns:
	Args:

=cut
sub failedBody{

	my ($html, $s_exp, $rh_locus_verdict, $rh_qc_verdict, $rh_qc_expected,  $rh_verified, $rh_observed, $ra_subject_pages, $ra_fails, $rh_consensus ) = @_;


	my $s_experiment = $s_exp;

	foreach my $s_ID (keys %$rh_verified){
		foreach my $s_loc (keys %{$$rh_verified{$s_ID}}){
			next unless defined $$ra_fails{$s_ID}{$s_loc};

			foreach my $s_typing (keys %{$$rh_verified{$s_ID}{$s_loc}}){

				print $html "\t\t<tr>\n";

				my ($n_cutoff1,$n_cutoff2) = $s_loc =~ /^[A|B|C]/ ? (5,16) : (4,13);
				
        my $s_checkbox = qq{
                  <td >
                    <div style="margin-left:4px;margin-top:10px;"><input type = "checkbox"
                       id = "myCheck"
                       value = "ham" />
                      </div>
                   </td>
                  };
        print $html $s_checkbox;

				my $id_link = "subjects/subject".$$ra_subject_pages{$s_ID}.".html";
				print $html  "\t\t\t<td><a href=\"$id_link\" class=\"id_links\">$s_ID</a></td>\n";

				print $html "\t\t\t<td>$s_loc</td>\n";
				if(defined $$rh_locus_verdict{$s_ID}{$s_loc} && $$rh_locus_verdict{$s_ID}{$s_loc} eq "FAIL"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>FAIL</b></td>\n";
				}elsif(!defined $$rh_locus_verdict{$s_ID}{$s_loc} || $$rh_locus_verdict{$s_ID}{$s_loc} eq "ERROR"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>DROPOUT</b></td>\n";
				}else{
					print $html  "\t\t\t<td>$$rh_locus_verdict{$s_ID}{$s_loc}</td>\n";
				}

				if($$rh_verified{$s_ID}{$s_loc}{$s_typing} eq "FAIL"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>FAIL</b></td>\n";
				}else{
					print $html  "\t\t\t<td>$$rh_verified{$s_ID}{$s_loc}{$s_typing}</td>\n";
				}

        if(defined $$rh_consensus{$s_ID}{$s_loc}){
          print $html "\t\t\t<td>$$rh_consensus{$s_ID}{$s_loc}</td>\n";
        }else{
          print $html "\t\t\t<td></td>\n";
        }

				print $html  "\t\t\t<td>$s_typing</td>\n";
				
				
				if(defined @{$$rh_observed{$s_ID}{$s_loc}}[0]){
	
          for(0..$#{$$rh_observed{$s_ID}{$s_loc}}){
            my $observed_gl = @{$$rh_observed{$s_ID}{$s_loc}}[$_];
            printObservedAlleles($html,$observed_gl,$s_typing,$n_cutoff2,$n_cutoff1);
          }

				}else{
					print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
					print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
				}

				print $html "</tr>\n";


			}
		}
	}

		my $end_body = qq{</tbody>
	            </table>
	          </div>
	        </div>
	      </div>
	    </div>
		};
		print $html $end_body,"\n";


}

################################################################################################################
=head2 drbxBody

  Title:     drbxBody
  Usage:
  Function:  Prints the body to the drbx.html file
  Returns:
  Args:

=cut
sub drbxBody{

  my ($html, $s_exp, $rh_locus_verdict, $rh_qc_verdict, $rh_qc_expected,  $rh_verified, $rh_observed, $ra_subject_pages, $ra_drbx, $rh_consensus ) = @_;


  my $s_experiment = $s_exp;

  foreach my $s_ID (keys %$rh_verified){
    foreach my $s_loc (keys %{$$rh_verified{$s_ID}}){
      next unless defined $$ra_drbx{$s_ID}{$s_loc};

      foreach my $s_typing (keys %{$$rh_verified{$s_ID}{$s_loc}}){

        print $html "\t\t<tr>\n";

        my ($n_cutoff1,$n_cutoff2) = $s_loc =~ /^[A|B|C]/ ? (5,16) : (4,13);
        
        my $s_checkbox = qq{
                  <td >
                    <div style="margin-left:4px;margin-top:10px;"><input type = "checkbox"
                       id = "myCheck"
                       value = "ham" />
                      </div>
                   </td>
                  };
        print $html $s_checkbox;

        my $id_link = "subjects/subject".$$ra_subject_pages{$s_ID}.".html";
        print $html  "\t\t\t<td><a href=\"$id_link\" class=\"id_links\">$s_ID</a></td>\n";

        print $html "\t\t\t<td>$s_loc</td>\n";
        if(defined $$rh_locus_verdict{$s_ID}{$s_loc} && $$rh_locus_verdict{$s_ID}{$s_loc} eq "FAIL"){
          print $html "\t\t\t<td style=\"color:#F7464A\"><b>FAIL</b></td>\n";
        }elsif(!defined $$rh_locus_verdict{$s_ID}{$s_loc} || $$rh_locus_verdict{$s_ID}{$s_loc} eq "ERROR"){
          print $html "\t\t\t<td style=\"color:#F7464A\"><b>DROPOUT</b></td>\n";
        }elsif(defined $$rh_locus_verdict{$s_ID}{$s_loc} && $$rh_locus_verdict{$s_ID}{$s_loc} eq "DRBX"){
          print $html "\t\t\t<td style=\"color:#F7464A\"><b>DRBX</b></td>\n";
        }else{
          print $html  "\t\t\t<td>$$rh_locus_verdict{$s_ID}{$s_loc}</td>\n";
        }

        if($$rh_verified{$s_ID}{$s_loc}{$s_typing} eq "FAIL"){
          print $html "\t\t\t<td style=\"color:#F7464A\"><b>FAIL</b></td>\n";
        }elsif($$rh_verified{$s_ID}{$s_loc}{$s_typing} eq "DRBX"){
          print $html "\t\t\t<td style=\"color:#F7464A\"><b>DRBX</b></td>\n";
        }else{
          print $html  "\t\t\t<td>$$rh_verified{$s_ID}{$s_loc}{$s_typing}</td>\n";
        }

        if(defined $$rh_consensus{$s_ID}{$s_loc}){
          print $html "\t\t\t<td>$$rh_consensus{$s_ID}{$s_loc}</td>\n";
        }else{
          print $html "\t\t\t<td></td>\n";
        }

        print $html  "\t\t\t<td>$s_typing</td>\n";
        
    
        if(defined @{$$rh_observed{$s_ID}{$s_loc}}[0]){
          for(0..$#{$$rh_observed{$s_ID}{$s_loc}}){
            my $observed_gl = @{$$rh_observed{$s_ID}{$s_loc}}[$_];
            printObservedAlleles($html,$observed_gl,$s_typing,$n_cutoff2,$n_cutoff1);
          }
        }else{
          print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
          print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
        }

        print $html "</tr>\n";


      }
    }
  }

    my $end_body = qq{</tbody>
              </table>
            </div>
          </div>
        </div>
      </div>
    };
    print $html $end_body,"\n";


}

################################################################################################################
=head2 failedHeader

	Title:     failedHeader
	Usage:
	Function:  Prints the header to the failed.html file
	Returns:
	Args:

=cut
sub qcHeader{
  my ( $html, $s_exp, $n_max_observed, $b_drbx, $b_error, $b_failed ) = @_;

  my $s_experiment = $s_exp;

  my $header = qq{
    <!DOCTYPE html>
  <html lang="en">
    <head>
         <meta charset="utf-8">
      <meta http-equiv="X-UA-Compatible" content="IE=edge">
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <meta name="description" content="">
      <meta name="author" content="">
      
      <title>Validation Report</title>
      <!-- Bootstrap core CSS -->
      <link href="http://netdna.bootstrapcdn.com/bootstrap/3.0.0/css/bootstrap.min.css" rel="stylesheet">
  
      <!-- Custom styles for this template -->
      <link   href="../css/dashboard.css" rel="stylesheet">
      <script type="text/javascript" src="../js/jquery-latest.min.js"></script>
      <script type="text/javascript" src="../js/jquery-ui.min.js"></script>
      <script type="text/javascript" src="../js/docs.min.js"></script>
      <script type="text/javascript" src="../js/Chart.js"></script>
 
      <!-- Tablesorter: theme -->
      <link class="theme default" rel="stylesheet" href="../css/theme.bootstrap.css">

      <!-- Tablesorter script: required -->
      <script src="http://mottie.github.io/tablesorter/js/jquery.tablesorter.js"></script>
      <script src="http://mottie.github.io/tablesorter/js/jquery.tablesorter.widgets.js"></script>
      <script src="http://mottie.github.io/tablesorter/addons/pager/jquery.tablesorter.pager.js"></script>

  <script id=\"js\">\$(function(){

    \$(\'#table1\').tablesorter({
      widthFixed : true,
      showProcessing: true,
      headerTemplate : '{content} {icon}', // Add icon for various themes

      widgets: [ 'stickyHeaders', 'filter' ],

      widgetOptions: {

        // extra class name added to the sticky header row
        stickyHeaders : '',
        // number or jquery selector targeting the position:fixed element
        stickyHeaders_offset : 50,
        // added to table ID, if it exists
        stickyHeaders_cloneId : '-sticky',
        // trigger "resize" event on headers
        stickyHeaders_addResizeEvent : false,
        // if false and a caption exist, it won't be included in the sticky header
        stickyHeaders_includeCaption : false,
        // The zIndex of the stickyHeaders, allows the user to adjust this to their needs
        stickyHeaders_zIndex : 2,
        // jQuery selector or object to attach sticky header to
        stickyHeaders_attachTo : null,
        // jQuery selector or object to monitor horizontal scroll position (defaults: xScroll > attachTo > window)
        stickyHeaders_xScroll : null,
        // jQuery selector or object to monitor vertical scroll position (defaults: yScroll > attachTo > window)
        stickyHeaders_yScroll : null,
        // scroll table top into view after filtering
        stickyHeaders_filteredToTop: true,

      }
    });
    \$(\'input[type=\"checkbox\"][value=\"change\"]\').change(function() {
         if(this.checked) {
            var checkboxes = \$(\'body\').find(\'input\');
            for(var i =0;i<checkboxes.length;i++){
              if(checkboxes[i].value !== "change"){
                checkboxes[i].checked = true;
              }
            }
         }
         if(!this.checked) {
            var checkboxes = \$(\'body\').find(\'input\');
            for(var i =0;i<checkboxes.length;i++){
              if(checkboxes[i].value !== "change"){
                checkboxes[i].checked = false;
              }
            }
         }

     });
    
});</script>
        <style type='text/css'>
    .my-legend .legend-scale ul {
      margin: 0;
      padding: 0;
      float: left;
      list-style: none;
      }
    .my-legend .legend-scale ul li {
      display: block;
      float: left;
      width: 50px;
      margin-bottom: 6px;
      text-align: center;
      font-size: 80%;
      list-style: none;
      }
    .my-legend ul.legend-labels li span {
      display: block;
      float: left;
      height: 15px;
      width: 50px;
      }
    .my-legend .legend-source {
      font-size: 70%;
      color: #999;
      clear: both;
      }
    .my-legend a {
      color: #777;
      }

      #expanding{
        display:none;
      }


      </style>
};
print $html $header;


my $header2 = qq{
    </head>

    <body>

      <nav class="navbar navbar-inverse navbar-fixed-top" role="navigation">
        <div class="container-fluid">
          <div class="navbar-header">
            <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
              <span class="sr-only">Toggle navigation</span>
              <span class="icon-bar"></span>
              <span class="icon-bar"></span>
              <span class="icon-bar"></span>
            </button>
            <a class="navbar-brand" href="../experiment.html">Validation Report - $s_experiment </a>
          </div>
          <div id="navbar" class="navbar-collapse collapse">
            <ul class="nav navbar-nav navbar-right">
              <li><a href="../experiment.html">Experiments</a></li>
              <li><a href="../log.html">Log</a></li>
              <li><a href="../help.html">Help</a></li>
            </ul>
          </div>
        </div>
      </nav>

      <div class="container-fluid">
        <div class="row">
          <div class="col-sm-3 col-md-2 sidebar">
            <ul class="nav nav-sidebar">
              <li><a href="index.html">Overview</a></li>
              <li><a href="results.html">Results</a></li>
    };

  print $html $header2;

  my $failed_inactive = qq{
    <li ><a href="fails.html">Failed</a></li>
  };
  print $html $failed_inactive if($b_failed);
  my $error_inactive  = qq{
    <li ><a href="errors.html">Dropout</a></li>
  };
  print $html $error_inactive if($b_error);
  my $drbx_inactive   = qq{
    <li ><a href="drbx.html">DRBX</a></li>
  };
  print $html $drbx_inactive if($b_drbx);
  my $qc_active = qq{
          <li  class="active"><a href="qc.html">QC<span class="sr-only">(current)</span></a></li>
  };
  print $html $qc_active;

  my $header3 = qq{
            </ul>
          </div>
          <div class="col-sm-9 col-sm-offset-3 col-md-10 col-md-offset-2 main">
          <div style="float:right;margin-top:20px;">
                <a class="hovering" >
              <span  onclick="Addfilter()" style="font-weight:bold;">Filter Rows</span>         
            </a>
            <span  style="font-weight:bold;"> | </span>
                <a class="hovering" >
                  <span  onclick="downloadcsv()" style="font-weight:bold;">Download to CSV</span>
                </a>
            </div>

            <h1 class="page-header">Validation Results</h1>
            <div class="table-responsive">
               <table class="table table-striped tablesorter-bootstrap" id="table1">
                <thead>
                  <tr>
                  <td style="border-color:#dddddd;border-style:solid;border-left-width:0px;border-right-width:0px;border-bottom-width:2px;" class="filter-false sorter-false">
                      <div style="margin-top:17px;"><input type = "checkbox" value="change"/></div>
                  </td>                  
                    <th>Sample ID</th>
                    <th>QC Verdict</th>
                    <th>Locus Verdict</th>
                    <th>Allele Verdict</th>
                    <th>QC Data</th>
                    <th>Expected Allele</th>
          };
          print $html $header3;
          for(1..$n_max_observed){
            print $html "<th>Observed Allele Call ",$_,"</th>\n";
          }
          my $header4 = qq{
                  </tr>
                </thead>
                <tbody>
          };
          print $html $header4;

}
################################################################################################################
=head2 failedBody

	Title:     failedBody
	Usage:
	Function:  Prints the body to the failed.html file
	Returns:
	Args:

=cut
sub qcBody{

	my ($html, $s_exp, $rh_locus_verdict, $rh_qc_verdict, $rh_qc_expected,  $rh_verified, $rh_observed, $ra_subject_pages, $ra_fails ) = @_;


	my $s_experiment = $s_exp;

	foreach my $s_ID (keys %$rh_verified){
		foreach my $s_loc (keys %{$$rh_verified{$s_ID}}){

			foreach my $s_typing (keys %{$$rh_verified{$s_ID}{$s_loc}}){

        next if !defined $$rh_qc_expected{$s_ID}{$s_loc} || $$rh_qc_expected{$s_ID}{$s_loc} !~ /\S/;

				print $html "\t\t<tr>\n";

                              
        my $s_checkbox = qq{
                  <td >
                    <div style="margin-left:4px;margin-top:10px;"><input type = "checkbox"
                       id = "myCheck"
                       value = "ham" />
                      </div>
                   </td>
                  };
        print $html $s_checkbox;

				my ($n_cutoff1,$n_cutoff2) = $s_loc =~ /^[A|B|C]/ ? (5,16) : (4,13);
				
				my $id_link = "subjects/subject".$$ra_subject_pages{$s_ID}.".html";
				print $html  "\t\t\t<td><a href=\"$id_link\" class=\"id_links\">$s_ID</a></td>\n";

				if(defined $$rh_qc_verdict{$s_ID}{$s_loc} && $$rh_qc_verdict{$s_ID}{$s_loc} eq "FAIL"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>FAIL</b></td>\n";
				}elsif(defined $$rh_qc_verdict{$s_ID}{$s_loc} && $$rh_qc_verdict{$s_ID}{$s_loc} eq "ERROR"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>ERROR</b></td>\n";
				}elsif(defined $$rh_qc_verdict{$s_ID}{$s_loc} && $$rh_qc_verdict{$s_ID}{$s_loc} eq "PASS"){
					print $html  "\t\t\t<td>PASS</td>\n";
				}else{
					print $html  "\t\t\t<td></td>\n";
				}

				if(defined $$rh_locus_verdict{$s_ID}{$s_loc} && $$rh_locus_verdict{$s_ID}{$s_loc} eq "FAIL"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>FAIL</b></td>\n";
				}elsif(!defined $$rh_locus_verdict{$s_ID}{$s_loc} || $$rh_locus_verdict{$s_ID}{$s_loc} eq "ERROR"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>ERROR</b></td>\n";
				}else{
					print $html  "\t\t\t<td>$$rh_locus_verdict{$s_ID}{$s_loc}</td>\n";
				}

				if($$rh_verified{$s_ID}{$s_loc}{$s_typing} eq "FAIL" && $$rh_locus_verdict{$s_ID}{$s_loc} ne "ERROR"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>FAIL</b></td>\n";
				}elsif($$rh_verified{$s_ID}{$s_loc}{$s_typing} eq "FAIL" && $$rh_locus_verdict{$s_ID}{$s_loc} eq "ERROR"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>ERROR</b></td>\n";
				}else{
					print $html  "\t\t\t<td>PASS</td>\n";
				}

				if(defined $$rh_qc_expected{$s_ID}{$s_loc}){

          if($$rh_qc_expected{$s_ID}{$s_loc} !~ /\|/){
              my ($h1,$h2) = split(/\+/,$$rh_qc_expected{$s_ID}{$s_loc});
              $h1 =~ s/ //g;$h2 =~ s/ //g;
              $h1 = g2p($h1) eq g2p($s_typing) ? "<b>".$h1."</b>" : $h1;
              $h2 = g2p($h2) eq g2p($s_typing) ? "<b>".$h2."</b>" : $h2;
              my $hap = join("+",$h1,$h2);
              print $html  "\t\t\t<td>$hap</td>\n";
          }else{

            my %h_haps   = split(/\|/,$$rh_qc_expected{$s_ID}{$s_loc});
            my $num_haps = keys %h_haps;
            my $cnt = 1;
            print $html "<td>";
            foreach my $hap (keys %h_haps){
              my ($h1,$h2) = split(/\+/,$hap);
              $h1 = g2p($h1) eq g2p($s_typing) ? "<b>".$h1."</b>" : $h1;
              $h2 = g2p($h2) eq g2p($s_typing) ? "<b>".$h2."</b>" : $h2;
              $hap = join("+",$h1,$h2);
              if($cnt % 2 == 0 && $cnt != $num_haps){
                $hap =~ s/ //g;
                print $html $hap."|<br>";
                $cnt++;
              }elsif($cnt != $num_haps){
                $hap =~ s/ //g;
                print $html $hap."|";
                $cnt++;
              }else{
                $hap =~ s/ //g;
                print $html $hap;
              }
            }
            print $html "</td>";
          }

				}else{
					print $html "\t\t\t<td></td>";
				}

				print $html  "\t\t\t<td>$s_typing</td>\n";

		
				if(defined @{$$rh_observed{$s_ID}{$s_loc}}[0]){
          for(0..$#{$$rh_observed{$s_ID}{$s_loc}}){
            my $observed_gl = @{$$rh_observed{$s_ID}{$s_loc}}[$_];
            printObservedAlleles($html,$observed_gl,$s_typing,$n_cutoff2,$n_cutoff1);
          }
				}else{
					print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
					print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
				}

				print $html "</tr>\n";


			}
		}
	}

		my $end_body = qq{</tbody>
	            </table>
	          </div>
	        </div>
	      </div>
	    </div>
		};
		print $html $end_body,"\n";


}
################################################################################################################
=head2 qcFooter

	Title:     qcFooter
	Usage:
	Function:
	Returns:
	Args:

=cut
sub qcFooter{

	my ( $html ) = @_;

	print $html "<script>\n";
	print $html "\$(document).ready(function(){ \n";
	print $html "\$(\"#myTable\").tablesorter(); \n";
	print $html "}\n";
	print $html ");\n";
	print $html "</script>\n";

	my $footer = qq{
			      <script>
	     // If the user clicks in the window, set the background color of <body> to yellow
	      function expand(obj) {
	        if(obj.getElementsByTagName("span")[0].style.display == "table"){
	          obj.getElementsByTagName("span")[0].style.display = "none";
	        }else{
	           obj.getElementsByTagName("span")[0].style.display = "table";
	        }
	      }
	    </script>
	    			<footer>
				<div class="footer">
					<img src="../img/bethematch.jpeg" alt="img02">
          <p><a href="https://github.com/nmdp-bioinformatics/pipeline/blob/master/LICENSE">Copyright ©</a> 2014 - 2015 National Marrow Donor Program. All Rights Reserved.</p>    </div>
				</div>
			</footer>	
	    <style type="text/css">
	        .bargraph {
	          list-style: none;
	          padding-top: 20px;
	          width:560px;
	        }
	        a.id_links {
				text-decoration: none;
				font:bold;
				color:black;
	        }
	        a.id_links:hover {
				color:#A8A8A8;
	        }
	        footer {
					color: #888;
					clear:both;
					position: relative; 
					bottom: 10px;
					left: 0; width: 100%; 
					text-align: center;
					font-size:9px;
					
					background: rgba(255, 255, 255, 0.6);
			}	        
	    </style>
	    <!--<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script> -->
	    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.0.0/js/bootstrap.min.js"></script>
	    <script src="../js/docs.min.js"></script>
	    <!-- IE10 viewport hack for Surface/desktop Windows 8 bug -->
	    <script src="../js/ie10-viewport-bug-workaround.js"></script>
	  
	  </body>
	</html>
	};
	print $html $footer;

}

################################################################################################################
=head2 failedBody

	Title:     failedBody
	Usage:
	Function:  Prints the body to the failed.html file
	Returns:
	Args:

=cut
sub errorsBody{

	my ($html, $s_exp, $rh_locus_verdict, $rh_qc_verdict, $rh_qc_expected, $rh_verified, $rh_observed, $ra_subject_pages, $ra_errors, $rh_consensus ) = @_;

	my $s_experiment = $s_exp;

	foreach my $s_ID (keys %$rh_verified){
		foreach my $s_loc (keys %{$$rh_verified{$s_ID}}){
			next unless defined $$ra_errors{$s_ID}{$s_loc};

			foreach my $s_typing (keys %{$$rh_verified{$s_ID}{$s_loc}}){

				print $html "\t\t<tr>\n";

				my ($n_cutoff1,$n_cutoff2) = $s_loc =~ /^[A|B|C]/ ? (5,16) : (4,13);
				
				next if(!defined $$ra_subject_pages{$s_ID});

				my $s_checkbox = qq{
                  <td >
                    <div style="margin-left:4px;margin-top:10px;"><input type = "checkbox"
                       id = "myCheck"
                       value = "ham" />
                      </div>
                   </td>
                  };
        print $html $s_checkbox;


        my $id_link = "subjects/subject".$$ra_subject_pages{$s_ID}.".html";
				print $html  "\t\t\t<td><a href=\"$id_link\" class=\"id_links\">$s_ID</a></td>\n";

				print $html "\t\t\t<td>$s_loc</td>";
				if($$rh_locus_verdict{$s_ID}{$s_loc} eq "FAIL"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>FAIL</b></td>\n";
				}elsif($$rh_locus_verdict{$s_ID}{$s_loc} eq "ERROR"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>DROPOUT</b></td>\n";
				}else{
					print $html  "\t\t\t<td>$$rh_locus_verdict{$s_ID}{$s_loc}</td>\n";
				}

				if($$rh_verified{$s_ID}{$s_loc}{$s_typing} eq "FAIL" && $$rh_locus_verdict{$s_ID}{$s_loc} eq "ERROR"){
					print $html "\t\t\t<td style=\"color:#F7464A\"><b>DROPOUT</b></td>\n";
				}elsif($$rh_verified{$s_ID}{$s_loc}{$s_typing} eq "FAIL"){
					print $html  "\t\t\t<td>$$rh_verified{$s_ID}{$s_loc}{$s_typing}</td>\n";
				}else{
					print $html  "\t\t\t<td>$$rh_verified{$s_ID}{$s_loc}{$s_typing}</td>\n";
				}

        if(defined $$rh_consensus{$s_ID}{$s_loc}){
          print $html "\t\t\t<td>$$rh_consensus{$s_ID}{$s_loc}</td>\n";
        }else{
          print $html "\t\t\t<td></td>\n";
        }
				
				print $html  "\t\t\t<td>$s_typing</td>\n";
				
				if(defined @{$$rh_observed{$s_ID}{$s_loc}}[0]){
          for(0..$#{$$rh_observed{$s_ID}{$s_loc}}){
            my $observed_gl = @{$$rh_observed{$s_ID}{$s_loc}}[$_];
            printObservedAlleles($html,$observed_gl,$s_typing,$n_cutoff2,$n_cutoff1);
          }
				}else{
					print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
					print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
				}

				print $html "</tr>\n";


			}
		}
	}

		my $end_body = qq{</tbody>
	            </table>
	          </div>
	        </div>
	      </div>
	    </div>
		};
		print $html $end_body,"\n";


}

################################################################################################################
=head2 subjectsBody

	Title: subjectsBody
	Usage:
	Function:
	Returns:
	Args:

=cut
sub subjectsBody{

	my ( $html, $s_exp, $s_ID, $s_blast_dir, $rh_locus_verdict, $rh_qc_verdict, $rh_qc_expected, $rh_expected, $rh_observed, $rh_sequence_data, $rh_seq_ids  ) = @_;

	my $s_experiment = $s_exp;
	foreach my $s_loc (sort keys %{$$rh_locus_verdict{$s_ID}}){
		print $html "\t\t<tr>\n";
		print $html "\t\t\t<th>$s_loc</th>\n";

		if($$rh_locus_verdict{$s_ID}{$s_loc} eq "ERROR"){
			print $html "\t\t\t<td style=\"color:#F7464A\"><b>ERROR</b></td>\n";
		}elsif($$rh_locus_verdict{$s_ID}{$s_loc} eq "FAIL"){
			print $html "\t\t\t<td style=\"color:#F7464A\"><b>FAIL</b></td>\n";
		}elsif($$rh_locus_verdict{$s_ID}{$s_loc} eq "DRBX"){
      print $html "\t\t\t<td style=\"color:#F7464A\"><b>DRBX</b></td>\n";
    }else{
			print $html "\t\t\t<td>PASS</td>\n";
		}

		print $html "\t\t\t<td>$$rh_expected{$s_ID}{$s_loc}</td>\n";
		my @a_expected_haps = split(/\|/,$$rh_expected{$s_ID}{$s_loc});

		my %h_expected_alleles;
		foreach(@a_expected_haps){
			my($x,$y) = split(/\+/,$_);
			$x = g2p($x);
			$y = g2p($y);
			$h_expected_alleles{$x}++;$h_expected_alleles{$y}++;
		}


		if(!defined @{$$rh_observed{$s_ID}{$s_loc}}[0]){
			print $html "\t\t\t<td style=\"color:#F7464A\"><b>ERROR</b></td>\n";
      print $html "\t\t\t<td style=\"color:#F7464A\"><b>ERROR</b></td>\n";
			next;
		}


		my ($n_cutoff1,$n_cutoff2) = $s_loc =~ /^[A|B|C]/ ? (4,13) : (3,10);

		my $printed = 0;my $printed2 = 0;
		if(defined @{$$rh_observed{$s_ID}{$s_loc}}[0]){

			my %h_alleles = map{ $_ => 1 } split(/\//,@{$$rh_observed{$s_ID}{$s_loc}}[0]);
			my $num_alleles = keys %h_alleles;

			if($num_alleles > 9){

				my @a_alleles = keys %h_alleles;
				for(my $i=1;$i<=$#a_alleles+1;$i++){
					print $html "<td><span  onclick=\"expand(this)\">\n" if $i == 1;
					my $s_end = $i == $#a_alleles+1 ? "" : "/";	

					if($i == $n_cutoff2){
						if(defined $h_expected_alleles{$a_alleles[$i-1]}){
							print $html "...<br><span id=\"expanding\"><b>$a_alleles[$i-1]".$s_end."</b>";
						}else{
							print $html "...<br><span id=\"expanding\">$a_alleles[$i-1]".$s_end;
						}
					}elsif($i % $n_cutoff1 == 0) {
						if(defined $h_expected_alleles{$a_alleles[$i-1]}){
							print $html "<b>$a_alleles[$i-1]".$s_end."</b><br>";
						}else{
							print $html "$a_alleles[$i-1]".$s_end."<br>\n";
						}
					}else{
						if(defined $h_expected_alleles{$a_alleles[$i-1]}){
							print $html "<b>$a_alleles[$i-1]".$s_end."</b>";
						}else{
							print $html "$a_alleles[$i-1]".$s_end;
						}
					}

				}

				print $html "</span></span></td>\n";
			}else{
				my @a_alleles = keys %h_alleles;
				print $html "<td>";
				for(my $i=1;$i<=$#a_alleles+1;$i++){
					my $s_end = $i == $#a_alleles+1 ? "" : "/";	

					if($i % $n_cutoff1 == 0) {
						if(defined $h_expected_alleles{$a_alleles[$i-1]}){
							print $html "<b>$a_alleles[$i-1]".$s_end."</b><br>\n";
						}else{
							print $html "$a_alleles[$i-1]".$s_end."<br>\n";
						}
					}else{
						if(defined $h_expected_alleles{$a_alleles[$i-1]}){
							print $html "<b>$a_alleles[$i-1]".$s_end."</b>";
						}else{
							print $html "$a_alleles[$i-1]".$s_end;
						}
					}

				}
				print $html "</td>\n";

			}

			
			my $gl2 = !defined @{$$rh_observed{$s_ID}{$s_loc}}[1] ? @{$$rh_observed{$s_ID}{$s_loc}}[0] : @{$$rh_observed{$s_ID}{$s_loc}}[1];

			my %h_alleles2 = map{ $_ => 1 } split(/\//,$gl2);
			my $num_alleles2 = keys %h_alleles2;

			if($num_alleles2 > 9){
				
				my @a_alleles = keys %h_alleles2;
				for(my $i=1;$i<=$#a_alleles+1;$i++){
					print $html "<td><span  onclick=\"expand(this)\">\n" if $i == 1;
					my $s_end = $i == $#a_alleles+1 ? "" : "/";	

					if($i == $n_cutoff2){
						if(defined $h_expected_alleles{$a_alleles[$i-1]}){
							print $html "...<br><span id=\"expanding\"><b>$a_alleles[$i-1]".$s_end."</b>";
						}else{
							print $html "...<br><span id=\"expanding\">$a_alleles[$i-1]".$s_end;
						}
					}elsif($i % $n_cutoff1 == 0) {
						if(defined $h_expected_alleles{$a_alleles[$i-1]}){
							print $html "<b>$a_alleles[$i-1]".$s_end."</b><br>\n";
						}else{
							print $html "$a_alleles[$i-1]".$s_end."<br>\n";
						}
					}else{
						if(defined $h_expected_alleles{$a_alleles[$i-1]}){
							print $html "<b>$a_alleles[$i-1]".$s_end."</b>";
						}else{
							print $html "$a_alleles[$i-1]".$s_end;
						}
					}

				}

				print $html "</span></span></td>\n";
			}else{

				my @a_alleles = keys %h_alleles2;
				print $html "<td>";
				for(my $i=1;$i<=$#a_alleles+1;$i++){
					my $s_end = $i == $#a_alleles+1 ? "" : "/";	

					if($i % $n_cutoff1 == 0) {
						if(defined $h_expected_alleles{$a_alleles[$i-1]}){
							print $html "<b>$a_alleles[$i-1]".$s_end."</b><br>\n";
						}else{
							print $html "$a_alleles[$i-1]".$s_end."<br>\n";
						}
					}else{
						if(defined $h_expected_alleles{$a_alleles[$i-1]}){
							print $html "<b>$a_alleles[$i-1]".$s_end."</b>";
						}else{
							print $html "$a_alleles[$i-1]".$s_end;
						}
					}

				}
				print $html "</td>\n";

			}


		}else{
			print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
			print $html "<td><b style=\"color:#F7464A\">NA</b></td>\n";
		}

		print $html "</tr>\n";
		print $html  "</tr>\n";


	}


	my $sequence_start = qq{
	          	<div style="margin-top:10px;" class="table-responsive">

		            <table class="table table-striped">
		            <h2>Consensus Sequences</h1>
		              <thead>
		                <tr>
		                  <th>Locus</th>
		                  <th>Zygousity</th>
		                  <th>Sequence 1</th>
		                  <th>Sequence 1</th>
		                </tr>
		              </thead>
		           	<tbody>
	};

	print $html $sequence_start;
	
	foreach my $s_loc (sort keys %{$$rh_sequence_data{$s_ID}}){
		print $html "<tr>\n";
		print $html "<th>$s_loc</th>\n";
		my $num_seqs = keys %{$$rh_sequence_data{$s_ID}{$s_loc}};

		print $html "<td>$num_seqs</td>\n";

		my @a_seqs = keys %{$$rh_sequence_data{$s_ID}{$s_loc}};

		for(my $i=0;$i<=$#a_seqs;$i++){
			my @a_nucls = split(//,$a_seqs[$i]);
			   print $html "<td>";
          my $pwd = `pwd`;chomp($pwd);
          my $full_blast_dir = $pwd."/report/blast/".$s_exp."/".$$rh_seq_ids{$s_ID}{$s_loc}{$a_seqs[$i]}.".html";
          if(defined $$rh_seq_ids{$s_ID}{$s_loc}{$a_seqs[$i]} && defined $s_blast_dir ){
            my $s_blast_file = "../../blast/".$s_exp."/".$$rh_seq_ids{$s_ID}{$s_loc}{$a_seqs[$i]}.".html";
            if(-e $full_blast_dir){
              print $html "<a href=\"".$s_blast_file."\" class=\"hovering\">";
            }
          }

			for(my $j=1;$j<=$#a_nucls+1;$j++){

				if($j % 50 == 0) {
					print $html "$a_nucls[$j-1]<br>\n";
				}else{
					print $html "$a_nucls[$j-1]";
				}
			}
      print $html "</a>" if(-e $full_blast_dir);
			print $html "</td>\n"
		}

		print $html "</tr>\n";

	}


	my $end_body = qq{
          		</tbody>
          	</table>
          </div>
          </div>
        </div>
      </div>
    </div>
	};
	print $html $end_body,"\n";

}

################################################################################################################
=head2 g2p

	Title:     g2p
	Usage:     my $s_allele = g2p($s_genomic_allele)
	Function:  This turns genomic level alleles into protein level alleles 
			   This will r
	Returns:   protein level alleles
	Args:      genomic level allele

=cut
sub g2p{
	my $allele    = shift;
	return $allele if $allele !~ /:/;
	my @a_alleles = split(/:/,$allele);
	return join(":",$a_alleles[0],$a_alleles[1]);
}


1;

