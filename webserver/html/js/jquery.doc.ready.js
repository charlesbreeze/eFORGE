var UCL=function(){};

$(document).ready(function(){
    var uclSetup = new UCL();
    fontResizer('14px','16px','18px');


    if ($("#access"))
    {
        var accessCookie = $.cookie('accessCookie');
        var accessTxt = "<p><span id=\"normalView\" style=\"cursor:pointer; font-weight:bold\">Normal view</span> | <span id=\"hcView\" style=\"cursor:pointer; font-weight:bold\">High contrast view</span></p>";
        $("#access").html(accessTxt);
        if (accessCookie == "true")
        {
            showContrast();
        }
        else
        {
            $("#hcView").click(function(){
                showContrast();
            });
        }
        $("#normalView").click(function(){
            hideContrast();
        });
    }
    
    function showContrast(){
       $("#contrastCSS").attr("href","//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/css/screen/ucl_contrast_view.css");
        $.cookie('accessCookie','true', { path: '/', expires: 30 });
    }
    function hideContrast(){
       $("#contrastCSS").attr("href","//static.ucl.ac.uk/silva/UCLDefaultLayoutV3/css/screen/ucl_normal_view.css");
        $.cookie('accessCookie','false', { path: '/', expires: 30 });
    }

//Start Google Analytics download tracking code
    var filetypes = /\.(zip|exe|pdf|doc*|xls*|ppt*|mp3)$/i;
    var baseHref = '';

    if ($('base').attr('href') != undefined)
        baseHref = $('base').attr('href');

    $('a').each(function() {
        var href = $(this).attr('href');
        if (href && (href.match(/^https?\:/i)) && (!href.match(document.domain))) {
            $(this).click(function() {
                var extLink = href.replace(/^https?\:\/\//i, '');
                _gaq.push(['_trackEvent', 'External', 'Click', extLink]);
                if ($(this).attr('target') != undefined && $(this).attr('target').toLowerCase() != '_blank') {
                    setTimeout(function() { location.href = href; }, 200);
                    return false;
                }
            });
        }
        else if (href && href.match(/^mailto\:/i)) {
            $(this).click(function() {
                var mailLink = href.replace(/^mailto\:/i, '');
                _gaq.push(['_trackEvent', 'Email', 'Click', mailLink]);
            });
        }
        else if (href && href.match(filetypes)) {
            $(this).click(function() {
                var extension = (/[.]/.exec(href)) ? /[^.]+$/.exec(href) : undefined;
                var filePath = href;
                _gaq.push(['_trackEvent', 'Download', 'Click-' + extension, filePath]);
                if ($(this).attr('target') != undefined && $(this).attr('target').toLowerCase() != '_blank') {
                    setTimeout(function() { location.href = href; }, 200);
                    return false;
                }
            });
        }
    });
//End Google Analytics download tracking code

});


