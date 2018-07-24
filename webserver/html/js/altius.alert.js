function positionOverlayEmbed() {
    x = $(window).width()/2;
    y = $(window).height()/2 + $(window).scrollTop();
    x -= $("#alert-message").width()/2;
    y -= $("#alert-message").height()/2;
    if (x < 0)
        x = 0;
    if (y < 0)
        y = 0;
    var cx = $("#alert-message").width()-19;
//    var offset = $("#alert-message").offset();
//    if (y > offset["top"] + 5 || y < offset["top"] -5 || x > offset["left"] + 5 || x < offset["left"] - 5)
//        $("#alert-message").offset({top:y,left:x});
    $("#alert-message").css("top",y + "px");
    $("#alert-message").css("left",x + "px");
    $("#alert-message").css("z-index","999");
    $("#alert-message").css("display","none");
}

$(document).ready(function() {
    if ($("#alert-message-overlay").is("div") && $.cookie("alertmessagecheck") != "false") {
//        if (!$("#alert-message-overlay").is("div"))
//            $("body").prepend("<div id=\"lightbox-screen\">&nbsp;</div>");
        // set the width and height to the window width and height
        $("#alert-message-overlay").css("width",$(window).width() + "px");
        $("#alert-message-overlay").css("height",$(window).height() + "px");
        $("#alert-message-overlay").fadeIn("fast");
        // now position the embedded window
        $("#alert-message").css("display","block");
        positionOverlayEmbed();
        $("#alert-message").css("opacity","1");
        $("#alert-message").css("display","none");
        $("#alert-message").fadeIn("slow");
        // now set up the closer click handlers
        // assume $uery.cookie is installed, as it should be
        $("#alert-message").click(function() {
            var expiry = new Date();
            expiry.setTime(expiry.getTime() + 1200000); 
//            var expiry = new Date() + 12000000;
            $.cookie("alertmessagecheck","false",{path:'/',expires:expiry});
            $("#alert-message-overlay").fadeOut("slow");
            $("#alert-message").fadeOut("slow");
        });
    }
/*    if (location.href.indexOf("www.ucl.ac.uk") == -1) {
        var servers = [
            {
                url:"localhost",
                msg:"This is your local machine"
            },
            {
                url:"127.0.0.1",
                msg:"This is your local machine"
            },
            {
                url:"ltms156.ltms-isd.ucl.ac.uk",
                msg:"This is ltms156"
            },
            {
                url:"www.silva-sandbox.ucl.ac.uk",
                msg:"This is the Silva Sandbox"
            },
            {
                url:"wwwcm-d.ucl.ac.uk/",
                msg:"This is wwwcm-d"
            }
        ];
        for (var i = 0; i < servers.length; i++) {
            if (location.href.indexOf(servers[i]["url"]) > -1) {
                if (!$("#sandbox-info").is("div"))
                    $("#header").prepend("<div id=\"sandbox-info\" style=\"position:absolute; top:81px; left:300px; height:16px; padding:6px 20px;background:#FFF; z-index:1\">" + servers[i]["msg"] + "</div>");
                else
                    $("#sandbox-info").html(servers[i]["msg"]);
                var w = ($("#nav > .hlist").width() -  $("#sandbox-info").width())/2;
                $("#sandbox-info").css("left",w + "px");
                break;
            }
        }*/
/*        var txt = $.get("/silva/sandbox-info", function(data) {
            if (!$("#sandbox-info").is("div"))
                $("#header").prepend("<div id=\"sandbox-info\" style=\"position:absolute; top:81px; left:300px; height:16px; padding:6px 20px;background:#FFF; z-index:1\">" + data + "</div>");
            else
                $("#sandbox-info").html(data);
            var w = ($("#nav > .hlist").width() -  $("#sandbox-info").width())/2;
            $("#sandbox-info").css("left",w + "px");
        });
*/
//    }
});
