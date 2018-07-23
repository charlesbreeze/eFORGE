function createAjaxObj(){
var httprequest=false
if (window.XMLHttpRequest){ // if Mozilla, Safari etc
httprequest=new XMLHttpRequest()
if (httprequest.overrideMimeType)
httprequest.overrideMimeType('text/xml')
}
else if (window.ActiveXObject){ // if IE
try {
httprequest=new ActiveXObject("Msxml2.XMLHTTP");
}
catch (e){
try{
httprequest=new ActiveXObject("Microsoft.XMLHTTP");
}
catch (e){}
}
}
return httprequest
}

// -------------------------------------------------------------------
// Main RSS Ticker Object function
// rss_ticker(RSS_id, cachetime, divId, divClass, delay, optionalswitch)
// -------------------------------------------------------------------

function rss_ticker(RSS_id, ticker, optionalswitch){
CurrentLength=0;
this.RSS_id=RSS_id //Array key indicating which RSS feed to display
this.ticker=ticker
this.cachetime=120 //Time to cache feed, in minutes. 0=no cache.
this.tickerid='newsbox' //ID of ticker div to display information
this.delay=3000 //Delay between msg change, in miliseconds.
this.logicswitch=(typeof optionalswitch!="undefined")? optionalswitch : -1
this.mouseoverBol=0 //Boolean to indicate whether mouse is currently over ticker (and pause it if it is)
this.pointer=0
divId='newsbox';
divClass='newsclass';
this.ajaxobj=createAjaxObj()
document.write('<div id=\"newsbox\" class=\"newsclass\"></div>')
//document.write('<div id="'+divId+'" class="'+divClass+'"></div>')
this.getAjaxcontent()
}

// -------------------------------------------------------------------
// getAjaxcontent()- Makes asynchronous GET request to "rssfetch.php" with the supplied parameters
// -------------------------------------------------------------------

rss_ticker.prototype.getAjaxcontent=function(){
if (this.ajaxobj){
var instanceOfTicker=this
var parameters="id="+encodeURIComponent(this.RSS_id)
this.ajaxobj.onreadystatechange=function(){instanceOfTicker.initialize()}
this.ajaxobj.open('GET', this.RSS_id, true)
this.ajaxobj.send(null)
}
}

// -------------------------------------------------------------------
// initialize()- Initialize ticker method.
// -Gets contents of RSS content and parse it using JavaScript DOM methods
// -------------------------------------------------------------------

rss_ticker.prototype.initialize=function(){
if (this.ajaxobj.readyState == 4){ //if request of file completed
if (this.ajaxobj.status==200){ //if request was successful
var xmldata=this.ajaxobj.responseXML
if (xmldata.getElementsByTagName("item").length==0){ //if no <item> elements found in returned content
document.getElementById(this.tickerid).innerHTML="<b>Error</b> fetching remote RSS feed!<br />"+this.ajaxobj.responseText
return
}
var instanceOfTicker=this
this.feeditems=xmldata.getElementsByTagName("item")

//Cycle through RSS XML object and store each piece of the item element as an attribute of the element
for (var i=0; i<this.feeditems.length; i++){
this.feeditems[i].setAttribute("ctitle", this.feeditems[i].getElementsByTagName("title")[0].firstChild.nodeValue)
this.feeditems[i].setAttribute("clink", this.feeditems[i].getElementsByTagName("link")[0].firstChild.nodeValue)
this.feeditems[i].setAttribute("cdescription", this.feeditems[i].getElementsByTagName("description")[0].firstChild.nodeValue)
}
document.getElementById(this.tickerid).onmouseover=function(){instanceOfTicker.mouseoverBol=1}
document.getElementById(this.tickerid).onmouseout=function(){instanceOfTicker.mouseoverBol=0}


this.rotatemsg()
}
}
}

// -------------------------------------------------------------------
// rotatemsg()- Rotate through RSS messages and displays them
// -------------------------------------------------------------------

rss_ticker.prototype.rotatemsg=function(){
var instanceOfTicker=this
if (this.mouseoverBol==1) //if mouse is currently over ticker, do nothing (pause it)
setTimeout(function(){instanceOfTicker.rotatemsg()}, 100)
else {



var tickerDiv=document.getElementById(this.tickerid);
//var tickercontent='<a href="'+this.feeditems[this.pointer].getAttribute("clink")+'">'+this.feeditems[this.pointer].getAttribute("ctitle")+'</a>'
var link='<p>Latest: <a href="'+this.feeditems[this.pointer].getAttribute("clink")+'">';
var title=this.feeditems[this.pointer].getAttribute("ctitle");
var endlink='</a></p>';

if (this.ticker=='True') {
tickercontent=link+title.substring(0, CurrentLength)+endlink;
if (this.logicswitch=="True"){
tickercontent+="<br />"+this.feeditems[this.pointer].getAttribute("cdescription");}
tickerDiv.innerHTML=tickercontent;

if (CurrentLength!=title.length) {
CurrentLength++;
this.delay=50;
}
else {
this.delay=3000;
this.pointer=(this.pointer<this.feeditems.length-1)? this.pointer+1 : 0
CurrentLength=0;
}
} else {

tickercontent=link+title+endlink;
if (this.logicswitch=="True") {
tickercontent+="<br />"+this.feeditems[this.pointer].getAttribute("cdescription");}
tickerDiv.innerHTML=tickercontent;
this.pointer=(this.pointer<this.feeditems.length-1)? this.pointer+1 : 0
}

//if (this.logicswitch=="showdescription")
//tickercontent+="<br />"+this.feeditems[this.pointer].getAttribute("cdescription")
//tickerDiv.innerHTML=tickercontent
setTimeout(function(){instanceOfTicker.rotatemsg()}, this.delay) //update container every second
}
}


