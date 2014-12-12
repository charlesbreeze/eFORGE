#!/usr/bin/env Rscript

args <- commandArgs(trailing=T);

outdir <- ".";
my.colour <- "#29A6C9";
input.size.unit <- "bytes";
sqlite.file <- NA;

i <- 1;
while (i < length(args)) {
    if (args[i] == "--outdir") {
        outdir <- args[i+1];
        i <- i +2;
    } else if (args[i] == "--colour") {
        my.colour <- args[i+1];
        i <- i +2;
    } else if (args[i] == "--input_size") {
        input.size.unit <- args[i+1];
        i <- i +2;
    } else {
        sqlite.file <- args[i];
        i <- i + 1;
    }
}

if (is.na(sqlite.file)) {
    stop(sprintf("Usage: log.R [--outdir '.'] [--colour #29A6C9] [--input_size bytes] sqlfile"))
}

if (file.access(sqlite.file) == -1) {
    stop(sprintf("Specified SQLite file ( %s ) does not exist", sqlite.file))
} else {
    suppressPackageStartupMessages(library(RSQLite))
    drv <- dbDriver("SQLite");
    con <- dbConnect(drv, sqlite.file)
}

if (!exists("con")) {
  stop("No connection")
}

statement <- "SELECT timestamp, runtime_sec, input_size, cmd_line, error from server_log";

res <- dbSendQuery(con, statement = statement)
data <- dbFetch(res)
dbClearResult(res)

data$days = as.numeric(difftime(data$timestamp, Sys.Date(), units="days")+6);

png(paste0(outdir, "/last_week.png"), width=680, height=300);
hist(subset(data, days>0)$days, breaks=0:28/4, xlab="days", ylab="Num of jobs", col=my.colour,
    main=paste0("Last week (n=", dim(subset(data, days>0))[1], ")"), axes=F)
axis(1, at=0:7, labels=F)
axis(1, at=(0:6)+0.5, tick=F, labels=weekdays(Sys.Date()-6:0, abbreviate = T))
axis(2)
dev.off()

data$days = as.numeric(difftime(data$timestamp, Sys.Date(), units="days")+30);

png(paste0(outdir, "/last_month.png"), width=680, height=300);
hist(subset(data, days>0)$days, breaks=0:31, xlab="date", ylab="Num of jobs", col=my.colour,
    main=paste0("Last month (n=", dim(subset(data, days>0))[1], ")"), axes=F)
axis(1, at=0:31, labels=F)
axis(1, at=0:4*7+1.5, tick=F, labels=format(Sys.Date()-(4:0)*7-1, format="%d/%m"))
axis(2)
dev.off()

png(paste0(outdir, "/scalability.png"), width=680, height=500);
smoothScatter(data$input_size, data$runtime_sec, colramp = colorRampPalette(c("white", my.colour)),
    xlab=paste0("Input size (",input.size.unit,")"), ylab="Runtime (s)", main="Scalability")
dev.off()
