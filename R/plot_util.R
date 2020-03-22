#################################################################################
##
## Author:  Nat Goodman
## Created: 20-03-22
##          from plot_nudge.R created 20-03-19
##          from misisg/R/plot.R created 19-01-09
##          uses code from repwr/R/plot.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Plot utility functions
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## ---- Plot Utility Functions ----
## display color palette - used to be in util.R
pal=function(col,border="light gray",...) {
 n=length(col)
 plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab="",ylab="",...)
 rect(0:(n-1)/n,0,1:n/n,1,col=col,border=border)
}
## auto-scale title
cex_title=function(title) {
  xyplt=par('plt');                     # dimensions of plot region
  xplt=xyplt[2]-xyplt[1];               # width of plot region
  min(1,xplt/strwidth(title,units='fig'));
}
## make colors from RColorBrewer palettes
## palettes - RColorBrewer palette names
## n - number of colors
## names - names to index colors. if set, overrides n
## skip - for seqeuential palettes, number of colors to skip  - 1st two usually too light
## dark.first - for seqeuential palettes, reverse so darker colors first
col_brew=function(palettes,n=0,names=NULL,skip=2,dark.first=TRUE) {
  if (is.null(palettes)||(missing(n)&&missing(names))) return(NULL);
  bad=palettes %notin% rownames(brewer.pal.info);
  if (any(bad)) stop(paste('Invalid palette name(s):',paste(collapse=', ',palettes[bad])));
  if (!is.null(names)) n=length(names);
  n.pal=length(palettes);
  ns=if(n.pal==1) n else as.integer(table(cut(1:n,n.pal)));
  ## col=do.call(c,lapply(1:length(palettes),
  ##                      function(i) col_brew_(palettes[i],ns[i],skip,dark.first)));
  col=do.call(c,col_brew_(palettes,ns,skip,dark.first));
  setNames(col,names);
}
col_brew_=Vectorize(function(pal,n,skip,dark.first) {
  if (n==0) return(NULL);
  if (brewer.pal.info[pal,'category']!='seq') m=n else m=n+skip;
  m=min(m,brewer.pal.info[pal,'maxcolors']);
  m=max(3,m);                         # all palettes have min 3 colors
  col=brewer.pal(m,pal);
  if (brewer.pal.info[pal,'category']=='seq') {
    if (skip>0) col=tail(col,-skip);
    if (dark.first) col=rev(col);
  }
  ## if want more colors than in palette, colorRampPalette will make more
  if (n>length(col)) col=colorRampPalette(col)(n);
  col;
},vectorize.args=cq(pal,n),SIMPLIFY=FALSE,USE.NAMES=FALSE);
