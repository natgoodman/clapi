#################################################################################
##
## Author:  Nat Goodman
## Created: 20-03-19
##          from misig/R/plot.R created 19-01-09
##          uses code from repwr/R/plot.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Plotting code for doc_nudge
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## ---- Plot Functions ----
## plot multiple lines - my adaptation of matplot - adapted from repwr/plotratm
## x is vector or matrix of x values
## y is vector or matrix of y values
##   like matplot, each line is column of x or y
##   unlike matplot, code currently assumes at most one of x,y is matrix - for smoothing to work
## col, lty, lwd are the usual line properties
## title, cex.title are title and cex for title
## xaxt, yaxt control labels on axes - NOT YET IMPLEMENTED
##   's' or NULL means let R do it
##    else list of axis params, eg, at, labels
##           title='',cex.title='auto'
## type is plot type - passed to matplot. 'n' also turns off extra lines, legend, grid
##   TODO: type='p' should cause legend to draw points instead of lines
## vline,hline are vectors of x or y positions for extra vertical or horizontal lines
## vhlty, vhcol, vhlwd are lty, col, lwd for these extra lines
## vlab, hlab contol writing vline, hline values along axes
## vhdigits is number of digits for these values
## smooth whether to smooth data to make plot prettier
##   aspline, spline, loess, linear, approx (sames as linear) none
##   default aspline
##   TRUE means aspline, FALSE means none
## smooth.xy tells which axis is domain of smoothing
##   default: 'x' if both x and y are vector-like, else whichever is vector-like
## smooth.args are additional args passed to smooth method.
##   note defaults for aspline, spline, loess
## legend tells whether and where to draw legend
##   TRUE or word like 'right' means draw legend, FALSE or NULL means no legend
## legend.title is legend title
## legend.args are further legend params
plotm=
  function(x,y,col='black',lty='solid',lwd=1,title='',cex.title='auto',type='l',
           xaxt='s',yaxt='s',
           vline=NULL,hline=NULL,vhlty='dashed',vhcol='grey50',
           vhlwd=1,vlab=TRUE,hlab=TRUE,vhdigits=2,
           ## NG 20-01-02: redo 'smooth' logic again to set smooth.args default
           smooth=cq(aspline,spline,loess,linear,approx,none),
           smoothx=if(is.logical(smooth)) {
                     if (smooth) 'aspline' else 'none';
                   } else match.arg(smooth),
           smooth.args=
             switch(smoothx,
                    aspline=list(method='improved'),spline=list(spar=0.5),loess=list(span=0.75),
                    list()),
           smooth.xy=cq(x,y),
           legend=if(is.vector(y)) FALSE else 'right',legend.title=NULL,
           legend.labels=if(is.vector(y)) NULL else colnames(y),
           legend.args=list(where=NULL,x=NULL,y=NULL,cex=0.8,
                            title=legend.title,labels=legend.labels,col=col,lty=lty,lwd=lwd),
           ...) {
    smooth=smoothx;                  # to avoid confusion later
    if (is.null(x)) stop("Nothing to plot: 'x' is NULL");
    if (is.null(y)) stop("Nothing to plot: 'y' is NULL");
    if (is.vector(x)) x=as.matrix(x)
    else if (length(dim(x))!=2) stop("'x' must be vector or 2-dimensional matrix-like object");
    if (is.vector(y)) y=as.matrix(y)
    else if (length(dim(y))!=2) stop("'y' must be vector or 2-dimensional matrix-like object");
    if (nrow(x)!=nrow(y)) stop("'x' and 'y' must have same number of rows");
    ## NG 19-12-31: allow both x,y to be matrix-like as in matplot
    ##   like matplot, if both have multiple columns, must be same number of columns
    if (ncol(x)>1&&ncol(y)>1&&ncol(x)!=ncol(y))
      stop("When 'x' and 'y' both have multiple columns, must have same number of columns");
    if (is.null(cex.title)|cex.title=='auto') cex.title=cex_title(title);
    ## NG 20-01-02: redo 'smooth' logic. start with TRUE. match choices if not logical
    if (smooth!='none') {
      if (missing(smooth.xy))
        smooth.xy=if(ncol(x)==1) 'x' else if (ncol(y)==1) 'y' else 'x';
      ## smooth
      if (smooth.xy=='x') {
        xy=smooth(x,y,method=smooth,method.args=smooth.args);
        x=xy$x; y=xy$y;
        ## x.smooth=seq(min(x),max(x),len=100);
        ## y=smooth(x,y,xout=x.smooth,method=smooth,spar=spar,span=span);
        ## x=x.smooth;
      } else {
        xy=smooth(y,x,method=smooth,method.args=smooth.args);
        x=xy$y; y=xy$x;
        ## y.smooth=seq(min(y),max(y),len=100);
        ## x=smooth(x=y,y=x,xout=y.smooth,method=smooth,spar=spar,span=span);
        ## y=y.smooth;
      }
    }
    matplot(x,y,main=title,cex.main=cex.title,col=col,lty=lty,lwd=lwd,type=type,
            xaxt=xaxt,yaxt=yaxt,...);
    if (type!='n') {
      grid();
      ## plot extra lines & values if desired. nop if vline, hline NULL
      vhline(vline=vline,hline=hline,vlab=vlab,hlab=hlab,vhdigits=vhdigits,
             lty=vhlty,col=vhcol,lwd=vhlwd);
      ## draw legend if desired
      if (is.null(legend)) legend=FALSE
      else if (!is.logical(legend)) {
        legend.args$where=legend;
        legend=TRUE;
      }
      if (legend) do.call(plotm_legend,legend.args);
    }
  }
## wrapper for plotm that pulls values from data frame
## x,y,z are column names
plotm_df=
  function(data,x,y,z=NULL,
           xlab=if(length(x)==1) x else 'x',
           ylab=if(length(y)==1) y else 'y',
           zlab=if(length(z)==1) z else 'z',
           ...) {
  force(xlab); force(ylab); force(zlab);
  if (!is.data.frame(data)) stop ('data must be data frame');
  ## NG 20-01-05: allow both x,y to have multiple columns as in matplot
  ##   like matplot, if both have multiple columns, must be same number of columns
  if (length(x)>1&&length(y)>1&&length(x)!=length(y))
    stop("When 'x' and 'y' both have multiple columns, must have same number of columns");
  if (!is_subset(x,colnames(data))) stop("'x' must contain column names");
  if (!is_subset(y,colnames(data))) stop("'y' must contain column names");
  if (!is.null(z)&&((length(z)!=1)||(z %notin% colnames(data))))
    stop('z must be NULL or contain exactly one column name');
  if (is.null(z)) {
    ## if no 'z', pass 'x', 'y' columns directly to plotm
    xdata=data[,x,drop=F];
    ydata=data[,y,drop=F];
  } else {
    ## split by 'z'
    by=split(data,data[,z]);
    xdata=do.call(cbind,lapply(by,function(data) data[,x,drop=F]));
    ydata=do.call(cbind,lapply(by,function(data) data[,y,drop=F]));
    ## simplify if all columns identical
    if (length(x)==1&&ncol(xdata)>1&&all(sapply(xdata,identical,xdata[,1])))
      xdata=xdata[,x,drop=F];
    if (length(y)==1&&ncol(ydata)>1&&all(sapply(ydata,identical,ydata[,1])))
      ydata=ydata[,y,drop=F];
  }
  plotm(xdata,ydata,xlab=xlab,ylab=ylab,legend.title=zlab,...);
  }
## plot multiple hyper, binom for multiple cases
## adapted from misig/plotmg, in turn adapted from repwr/plot_ragm and others
## what - 'hyper', 'binom' or both
## distribution params - can be vectors
##   g1,g2 - number of elements in the two groups
##   n - number of elements selected
## col, lty, lwd - the usual line properties
## palette - RColorBrewer sequential palette name
## title, cex.title - title and cex for title
## type - plot type - passed to matplot. 'n' also turns off extra lines, legend, grid
##   TODO: type='p' should cause legend to draw points instead of lines
## extend.x - extend xlim to g1 - CAUTION: often scrunches informative part of plot
## vline,hline are vectors of x or y positions for extra vertical or horizontal lines
## vhlty, vhcol, vhlwd are lty, col, lwd for these extra lines
## vlab, hlab contol writing vline, hline values along axes
## vhdigits is number of digits for these values

## legend tells whether and where to draw legend
##   TRUE or word like 'right' means draw legend, FALSE or NULL means no legend
## legend.title is legend title
## legend.args are further legend params

## legends is list of legend.args passed to multi_legend

plotmg=
  function(what=cq(hyper,binom),g1,g2,n,extend.x=FALSE,
           col=NULL,lty=setNames(cq(solid,dotted),cq(hyper,binom)),lwd=1,type='l',
           palette=cq(),
           title='',cex.title='auto',
           xaxt='s',yaxt='s',
           vline=NULL,hline=NULL,vhlty='dashed',vhcol='grey50',
           vhlwd=1,vlab=TRUE,hlab=TRUE,vhdigits=2,
           legend='right',legend.title=NULL,legend.labels=NULL,
           legend.args=list(where=NULL,x=NULL,y=NULL,cex=0.8,
                            title=legend.title,labels=legend.labels,col=col,lty=lty,lwd=lwd)

           
           ...) {
    smooth=smoothx;                  # to avoid confusion later
    default.args=list(x=x,y=y,col=col,lty=lty,lwd=lwd,
                      smooth=smooth,smooth.args=smooth.args,legend=FALSE,...);
    if (!is.null(plotms)) 
      if (is_list(x)||is_list(y))
        stop("When plotms is set, 'x' and 'y' cannot have multiple groups");
    if (is.null(plotms)) {
      if (is_list(x)&&is_list(y))
        stop("At most one of 'x' or 'y' can have multiple groups");
      ## create plotms from x, y
      if (is_list(x)) plotms=lapply(x,function(x) list(x=x,y=y))
      else if (is_list(y)) plotms=lapply(y,function(y) list(x=x,y=y))
      else plotms=list(list(x=x,y=y));
    }
    ## fill defaults
    plotms=lapply(plotms,function(plotm.args) fill_defaults(default.args,plotm.args));
    ## if plotms is singleton, just call plotm
    if (length(plotms)==1) {
      do.call(plotm,plotms[[1]]);
    } else {
      ## general case: plotms has multiple groups
      ## setup plot: calculate range, then call plotm with type='n'
      x.range=range(sapply(plotms,function(plotm.args) range(plotm.args$x)));
      y.range=range(sapply(plotms,function(plotm.args) range(plotm.args$y)));
      plotm(x=x.range,y=y.range,type='n',smooth='none',...);
      ## do plot! w/ add=T, w/o legend
      lapply(plotms,function(plotm.args) {
        plotm.args$add=TRUE;
        do.call(plotm,plotm.args);
      });
    }
    ## draw legend if desired
    if (is.null(legend)) legend=FALSE;
    if (is.logical(legend)&&legend) legend='right';
    if (!is.logical(legend)) multi_legend(legends=legends,legend=legend);
    return();
  }
## wrapper for plotmg that pulls values from data frames
## interface optimized for common case:
##   all groups use same colors; each group may use different lty, lwd
## data is data frame
## x,y,z are column names
##   x,y must be singleton or have same length as list of data frames
##   z must be singleton
## col, lty, lwd are vectors
##   col applies to all groups
##   elements of lty, lwd apply to each group in turn
## legend TBD
plotmg_df=
  function(data,x,y,z=NULL,
           xlab=if(length(x)==1) x else 'x',
           ylab=if(length(y)==1) y else 'y',
           zlab=if(length(z)==1) z else 'z',
           col='black',lty='solid',lwd=1,
           ...) {
    force(xlab); force(ylab); force(zlab);
    if (!is.data.frame(data)) stop ('data must be data frame');
    ##   like matplot, if both have multiple columns, must be same number of columns
    if (length(x)>1&&length(y)>1&&length(x)!=length(y))
      stop("When 'x' and 'y' both have multiple columns, must have same number of columns");

    if (!is_subset(x,colnames(data))) stop("'x' must contain column names");
    if (!is_subset(y,colnames(data))) stop("'y' must contain column names");
    if (!is.null(z)&&((length(z)!=1)||(z %notin% colnames(data))))
      stop('z must be NULL or contain exactly one column name');

    ngroup=max(length(x),length(y));
    x=rep(x,length=ngroup);
    y=rep(y,length=ngroup);
    lty=rep(lty,length=ngroup);
    lwd=rep(lwd,length=ngroup);

    plotms=lapply(1:ngroup,function(i) {
      x=x[i]; y=y[i]; lty=lty[i]; lwd=lwd[i]; 
      if (is.null(z)) {
        ## if no 'z', construct plotms from 'x', 'y' columns
        xdata=data[,x,drop=F];
        ydata=data[,y,drop=F];
        col=col[1];
      } else {
        ## split by 'z'
        by=split(data,data[,z]);
        xdata=do.call(cbind,lapply(by,function(data) data[,x,drop=F]));
        ydata=do.call(cbind,lapply(by,function(data) data[,y,drop=F]));
        col=rep(col,length=length(by));
        ## simplify if all columns identical
        if (length(x)==1&&ncol(xdata)>1&&all(sapply(xdata,identical,xdata[,1])))
          xdata=xdata[,x,drop=F];
        if (length(y)==1&&ncol(ydata)>1&&all(sapply(ydata,identical,ydata[,1])))
          ydata=ydata[,y,drop=F];
      }
      list(x=xdata,y=ydata,col=col,lwd=lwd,lty=lty);
      });
    plotmg(plotms,xlab=xlab,ylab=ylab,col=col,lty=lty,lwd=lwd,...);
    return()
  }
## empty plot - just title & axes
plotempty=
  function(title='',cex.title='auto',xlab='x',ylab='y',xlim=c(0,1),ylim=c(0,1),
           xaxp=c(xlim,1),yaxp=c(ylim,1),...) {
    if (is.null(cex.title)|cex.title=='auto') cex.title=cex_title(title);
    plot(x=NULL,y=NULL,type='n',main=title,cex.main=cex.title,xlab=xlab,ylab=ylab,
         xlim=xlim,ylim=ylim,xaxp=xaxp,yaxp=yaxp,...);
    xylim=par('usr');                   # limits of disply region
    xmid=mean(xylim[1:2]);
    ymid=mean(xylim[3:4]);
    text(xmid,ymid,'PLOT DELIBERATELY LEFT BLANK',adj=c(0.5,0.5));
    invisible();
}

## helper functions to plot horizontal and vertical line segments
vhline=function(vline=NULL,hline=NULL,vlab=TRUE,hlab=TRUE,vhdigits=2,col=NA,cex=0.75,...) {
  xylim=par('usr');
  vline=vline[which(between(vline,xylim[1],xylim[2]))];
  hline=hline[which(between(hline,xylim[3],xylim[4]))];
  abline(v=vline,h=hline,col=col,...);
  ## write vhline values along axes
  vline=vline[vlab];
  if (length(vline)>0)
    mtext(round(vline,vhdigits),side=1,at=vline,col=col,line=0.25,cex=cex*par('cex'));
  hline=hline[hlab];
  if (length(hline)>0)
    mtext(round(hline,vhdigits),side=2,at=hline,col=col,line=0.25,cex=cex*par('cex'));
}
hline=
  function(y,x0=0,x,col='black',lty='solid',lwd=1,cex=0.75,text=NULL,
           label=list(text=text,side=2,at=y,col=col,line=0.25,cex=cex*par('cex'),las=1)) {
    segments(x0=x0,x1=x,y0=y,y1=y,col=col,lty=lty,lwd=lwd);
    if (!is.null(text)) do.call(mtext,label);
  }
vline=
  function(x,y0=0,y,col='black',lty='solid',lwd=1,cex=0.75,text=NULL,
           label=list(text=text,side=1,at=x,col=col,line=0.25,cex=cex*par('cex'),las=1)) {
    segments(x0=x,x1=x,y0=y0,y1=y,col=col,lty=lty,lwd=lwd);
    if (!is.null(text)) do.call(mtext,label);
  }

## draw pval legend. works for big picture figure and probability plots
pval_legend=function(x.scale,y.scale,x0=NULL,cex,label='p-value') {
  param(brk.pval,col.pval,steps.pvcol,sig.level);
  ## plt=par('usr');                       # plot region in user coordinates
  ## names(plt)=cq(left,right,bottom,top);
  xtkl=par('xaxp');                    # x tick locations
  ytkl=par('yaxp');                    # y tick locations
  names(xtkl)=names(ytkl)=cq(lo,hi,num);
  if (is.null(x0)) x0=xtkl['lo'];
  width=x.scale*(xtkl['hi']-x0);
  x1=x0+width;
  y1=ytkl['hi'];
  height=y.scale*(y1-ytkl['lo']);
  y0=y1-height;
  ## image sometimes leaves blank space when y0 is between tick marks. roundoff problem, I think
  ## works better to stretch y0 to next lower tick
  tkl=seq(ytkl['lo'],y1,len=ytkl['num']+1);
  y0=tkl[findInterval(y0,tkl)];
  height=y1-y0;                         # adjust height for new y0
  x=c(x0,x1);
  y=seq(y0,y1,length.out=2*steps.pvcol+1)[1:(2*steps.pvcol)]
  z=t(as.matrix(rev(head(brk.pval,-1))));
  image(x,y,z,add=TRUE,breaks=brk.pval,col=col.pval);
  ## add legend text
  x1=x1+strwidth(' ',cex=cex);
  text(x1,y0,0,adj=c(0,0),cex=cex)
  text(x1,y0+height/2,sig.level,adj=c(0,0.5),cex=cex)
  text(x1,y0+height,1,adj=c(0,1),cex=cex)
  y1=y1+strheight('p-value',cex=cex);
  text(x0+width/2,y1,label,adj=c(0.5,0.5),cex=cex);
}
## draw plotm legend. adapted from repwr/mesr_legend
## labels and legend are synonyms
plotm_legend=
  function(where=NULL,x=NULL,y=NULL,cex=0.8,bty='n',
           title=NULL,title.col='black',
           col='black',lty='solid',lwd=1,labels=NULL,legend=labels,...) {
    if (is.null(legend)) return();      # nothing to draw
    if (is.null(x)) x=where;
    legend(x,y,bty=bty,legend=legend,cex=cex,col=col,lwd=lwd,lty=lty,
          title=title,title.col=title.col,...);
  }
## plot multiple legends. adapted from repwr/ragm_legend
## legends is list of legend.args - arguments to base::legend or plotm_legend
## where, x, y are starting position
## others used as defaults in each legend
multi_legend=
  function(legends,legend='right',where=legend,x=NULL,y=NULL,cex=0.8,bty='n',
           title=NULL,title.col='black',col='black',lty='solid',lwd=1,
           ## include legend-related plotm args (except legend defined ablove)
           legend.title=NULL,legend.labels=NULL,labels=NULL,
           ...) {
    default.args=
      list(cex=cex,bty=bty,title=title,title.col=title.col,col=col,lty=lty,lwd=lwd,
           ## legend-related plotm args
           legend.title=legend.title,legend.labels=legend.labels,labels=labels,legend=NULL,
           ...);
    if (is.null(x)) x=where;
    sapply(legends,function(legend.args) {
      if (is.null(legend.args)) return();
      legend.args$x=x;
      legend.args$y=y;
      legend.args=fill_defaults(default.args,legend.args);
    ## handle legend-related plotm args
      ##   legend.title, legend.labels, labels, legend (when legend or labels not set)
      legend.args=within(legend.args,{
        if (is.null(title)) title=legend.title;
        if (is.null(legend)) {
          legend=if(is.null(labels)) legend.labels else labels;
        }
        rm(legend.title,legend.labels,labels);
      });
      if (is.null(legend.args$legend)) return();
      ## use 'graphics::legend' to avoid collision with legend arg
      where.next=do.call(graphics::legend,legend.args);
      ## <<- assigns to variables in outer scope, ie, function scope
      ##   from stackoverflow.com/questions/13640157. Thanks!
      ## could also use, eg, assign('x',where.next$rect$left,,envir=parent.frame(n=3))
      x<<-where.next$rect$left;
      y<<-where.next$rect$top-where.next$rect$h;
    })
    return();
  }

## fill tail of probability density
## adapted from https://stackoverflow.com/questions/45527077. Thx!
fill_tail=function(tail=cq(upper,lower),n,d,d0,d.crit) {
  tail=match.arg(tail);
  x=if(tail=='upper') x=d[d>d.crit] else x=d[d<(-d.crit)]
  y=d_d2t(n,x,d0);
  ## toss in extra x, y values to close polygon
  x=c(min(x),x,max(x));
  y=c(0,y,0);
  ## looks nicer if y stops just short of curve
  y=sapply(y,function(y) max(0,y-1e-3));
  ## do it!
  polygon(x=x,y=y,col='grey',border=NA);
}
## TODO: y.gap doesn't really do what I want...
fill_area=function(n,d0,d,col='grey',y.gap=1e-3) {
  y=d_d2t(n,d,d0);
  ## toss in extra x, y values to close polygon
  x=c(min(d),d,max(d));
  y=c(0,y,0);
  ## looks nicer if y stops just short of curve
  y=sapply(y,function(y) max(0,y-y.gap));
  ## do it!
  polygon(x=x,y=y,col=col,border=NA);
}

pval2col=function(pval) {
  param(col.pval,brk.pval,min.pvcol);
  col.pval[findInterval(-log10(clamp_pval(pval,min.pvcol)),brk.pval,all.inside=TRUE)];
}
d2col=function(n,sd.het,distribution,d) {
  pval=if(distribution=='d2t') d2pval(n,d) else d2htpval(n,sd.het,d);
  pval2col(pval);
}
clamp_pval=function(pval,min.pvcol) sapply(pval,function(pval) max(min(pval,1),min.pvcol));
