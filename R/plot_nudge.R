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
## plotm_distr. plot hyper, binom for multiple cases
## plotm_nudge (below) is wrapper that pulls data from nudge data frame
## adapted from misig/plotmg, in turn adapted from repwr/plot_ragm and others
## distr - 'hyper', 'binom' or both
## distribution params - can be vectors
##   g1,g2 - number of elements in the two groups
##   n - number of elements selected
## extend.x - extend xlim to g1 - CAUTION: often scrunches informative part of plot
## qlim - quantile bounds for distributions - used to set xlim
## col, lty, lwd - the usual line properties
## palette - RColorBrewer sequential palette name
## title, cex.title - title and cex for title
## type - plot type - passed to matplot. 'n' also turns off extra lines, legend, grid
##   TODO: type='p' should cause legend to draw points instead of lines
## vline,hline are vectors of x or y positions for extra vertical or horizontal lines
## vhlty, vhcol, vhlwd are lty, col, lwd for these extra lines
## vlab, hlab contol writing vline, hline values along axes
## vhdigits is number of digits for these values
## legend tells whether and where to draw legend
##   TRUE or word like 'right' means draw legend, FALSE or NULL means no legend
## legend.title is legend title
## more TBD
plotm_distr=
  function(g1,g2,n,casenames=NULL,distr=cq(hyper,binom),extend.x=FALSE,
           qlim=c(0.001,0.999),xlim=NULL,
           col=NULL,lty=cq(solid,dotted),lwd=1,type='l',
           palettes='Set1',
           title='Cumulative probability',title.distr='auto',title.vars='auto',cex.title='auto',
           xlab='hits',ylab='cum prob',
           vline=g1,hline=NULL,vhlty='dashed',vhcol='grey50',
           vhlwd=1,vlab=TRUE,hlab=TRUE,vhdigits=2,
           legend='right',legend.title=NULL,legend.labels=NULL,legend.vars='auto',
           legend.lty='solid',legend.lwd=1,legend.cex=0.8,
           ...) {
    distr=match.arg(distr,several.ok=TRUE);
    cases=data.frame(g1=g1,g2=g2,n=n);
    col=col_brew(palettes,n=nrow(cases));
    lty=setNL(lty,distr);
    lwd=setNL(lwd,distr);
    if (is.null(xlim)) {
      ## use merge to expand data frames. from stackoverflow.com/questions/11693599. Thanks!
      cases.xl=merge(cases,expand.grid(q=qlim,distr=distr));
      xlim=with(cases.xl,range(ifelse(distr=='hyper',qhyper(q,g1,g2,n),qbinom(q,n,g1/(g1+g2)))));
      if (extend.x) xlim=range(xlim,g1);
    }
    x=xlim[1]:xlim[2];
    y=do.call(cbind,lapply(distr,function(distr)
      do.call(cbind,
              if (distr=='hyper') withrows(cases,case,phyper(x,g1,g2,n))
              else withrows(cases,case,pbinom(x,n,g1/(g1+g2))))));
    col=rep(col,length(distr))
    lty=rep(lty[distr],each=nrow(cases));
    lwd=rep(lwd[distr],each=nrow(cases));
    ## add 'distr' to title if desired
    if (title.distr=='auto') title=paste(title,'for',paste(collapse=', ',distr));
    ## put constant vars in title, varying ones in legend
    count=apply(cases,2,function(x) length(unique(x)));
    if (title.vars=='auto') {
      if (length(casenames)==1) title=paste0(title,': case=',casenames)
      else {
        title.vars=names(count)[count==1];
        if (length(title.vars)!=0) title=paste0(title,': ',nv(names=title.vars,SEP=', '));
      }
    }
    if (legend.vars=='auto') {
      if (length(casenames)>1) legend.labels=casenames
      else {
        legend.vars=names(count)[count>1];
        legend.labels=unlist(
          withrows(cases,case,
                   paste(collapse=' ',c(legend.labels,nv(names=legend.vars,SEP=', ')))));
      }
    }
    if (cex.title=='auto') cex.title=cex_title(title);
    matplot(x,y,lty=lty,lwd=lwd,col=col,xlim=xlim,ylim=c(0,1),
            main=title,cex.main=cex.title,xlab=xlab,ylab=ylab,type=type);
    if (type!='n') {
      grid();
      ## plot extra lines & values if desired. nop if vline, hline NULL
      vhline(vline=vline,hline=hline,vlab=vlab,hlab=hlab,vhdigits=vhdigits,
             lty=vhlty,col=vhcol,lwd=vhlwd);
      ## draw legend if desired
      if (is.null(legend)) legend=FALSE
      else if (!is.logical(legend)) {
        where=legend;
        legend=TRUE;
      }
      if (legend) {
        legends=NULL;
        legend.cases=
          if (nrow(cases)<=1) NULL
          else list(title='case',labels=legend.labels,col=col,lty=legend.lty,lwd=legend.lwd);
        legend.distrs=
          if (length(distr)==1) NULL
          else list(title='distr',labels=distr,lty=lty[distr],lwd=lwd[distr]);
        add_legend(list(legend.cases,legend.distrs),where=where,cex=legend.cex);
      }
    }
    return();
  }
## plotm_nudge. wrapper for plotm_distr that pulls data from nudge data frame
plotm_nudge=function(nudge,...) 
  with(nudge,plotm_distr(g1=n.treated,g2=n.control,n=cvr.all,casenames=rownames(nudge),...))

## plotq_distr. QQ plot of hyper vs binom (or vice versa) for multiple cases
## plotq_nudge (below) is wrapper that pulls data from nudge data frame
## distr - 'hyper', 'binom' or both
## distribution params - can be vectors
##   g1,g2 - number of elements in the two groups
##   n - number of elements selected
## qlim,qby,q - probs for which function computes quantiles
## col, pch, cex - the usual line properties
## palette - RColorBrewer sequential palette name
## title, cex.title - title and cex for title
## type - plot type - passed to matplot. 'n' also turns off extra lines, legend, grid
## vline,hline are vectors of x or y positions for extra vertical or horizontal lines
## vhlty, vhcol, vhlwd are lty, col, lwd for these extra lines
## vlab, hlab contol writing vline, hline values along axes
## vhdigits is number of digits for these values
## legend tells whether and where to draw legend
##   TRUE or word like 'right' means draw legend, FALSE or NULL means no legend
## legend.title is legend title
## more TBD
plotq_distr=
  function(g1,g2,n,casenames=NULL,xdistr='hyper',ydistr='binom',
           ## qlim=c(0.001,0.999),qby=qlim[1],q=seq(qlim[1],qlim[2],by=qby),
           qlim=c(0.01,0.99),qby=qlim[1],q=seq(qlim[1],qlim[2],by=qby),
           col=NULL,pch=19,cex=0.8,type='p',
           palettes='Set1',
           title='QQ plot',title.distr='auto',title.vars='auto',cex.title='auto',
           xlab=xdistr,ylab=ydistr,
           vline=NULL,hline=NULL,vhlty='dashed',vhcol='grey50',
           vhlwd=1,vlab=TRUE,hlab=TRUE,vhdigits=2,
           legend='right',legend.title=NULL,legend.labels=NULL,legend.vars='auto',
           legend.pch=19,legend.cex=0.8,
           ...) {
    distr=cq(hyper,binom);
    distr.mess=paste(collapse=', ',distr);
    if (xdistr %notin% distr) stop(paste('xdistr must be one of',distr.mess,'not',xdistr));
    if (ydistr %notin% distr) stop(paste('ydistr must be one of',distr.mess,'not',ydistr));
    cases=data.frame(g1=g1,g2=g2,n=n);
    col=col_brew(palettes,n=nrow(cases));
    x=qdistr(xdistr,q,g1,g2,n);
    y=do.call(cbind,withrows(cases,case,qdistr(ydistr,q,g1,g2,n)));
    ## add 'distr' to title if desired
    if (title.distr=='auto') title=paste(title,'of',ydistr,'vs',xdistr);
    ## put constant vars in title, varying ones in legend
    count=apply(cases,2,function(x) length(unique(x)));
    if (title.vars=='auto') {
      if (length(casenames)==1) title=paste0(title,': case=',casenames)
      else {
        title.vars=names(count)[count==1];
        if (length(title.vars)!=0) title=paste0(title,': ',nv(names=title.vars,SEP=', '));
      }
    }
    if (legend.vars=='auto') {
      if (length(casenames)>1) legend.labels=casenames
      else {
        legend.vars=names(count)[count>1];
        legend.labels=unlist(
          withrows(cases,case,
                   paste(collapse=' ',c(legend.labels,nv(names=legend.vars,SEP=', ')))));
      }
    }
    if (cex.title=='auto') cex.title=cex_title(title);
    matplot(x,y,col=col,pch=pch,cex=cex,
            main=title,cex.main=cex.title,xlab=xlab,ylab=ylab,type=type);
    if (type!='n') {
      grid();
      abline(a=0,b=1,lty='dotted',col='grey');
      ## plot extra lines & values if desired. nop if vline, hline NULL
      vhline(vline=vline,hline=hline,vlab=vlab,hlab=hlab,vhdigits=vhdigits,
             lty=vhlty,col=vhcol,lwd=vhlwd);
      ## draw legend if desired
      if (is.null(legend)) legend=FALSE
      else if (!is.logical(legend)) {
        where=legend;
        legend=TRUE;
      }
      if (legend&&nrow(cases)>1) 
        add_legend(where=where,title='case',labels=legend.labels,col=col,
          pch=legend.pch,cex=legend.cex);
    }
    return();
  }
## plotq_nudge. wrapper for plotq_distr that pulls data from nudge data frame
plotq_nudge=function(nudge,...) 
  with(nudge,plotq_distr(g1=n.treated,g2=n.control,n=cvr.all,casenames=rownames(nudge),...))

## TODO: maybe move this to general 'stats_nudge' along with wrappers for other distr functions
qdistr=Vectorize(function(distr=cq(hyper,binom),q,g1,g2,n,...) {
  distr=match.arg(distr);
  if (distr=='hyper') qhyper(q,g1,g2,n,...) else qbinom(q,n,g1/(g1+g2),...);
},vectorize.args=cq(g1,g2,n))
