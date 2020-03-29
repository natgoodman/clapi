#################################################################################
##
## Author:  Nat Goodman
## Created: 20-02-27
##
## Copyright (C) 2020 Nat Goodman.
## 
## Analyze data extracted by dat_nudge.R from
## "Health Insurance and Mortality: Experimental Evidence from Taxpayer Outreach"
## by Jacob Goldin, Ithai Z. Lurie, and Janet McCubbin, December 201
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
library(BiasedUrn);
## -- Analyze data from taxpayer outreach paper --
doc_nudge=function(sect=NULL,TRN=FALSE,plaus.cutoff=0.05) {
  sect.all=cq(hyper,binom,plotm,plotq);
  if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
  sapply(sect,function(sect) {
    sect_start(sect,sect.all);
##### hyper
    if (sect=='hyper') {
      load_data(nudge);
      ## stage 1. show that NULL hypothesis is implausible
      ##   NULL hypothesis is that treatment has no effect: the differences
      ##   between control and treatment groups is sampling artifcact.
      ## Analyze as selection without replacement: from an urn with
      ##   'treament' and 'control' balls, we select 'cvr.all' subjects;
      ##   what's the probability of getting the observed split?
      ## This is a classic hypergeometric distribution problem. 
      nudge$pval=with(nudge,phyper(cvr.treated,n.treated,n.control,cvr.all,lower=FALSE));
      ## The pvalues are impressively low (4e-141 to 2e-03) indicating the NULL
      ##   hypothesis is quite implausible. The m11 cases has the biggest pvalue;
      ##   arguably, for this case, the NULL is almost plausible
      pval=nudge;
      dotbl(pval);
      ## stage 2. compute selection odds and plausibility limits using
      ##  biased urn, aka noncentral hypergeometric distribution.
      ##  convert odds plaus limits to coverage plaus limits
      nudge$odds=with(nudge,odds_hyper(mu=cvr.treated,m1=n.treated,m2=n.control,n=cvr.all));
      odds=nudge;
      dotbl(odds);
      ## stage 3. compute plausibility limits on odds and coverage
      plodds=with(nudge,plodds_hyper(mu=cvr.treated,m1=n.treated,m2=n.control,n=cvr.all,
                                     plaus.cutoff=plaus.cutoff));
      nudge$odds.lo=plodds[1,]; nudge$odds.hi=plodds[2,];
      cvr.treated.lo=with(nudge,mean_hyper(m1=n.treated,m2=n.control,n=cvr.all,odds=odds.lo));
      cvr.treated.hi=with(nudge,mean_hyper(m1=n.treated,m2=n.control,n=cvr.all,odds=odds.hi));
      nudge$cvr.treated.lo=round(cvr.treated.lo);
      nudge$cvr.treated.hi=round(cvr.treated.hi);
      plaus=nudge;
      dotbl(plaus);
      ## stage 4. try to explain difference between groups by imbalance of
      ##   some unknown confound. compute necessary imbalance assuming
      ##   most extreme case where every subject with confound gets
      ##   coverage
      t1.rate=seq(0.1,1,by=0.1);
      c1.rate=sapply(t1.rate,function(t1.rate)
        with(nudge,
             c1_rate(n.all,n.treated,n.control,cvr.all,cvr.treated,cvr.control,t1.rate)));
      colnames(c1.rate)=t1.rate;
      rownames(c1.rate)=rownames(nudge);
      c1t1.imbalance=sapply(t1.rate,function(t1.rate) c1.rate[,as.character(t1.rate)]/t1.rate);
      colnames(c1t1.imbalance)=t1.rate;
      dotbl(c1t1.imbalance);
      ## TODO. really should use same c1.rate for all nudge-rows
      ## compute pvals for imbalance
      t1.ratem=repr(rbind(t1.rate),nrow(nudge)); # extend t1.rate into matrix
      t1=t1.ratem*nudge$n.treated;               # size of each t1 cell
      c1=c1.rate*nudge$n.control;                # size of each c1 cell
      n1=t1+c1;                                  # number of 'balls' selected
      pval.imbalance=with(nudge,phyper(t1,n.treated,n.control,n1,lower=FALSE));
      colnames(pval.imbalance)=t1.rate;
      dotbl(pval.imbalance);
    }
##### binom
    if (sect=='binom') {
      load_data(nudge);
      ## stage 1. show that NULL hypothesis is implausible
      ##   NULL hypothesis is that treatment has no effect: the differences
      ##   between control and treatment groups is sampling artifcact.
      ## Analyze as selection with replacement: for each subject
      ##   flip coin to decide if 'treated' or 'control'
      ##   do this 'cvr.all' times; what's the probability of getting the observed split?
      ## This is a classic binomial  distribution problem. 
      nudge$pval=with(nudge,pbinom(cvr.treated,cvr.all,n.treated/n.all,lower=FALSE));
      ## Some pvalues are impressively low (4e-141 to 2e-03) indicating the NULL
      ##   hypothesis is quite implausible.
      ## The m6_10 and m11 cases have non-sig p-values (0.1,0.3); arguably, for these cases,
      ##   the NULL is almost plausible. For these cases 'cvr.all' is nearly everyone
      ##   which perforce drives the p-value up
      pval=nudge;
      dotbl(pval);
      ## stage 2. compute selection odds and plausibility limits 
      ##  convert odds plaus limits to coverage plaus limits
      nudge$odds=with(nudge,odds_binom(mu=cvr.treated,m1=n.treated,m2=n.control,n=cvr.all));
      odds=nudge;
      dotbl(odds);
      ## stage 3. compute plausibility limits on odds and coverage
      plodds=with(nudge,plodds_binom(mu=cvr.treated,m1=n.treated,m2=n.control,n=cvr.all,
                                     plaus.cutoff=plaus.cutoff));
      nudge$odds.lo=plodds[1,]; nudge$odds.hi=plodds[2,];
      cvr.treated.lo=with(nudge,mean_binom(m1=n.treated,m2=n.control,n=cvr.all,odds=odds.lo));
      cvr.treated.hi=with(nudge,mean_binom(m1=n.treated,m2=n.control,n=cvr.all,odds=odds.hi));
      nudge$cvr.treated.lo=round(cvr.treated.lo);
      nudge$cvr.treated.hi=round(cvr.treated.hi);
      plaus=nudge;
      dotbl(plaus);
      ## stage 4. try to explain difference between groups by imbalance of
      ##   some unknown confound. compute necessary imbalance assuming
      ##   most extreme case where every subject with confound gets
      ##   coverage
      ## CAUTION: I don't know how to do this right, ie, with multinomial distr
      ##   Instead, approximate with hyper with many more balls in urn
      t1.rate=seq(0.1,1,by=0.1);
      c1.rate=sapply(t1.rate,function(t1.rate)
        with(nudge,
             c1_rate(n.all*100,n.treated*100,n.control*100,
                     cvr.all,cvr.treated,cvr.control,t1.rate)));
      colnames(c1.rate)=t1.rate;
      rownames(c1.rate)=rownames(nudge);
      c1t1.imbalance=sapply(t1.rate,function(t1.rate) c1.rate[,as.character(t1.rate)]/t1.rate);
      colnames(c1t1.imbalance)=t1.rate;
      dotbl(c1t1.imbalance);
      ## TODO. really should use same c1.rate for all nudge-rows
      ## compute pvals for imbalance
      t1.ratem=repr(rbind(t1.rate),nrow(nudge)); # extend t1.rate into matrix
      t1=t1.ratem*nudge$n.treated;               # size of each t1 cell
      c1=c1.rate*nudge$n.control;                # size of each c1 cell
      n1=t1+c1;                                  # number of 'balls' selected
      pval.imbalance=with(nudge,pbinom(round(t1),round(n1),n.treated/n.all,lower=FALSE));
      colnames(pval.imbalance)=t1.rate;
      dotbl(pval.imbalance);
    }
##### plotm. line graphs comparing the distributions
   if (sect=='plotm') {
     load_data(nudge);
     sapply(rownames(nudge),function(case)
       dofig(plotm_nudge,case,nudge=nudge[case,],title='Cumulative probability'));
   }
##### plotq. QQ plots comparing the distributions
    if (sect=='plotq') {
      load_data(nudge);
      sapply(rownames(nudge),function(case)
        dofig(plotq_nudge,case,nudge=nudge[case,],title='Cumulative probability'));
    }
  });
  sect;
}

## TODO (maybe): move into stats.R or stats_nudge.R
## wrappers for BiasedUrn functions. only need Fisher versions
## univariate: fully vectorized odds and mean
odds_hyper=Vectorize(oddsFNCHypergeo);
mean_hyper=Vectorize(meanFNCHypergeo);
## plausibility limit (aka confidence interval) on adds
plodds_hyper=function(mu,m1,m2,n,simplify=TRUE,plaus.cutoff=0.05,interval=c(0.5,1.5)) {
  pl=plodds_hyper_(mu,m1,m2,n,plaus.cutoff,interval);
  if (simplify&ncol(pl)==1) pl=as.vector(pl);
  pl;
}
plodds_hyper_=Vectorize(function(mu,m1,m2,n,plaus.cutoff=0.05,interval=c(0.5,1.5)) {
  p0=plaus.cutoff/2; p1=1-p0;
  lo=suppressWarnings(
    uniroot(function(odds) pFNCHypergeo(mu,m1,m2,n,odds=odds,lower.tail=FALSE)-p0,
            interval=interval)$root);
  hi=suppressWarnings(
    uniroot(function(odds) pFNCHypergeo(mu,m1,m2,n,odds=odds,lower.tail=FALSE)-p1,
            interval=interval)$root);
  c(lo,hi);
},vectorize.args=cq(mu,m1,m2,n,plaus.cutoff))

c1_rate=
  Vectorize(function(n.all,n.treated,n.control,cvr.all,cvr.treated,cvr.control,t1.rate) {
    odds=c(1,100,1,100);
    uniroot(function(c1.rate) {
      c1=round(c1.rate*n.control);
      c0=n.control-c1;
      t1=round(t1.rate*n.treated);
      t0=n.treated-t1;
      mu=meanMFNCHypergeo(m=c(c0,c1,t0,t1),n=cvr.all,odds=odds);
      mu[1]+mu[2]-cvr.control;
    },interval=c(0,1))$root;
  })

##### binom versions of the above. not terribly deep, but... 
odds_binom=function(mu,m1,m2,n) {
  prob=mu/n;
  m2*prob/(m1*(1-prob));
}
mean_binom=function(m1,m2,n,odds=1) {
  m1=odds*m1
  prob=m1/(m1+m2);
  ## round(n*prob);
  n*prob;
}
## plausibility limit (aka confidence interval) on adds
plodds_binom=function(mu,m1,m2,n,simplify=TRUE,plaus.cutoff=0.05,interval=c(0.01,10)) {
  pl=plodds_binom_(mu,m1,m2,n,plaus.cutoff,interval=interval);
  if (simplify&ncol(pl)==1) pl=as.vector(pl);
  pl;
}
plodds_binom_=Vectorize(function(mu,m1,m2,n,plaus.cutoff=0.05,interval=c(0.01,10)) {
  p0=plaus.cutoff/2; p1=1-p0;
  lo=suppressWarnings(
    uniroot(function(odds) {
      m1=m1*odds;
      pbinom(mu,n,m1/(m1+m2),lower.tail=FALSE)-p0;
    },interval=interval)$root);
  hi=suppressWarnings(
    uniroot(function(odds) {
      m1=m1*odds;
      pbinom(mu,n,m1/(m1+m2),lower.tail=FALSE)-p1;
    },interval=interval)$root);
  c(lo,hi);
},vectorize.args=cq(mu,m1,m2,n,plaus.cutoff))

c1_rate_binom=
  Vectorize(function(n.all,n.treated,n.control,cvr.all,cvr.treated,cvr.control,t1.rate) {
    odds=c(1,100,1,100);
    uniroot(function(c1.rate) {
      c1=round(c1.rate*n.control);
      c0=n.control-c1;
      t1=round(t1.rate*n.treated);
      t0=n.treated-t1;
      ## WRONG! has to be multinomial
      mu=mean_binom(m=c(c0,c1,t0,t1),n=cvr.all,odds=odds);
      mu[1]+mu[2]-cvr.control;
    },interval=c(0,1))$root;
  })


## multivariate mean. NOT USED
##   m, odds vectors or matrices - if matrices, each column is a case
##   n single number or vector
mean_hyperm=function(m,n,odds,simplify=TRUE,round=TRUE) {
  ## block of code adapted from misig/R/plot.R:plotm
  if (is.null(m)||is.null(odds)) return(NULL);
  if (is.vector(m)) m=as.matrix(m)
  else if (length(dim(m))!=2) stop("'m' must be vector or 2-dimensional matrix-like object");
  if (is.vector(odds)) odds=as.matrix(odds)
  else if (length(dim(odds))!=2) stop("'odds' must be vector or 2-dimensional matrix-like object");
  if (nrow(m)!=nrow(odds)) stop("'m' and 'odds' must have same number of rows");
  ## like matplot, if both have multiple columns, must be same number of columns
  if (ncol(m)>1&&ncol(odds)>1&&ncol(m)!=ncol(odds))
    stop("When 'm' and 'odds' have multiple columns, must have same number of columns");
  ## recycle smaller to width of larger
  if (ncol(m)>ncol(odds)) odds=repc(odds,length=ncol(m))
  else if (ncol(m)<ncol(odds)) m=repc(m,length=ncol(odds));
  n=rep(n,length=ncol(m));
  mu=do.call(cbind,lapply(seq_len(ncol(m)),function(j) meanMFNCHypergeo(m[,j],n[j],odds[,j])));
  if (round) mu=round(mu);
  if (simplify&ncol(mu)==1) mu=as.vector(mu);
  mu;
}

## try to understand odd behavior of hypergeometric distribution. NOT USED in running code
foo_hyper=function(g1,g2,n,m=param(m.sim)) foo_distr('hyper',g1,g2,n,m);
foo_binom=function(g1,g2,n,m=param(m.sim)) foo_distr('binom',g1,g2,n,m);
foo_distr=function(what=cq(hyper,binom),g1,g2,n,m=param(m.sim)) {
  what=match.arg(what,several.ok=TRUE);
  hits=lapply(what,function(what) {
    file=filename_sim(what,g1,g2,n,m);
    if (!file.exists(file)) dosim(what,g1,g2,n,m);
    sim=load_sim(g1,g2,n,m,what=what,file=file);
    hits=apply(sim,2,function(x) length(which(x<=g1)));
  })
  xlim=range(g1,do.call(c,hits));
  plot(x=NULL,y=NULL,type='n',xlim=xlim,ylim=c(0,1),xlab='hits',ylab='cum prob',
       main=paste(paste(collapse=',',what),'sim and theory:',nvq(g1,g2,n,m)));
  col=cq(blue,green);
  sapply(1:length(hits),function(i) plot(ecdf(hits[[i]]),add=T,col=col[i],pch=19));
  x=xlim[1]:xlim[2];
  y=do.call(cbind,lapply(what,function(what)
    if(what=='hyper') phyper(x,g1,g2,n) else pbinom(x,n,g1/(g1+g2))));
  matlines(x,y,type='l',lty='solid',col=col);
  grid();
  abline(v=g1,lty='dotted',col='red');
  legend('right',legend=what,title='distribution',col=col,lty='solid',pch=19,seg.len=4,bty='n');
}
theo_distr=function(what=cq(hyper,binom),g1,g2,n,m=NULL,extend.x=TRUE) {
  what=match.arg(what,several.ok=TRUE);
  col=cq(blue,green);
  q=c(0.001,0.999)                      # quantile bounds for computing probs 
  xlim=range(do.call(c,lapply(what,function(what)
    if(what=='hyper') qhyper(q,g1,g2,n) else qbinom(q,n,g1/(g1+g2)))));
  if (extend.x) xlim=range(xlim,g1);
  x=xlim[1]:xlim[2];
  y=do.call(cbind,lapply(what,function(what)
    if(what=='hyper') phyper(x,g1,g2,n) else pbinom(x,n,g1/(g1+g2))));
  matplot(x,y,type='l',lty='solid',col=col,xlim=xlim,ylim=c(0,1),xlab='hits',ylab='cum prob',
          main=paste(paste(collapse=',',what),nvq(g1,g2,n)));
  grid();
  abline(v=g1,lty='dotted',col='red');
  legend('right',legend=what,title='distribution',col=col,lty='solid',seg.len=4,bty='n');
}
## for now, just one distribution, multiple n
theo_multi=function(what=cq(hyper,binom),g1,g2,n,m=NULL,extend.x=TRUE) {
  ## what=match.arg(what,several.ok=TRUE);
  what=match.arg(what);
  col=cq(blue,green);
  col=setNames(RColorBrewer::brewer.pal(max(3,length(n)),'Set1'),n);
  q=c(0.001,0.999)                      # quantile bounds for computing probs 
  ## xlim=range(do.call(c,lapply(what,function(what)
  ##   if(what=='hyper') qhyper(q,g1,g2,n) else qbinom(q,n,g1/(g1+g2)))));
  xlim=range(do.call(c,lapply(n,function(n)
    if(what=='hyper') qhyper(q,g1,g2,n) else qbinom(q,n,g1/(g1+g2)))));
  if (extend.x) xlim=range(xlim,g1);
  x=xlim[1]:xlim[2];
  ## y=do.call(cbind,lapply(what,function(what)
  ##   if(what=='hyper') phyper(x,g1,g2,n) else pbinom(x,n,g1/(g1+g2))));
  y=do.call(cbind,lapply(n,function(n)
    if(what=='hyper') phyper(x,g1,g2,n) else pbinom(x,n,g1/(g1+g2))));
  matplot(x,y,type='l',lty='solid',col=col,xlim=xlim,ylim=c(0,1),xlab='hits',ylab='cum prob',
          main=paste(paste(collapse=',',what),nvq(g1,g2)));
  grid();
  abline(v=g1,lty='dotted',col='red');
  ## legend('right',legend=what,title='distribution',col=col,lty='solid',seg.len=4,bty='n');
  legend('right',legend=n,title='n',col=col,lty='solid',seg.len=4,bty='n');
}

## NOT USED
foo_=Vectorize(function(n.all,n.treated,n.control,cvr.all,cvr.treated,cvr.control,t1.rate) {
  odds=c(1,100,1,100);
  uniroot(function(c1.rate) {
    c1=round(c1.rate*n.control);
    c0=n.control-c1;
    t1=round(t1.rate*n.treated);
    t0=n.treated-t1;
    mu=meanMFNCHypergeo(m=c(c0,c1,t0,t1),n=cvr.all,odds=odds);
    mu[1]+mu[2]-cvr.control;
  },interval=c(0,1))$root;
})
## NOT USED
foo=function(t1.rate=seq(0.1,1,by=0.1),odds=c(1,100,1,100)) {
  m=lapply(t1.rate,function(t1.rate) {
    m=do.call(cbind,withrows(nudge,case,{
      c1.rate=foo_(n.all,n.treated,n.control,cvr.all,cvr.treated,cvr.control,t1.rate);
      c1=round(c1.rate*n.control);
      c0=n.control-c1;
      t1=round(t1.rate*n.treated);
      t0=n.treated-t1;
      m=c(c0,c1,t0,t1);
    }))
    rownames(m)=cq(c0,c1,t0,t1);
    colnames(m)=rownames(nudge);
    m;
    ## mean_hyperm(m=c(c0,c1,t0,t1),n=cvr.all,odds=odds,simplify=FALSE);
  })
  names(m)=t1.rate;
}
