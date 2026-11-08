# =============================================================================
# old_linkMap_functions.R
# Archive of the original base-graphics linkMap.cross(), linkMap.wgaim() and
# linkMap.default() removed from wgaim16.R on 2025-05-09.
# Superseded by:
#   R/linkMap.cross.R   -- ggplot2 linkMap.cross()
#   R/linkMap.QTLAim.R  -- ggplot2 linkMap.QTLAim()
# Kept here for reference only. Do not source().
# =============================================================================


linkMap <- function(object, ...)
    UseMethod("linkMap")

linkMap.cross <- function(object, chr, chr.dist, marker.names = "markers", tick = FALSE, squash = TRUE, m.cex = 0.6, ...){
  circ <- function(x, y, shiftx = 0, shifty = 0, ely = 1, elx = 1)
    ((x - shiftx)^2)/elx + ((y - shifty)^2)/ely
  dots <- list(...)
  old.xpd <- par("xpd")
  par(xpd = TRUE)
  on.exit(par(xpd = old.xpd))
  map <- pull.map(object)
  if(!missing(chr)) {
    if(any(is.na(pmatch(chr, names(map)))))
      stop("Some names of chromosome(s) subset do not match names of map.")
    map <- map[chr]
  }
  n.chr <- length(map)
  mt <- list()
  if(!missing(chr.dist)) {
    if(!all(names(chr.dist) %in% c("start","end")))
        stop("names of chr.dist must be \"start\" and/or \"end\"")
    dname <- names(chr.dist)
     if(!is.null(chr.dist$start)) {
        if(length(chr.dist$start) == 1)
            chr.dist$start <- rep(chr.dist$start, length(map))
        if(length(chr.dist$start) != length(map))
            stop("Length of user specified chromosome starting distances need to be 1 or equal to the number of chromosomes")
    } else chr.dist$start <- 0
    if(!is.null(chr.dist$end)) {
        if(length(chr.dist$end) == 1)
            chr.dist$end <- rep(chr.dist$end, length(map))
        if(length(chr.dist$end) != length(map))
            stop("Length of user specified chromosome ending distances need to be 1 or equal to the number of chromosomes")
    } else chr.dist$end <- unlist(lapply(map, max))
    chr.dist <- do.call("cbind.data.frame", chr.dist)
    nm <- names(map)
    map <- lapply(1:nrow(chr.dist), function(el, map, chr.dist){
        tmap <- map[[el]]
        if(max(tmap) < chr.dist$start[el]) NULL
        else tmap[(tmap >= chr.dist$start[el]) & (tmap <= chr.dist$end[el])]
    }, map, chr.dist)
    names(map) <- nm
    map <- map[!sapply(map, is.null)]
    n.chr <- length(map)
  }
  maxlen <- max(unlist(lapply(map, max)))
  minlen <- min(unlist(lapply(map, min)))
  omap <- map
  if(is.null(marker.names)) {
    chrpos <- 1:n.chr
    thelim <- range(chrpos) + c(-0.5, 0.5)
  }
  else {
    if(!is.null(fmark <- attr(object, "flanking"))){
        map <- lapply(map, function(el, fmark){
            el <- el[names(el) %in% fmark]
            if(length(el) == 0) NULL
            else el
        }, fmark)
        omap <- omap[!sapply(map, is.null)]
        map <- map[!sapply(map, is.null)]
        n.chr <- length(map)
    }
    if(all(is.na(pmatch(marker.names, c("markers","dist")))))
      stop("marker.names argument must be either \"dist\", or \"markers\".")
    if(!is.na(pmatch("cex", names(dots))))
      dots$cex <- NULL
#    else cex <- par("cex")
    if(!squash)
      chrpos <- seq(1, n.chr * 3, by = 3)
    else
      chrpos <- seq(1, n.chr * 2, by = 2)
    thelim <- range(chrpos) + c(-1.6, 1.35)
    for(i in 1:n.chr){
        mt[[i]] <- map[[i]]
        if(length(mt[[i]]) > 1){
            conv <- par("pin")[2]/maxlen
            for(j in 1:(length(mt[[i]]) - 1)){
                ch <- mt[[i]][j + 1]*conv - (mt[[i]][j]*conv + 10*par("csi")*m.cex/9)
                if(ch < 0){
                    temp <- mt[[i]][j + 1]*conv + abs(ch)
                    mt[[i]][j + 1] <- temp/conv
                }
            }
        }
    }
    maxlen <- max(c(unlist(lapply(omap, max)),unlist(lapply(mt, max))))
    names(mt) <- names(map)
  }
  plot(0, 0, type = "n", ylim = c(maxlen, minlen), xlim = thelim,
       xaxs = "i", ylab = "Location (cM)", xlab = "Chromosome",
       axes = FALSE, ...)
  axis(side = 2,  ylim = c(maxlen, minlen))
  pins <- par()$plt
  for(i in 1:n.chr) {
      if(!is.null(marker.names)) {
          if(marker.names == "dist")
              alis <- list(x = chrpos[i] + 0.50, y =  mt[[i]], labels = as.character(round(map[[i]], 2)),
                           adj = c(0, 0.5), cex = m.cex)
          else
              alis <- list(x = chrpos[i] + 0.50, y = mt[[i]], labels = names(map[[i]]), adj = c(0, 0.5), cex = m.cex)
          do.call("text", c(alis, dots))
          segments(chrpos[i] + 0.25, map[[i]], chrpos[i] + 0.3, map[[i]])
          segments(chrpos[i] + 0.3, map[[i]], chrpos[i] + 0.4, mt[[i]])
          segments(chrpos[i] + 0.40, mt[[i]], chrpos[i] + 0.45, mt[[i]])
    }
    map[[i]] <- omap[[i]]
    barl <- chrpos[i]- 0.03
    barr <- chrpos[i]+ 0.03
    segments(barl, min(map[[i]]), barl, max(map[[i]]), lwd = 1)
    segments(barr, min(map[[i]]), barr, max(map[[i]]), lwd = 1)
    segments(barl - 0.17, map[[i]], barr + 0.17, map[[i]])
    # attempt to put curves at ends of chromosomes
    xseq <- seq(barl, barr, length = 20) - chrpos[i]
    yseq <- circ(xseq, xseq, ely = 1, elx = 0.07/maxlen)
    yseq <- yseq - max(yseq)
    lines(xseq + chrpos[i], min(map[[i]]) + yseq)
    lines(xseq + chrpos[i], max(map[[i]]) - yseq)
  }
  axis(side = 1, at = chrpos, labels = names(map), tick = FALSE)
  if(is.na(pmatch("main", names(dots))) & !as.logical(sys.parent()))
    title("Genetic Map")
  invisible(list(mt = mt, map = map, chrpos = chrpos))
}

linkMap.wgaim <- function (object, intervalObj, chr, chr.dist, marker.names = "markers", flanking = TRUE, list.col = list(q.col = "light blue", m.col = "red", t.col = "light blue"), list.cex = list(t.cex = 0.6, m.cex = 0.6), trait.labels = NULL, tick = FALSE, ...)
{
    dots <- list(...)
    if (missing(intervalObj))
        stop("intervalObj is a required argument")
    if (!inherits(intervalObj, "cross"))
        stop("intervalObj is not of class \"cross\"")
    if (!length(wchr <- object$QTL$effects)) {
        warning("There are no significant QTL's. Plotting map only...")
        linkMap(intervalObj, chr, chr.dist, marker.names = marker.names,
            tick = tick, squash = FALSE, ...)
        return(invisible())
    }
    qtlm <- getQTL(object, intervalObj)
    wchr <- qtlm[,1]
    if(object$QTL$type == "interval")
        qtlm <- qtlm[,3:6]
    else
        qtlm <- cbind(qtlm[,3:4],qtlm[,3:4])
    if (is.null(list.col$q.col))
        list.col$q.col <- "light blue"
    if (is.null(list.col$t.col))
        list.col$t.col <- list.col$q.col
    if (missing(chr))
        chr <- unique(wchr)[order(unique(wchr))]
    if(is.null(list.cex$m.cex))
        list.cex$m.cex <- 0.6
    if(is.null(list.cex$t.cex))
        list.cex$t.cex <- 0.6
    if(flanking)
        attr(intervalObj, "flanking") <- unique(c(qtlm[,1],qtlm[,3]))
    lmap <- linkMap(intervalObj, chr, chr.dist, marker.names = marker.names,
        tick = tick, squash = TRUE, m.cex = list.cex$m.cex, ...)
    map <- lmap$map
    if (is.null(trait <- trait.labels))
        trait <- rep(as.character(object$call$fixed[[2]]), length(wchr))
    if (length(trait.labels) == 1)
        trait <- rep(trait.labels, length(wchr))
    qtrait <- unique(trait)
    qtlm <- cbind.data.frame(qtlm, trait = factor(trait, levels = unique(trait)))
    qtlm[, 2] <- as.numeric(as.character(qtlm[, 2]))
    qtlm[, 4] <- as.numeric(as.character(qtlm[, 4]))
    qtlList <- lapply(split(qtlm, wchr), function(el) el[order(el[,
        2]), ])
    qtlm <- do.call("rbind", qtlList)
    wchr <- wchr[order(wchr)]
    chr <- names(map)
    if (!missing(chr)) {
        oc <- wchr %in% chr
        if(!all(oc)) {
            warning("Some QTL's exist outside chromosome(s) subset, Omitting QTL's....")
            qtlm <- qtlm[oc,]
            wchr <- wchr[oc]
        }
    }
    if (!missing(chr.dist)) {
        mins <- unlist(lapply(map, min))
        maxs <- unlist(lapply(map, max))
        mins <- mins[wchr]; maxs <- maxs[wchr]
        om <- qtlm[, 4] < maxs & qtlm[,2] > mins
        if (!all(om)) {
            warning("Some QTL regions outside distances specified. Omitting QTL's....")
            qtlm <- qtlm[om, ]
            wchr <- wchr[om]
        }
    }
    n.chr <- length(map)
    maxlen <- max(unlist(lapply(map, max)))
    chrpos <- lmap$chrpos
    mt <- lmap$mt
    if (!is.na(cind <- pmatch("col", names(dots))))
        dots <- dots[-cind]
    if (is.null(dim(qtlm)))
        qtlm <- matrix(qtlm, nrow = 1, byrow = FALSE)
    tlis <- list()
    if (!is.na(pmatch("cex", names(dots))))
        p.cex <- dots$cex
    else p.cex <- par("cex")
    dots$cex <- NULL
    for (i in 1:n.chr) {
        conv <- par("pin")[2]/maxlen
        ind <- wchr %in% names(map)[i]
#        if (as.logical(length(ind <- wchr %in% names(map)[i]))) {
         if(any(ind)){
            tlis[[i]] <- as.vector(qtlm[ind, 2] + qtlm[ind, 4])/2
            names(tlis[[i]]) <- as.character(qtlm[ind, 5])
            if (length(tlis[[i]]) > 1) {
                for (j in 1:(length(tlis[[i]]) - 1)) {
                  ch <- tlis[[i]][j + 1] * conv - (tlis[[i]][j] *
                    conv + 10 * par("csi") * list.cex$t.cex/9)
                  if (ch < 0) {
                    temp <- tlis[[i]][j + 1] * conv + abs(ch)
                    tlis[[i]][j + 1] <- temp/conv
                  }
                }
            }
        }
    }
    qtld <- qtlm[, 1:4]
    nodup <- !duplicated(do.call("paste", qtld))
    qtls <- qtld[nodup, ]
    whd <- pmatch(do.call("paste", qtld), do.call("paste", qtls),
        duplicates.ok = TRUE)
    dlis <- split(as.character(qtlm[, 5]), whd)
    qtlm <- qtls
    wchr <- wchr[nodup]
    for (i in 1:n.chr) {
        ind <- wchr %in% names(map)[i]
        if(any(ind)){
            ind <- (1:length(wchr))[ind]
            for (j in ind) {
                if (!is.null(marker.names)) {
                    wh <- mt[[i]][c(as.character(qtlm[j, 1]), as.character(qtlm[j, 3]))]
                  if (marker.names == "dist") {
                    dist <- map[[i]][c(as.character(qtlm[j, 1]), as.character(qtlm[j, 3]))]
                    dlabs <- as.character(round(as.numeric(unlist(dist)), 2))
                    alis <- list(x = chrpos[i] + 0.5, y = wh,
                      labels = dlabs, adj = c(0, 0.5), col = list.col$m.col, cex = list.cex$m.cex)
                  }
                  else alis <- list(x = chrpos[i] + 0.5, y = wh,
                    labels = names(wh), adj = c(0, 0.5), col = list.col$m.col, cex = list.cex$m.cex)
                  do.call("text", c(alis, dots))
                }
                yv <- c(qtlm[j, 2], qtlm[j, 4])
                yv <- c(yv, rev(yv))
                dind <- dlis[[j]]
                q.cols <- list.col$q.col[pmatch(dind, qtrait)]
                qind <- 1:length(dind)
                if (length(dlis[[j]]) > 1) {
                  int <- seq(chrpos[i] - 0.2, chrpos[i] + 0.2,
                    length = length(dind) + 1)
                  for (k in 1:length(dind)) {
                    xv <- c(rep(int[k], 2), rep(int[k + 1], 2))
                    if(object$QTL$type == "interval")
                        polygon(xv, y = yv, border = NA, col = q.cols[k])
                    else {
                        plis <- list(x = (xv[1] + xv[3])/2,y = yv[1], col = q.cols[k], cex = p.cex)
                        do.call("points", c(plis, dots))
                    }
                  }
                }
                else {
                  xv <- c(rep(chrpos[i] - 0.2, 2), rep(chrpos[i] +
                    0.2, 2))
                  if(object$QTL$type == "interval")
                      polygon(xv, y = yv, border = NA, col = q.cols)
                  else {
                      plis <- list(x = chrpos[i], y = yv[1], col = q.cols, cex = p.cex)
                      do.call("points", c(plis,dots))
                  }
                }
                segments(chrpos[i] - 0.25, yv[1], chrpos[i] -
                  0.25, yv[2])
                segments(chrpos[i] - 0.25, sum(yv[1:2])/2, chrpos[i] -
                  0.3, sum(yv[1:2])/2)
                segments(chrpos[i] - 0.3, sum(yv[1:2])/2, chrpos[i] -
                  0.4, tlis[[i]][qind])
                segments(chrpos[i] - 0.4, tlis[[i]][qind], chrpos[i] -
                  0.45, tlis[[i]][qind])
                if (length(list.col$t.col) > 1)
                  t.cols <- list.col$t.col[pmatch(dind, qtrait)]
                else t.cols <- list.col$t.col
                text(chrpos[i] - 0.5, tlis[[i]][qind], names(tlis[[i]][qind]),
                  adj = c(1, 0.3), col = t.cols, cex = list.cex$t.cex)
                tlis[[i]] <- tlis[[i]][-qind]
            }
        }
        segments(chrpos[i] - 0.2, map[[i]], chrpos[i] + 0.2,
            map[[i]])
    }
    if (is.na(pmatch("main", names(dots))))
        title("Genetic Map with QTL")
}

linkMap.default <- function (object, intervalObj, chr, chr.dist, marker.names = "markers", flanking = TRUE, list.col = list(q.col = rainbow(length(object)), m.col = "red", t.col = rainbow(length(object))),
    list.cex = list(m.cex = 0.6, t.cex = 0.6), trait.labels = NULL, tick = FALSE, ...)
{
    old.par <- par(no.readonly = TRUE)
    par(mar = c(5, 4, 5, 2) + 0.1)
    par(xpd = TRUE)
    on.exit(par(old.par))
    dlist <- list()
    lclass <- lapply(object, function(el) {
        if (!inherits(el, "wgaim"))
            stop("list objects need to inherit from class \"wgaim\"")
        class(el)[1]
    })
    lclass <- unique(unlist(lclass))
    if (any(is.na(pmatch(lclass, "wgaim"))))
        stop("link.map method is only for list objects of class \"wgaim\"")
    type <- unique(unlist(lapply(object, function(el) el$QTL$type)))
    if(length(type) > 1)
        stop("Models need to have same analyses \"type\".")
    dlist$QTL$type <- type
    effects <- lapply(object, function(el) el$QTL$effects)
    len <- lapply(effects, length)
    if (!is.null(trait.labels)) {
        if (length(trait.labels) != length(object))
            stop("Length of trait labels does not equal number of models specified.")
    }
    else trait.labels <- unlist(lapply(object, function(el) as.character(el$call$fixed[[2]])))
    tt <- table(trait.labels)
    tt <- tt[tt > 1]
    if(length(tt)){
        tn <- names(tt)
        for(i in length(tt)) {
            whl <- trait.labels %in% tn[i]
            trait.labels[whl] <- paste(trait.labels[whl], 1:tt[i], sep = "")
        }
    }
    trait <- rep(trait.labels, times = len)
    dlist$QTL$effects <- unlist(effects)
    class(dlist) <- "wgaim"
    if (length(list.col$q.col) != length(object)) {
        warning("QTL colours not of the same length as the number of traits,\n Choosing \"q.col = rainbow(",
            length(object), ")\".")
        list.col$q.col <- rainbow(length(object))
    }
    if (length(list.col$t.col) != length(object)) {
        if (length(list.col$t.col) != 1) {
            warning("Inappropriate length for trait name colours, using QTL colours")
            list.col$t.col <- list.col$q.col
        }
    }
    linkMap(dlist, intervalObj, chr, chr.dist, marker.names = marker.names,
        list.col = list.col, list.cex = list.cex, trait.labels = trait, tick = tick, flanking = flanking, ...)
}
