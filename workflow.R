# Analysis outline for "Systems-Level Annotation..."


# Load Data
library(mzR)
ramp = mzR:::rampOpen(file)
data = mzR:::rampRawData(ramp)
mzR:::rampClose(ramp)


# Detect ROIs
dengroup.ppm = function (features, ppm, rtwid, minlength) {
  l = nrow(features)
  
  groups = vector(mode = "list", length = l)
  action = vector(mode = 'numeric', length = l)
  unchanged = vector(mode = 'numeric', length = l)
  
  groups[[1]] = seq(l)
  action[1] = 0
  unchanged[1] = 0
  tg = 0
  i = 1
  groupmaxindex = 1
  repeat {
    if (i %% 100 == 0) { cat("\rProgress indicator:", groupmaxindex - i, "groups remaining to analyze. (This will grow before starting to decrease.)            ") }
    
    tg = groups[[i]]
    tg <<- tg
    
    if (length(tg) >= minlength - 1) {
      
      if (action[i] == 0) {
        ng = kernelsplitmass(features[tg,,drop=F], ppm)
        action[i] = 1
      } else {
        ng = kernelsplitrt(features[tg,,drop=F], rtwid)
        action[i] = 0
      }
      
      if (max(ng) > 1) {
        newgroups = split(tg, ng) %>% unname
        ngroups = length(newgroups)
        
        indices_to_replace = c(i,(groupmaxindex+1):(groupmaxindex+ngroups-1))
        
        groups[indices_to_replace] = unname(newgroups)
        action[indices_to_replace] = action[i]
        unchanged[indices_to_replace] = 0
        groupmaxindex = groupmaxindex + ngroups - 1
        
      } else if (unchanged[i] == 3) {
        i = i+1
        
      } else {
        unchanged[i] = unchanged[i] + 1
      }
    } else { i = i+1 }
    
    if ((i == groupmaxindex & i > 1) | (i == groupmaxindex & unchanged[i] > 2)) {
      cat("\rProgress indicator: 0 groups remaining to analyze. (This will grow before starting to decrease.)            ")
      return(groups[1:groupmaxindex])
    }
  }
}


# Calculate Baselines per ROI
library(baseline)
for (roi in rois) {
  baseline(roi, lambda1, lambda2, method = "irls")
}


# Initialize Peaks
# Implementation of Scale-Space Peak Picking, Antoine Liutkus, https://hal.inria.fr/hal-01103123/file/SSPP.pdf
sspp = function(v, S = 1:20, maxdist = 5, do.plot = T) {
  
  dC = C = vector(mode = "numeric", length=length(v))
  O = NULL
  N = length(v)
  v.original = v.prev = v
  
  #dC.m Might ofer additional diagnostic information.  Could exclude peaks that moved around.  Probably just unnecesarily complex
  #dC.m = matrix(ncol = length(dC), nrow = length(S))
  
  for (i in seq_along(S)) {
    s = S[i]
    #s = round(max(c(1,i*N / (10 * S))))
    
    W = signal::hamming(s)
    W = W/sum(W)
    
    v.prev = v
    
    v = as.numeric(stats::filter(v,W))
    v[is.na(v)] = v.prev[is.na(v)]
    
    P = localMaxima(v)
    
    if (length(P) == 0) break
    
    if (i == 1) {
      dC[P] = v[P]
      O = P
    } else {
      dC[] = 0
      
      dist = abs(outer(O, P, "-")) < maxdist
      
      ints = matrix(v.original[O], ncol = ncol(dist), nrow = nrow(dist))
      ints[!dist] = NA
      maxints = matrix(matrixStats::colMaxs(ints, na.rm=T), ncol = ncol(ints), nrow = nrow(ints), byrow = T)
      
      neighbors = which(ints == maxints, arr.ind =T)
      
      dC[O[neighbors[,1]]] = v[P[neighbors[,2]]] * s^1.4
    }
    C = C + dC
    #dC.m[i,] = dC
    #O = union(O, P)
  }
  
  if (do.plot) {
    scores = cbind(which(C>1), C[C>1])
    scores[,2] = scores[,2] * max(v.original)/max(scores[,2])
    
    par(mfcol=c(1,1))
    plot(v.original, type="l")
    lines(scores, type="h", col = "red")
    lines(v, type="l", col="grey")
    points(y=v.original[O], x=O, col="red")
  }
  
  C
}


# Fit Model to Seeded Peaks
library(sn)
fit = function(p, x, y) {
  yhat = p[4] * suppressWarnings( dsn(rep(x, each = ncol(parmat)), xi = p[1], omega = p[2], alpha = p[3]) )
  sum((y-yhat)^2)
  }

for (seed in seeds) {
  optim(parvec, fit(p, scans, roi), method = "L-BFGS-B", control = list(factr = 1E-3, ndeps = c(0.01, .01, .01, .01,.01)))
  }


# Group Features
dengroup.ppm(features, ppm, rtwid, minlength)


# Warpgroup To Find Consensus Features
library(warpgroup)
for (group in groups) {
  warpgroup(ps, eic.mat.s, sc.aligned.lim = 5, min.peaks = 1, tw = "dtw", pct.pad = 0.1)
  }
