brainsync <- function(

  data_wd = "",

  filepattern = "",
  sep = ",",

  channel.lim = c(),

  freq.lim = c(),
  freq.steps,

  freq.range.min = NA,
  freq.range.max = NA,

  itermax = 10000,
  nstart = 1000,

  alpha_chitest = 0.05,
  pval_adjustment = FALSE,
  method = p.adjust.methods,

  detailed_running_output = FALSE,

  writechisquaretable = TRUE,
  comparisonname = ""

) {

  required_packages()


  # Input split

  freq.min <- freq.lim[1]; freq.max <-freq.lim[2]
  a <- channel.lim[1]; b <- channel.lim[2]
  # Detailed Output

  if (detailed_running_output == TRUE) {
    print_clusters = TRUE
    print_fx = TRUE
    printfilenumber = TRUE
    } else {print_clusters = FALSE; print_fx = FALSE;   printfilenumber = FALSE}

  # Fx Generation and general properties
  freq.stepsb = 1/freq.steps
  Fx <- seq(from = freq.min, to = freq.max, by = freq.stepsb)

  d1 = length(Fx)



  if (is.na(freq.range.min) == TRUE) {
    freq.range.min <- freq.min
 }

  if (is.na(freq.range.max) == TRUE) {
    freq.range.max <- freq.max
  }

  freq.range.min.calc <- which(Fx %in% c(freq.range.min))
  freq.range.max.calc <- which(Fx %in% c(freq.range.max))
  if (print_fx == TRUE) {
    print("Minimal and Maximal Frequencies of interest")
    print(paste(freq.range.min.calc,freq.range.max.calc))

  }
  Fx <- Fx[freq.range.min.calc:freq.range.max.calc]

  if (print_fx == TRUE) {
    print("Data Frequency Table")
    print(Fx)
  }

  # Read data

  data_files <- list.files(path = data_wd, pattern = filepattern, full.names = TRUE)
  lengthdf <- length(data_files)

  # First Loop

  clusters <- list();
  distrib <- vector()

  for (k in 1:lengthdf) {

    data <- read.csv(data_files[k], sep = sep, header = FALSE)
    data <- data[freq.range.min.calc:freq.range.max.calc,a:b]
    res <- cor(data);

    cor(data$V1, data$V2)

    # dendogram fir dissimilarity
    corloads = cor(data)
    dissimilarity = 1 - corloads
    distance = as.dist(dissimilarity)

    # find optimal number of clusters
    clus <- hclust(distance)
    op_k <- kgs(clus, distance, maxclus = (dim(data)[2]-1))

    rows = as.integer(names(op_k[which(op_k == min(op_k))]));

    kmeans_class <- kmeans(res, rows, iter.max = itermax, nstart = nstart);
    kmc <- kmeans_class$cluster; kmx <- unname(kmc);

    clusters[[k]] = kmx

    if (printfilenumber == TRUE) {
      print("File Number")
      print(k)}


    progressbar(i = k, n = lengthdf, input = "Clustered Data Samples")
  }

  if (print_clusters == TRUE) {
    print("Clusters")
    print(clusters)}

  channels <- colnames(data);

  df <- data.frame(channels = channels, xx = clusters);
  dimdf <- dim(df)[2]
  for (i in 2:dimdf) {j = i-1; name = paste("P",j, sep = "");
  colnames(df)[i] = name}

  dat <- df
  dimy2 <- dim(dat)[2];
  dimy1 <- dim(dat)[1]
  dat <- dat[2:dimy2]

  dimy2 <- dimy2-1

  vec <- vector()
  full <- list()

  #

  for (k in 1:dimy2) {          # for all participants

    o = 1;

    for (n in 1:dimy1) {        # channel in position 1

      for (m in 1:dimy1) {    # channel in position 2

        if (dat[n,k] == dat[m,k]) {vec[o] = 1} else {vec[o] = 0}

        o = o + 1;

      }

      full[[k]] = vec;

    }

    progressbar(i = k, n = dimy2, input = "Participant-level Channel-Grouping")
  }

  naming <- vector();
  o = 1; for (n in 1:dimy1) {        # channel in position 1

    for (m in 1:dimy1) {    # channel in position 2

      naming[o] = paste("V",n,"_V",m, sep="")
      o = o + 1;
    }

  }; head(naming)
  l1 <- length(naming) #-1

  chtest <- data.frame(
    group = naming[1:l1],
    indic = full
  )

  dimtest <- dim(chtest)[2]
  for (i in 2:dimtest) {j = i-1; name = paste("P",j, sep = "");colnames(chtest)[i] = name}


  ctest_red <- data.frame(indic = full)
  dimtest <- dim(ctest_red)[2]
  for (i in 1:dimtest) {j = i; name = paste("P",j, sep = "");colnames(ctest_red)[i] = name}
    xy <- ctest_red

  dimxy1 <- dim(xy)[1]; dimxy2 <- dim(xy)[2]

  vecxy0 <- vector()
  vecxy1 <- vector()


  o = 1; for (i in 1:dimxy1) {

    count0 = 0; count1 = 0;

    for (j in 1:dimxy2) {

      if (xy[i,j] == 1) {count1 = count1 + 1} else {count0 = count0 + 1}

    }

    vecxy0[o] = count0;
    vecxy1[o] = count1;

    o = o + 1

    progressbar(i = i, n = dimxy1, input = "Between-Participant Channel-Grouping")

  }

  # Chi Square test

  pvals <- vector()

  for (i in 1:length(vecxy0)) {
    table <- data.frame(yes = vecxy0[i], no = vecxy1[i]);
    chi <- chisq.test(table) #here for two sided chi square: sign. indicates unequal group
    # unequal group means in one of the two groups there is sign more of either yes or no
    pvals[i] = chi$p.value
  };


  alpha <- alpha_chitest

  names2 <- naming[1:l1]
  df <- data.frame(comp = names2, pval = pvals); df

  df.pval.1 <- df

  names3 <- vector(); pvals2 <- vector(); yes = vector(); no <- vector()
  #print(df$pval)
  if (pval_adjustment == TRUE) {p2 <- p.adjust(p = df$pval, method = method)} else {
  }; #print(p2)

  o = 1; for (i in 1:length(names2)) {

    p <- p2





    if (p[i] <= alpha) {
      names3[o] = names2[i];
      pvals2[o] = df$pval[i];
      yes[o] = vecxy1[i];
      no[o] = vecxy0[i];

      o = o+1;
    } else {}

  }

  df2 <- data.frame(
    col <- names3,
    pvals3 <- pvals2,
    yes = yes,
    no = no
  );


  colnames(df2)[1] = "comp"
  colnames(df2)[2] = "pval4"
  head(df2)

  dimcomp <- length(df2$comp)
  namesx <- vector(); pvalsx <- vector(); yesx = vector(); nox <- vector()

  p = 1; for (i in 1:dimcomp) {

    if (df2$yes[i] >= df2$no[i]) {
      namesx[p] = df2$comp[i];
      pvalsx[p] = df2$pval4[i];
      yesx[p] = df2$yes[i];
      nox[p] = df2$no[i];
      p = p+1;
    }

  }

  df3 <- data.frame(
    col <- namesx,
    pval <- pvalsx,
    yes = yesx,
    no = nox
  );

  colnames(df3)[1] = "Comp"
  colnames(df3)[2] = "pVal"
 # head(df3)

  o = 1; for (i in 1:dim(df3)[1]) {
    if (df3$yes[i] == 86) {o = o+1}
  };o-1

  #print(df3)


  namesx4 <- vector(); pvalsx4 <- vector(); yesx4 = vector(); nox4 <- vector()
  dim4 <- dim(df3)[1]
  p = 1; for (i in 1:dim4) {

    if (df3$yes[i] == 86) {
    } else {
      namesx4[p] = df3$Comp[i];
      pvalsx4[p] = df3$pVal[i];
      yesx4[p] = df3$yes[i];
      nox4[p] = df3$no[i];
      p = p+1;
    }

  }


  df4 <- data.frame(
    col <- namesx4,
    pval <- pvalsx4,
    yes = yesx4,
    no = nox4
  );

  colnames(df4)[1] = "Comp"
  colnames(df4)[2] = "pVal"



data <- df4

# Graph Building

i = 28; x1 <- vector(); x2 <- vector()

for (i in 1:dim(data)[1]) {

  test <- data[i,]; test
  x <- nchar(test[1])

  if (substring(test[1], 3, 3) == "_") {
    # for 1 - 9
    x1[i] <- substring(test[1], 2, 2)
    # 2nd V at position 4
    x2[i] <- substring(test[1], 5)
  } else {
    if (substring(test[1], 4, 4) == "_") {
      # for 10 - 99
      x1[i] <- substring(test[1], 2, 3)
      x2[i] <- substring(test[1], 6)
    }  else {
      if (substring(test[1], 5, 5) == "_") {
        # for 100+
        x1[i] <- substring(test[1], 2, 4)
        x2[i] <- substring(test[1], 7)
      }
    }
  }

}

df5 <- data.frame(
  Comp <- df4$Comp,
  x1 <- x1,
  x2 <- x2,
  pVal <- df4$pVal,
  yes = df4$yes,
  no = df4$no
);

colnames(df5)[1] = "Comp"
colnames(df5)[4] = "pVal"
colnames(df5)[2] = "C1"
colnames(df5)[3] = "C2"



df6 <- data.frame(
  x1 <- x1,
  x2 <- x2,
  pVal <- df4$pVal,
  yes = df4$yes,
  no = df4$no
);

colnames(df6)[3] = "pVal"
colnames(df6)[1] = "C1"
colnames(df6)[2] = "C2"

df6$C1 <- as.numeric(df6$C1)
df6$C2 <- as.numeric(df6$C2)

df7 = df6

if (writechisquaretable == TRUE) {
  write.table(df6, paste(comparisonname, ".csv", sep = ""),row.names = FALSE,quote = FALSE,sep = ';')
}

c1x <- vector();  c2x <- vector(); pvalx <- vector();

o = 1; for (i in 1:dim(df7)[1]) {

  if (df7$C1[i] >= df7$C2[i]) {} else {

    c1x[o] <- df7$C1[i];
    c2x[o] <- df7$C2[i];
    pvalx[o] <- df7$pVal[i];

    o = o + 1;
  }
};

df8 <- data.frame(
  C1 = c1x,
  C2 = c2x,
  pVal = pvalx
)


c1u <- unique(df8$C1)
c2u <- unique(df8$C2)
c12u <- c(c1u, c2u)
c12u2 <- unique(c12u)

nodes = cbind('id' = c12u2)

links = cbind('from'=df8$C1, 'to'=df8$C2, 'weight'=rep(1,dim(df8)[1]))

net = graph_from_data_frame(links,vertices = nodes,directed = F)

plot(net, vertex.color="gold", vertex.label.color="black", vertex.label.cex=0.8,
     vertex.label.dist=0, edge.curved=0)

# save plot
pdf(paste(comparisonname, ".pdf", sep =""), height = 8, width = 8)
plot(net, vertex.color="gold", vertex.label.color="black", vertex.label.cex=0.8,
     vertex.label.dist=0, edge.curved=0)
dev.off()

graphl <- split(names(V(net)), components(net)$membership)


print("Function successfully run.")

outputlist <- list("ChiSquareChannelComparisons" = df4, "ClusteredGraphs" = graphl, "UnadjustedPValues" = df.pval.1)



} # Overall Function
