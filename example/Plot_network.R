library(igraph)

inf_total <- read.table("termList_to_keep.txt", header = FALSE, sep="\t")
term_total <- as.character(inf_total[, 1])
term_total <- c(term_total, "Simvastatin")


inf <- read.table("matrix.txt",header = TRUE,row.names="Substances", sep="\t")
termDocMatrix <- as.matrix(inf)
# change it to a Boolean matrix
termDocMatrix[termDocMatrix>=1] <- 1
# transform into a term-term adjacency matrix
termMatrix <- termDocMatrix %*% t(termDocMatrix)

# build a graph from the above matrix
g <- graph.adjacency(termMatrix, weighted=T, mode = "undirected")
# remove loops
g <- simplify(g)
# set labels and degrees of vertices
V(g)$label <- V(g)$name
V(g)$degree <- degree(g)
inf1 <- read.table("vertex_attribute.txt",header = TRUE,sep="\t")
targetName <- inf1$Substances[inf1$Type == "Target"]
V(g)$Type <- as.character(inf1$Type[match(V(g)$label,inf1$Substances)])
V(g)$color <- V(g)$Type
V(g)$color <- gsub("Compounds","yellow",V(g)$color)
V(g)$color <- gsub("Proteins","red",V(g)$color)
V(g)$color <- gsub("Target","blue",V(g)$color)

# Check whether all the terms in term_total exit in V(g)$name
for (termTmp in term_total) {
  if (!(termTmp %in% V(g)$name)) {
    #print(termTmp)
    invisible()
  }
}

unmapped_list <- c()
for (termTmp in V(g)$name) {
  if (!(termTmp %in% term_total)) {
    unmapped_list <- c(unmapped_list, termTmp)
  }
}

# Remove the vertices in unmapped_list from g
g <- g - unmapped_list
inf1 <- inf1[!(inf1$Substances %in% unmapped_list), ]

for (i in inf1$Substances[inf1$Type == "Compounds"]) {
  for (k in inf1$Substances[inf1$Type == "Compounds"]) {
    g[i,k] <- NULL
  }  
}

#Remove edges linking compounds and Ibuprofen
for (i in inf1$Substances[inf1$Type == "Compounds"]) {
  for (k in inf1$Substances[inf1$Type == "Target"]) {
    g[i,k] <- NULL
  }  
}

#Delete vertex whose degree = 0
delet_vertex <- V(g)[degree(g)<1]
#write.csv(names(delet_vertex), "removed_term.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
removed_from_paperList = c()
for (termTmp in names(delet_vertex)) {
  if (termTmp %in%term_total) {
    removed_from_paperList <- c(removed_from_paperList, termTmp)
  }
}

#write.csv(removed_from_paperList, "removed_term_from_paper.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

g <- delete.vertices(g, delet_vertex)

removeList <- c()

# Rmove the vertices that are not connected with target
for (i in V(g)[V(g)$color=='yellow']) {
  for (k in V(g)[V(g)$color=='blue']) {
    if (is.infinite(distances(g, v=i, to=k))) {
      removeList <- c(removeList,i)
    }
  }
}
for (i in V(g)[V(g)$color=='red']) {
  for (k in V(g)[V(g)$color=='blue']) {
    if (is.infinite(distances(g, v=i, to=k))) {
      removeList <- c(removeList,i)
    }
  }
}

if (!(is.null(removeList))) {
  remainList <- V(g)$name[-removeList]
  g <- induced.subgraph(g,V(g)$name%in%remainList)
}

#write.csv(V(g)$name, "kept_term.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

for (termTmp in term_total) {
  if (!(termTmp %in% V(g)$name)) {
    #print(termTmp)
    invisible()
  }
}

# set seed to make the layout reproducible
set.seed(3952)
pdf(paste(targetName,"_network.pdf"))
plot(g, layout=layout.kamada.kawai,edge.color="black")
#tkplot(g, layout=layout.kamada.kawai,edge.color="black")
dev.off()